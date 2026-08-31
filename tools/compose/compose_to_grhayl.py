#!/usr/bin/env python3
"""Convert one fixed CompOSE HDF5 schema into a regularized GRHayL table."""

import argparse
from collections import namedtuple
import itertools
import json
import math
import os
from pathlib import Path
import secrets
import stat
import sys

import h5py
import numpy as np


MEV_TO_ERG = 1.602176634e-6
C_SQUARED = 8.987551787368176e20
MEV_TO_GRAM = MEV_TO_ERG / C_SQUARED
INT_MAX = 2_147_483_647
MAX_PLANE_CELLS = 65_536
MAX_CONTROL_BYTES = 4_096
MAX_CONTROL_TOKENS = 128
ROOT_RELATIVE_TOLERANCE = 1.0e-10
LOW_BULK_FRACTION = 1.0e-14
HIGH_BULK_FRACTION = 1.0 - 1.0e-12
NODE_CLOSURE_TOLERANCE = 1.0e-12

THERMO_IDS = (1, 2, 3, 4, 5, 10, 11, 12, 13, 21)
PAIR_IDS = (10, 11, 4002)
QUADRUPLE_IDS = (999,)
INPUT_FILES = ("eoscompose.h5", "eos.parameters", "eos.quantities", "eos.thermo")
OUTPUT_FIELDS = (
    "Abar", "Xa", "Xh", "Xn", "Xp", "Zbar", "cs2", "dedt",
    "dpderho", "dpdrhoe", "entropy", "gamma", "logenergy",
    "logpress", "mu_e", "mu_n", "mu_p", "muhat", "munu",
)


class ConversionError(Exception):
    """Raised when input, regularization, or publication fails."""


def _remove_owned_file(path):
    try:
        os.unlink(path)
    except FileNotFoundError:
        pass


def _cleanup(actions, primary_error=None):
    failures = []
    for description, action in actions:
        try:
            action()
        except Exception as error:
            failures.append(f"{description} failed: {error}")
    if not failures:
        return
    if primary_error is None:
        raise ConversionError("; ".join(failures))
    previous = getattr(primary_error, "cleanup_failures", ())
    primary_error.cleanup_failures = (*previous, *failures)


Profile = namedtuple(
    "Profile",
    (
        "name program_banner manual_version table_id nrho nt nye nb_min nb_max "
        "temperature_min temperature_max ye_min ye_max thermo_ids pair_ids "
        "quadruple_ids thermo_index_name m_n_mev m_p_mev i_l "
        "m_ref_rho_mev m_ref_mu_mev rho_axis temperature_axis ye_axis "
        "beta_equilibrium fixed_entropy nn"
    ),
)


PROFILES = {
    "sro-sly4-sna-141-regularized-v1": Profile(
        name="sro-sly4-sna-141-regularized-v1",
        program_banner="2.17",
        manual_version="3.01",
        table_id=141,
        nrho=391,
        nt=163,
        nye=66,
        nb_min=6.3095733426599830e-13,
        nb_max=6.3095733426599834,
        temperature_min=1.0e-3,
        temperature_max=251.18864315095823,
        ye_min=5.0e-3,
        ye_max=0.655,
        thermo_ids=THERMO_IDS,
        pair_ids=PAIR_IDS,
        quadruple_ids=QUADRUPLE_IDS,
        thermo_index_name="index_thermo",
        m_n_mev=939.5654,
        m_p_mev=938.2754,
        i_l=1,
        m_ref_rho_mev=939.5654,
        m_ref_mu_mev=939.5654,
        rho_axis="log",
        temperature_axis="log",
        ye_axis="linear",
        beta_equilibrium=False,
        fixed_entropy=False,
        nn=None,
    ),
    "test-analytic": Profile(
        name="test-analytic",
        program_banner="2.17-test",
        manual_version="3.01-test",
        table_id=0,
        nrho=5,
        nt=4,
        nye=3,
        nb_min=1.0e-5,
        nb_max=1.0e-1,
        temperature_min=1.0e-2,
        temperature_max=1.0e1,
        ye_min=0.2,
        ye_max=0.4,
        thermo_ids=THERMO_IDS,
        pair_ids=PAIR_IDS,
        quadruple_ids=QUADRUPLE_IDS,
        thermo_index_name="index",
        m_n_mev=939.5654,
        m_p_mev=938.2754,
        i_l=1,
        m_ref_rho_mev=939.5654,
        m_ref_mu_mev=939.5654,
        rho_axis="log",
        temperature_axis="log",
        ye_axis="linear",
        beta_equilibrium=False,
        fixed_entropy=False,
        nn=None,
    ),
}


REGULARIZATION_POLICY = {
    "absent_heavy_nucleus": "A=1,Z=0,Xh=0 with three-species L2 projection",
    "bulk_fraction_bounds": [LOW_BULK_FRACTION, HIGH_BULK_FRACTION],
    "composition_projection": "minimum-L2 active-set",
    "density_derivative": "second-order nonuniform",
    "energy_shift_floor_mev": 20.0,
    "entropy_anchor": "max(raw_S(Tmin),0)",
    "entropy_integrator": "logarithmic-temperature mean",
    "isotonic_weights": "equal",
    "heavy_charge_roundoff": "clip Z to [0,A] within 1e-8*max(1,A)",
    "log_step": "max(64*float64-epsilon*scale,8e-10/ln(10))",
    "node_closure_tolerance": NODE_CLOSURE_TOLERANCE,
    "offgrid_charge_closure": "diagnostic only; current consumers use Xn/Xp",
    "temperature_derivative": "PCHIP weighted harmonic",
}


def _checked_cell_count(nrho, nt, nye):
    dimensions = (nrho, nt, nye)
    if any(not isinstance(value, (int, np.integer)) or value <= 0 for value in dimensions):
        raise ConversionError("every grid dimension must be a positive integer")
    cells = int(nrho) * int(nt) * int(nye)
    if 20 * cells > INT_MAX:
        raise ConversionError("20*Ncells exceeds INT_MAX")
    return cells


def _checked_plane_cells(nrho, nt):
    cells = int(nrho) * int(nt)
    if nrho <= 0 or nt <= 0 or cells > MAX_PLANE_CELLS:
        raise ConversionError("one Y_e plane exceeds the 65,536-cell plane bound")
    return cells


def _regularization_log_step(values):
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0 or not np.all(np.isfinite(values)):
        raise ConversionError("regularization values must be finite")
    scale = max(1.0, float(np.max(np.abs(values))))
    floating_margin = 64.0 * np.finfo(np.float64).eps * scale
    inverse_margin = 8.0 * ROOT_RELATIVE_TOLERANCE / math.log(10.0)
    return max(floating_margin, inverse_margin)


def _strict_isotonic(values, minimum_step):
    values = np.asarray(values, dtype=np.float64)
    if values.ndim != 1 or values.size < 2 or not np.all(np.isfinite(values)):
        raise ConversionError("strict isotonic input must be a finite one-dimensional array")
    if not math.isfinite(minimum_step) or minimum_step <= 0.0:
        raise ConversionError("strict isotonic minimum step must be positive and finite")

    shifted = values - minimum_step * np.arange(values.size, dtype=np.float64)
    means = []
    weights = []
    for value in shifted:
        means.append(float(value))
        weights.append(1)
        while len(means) > 1 and means[-2] > means[-1]:
            weight = weights[-2] + weights[-1]
            mean = (weights[-2] * means[-2] + weights[-1] * means[-1]) / weight
            means[-2:] = [mean]
            weights[-2:] = [weight]

    result = np.empty_like(values)
    offset = 0
    for mean, weight in zip(means, weights):
        result[offset:offset + weight] = mean
        offset += weight
    result += minimum_step * np.arange(values.size, dtype=np.float64)
    if not np.all(np.diff(result) > 0.0):
        raise ConversionError("strict isotonic projection did not produce increasing values")
    return result


def _monotone_derivative(axis, values):
    axis = np.asarray(axis, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    if axis.ndim != 1 or values.ndim != 1 or axis.shape != values.shape or axis.size < 3:
        raise ConversionError("monotone derivative requires matching one-dimensional arrays")
    if not np.all(np.isfinite(axis)) or not np.all(np.diff(axis) > 0.0):
        raise ConversionError("derivative axis must be finite and strictly increasing")
    intervals = np.diff(axis)
    secants = np.diff(values) / intervals
    if not np.all(np.isfinite(secants)) or not np.all(secants > 0.0):
        raise ConversionError("derivative values must be finite and strictly increasing")

    derivative = np.empty_like(values)
    derivative[0] = secants[0]
    derivative[-1] = secants[-1]
    for index in range(1, values.size - 1):
        left = intervals[index - 1]
        right = intervals[index]
        first_weight = 2.0 * right + left
        second_weight = right + 2.0 * left
        derivative[index] = (first_weight + second_weight) / (
            first_weight / secants[index - 1] + second_weight / secants[index]
        )
    return derivative


def _energy_shift_from_extrema(minimum, maximum_absolute):
    if not math.isfinite(minimum) or not math.isfinite(maximum_absolute):
        raise ConversionError("energy extrema must be finite")
    scale = max(1.0, maximum_absolute)
    margin = 64.0 * np.finfo(np.float64).eps * scale
    shift = 20.0 if minimum + 20.0 > margin else -minimum + margin
    if not math.isfinite(shift):
        raise ConversionError("energy shift is not finite")
    return shift


def _project_bulk(first_contribution, second_contribution, enthalpy_density):
    first_contribution = np.asarray(first_contribution, dtype=np.float64)
    second_contribution = np.asarray(second_contribution, dtype=np.float64)
    enthalpy_density = np.asarray(enthalpy_density, dtype=np.float64)
    if first_contribution.shape != second_contribution.shape or first_contribution.shape != enthalpy_density.shape:
        raise ConversionError("bulk projection arrays must have matching shapes")
    if not np.all(np.isfinite(enthalpy_density)) or not np.all(enthalpy_density > 0.0):
        raise ConversionError("bulk projection requires positive finite enthalpy density")

    raw_bulk = first_contribution + second_contribution
    lower = LOW_BULK_FRACTION * enthalpy_density
    upper = HIGH_BULK_FRACTION * enthalpy_density
    low_count = int(np.count_nonzero(raw_bulk < lower))
    high_count = int(np.count_nonzero(raw_bulk > upper))
    bulk = np.clip(raw_bulk, lower, upper)
    first = np.clip(0.5 * (bulk + first_contribution - second_contribution), 0.0, bulk)
    second = bulk - first
    return first, second, low_count, high_count


def _project_composition(raw, heavy_charge_fraction, ye):
    raw = np.asarray(raw, dtype=np.float64)
    heavy_charge_fraction = np.asarray(heavy_charge_fraction, dtype=np.float64)
    ye = np.asarray(ye, dtype=np.float64)
    if raw.ndim < 2 or raw.shape[0] != 4 or raw.shape[1:] != heavy_charge_fraction.shape or ye.shape != heavy_charge_fraction.shape:
        raise ConversionError("composition projection arrays have incompatible shapes")
    if not np.all(np.isfinite(raw)) or not np.all(np.isfinite(heavy_charge_fraction)) or not np.all(np.isfinite(ye)):
        raise ConversionError("composition projection requires finite values")
    if np.any(heavy_charge_fraction < 0.0) or np.any(heavy_charge_fraction > 1.0) or np.any(ye < 0.0) or np.any(ye > 1.0):
        raise ConversionError("composition constraints are infeasible")

    original_shape = heavy_charge_fraction.shape
    raw_flat = raw.reshape(4, -1)
    q_flat = heavy_charge_fraction.reshape(-1)
    ye_flat = ye.reshape(-1)
    point_count = q_flat.size
    best = np.zeros_like(raw_flat)
    best_cost = np.full(point_count, np.inf)
    charge_coefficients = (
        np.zeros(point_count), np.ones(point_count),
        np.full(point_count, 0.5), q_flat,
    )

    for support_size in range(2, 5):
        for support in itertools.combinations(range(4), support_size):
            sum_charge = sum(charge_coefficients[index] for index in support)
            sum_charge_squared = sum(charge_coefficients[index] ** 2 for index in support)
            determinant = support_size * sum_charge_squared - sum_charge ** 2
            nonsingular = determinant > 64.0 * np.finfo(np.float64).eps
            if not np.any(nonsingular):
                continue

            raw_mass = sum(raw_flat[index] for index in support)
            raw_charge = sum(charge_coefficients[index] * raw_flat[index] for index in support)
            mass_residual = raw_mass - 1.0
            charge_residual = raw_charge - ye_flat
            lambda_mass = np.zeros(point_count)
            lambda_charge = np.zeros(point_count)
            lambda_mass[nonsingular] = (
                mass_residual[nonsingular] * sum_charge_squared[nonsingular]
                - charge_residual[nonsingular] * sum_charge[nonsingular]
            ) / determinant[nonsingular]
            lambda_charge[nonsingular] = (
                support_size * charge_residual[nonsingular]
                - sum_charge[nonsingular] * mass_residual[nonsingular]
            ) / determinant[nonsingular]

            candidate = np.zeros_like(raw_flat)
            feasible = nonsingular.copy()
            for index in support:
                values = raw_flat[index] - lambda_mass - charge_coefficients[index] * lambda_charge
                candidate[index] = values
                feasible &= values >= 0.0
            candidate[candidate < 0.0] = 0.0
            cost = np.sum((candidate - raw_flat) ** 2, axis=0)
            update = feasible & (cost < best_cost)
            best[:, update] = candidate[:, update]
            best_cost[update] = cost[update]

    if np.any(~np.isfinite(best_cost)):
        raise ConversionError("composition constraints are infeasible")
    projected = best.reshape((4,) + original_shape)
    mass_residual = np.sum(projected, axis=0) - 1.0
    charge_residual = projected[1] + 0.5 * projected[2] + heavy_charge_fraction * projected[3] - ye
    if np.max(np.abs(mass_residual)) > NODE_CLOSURE_TOLERANCE or np.max(np.abs(charge_residual)) > NODE_CLOSURE_TOLERANCE:
        raise ConversionError("composition projection did not meet node closure tolerance")
    return projected


def _read_control_tokens(path):
    try:
        with path.open("rb") as stream:
            payload = stream.read(MAX_CONTROL_BYTES + 1)
    except (OSError, UnicodeError) as error:
        raise ConversionError(f"control file {path.name} could not be read: {error}") from error
    if len(payload) > MAX_CONTROL_BYTES:
        raise ConversionError(f"control file {path.name} exceeds byte limit")
    try:
        text = payload.decode("ascii")
    except UnicodeError as error:
        raise ConversionError(f"control file {path.name} could not be read: {error}") from error
    tokens = []
    for line in text.splitlines():
        tokens.extend(line.split("#", 1)[0].split())
        if len(tokens) > MAX_CONTROL_TOKENS:
            raise ConversionError(f"control file {path.name} exceeds token limit")
    return tokens


def _require_regular_file(path):
    try:
        mode = path.lstat().st_mode
    except OSError as error:
        raise ConversionError(f"input file {path.name} is missing: {error}") from error
    if not stat.S_ISREG(mode):
        raise ConversionError(f"input file {path.name} must be a regular file, not a symlink")


def _validate_controls(input_dir, profile):
    for filename in INPUT_FILES:
        _require_regular_file(input_dir / filename)

    parameter_tokens = _read_control_tokens(input_dir / "eos.parameters")
    expected_parameters = (
        1.0, 1.0, 1.0,
        float(profile.beta_equilibrium), float(profile.fixed_entropy), 1.0,
        profile.temperature_min, profile.nb_min, profile.ye_min,
        profile.temperature_max, profile.nb_max, profile.ye_max,
        float(profile.nt), float(profile.nrho), float(profile.nye),
        1.0, 1.0, 0.0,
    )
    try:
        parameter_values = tuple(float(token) for token in parameter_tokens)
    except ValueError as error:
        raise ConversionError("eos.parameters control schema is not numeric") from error
    if len(parameter_values) != len(expected_parameters) or any(
        not math.isclose(actual, expected, rel_tol=2.0e-15, abs_tol=0.0)
        for actual, expected in zip(parameter_values, expected_parameters)
    ):
        raise ConversionError("eos.parameters control schema does not match fixed profile")

    quantity_tokens = _read_control_tokens(input_dir / "eos.quantities")
    expected_quantities = (
        10, 0, 0, *profile.thermo_ids,
        3, 1, *profile.pair_ids, *profile.quadruple_ids,
        0, 0, 2,
    )
    try:
        quantity_values = tuple(int(token) for token in quantity_tokens)
    except ValueError as error:
        raise ConversionError("eos.quantities control schema is not integral") from error
    if quantity_values != expected_quantities:
        raise ConversionError("eos.quantities control schema does not match fixed profile")

    try:
        with (input_dir / "eos.thermo").open("r", encoding="ascii") as stream:
            header_line = stream.readline(513)
    except (OSError, UnicodeError) as error:
        raise ConversionError(f"control file eos.thermo could not be read: {error}") from error
    if len(header_line) == 513 and not header_line.endswith("\n"):
        raise ConversionError("eos.thermo control header is overlong")
    thermo_tokens = header_line.split("#", 1)[0].split()
    if len(thermo_tokens) != 3:
        raise ConversionError("eos.thermo control header is incomplete")
    try:
        thermo_header = (float(thermo_tokens[0]), float(thermo_tokens[1]), int(thermo_tokens[2]))
    except ValueError as error:
        raise ConversionError("eos.thermo control header is malformed") from error
    expected_header = (profile.m_n_mev, profile.m_p_mev, profile.i_l)
    if not np.allclose(thermo_header[:2], expected_header[:2], rtol=2.0e-15, atol=2.0e-15) or thermo_header[2] != expected_header[2]:
        raise ConversionError("eos.thermo control header does not match fixed profile")


def _require_hard_link(group, name):
    try:
        link = group.get(name, getlink=True)
    except Exception as error:
        raise ConversionError(f"HDF5 link {group.name}/{name} could not be inspected") from error
    if not isinstance(link, h5py.HardLink):
        raise ConversionError(f"HDF5 link {group.name}/{name} must be a hard link")


def _validate_dataset(group, name, shape, dtype):
    _require_hard_link(group, name)
    value = group[name]
    if not isinstance(value, h5py.Dataset):
        raise ConversionError(f"HDF5 dataset {group.name}/{name} is not a dataset")
    if value.is_virtual:
        raise ConversionError(f"HDF5 dataset {value.name} uses virtual storage")
    creation = value.id.get_create_plist()
    if creation.get_external_count() != 0:
        raise ConversionError(f"HDF5 dataset {value.name} uses external storage")
    if value.chunks is not None or value.compression is not None or value.shuffle or value.fletcher32 or value.scaleoffset is not None:
        raise ConversionError(f"HDF5 dataset {value.name} has unsupported layout or filters")
    if creation.get_layout() != h5py.h5d.CONTIGUOUS:
        raise ConversionError(f"HDF5 dataset {value.name} has unsupported layout")
    if value.shape != shape:
        raise ConversionError(f"HDF5 dataset {value.name} has wrong shape or dimension order")
    if value.dtype != np.dtype(dtype):
        raise ConversionError(f"HDF5 dataset {value.name} has wrong dtype")
    return value


def _validate_attribute(group, name, expected, identifier_kind=None):
    if name not in group.attrs:
        raise ConversionError(f"HDF5 attribute {group.name}/{name} is missing")
    value = group.attrs[name]
    if not isinstance(value, np.ndarray) or value.shape != (1,) or value.dtype != np.dtype("int32"):
        raise ConversionError(f"HDF5 attribute {group.name}/{name} is malformed")
    if int(value[0]) != expected:
        if identifier_kind is not None:
            raise ConversionError(f"{identifier_kind} identifier count does not match fixed profile")
        raise ConversionError(f"HDF5 attribute {group.name}/{name} has wrong value")


def _reject_hard_link_aliases(objects):
    addresses = {}
    for value in objects:
        address = h5py.h5o.get_info(value.id).addr
        if address in addresses:
            raise ConversionError(
                f"HDF5 hard-link aliases {addresses[address]} and {value.name} are unsupported"
            )
        addresses[address] = value.name


def _axis_tolerance(values):
    return 256.0 * np.finfo(np.float64).eps * max(1.0, float(np.max(np.abs(values))))


def _validate_axis(values, mode, expected_min, expected_max, name):
    values = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(values)) or not np.all(values > 0.0 if mode == "log" else values >= 0.0):
        raise ConversionError(f"{name} axis must contain finite supported values")
    transformed = np.log10(values) if mode == "log" else values
    differences = np.diff(transformed)
    tolerance = _axis_tolerance(transformed)
    if not np.all(differences > 0.0) or not np.allclose(differences, differences[0], rtol=0.0, atol=tolerance):
        raise ConversionError(f"{name} axis must be strictly increasing and uniform")
    expected_lower = math.log10(expected_min) if mode == "log" else expected_min
    expected_upper = math.log10(expected_max) if mode == "log" else expected_max
    endpoint_tolerance = 512.0 * np.finfo(np.float64).eps * max(1.0, abs(expected_lower), abs(expected_upper))
    if not math.isclose(float(transformed[0]), expected_lower, rel_tol=0.0, abs_tol=endpoint_tolerance):
        raise ConversionError(f"{name} axis lower endpoint does not match profile")
    if not math.isclose(float(transformed[-1]), expected_upper, rel_tol=2.0e-14, abs_tol=endpoint_tolerance):
        raise ConversionError(f"{name} axis upper endpoint does not match profile")


def _open_validated_input(path, profile):
    try:
        h5 = h5py.File(path, "r")
    except (OSError, ValueError) as error:
        raise ConversionError(f"input HDF5 file could not be opened: {error}") from error
    try:
        expected_root = {
            "Composition_pairs", "Composition_quadrupels", "Parameters",
            "Thermo_qty", "metadata",
        }
        if set(h5.keys()) != expected_root:
            raise ConversionError("HDF5 root groups do not match fixed profile")
        for name in expected_root:
            _require_hard_link(h5, name)
            if not isinstance(h5[name], h5py.Group):
                raise ConversionError(f"HDF5 root object {name} is not a group")

        parameters = h5["Parameters"]
        quantities = h5["Thermo_qty"]
        pairs = h5["Composition_pairs"]
        quadruples = h5["Composition_quadrupels"]
        metadata = h5["metadata"]
        if set(parameters.keys()) != {"nb", "t", "yq"}:
            raise ConversionError("Parameters dataset set does not match fixed schema")
        if set(quantities.keys()) != {profile.thermo_index_name, "thermo"}:
            raise ConversionError("Thermo_qty dataset set does not match fixed schema")
        if set(pairs.keys()) != {"index_yi", "yi"}:
            raise ConversionError("Composition_pairs dataset set does not match fixed schema")
        if set(quadruples.keys()) != {"index_av", "aav", "nav", "yav", "zav"}:
            raise ConversionError("Composition_quadrupels dataset set does not match fixed schema")
        if set(metadata.keys()) or set(metadata.attrs.keys()) != {"date", "time"}:
            raise ConversionError("metadata group does not match official writer schema")
        objects = [parameters, quantities, pairs, quadruples, metadata]
        for group in objects.copy():
            for name in group:
                _require_hard_link(group, name)
                objects.append(group[name])
        _reject_hard_link_aliases(objects)
        for name in ("date", "time"):
            value = metadata.attrs[name]
            if not isinstance(value, (bytes, np.bytes_)):
                raise ConversionError(f"metadata attribute {name} is malformed")

        if set(parameters.attrs.keys()) != {"pointsnb", "pointst", "pointsyq", "tabulation_scheme"}:
            raise ConversionError("Parameters attribute set does not match fixed schema")
        if set(quantities.attrs.keys()) != {"pointsqty"} or set(pairs.attrs.keys()) != {"pointspairs"} or set(quadruples.attrs.keys()) != {"pointsav"}:
            raise ConversionError("HDF5 group attribute set does not match fixed schema")
        _validate_attribute(parameters, "pointsnb", profile.nrho)
        _validate_attribute(parameters, "pointst", profile.nt)
        _validate_attribute(parameters, "pointsyq", profile.nye)
        _validate_attribute(parameters, "tabulation_scheme", 1)
        _validate_attribute(quantities, "pointsqty", len(profile.thermo_ids), "thermodynamic")
        _validate_attribute(pairs, "pointspairs", len(profile.pair_ids), "pair")
        _validate_attribute(quadruples, "pointsav", len(profile.quadruple_ids), "quadruple")

        nb = _validate_dataset(parameters, "nb", (profile.nrho,), "float64")[...]
        temperature = _validate_dataset(parameters, "t", (profile.nt,), "float64")[...]
        ye = _validate_dataset(parameters, "yq", (profile.nye,), "float64")[...]
        thermo_index = _validate_dataset(quantities, profile.thermo_index_name, (len(profile.thermo_ids),), "int32")[...]
        pair_index = _validate_dataset(pairs, "index_yi", (len(profile.pair_ids),), "int32")[...]
        quadruple_index = _validate_dataset(quadruples, "index_av", (len(profile.quadruple_ids),), "int32")[...]
        if len(np.unique(thermo_index)) != len(thermo_index) or tuple(thermo_index) != profile.thermo_ids:
            raise ConversionError("thermodynamic identifier array is duplicate, missing, extra, or reordered")
        if len(np.unique(pair_index)) != len(pair_index) or tuple(pair_index) != profile.pair_ids:
            raise ConversionError("pair identifier array is duplicate, missing, extra, or reordered")
        if len(np.unique(quadruple_index)) != len(quadruple_index) or tuple(quadruple_index) != profile.quadruple_ids:
            raise ConversionError("quadruple identifier array is duplicate, missing, extra, or reordered")

        full_shape = (profile.nye, profile.nt, profile.nrho)
        _validate_dataset(quantities, "thermo", (len(profile.thermo_ids),) + full_shape, "float64")
        _validate_dataset(pairs, "yi", (len(profile.pair_ids),) + full_shape, "float64")
        for name in ("aav", "nav", "yav", "zav"):
            _validate_dataset(quadruples, name, (len(profile.quadruple_ids),) + full_shape, "float64")

        _validate_axis(nb, profile.rho_axis, profile.nb_min, profile.nb_max, "n_b")
        _validate_axis(temperature, profile.temperature_axis, profile.temperature_min, profile.temperature_max, "temperature")
        _validate_axis(ye, profile.ye_axis, profile.ye_min, profile.ye_max, "Y_q")
        return h5, nb, temperature, ye
    except Exception as error:
        _cleanup((("input HDF5 close", h5.close),), error)
        raise


def _scan_energy_extrema(thermo, profile):
    energy_index = profile.thermo_ids.index(21)
    minimum = math.inf
    maximum_absolute = 0.0
    for ye_index in range(profile.nye):
        values = np.asarray(thermo[energy_index, ye_index, :, :], dtype=np.float64) - profile.m_ref_rho_mev
        if not np.all(np.isfinite(values)):
            raise ConversionError("thermodynamic energy contains nonfinite values")
        minimum = min(minimum, float(np.min(values)))
        maximum_absolute = max(maximum_absolute, float(np.max(np.abs(values))))
    return minimum, maximum_absolute


def _regularize_plane(raw, rho, temperature, ye_value, energy_shift, profile):
    pressure_raw = raw["pressure"]
    energy_raw = raw["energy"]
    entropy_raw = raw["entropy"]
    shape = pressure_raw.shape
    if shape != (profile.nt, profile.nrho):
        raise ConversionError("plane has unexpected shape")
    for name, values in raw.items():
        if np.asarray(values).shape != shape or not np.all(np.isfinite(values)):
            raise ConversionError(f"plane field {name} has wrong shape or nonfinite values")
    if not np.all(pressure_raw > 0.0):
        raise ConversionError("raw pressure must be positive")
    if not np.all(energy_raw + energy_shift > 0.0):
        raise ConversionError("shifted raw energy must be positive")

    raw_logenergy = np.log10(energy_raw + energy_shift)
    raw_logpress = np.log10(pressure_raw)
    logenergy = np.empty_like(raw_logenergy)
    logpress = np.empty_like(raw_logpress)
    for rho_index in range(profile.nrho):
        energy_ray = raw_logenergy[:, rho_index]
        pressure_ray = raw_logpress[:, rho_index]
        logenergy[:, rho_index] = _strict_isotonic(energy_ray, _regularization_log_step(energy_ray))
        logpress[:, rho_index] = _strict_isotonic(pressure_ray, _regularization_log_step(pressure_ray))
    energy = np.power(10.0, logenergy) - energy_shift
    pressure = np.power(10.0, logpress)

    energy_per_baryon = energy / (MEV_TO_ERG / (profile.m_ref_rho_mev * MEV_TO_GRAM))
    entropy = np.empty_like(entropy_raw)
    entropy[0, :] = np.maximum(entropy_raw[0, :], 0.0)
    for temperature_index in range(1, profile.nt):
        left_temperature = temperature[temperature_index - 1]
        right_temperature = temperature[temperature_index]
        logarithmic_mean = (right_temperature - left_temperature) / math.log(right_temperature / left_temperature)
        entropy[temperature_index, :] = (
            entropy[temperature_index - 1, :]
            + (energy_per_baryon[temperature_index, :] - energy_per_baryon[temperature_index - 1, :]) / logarithmic_mean
        )

    energy_temperature_derivative = np.empty_like(energy)
    pressure_temperature_derivative = np.empty_like(pressure)
    for rho_index in range(profile.nrho):
        energy_temperature_derivative[:, rho_index] = _monotone_derivative(temperature, energy[:, rho_index])
        pressure_temperature_derivative[:, rho_index] = _monotone_derivative(temperature, pressure[:, rho_index])
    energy_density_derivative = np.gradient(energy, rho, axis=1, edge_order=2)
    pressure_density_derivative = np.gradient(pressure, rho, axis=1, edge_order=2)
    dpderho_raw = pressure_temperature_derivative / energy_temperature_derivative
    dpdrhoe_raw = pressure_density_derivative - dpderho_raw * energy_density_derivative

    density = rho[None, :]
    enthalpy_density = density * (C_SQUARED + energy) + pressure
    if not np.all(np.isfinite(enthalpy_density)) or not np.all(enthalpy_density > 0.0):
        raise ConversionError("regularized enthalpy must be positive and finite")
    first_raw = density * dpdrhoe_raw
    second_raw = pressure / density * dpderho_raw
    first, second, low_count, high_count = _project_bulk(first_raw, second_raw, enthalpy_density)
    dpdrhoe = first / density
    dpderho = second * density / pressure
    bulk = first + second
    cs2 = np.clip(
        C_SQUARED * bulk / enthalpy_density,
        LOW_BULK_FRACTION * C_SQUARED,
        HIGH_BULK_FRACTION * C_SQUARED,
    )
    gamma = bulk / pressure

    abar_raw = raw["abar"]
    zbar_raw = raw["zbar"]
    neutron_number = raw["neutron_number"]
    heavy_number_fraction = raw["heavy_number_fraction"]
    absent = (
        (abar_raw == 0.0) & (zbar_raw == 0.0)
        & (neutron_number == 0.0) & (heavy_number_fraction == 0.0)
    )
    neutron_residual = neutron_number - (abar_raw - zbar_raw)
    neutron_tolerance = 1.0e-8 * np.maximum(1.0, np.abs(abar_raw))
    if np.any(np.abs(neutron_residual) > neutron_tolerance):
        raise ConversionError("N999 is inconsistent with A999-Z999")
    if np.any((abar_raw <= 0.0) & ~absent):
        raise ConversionError("present heavy nuclei require positive A999")
    if np.any((zbar_raw < -neutron_tolerance) | (zbar_raw > abar_raw + neutron_tolerance)):
        raise ConversionError("A999/Z999 values are outside physical bounds")
    zbar_clipped = np.clip(zbar_raw, 0.0, abar_raw)
    charge_clamp = np.abs(zbar_clipped - zbar_raw)
    clamp_mask = (~absent) & (charge_clamp > 0.0)
    abar = np.where(absent, 1.0, abar_raw)
    zbar = np.where(absent, 0.0, zbar_clipped)
    xh = abar_raw * heavy_number_fraction
    composition_raw = np.stack((raw["xn"], raw["xp"], 4.0 * raw["alpha_number_fraction"], xh), axis=0)
    projected = _project_composition(
        composition_raw,
        zbar / abar,
        np.full(shape, ye_value, dtype=np.float64),
    )
    if np.any(absent):
        alpha_absent = np.clip(
            (
                1.0 - composition_raw[0] - composition_raw[1]
                + 2.0 * composition_raw[2]
            ) / 3.0,
            0.0,
            2.0 * np.minimum(ye_value, 1.0 - ye_value),
        )
        projected[2, absent] = alpha_absent[absent]
        projected[1, absent] = ye_value - 0.5 * alpha_absent[absent]
        projected[0, absent] = 1.0 - ye_value - 0.5 * alpha_absent[absent]
        projected[3, absent] = 0.0

    mu_n = raw["mu_b"] + profile.m_n_mev - profile.m_ref_mu_mev
    muhat = -raw["mu_q"]
    mu_p = mu_n - muhat
    munu = raw["mu_l"]
    mu_e = munu + muhat

    fields = {
        "Abar": abar,
        "Xa": projected[2],
        "Xh": projected[3],
        "Xn": projected[0],
        "Xp": projected[1],
        "Zbar": zbar,
        "cs2": cs2,
        "dedt": energy_temperature_derivative,
        "dpderho": dpderho,
        "dpdrhoe": dpdrhoe,
        "entropy": entropy,
        "gamma": gamma,
        "logenergy": logenergy,
        "logpress": logpress,
        "mu_e": mu_e,
        "mu_n": mu_n,
        "mu_p": mu_p,
        "muhat": muhat,
        "munu": munu,
    }
    for name, values in fields.items():
        if not np.all(np.isfinite(values)):
            raise ConversionError(f"regularized field {name} contains nonfinite values")
    distortion = {
        "max_logenergy_change": float(np.max(np.abs(logenergy - raw_logenergy))),
        "max_logpress_change": float(np.max(np.abs(logpress - raw_logpress))),
        "max_entropy_change": float(np.max(np.abs(entropy - entropy_raw))),
        "composition_max_change": float(np.max(np.abs(projected - composition_raw))),
        "bulk_clamped_low": low_count,
        "bulk_clamped_high": high_count,
        "heavy_nucleus_absent_sentinels": int(np.count_nonzero(absent)),
        "heavy_charge_roundoff_clamps": int(np.count_nonzero(clamp_mask)),
        "heavy_charge_roundoff_max": float(np.max(charge_clamp)),
    }
    return fields, distortion


def _read_raw_plane(h5, ye_index, profile):
    thermo = h5["Thermo_qty/thermo"]
    pairs = h5["Composition_pairs/yi"]
    quadruples = h5["Composition_quadrupels"]
    thermo_position = {identifier: position for position, identifier in enumerate(profile.thermo_ids)}
    pair_position = {identifier: position for position, identifier in enumerate(profile.pair_ids)}
    density_mass = profile.m_ref_rho_mev * MEV_TO_GRAM
    energy_unit = MEV_TO_ERG / density_mass
    return {
        "pressure": np.asarray(thermo[thermo_position[1], ye_index, :, :], dtype=np.float64) * MEV_TO_ERG * 1.0e39,
        "entropy": np.asarray(thermo[thermo_position[2], ye_index, :, :], dtype=np.float64),
        "mu_b": np.asarray(thermo[thermo_position[3], ye_index, :, :], dtype=np.float64),
        "mu_q": np.asarray(thermo[thermo_position[4], ye_index, :, :], dtype=np.float64),
        "mu_l": np.asarray(thermo[thermo_position[5], ye_index, :, :], dtype=np.float64),
        "energy": (np.asarray(thermo[thermo_position[21], ye_index, :, :], dtype=np.float64) - profile.m_ref_rho_mev) * energy_unit,
        "xn": np.asarray(pairs[pair_position[10], ye_index, :, :], dtype=np.float64),
        "xp": np.asarray(pairs[pair_position[11], ye_index, :, :], dtype=np.float64),
        "alpha_number_fraction": np.asarray(pairs[pair_position[4002], ye_index, :, :], dtype=np.float64),
        "abar": np.asarray(quadruples["aav"][0, ye_index, :, :], dtype=np.float64),
        "neutron_number": np.asarray(quadruples["nav"][0, ye_index, :, :], dtype=np.float64),
        "heavy_number_fraction": np.asarray(quadruples["yav"][0, ye_index, :, :], dtype=np.float64),
        "zbar": np.asarray(quadruples["zav"][0, ye_index, :, :], dtype=np.float64),
    }


def _create_output(path, profile, rho, temperature, ye, energy_shift):
    h5 = None
    owned = False
    try:
        h5 = h5py.File(path, "x")
        owned = True
        h5.create_dataset("pointsrho", data=np.array([profile.nrho], dtype=np.int32))
        h5.create_dataset("pointstemp", data=np.array([profile.nt], dtype=np.int32))
        h5.create_dataset("pointsye", data=np.array([profile.nye], dtype=np.int32))
        h5.create_dataset("have_rel_cs2", data=np.array([1], dtype=np.int32))
        h5.create_dataset("energy_shift", data=np.array([energy_shift], dtype=np.float64))
        h5.create_dataset("logrho", data=np.log10(rho).astype(np.float64))
        h5.create_dataset("logtemp", data=np.log10(temperature).astype(np.float64))
        h5.create_dataset("ye", data=np.asarray(ye, dtype=np.float64))
        shape = (profile.nye, profile.nt, profile.nrho)
        for name in OUTPUT_FIELDS:
            h5.create_dataset(name, shape=shape, dtype=np.float64)
        return h5
    except Exception as error:
        actions = []
        if h5 is not None:
            actions.append(("partial output close", h5.close))
        if owned:
            actions.append(("partial output removal", lambda: _remove_owned_file(path)))
        _cleanup(actions, error)
        raise


def _write_plane(output, ye_index, fields):
    for name in OUTPUT_FIELDS:
        output[name][ye_index, :, :] = fields[name]


def _cell_midpoints(values):
    return (
        values[:-1, :-1] + values[1:, :-1]
        + values[:-1, 1:] + values[1:, 1:]
    ) / 4.0


def _offgrid_charge_midpoint_max(output, profile):
    maximum = 0.0
    previous = {
        name: np.asarray(output[name][0, :, :], dtype=np.float64)
        for name in ("Abar", "Zbar", "Xp", "Xa", "Xh")
    }
    for ye_index in range(1, profile.nye):
        current = {
            name: np.asarray(output[name][ye_index, :, :], dtype=np.float64)
            for name in ("Abar", "Zbar", "Xp", "Xa", "Xh")
        }
        averaged = {}
        for name in ("Abar", "Zbar", "Xp", "Xa", "Xh"):
            ye_average = 0.5 * (previous[name] + current[name])
            averaged[name] = _cell_midpoints(ye_average)
        ye_midpoint = 0.5 * (output["ye"][ye_index - 1] + output["ye"][ye_index])
        residual = (
            averaged["Xp"] + 0.5 * averaged["Xa"]
            + averaged["Zbar"] / averaged["Abar"] * averaged["Xh"]
            - ye_midpoint
        )
        maximum = max(maximum, float(np.max(np.abs(residual))))
        previous = current
    return maximum


def _manifest(profile, distortion):
    return {
        "profile": profile.name,
        "program_banner": profile.program_banner,
        "manual_version": profile.manual_version,
        "table_id": profile.table_id,
        "thermo_ids": list(profile.thermo_ids),
        "pair_ids": list(profile.pair_ids),
        "quadruple_ids": list(profile.quadruple_ids),
        "m_ref_rho_mev": profile.m_ref_rho_mev,
        "m_ref_mu_mev": profile.m_ref_mu_mev,
        "ye_mapping": "Y_e=Y_q; electron/positron-only local charge neutrality",
        "regularization": "strict-isotonic-derivative-projection-v1",
        "regularization_policy": REGULARIZATION_POLICY,
        "distortion": distortion,
        "nn": profile.nn,
    }


def _validate_output(path, profile, expected_distortion=None):
    try:
        output = h5py.File(path, "r")
    except OSError as error:
        raise ConversionError(f"published candidate could not be reopened: {error}") from error
    try:
        expected_root = set(OUTPUT_FIELDS) | {
            "energy_shift", "grhayl_compose", "have_rel_cs2", "logrho",
            "logtemp", "pointsrho", "pointstemp", "pointsye", "ye",
        }
        if set(output.keys()) != expected_root:
            raise ConversionError("reopened output root does not match StellarCollapse schema")
        for name in expected_root:
            _require_hard_link(output, name)
        for name, value in (("pointsrho", profile.nrho), ("pointstemp", profile.nt), ("pointsye", profile.nye), ("have_rel_cs2", 1)):
            dataset = _validate_dataset(output, name, (1,), "int32")
            if int(dataset[0]) != value:
                raise ConversionError(f"reopened output scalar {name} has wrong value")
        energy_shift = float(_validate_dataset(output, "energy_shift", (1,), "float64")[0])
        logrho = _validate_dataset(output, "logrho", (profile.nrho,), "float64")[...]
        logtemp = _validate_dataset(output, "logtemp", (profile.nt,), "float64")[...]
        ye = _validate_dataset(output, "ye", (profile.nye,), "float64")[...]
        _validate_axis(np.power(10.0, logrho), profile.rho_axis, profile.nb_min * 1.0e39 * profile.m_ref_rho_mev * MEV_TO_GRAM, profile.nb_max * 1.0e39 * profile.m_ref_rho_mev * MEV_TO_GRAM, "reopened density")
        _validate_axis(np.power(10.0, logtemp), profile.temperature_axis, profile.temperature_min, profile.temperature_max, "reopened temperature")
        _validate_axis(ye, profile.ye_axis, profile.ye_min, profile.ye_max, "reopened Y_e")
        shape = (profile.nye, profile.nt, profile.nrho)
        for name in OUTPUT_FIELDS:
            _validate_dataset(output, name, shape, "float64")

        density = np.power(10.0, logrho)[None, :]
        sentinel_count = 0
        for ye_index in range(profile.nye):
            fields = {name: np.asarray(output[name][ye_index, :, :], dtype=np.float64) for name in OUTPUT_FIELDS}
            if any(not np.all(np.isfinite(values)) for values in fields.values()):
                raise ConversionError("reopened output field contains nonfinite values")
            pressure = np.power(10.0, fields["logpress"])
            shifted_energy = np.power(10.0, fields["logenergy"])
            energy = shifted_energy - energy_shift
            enthalpy = C_SQUARED + energy + pressure / density
            if not np.all(pressure > 0.0) or not np.all(shifted_energy > 0.0) or not np.all(enthalpy > 0.0) or not np.all(fields["dedt"] > 0.0):
                raise ConversionError("reopened output violates positive thermodynamic gates")
            for values in (energy, pressure, fields["entropy"], enthalpy):
                if not np.all(np.diff(values, axis=0) > 0.0):
                    raise ConversionError("reopened output temperature ray is not strictly invertible")

            enthalpy_density = density * (C_SQUARED + energy) + pressure
            bulk = density * fields["dpdrhoe"] + pressure / density * fields["dpderho"]
            expected_cs2 = C_SQUARED * bulk / enthalpy_density
            expected_gamma = bulk / pressure
            if np.any(fields["cs2"] < LOW_BULK_FRACTION * C_SQUARED) or np.any(fields["cs2"] > HIGH_BULK_FRACTION * C_SQUARED):
                raise ConversionError("reopened output sound speed is outside fixed causal bounds")
            if np.any(density * fields["dpdrhoe"] < 0.0) or np.any(pressure / density * fields["dpderho"] < 0.0):
                raise ConversionError("reopened derivative contributions are negative")
            if not np.allclose(fields["cs2"], expected_cs2, rtol=3.0e-13, atol=0.0) or not np.allclose(fields["gamma"], expected_gamma, rtol=3.0e-13, atol=0.0):
                raise ConversionError("reopened derivative, sound-speed, or Gamma identity failed")
            if not np.allclose(fields["muhat"], fields["mu_n"] - fields["mu_p"], rtol=0.0, atol=2.0e-12):
                raise ConversionError("reopened nucleon chemical identity failed")
            if not np.allclose(fields["munu"], fields["mu_e"] - fields["muhat"], rtol=0.0, atol=2.0e-12):
                raise ConversionError("reopened lepton chemical identity failed")

            fractions = np.stack((fields["Xn"], fields["Xp"], fields["Xa"], fields["Xh"]), axis=0)
            if np.any(fractions < 0.0) or np.any(fractions > 1.0) or np.any(fields["Abar"] <= 0.0) or np.any(fields["Zbar"] < 0.0) or np.any(fields["Zbar"] > fields["Abar"]):
                raise ConversionError("reopened composition is outside physical bounds")
            mass = np.sum(fractions, axis=0)
            charge = fields["Xp"] + 0.5 * fields["Xa"] + fields["Zbar"] / fields["Abar"] * fields["Xh"]
            if not np.allclose(mass, 1.0, rtol=0.0, atol=NODE_CLOSURE_TOLERANCE) or not np.allclose(charge, ye[ye_index], rtol=0.0, atol=NODE_CLOSURE_TOLERANCE):
                raise ConversionError("reopened node composition closure failed")
            sentinels = (fields["Abar"] == 1.0) & (fields["Zbar"] == 0.0)
            if np.any(fields["Xh"][sentinels] != 0.0):
                raise ConversionError("reopened absent-heavy sentinel is not inert")
            sentinel_count += int(np.count_nonzero(sentinels))

        manifest_group = output["grhayl_compose"]
        if not isinstance(manifest_group, h5py.Group) or set(manifest_group.keys()) != {"manifest_json"} or set(manifest_group.attrs.keys()):
            raise ConversionError("reopened manifest group has wrong shape")
        manifest_dataset = _validate_dataset(manifest_group, "manifest_json", (), h5py.string_dtype("utf-8"))
        try:
            manifest = json.loads(manifest_dataset.asstr()[()])
        except (ValueError, TypeError) as error:
            raise ConversionError("reopened manifest is malformed") from error
        required_manifest_keys = {
            "distortion", "m_ref_mu_mev", "m_ref_rho_mev", "manual_version",
            "nn", "pair_ids", "profile", "program_banner", "quadruple_ids",
            "regularization", "regularization_policy", "table_id", "thermo_ids",
            "ye_mapping",
        }
        if set(manifest) != required_manifest_keys or manifest["profile"] != profile.name or manifest["program_banner"] != profile.program_banner or manifest["manual_version"] != profile.manual_version or manifest["table_id"] != profile.table_id:
            raise ConversionError("reopened manifest facts do not match fixed profile")
        if manifest["m_ref_rho_mev"] != profile.m_ref_rho_mev or manifest["m_ref_mu_mev"] != profile.m_ref_mu_mev or manifest["regularization_policy"] != REGULARIZATION_POLICY:
            raise ConversionError("reopened manifest reference or regularization facts are wrong")
        if manifest["thermo_ids"] != list(profile.thermo_ids) or manifest["pair_ids"] != list(profile.pair_ids) or manifest["quadruple_ids"] != list(profile.quadruple_ids):
            raise ConversionError("reopened manifest identifier facts are wrong")
        if manifest["ye_mapping"] != "Y_e=Y_q; electron/positron-only local charge neutrality" or manifest["regularization"] != "strict-isotonic-derivative-projection-v1" or manifest["nn"] is not profile.nn:
            raise ConversionError("reopened manifest semantic facts are wrong")
        distortion = manifest["distortion"]
        expected_distortion_keys = {
            "bulk_clamped_high", "bulk_clamped_low", "composition_max_change",
            "energy_shift_mev", "heavy_charge_roundoff_clamps",
            "heavy_charge_roundoff_max", "heavy_nucleus_absent_sentinels",
            "max_entropy_change", "max_logenergy_change",
            "max_logpress_change", "offgrid_charge_midpoint_max",
        }
        if not isinstance(distortion, dict) or set(distortion) != expected_distortion_keys:
            raise ConversionError("reopened manifest distortion facts are malformed")
        for name, value in distortion.items():
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                raise ConversionError(f"reopened manifest distortion value {name} is invalid")
        if not isinstance(distortion["bulk_clamped_low"], int) or not isinstance(distortion["bulk_clamped_high"], int):
            raise ConversionError("reopened manifest bulk-clamp counts are invalid")
        for name in ("heavy_charge_roundoff_clamps", "heavy_nucleus_absent_sentinels"):
            if not isinstance(distortion[name], int):
                raise ConversionError("reopened manifest heavy-nucleus counts are invalid")
        if distortion["heavy_nucleus_absent_sentinels"] != sentinel_count:
            raise ConversionError("reopened absent-heavy sentinel count does not match manifest")
        shift_mev = energy_shift * profile.m_ref_rho_mev * MEV_TO_GRAM / MEV_TO_ERG
        if not math.isclose(distortion["energy_shift_mev"], shift_mev, rel_tol=2.0e-15, abs_tol=0.0):
            raise ConversionError("reopened manifest energy shift does not match scalar")
        if expected_distortion is not None and distortion != expected_distortion:
            raise ConversionError("reopened manifest distortion does not match converter results")
    finally:
        _cleanup((("validated output close", output.close),), sys.exc_info()[1])


def _create_exclusive_output(output, profile, rho, temperature, ye, energy_shift):
    if not output.parent.is_dir():
        raise ConversionError("output parent directory does not exist")
    for _ in range(128):
        candidate = output.parent / f".{output.name}.{secrets.token_hex(12)}.tmp"
        try:
            return candidate, _create_output(candidate, profile, rho, temperature, ye, energy_shift)
        except (FileExistsError, OSError) as error:
            if isinstance(error, FileExistsError) or getattr(error, "errno", None) == 17:
                continue
            raise
    raise ConversionError("could not allocate an exclusive converter temporary file")


def convert(input_dir, profile_name, output):
    input_dir = Path(input_dir)
    output = Path(output)
    if profile_name not in PROFILES:
        raise ConversionError(f"unknown profile {profile_name!r}")
    profile = PROFILES[profile_name]
    _checked_cell_count(profile.nrho, profile.nt, profile.nye)
    _checked_plane_cells(profile.nrho, profile.nt)
    if not input_dir.is_dir():
        raise ConversionError("input directory does not exist or is not a directory")
    if os.path.lexists(output):
        raise ConversionError(f"output {output} already exists")
    _validate_controls(input_dir, profile)

    input_h5 = None
    output_h5 = None
    temporary = None
    try:
        input_h5, nb, temperature, ye = _open_validated_input(input_dir / "eoscompose.h5", profile)
        minimum_energy, maximum_absolute_energy = _scan_energy_extrema(input_h5["Thermo_qty/thermo"], profile)
        shift_mev = _energy_shift_from_extrema(minimum_energy, maximum_absolute_energy)
        density_mass = profile.m_ref_rho_mev * MEV_TO_GRAM
        rho = nb * 1.0e39 * density_mass
        energy_shift = shift_mev * MEV_TO_ERG / density_mass
        temporary, output_h5 = _create_exclusive_output(
            output, profile, rho, temperature, ye, energy_shift
        )

        aggregate = {
            "max_logenergy_change": 0.0,
            "max_logpress_change": 0.0,
            "max_entropy_change": 0.0,
            "composition_max_change": 0.0,
            "bulk_clamped_low": 0,
            "bulk_clamped_high": 0,
            "heavy_nucleus_absent_sentinels": 0,
            "heavy_charge_roundoff_clamps": 0,
            "heavy_charge_roundoff_max": 0.0,
        }
        for ye_index in range(profile.nye):
            raw = _read_raw_plane(input_h5, ye_index, profile)
            fields, plane_distortion = _regularize_plane(
                raw, rho, temperature, float(ye[ye_index]), energy_shift, profile
            )
            _write_plane(output_h5, ye_index, fields)
            for name in (
                "max_logenergy_change", "max_logpress_change",
                "max_entropy_change", "composition_max_change",
            ):
                aggregate[name] = max(aggregate[name], plane_distortion[name])
            aggregate["bulk_clamped_low"] += plane_distortion["bulk_clamped_low"]
            aggregate["bulk_clamped_high"] += plane_distortion["bulk_clamped_high"]
            aggregate["heavy_nucleus_absent_sentinels"] += plane_distortion["heavy_nucleus_absent_sentinels"]
            aggregate["heavy_charge_roundoff_clamps"] += plane_distortion["heavy_charge_roundoff_clamps"]
            aggregate["heavy_charge_roundoff_max"] = max(
                aggregate["heavy_charge_roundoff_max"],
                plane_distortion["heavy_charge_roundoff_max"],
            )

        aggregate["energy_shift_mev"] = shift_mev
        aggregate["offgrid_charge_midpoint_max"] = _offgrid_charge_midpoint_max(output_h5, profile)
        ordered_distortion = {
            "bulk_clamped_high": int(aggregate["bulk_clamped_high"]),
            "bulk_clamped_low": int(aggregate["bulk_clamped_low"]),
            "composition_max_change": float(aggregate["composition_max_change"]),
            "energy_shift_mev": float(aggregate["energy_shift_mev"]),
            "heavy_charge_roundoff_clamps": int(aggregate["heavy_charge_roundoff_clamps"]),
            "heavy_charge_roundoff_max": float(aggregate["heavy_charge_roundoff_max"]),
            "heavy_nucleus_absent_sentinels": int(aggregate["heavy_nucleus_absent_sentinels"]),
            "max_entropy_change": float(aggregate["max_entropy_change"]),
            "max_logenergy_change": float(aggregate["max_logenergy_change"]),
            "max_logpress_change": float(aggregate["max_logpress_change"]),
            "offgrid_charge_midpoint_max": float(aggregate["offgrid_charge_midpoint_max"]),
        }
        manifest_group = output_h5.create_group("grhayl_compose")
        manifest_group.create_dataset(
            "manifest_json",
            data=json.dumps(_manifest(profile, ordered_distortion), sort_keys=True, separators=(",", ":")),
            dtype=h5py.string_dtype("utf-8"),
        )
        output_h5.flush()
        output_h5.close()
        output_h5 = None
        input_h5.close()
        input_h5 = None

        _validate_output(temporary, profile, ordered_distortion)
        try:
            os.link(temporary, output)
        except FileExistsError as error:
            raise ConversionError(f"output {output} already exists") from error
        except OSError as error:
            raise ConversionError(f"atomic publication failed: {error}") from error
        os.unlink(temporary)
        temporary = None
        return {"profile": profile_name, **ordered_distortion}
    except ConversionError:
        raise
    except Exception as error:
        converted = ConversionError(str(error))
        converted.cleanup_failures = getattr(error, "cleanup_failures", ())
        raise converted from error
    finally:
        actions = []
        if output_h5 is not None:
            actions.append(("temporary output close", output_h5.close))
        if input_h5 is not None:
            actions.append(("input HDF5 close", input_h5.close))
        if temporary is not None:
            actions.append(
                ("temporary output removal", lambda: _remove_owned_file(temporary))
            )
        _cleanup(actions, sys.exc_info()[1])


class _ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        raise ValueError(message)


def _parse_arguments(argv):
    parser = _ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--output", required=True)
    try:
        arguments = parser.parse_args(argv)
    except ValueError as error:
        print(f"compose_to_grhayl: CLI error: {error}", file=sys.stderr)
        return None
    if arguments.profile not in PROFILES:
        print(f"compose_to_grhayl: CLI error: unknown profile {arguments.profile!r}", file=sys.stderr)
        return None
    return arguments


def main(argv=None):
    arguments = _parse_arguments(argv)
    if arguments is None:
        return 2
    try:
        report = convert(arguments.input_dir, arguments.profile, arguments.output)
    except ConversionError as error:
        print(f"compose_to_grhayl: {error}", file=sys.stderr)
        for failure in getattr(error, "cleanup_failures", ()):
            print(f"compose_to_grhayl: cleanup: {failure}", file=sys.stderr)
        return 1
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
