import contextlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest import mock

import h5py
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = ROOT / "tools" / "compose" / "compose_to_grhayl.py"
SPEC = importlib.util.spec_from_file_location("compose_to_grhayl", MODULE_PATH)
compose = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(compose)

MEV_TO_ERG = 1.602176634e-6
C_SQUARED = 8.987551787368176e20
MEV_TO_GRAM = MEV_TO_ERG / C_SQUARED
M_REF_MEV = 939.5654
M_REF_GRAM = M_REF_MEV * MEV_TO_GRAM
CGS_TO_CODE_PRESSURE = 1.80123683248503e-39

THERMO_IDS = np.array([1, 2, 3, 4, 5, 10, 11, 12, 13, 21], dtype=np.int32)
PAIR_IDS = np.array([10, 11, 4002], dtype=np.int32)
OUTPUT_FIELDS = (
    "Abar", "Xa", "Xh", "Xn", "Xp", "Zbar", "cs2", "dedt",
    "dpderho", "dpdrhoe", "entropy", "gamma", "logenergy",
    "logpress", "mu_e", "mu_n", "mu_p", "muhat", "munu",
)


def _write_controls(directory):
    (directory / "eos.parameters").write_text(
        """1 1 1
0 0
1
1.0e-2 1.0e-5 2.0e-1
1.0e1 1.0e-1 4.0e-1
4 5 3
1 1 0
""",
        encoding="ascii",
    )
    (directory / "eos.quantities").write_text(
        """10 0 0
1 2 3 4 5 10 11 12 13 21
3 1
10 11 4002 999
0

0

2
""",
        encoding="ascii",
    )
    (directory / "eos.thermo").write_text(
        "939.5654 938.2754 1\n", encoding="ascii"
    )


def _fixture_arrays():
    nb = np.logspace(-5.0, -1.0, 5, dtype=np.float64)
    temperature = np.logspace(-2.0, 1.0, 4, dtype=np.float64)
    ye = np.array([0.2, 0.3, 0.4], dtype=np.float64)

    iy, it, ir = np.indices((3, 4, 5), dtype=np.float64)
    density_coordinate = nb[None, None, :] / nb[0]
    pressure = 1.0e-5 * (1.0 + density_coordinate) * (1.0 + temperature[None, :, None]) + 0.0 * iy
    energy = 2.0 + 0.5 * temperature[None, :, None] + 0.1 * iy + 0.0 * ir

    # One asymmetric ray exercises strict regression in full conversion.
    pressure[1, :, 2] = 1.0e-2 * np.array([1.0, 4.0, 2.0, 3.0])
    energy[1, :, 2] = np.array([2.0, 5.0, 3.0, 4.0])

    thermo = np.empty((10, 3, 4, 5), dtype=np.float64)
    thermo[0] = pressure
    thermo[1] = 0.1 + 0.01 * it + 0.001 * ir
    thermo[2] = 4.0 + iy + 0.1 * it + 0.01 * ir
    thermo[3] = -2.0 + 0.1 * iy - 0.01 * ir
    thermo[4] = 6.0 + 0.2 * iy + 0.01 * it
    thermo[5] = -1000.0 - ir
    thermo[6] = -2000.0 - it
    thermo[7] = 9.0
    thermo[8] = -7.0
    thermo[9] = M_REF_MEV + energy

    abar = 80.0 + ir + it
    qheavy = 0.2 + 0.01 * iy
    zbar = qheavy * abar
    xa = np.full_like(abar, 0.04)
    xh = np.full_like(abar, 0.20)
    xp = ye[:, None, None] - 0.5 * xa - qheavy * xh
    xn = 1.0 - xa - xh - xp
    xp = xp + 1.0e-3 * (1.0 + ir)  # Deliberate closure error.

    pairs = np.stack((xn, xp, xa / 4.0), axis=0)

    # Official table convention: all-zero quadruple means no heavy nucleus.
    abar[2, 3, 4] = 0.0
    zbar[2, 3, 4] = 0.0
    xh[2, 3, 4] = 0.0
    # Present-node neutron-number roundoff makes Z exceed A slightly.
    zbar[2, 2, 3] = abar[2, 2, 3] + 1.0e-12
    heavy_number_fraction = np.divide(
        xh, abar, out=np.zeros_like(xh), where=abar != 0.0
    )
    quadruples = {
        "aav": abar[None, ...],
        "nav": (abar - zbar)[None, ...],
        "yav": heavy_number_fraction[None, ...],
        "zav": zbar[None, ...],
    }
    return nb, temperature, ye, thermo, pairs, quadruples


def _write_fixture(directory):
    _write_controls(directory)
    nb, temperature, ye, thermo, pairs, quadruples = _fixture_arrays()
    with h5py.File(directory / "eoscompose.h5", "w") as h5:
        parameters = h5.create_group("Parameters")
        parameters.attrs.create("pointsnb", np.array([5], dtype=np.int32))
        parameters.attrs.create("pointst", np.array([4], dtype=np.int32))
        parameters.attrs.create("pointsyq", np.array([3], dtype=np.int32))
        parameters.attrs.create("tabulation_scheme", np.array([1], dtype=np.int32))
        parameters.create_dataset("nb", data=nb)
        parameters.create_dataset("t", data=temperature)
        parameters.create_dataset("yq", data=ye)

        quantities = h5.create_group("Thermo_qty")
        quantities.attrs.create("pointsqty", np.array([10], dtype=np.int32))
        quantities.create_dataset("index", data=THERMO_IDS)
        quantities.create_dataset("thermo", data=thermo)

        pairs_group = h5.create_group("Composition_pairs")
        pairs_group.attrs.create("pointspairs", np.array([3], dtype=np.int32))
        pairs_group.create_dataset("index_yi", data=PAIR_IDS)
        pairs_group.create_dataset("yi", data=pairs)

        quadruples_group = h5.create_group("Composition_quadrupels")
        quadruples_group.attrs.create("pointsav", np.array([1], dtype=np.int32))
        quadruples_group.create_dataset("index_av", data=np.array([999], dtype=np.int32))
        for name, values in quadruples.items():
            quadruples_group.create_dataset(name, data=values)

        metadata = h5.create_group("metadata")
        metadata.attrs["date"] = np.bytes_("discarded")
        metadata.attrs["time"] = np.bytes_("discarded")


class ComposeConverterTest(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)
        self.input_dir = self.root / "input"
        self.input_dir.mkdir()
        _write_fixture(self.input_dir)
        self.output = self.root / "table.h5"

    def tearDown(self):
        self.tempdir.cleanup()

    def test_strict_isotonic_regression_and_minimum_slope(self):
        raw = np.array([0.0, 2.0, 1.0, 1.0, 4.0])
        cleaned = compose._strict_isotonic(raw, 0.125)
        self.assertTrue(np.all(np.diff(cleaned) >= 0.125 - 4 * np.finfo(float).eps))
        np.testing.assert_allclose(cleaned, [0.0, 29.0 / 24.0, 4.0 / 3.0, 35.0 / 24.0, 4.0])
        np.testing.assert_allclose(compose._strict_isotonic(np.arange(5.0), 0.125), np.arange(5.0))
        with self.assertRaisesRegex(compose.ConversionError, "finite"):
            compose._strict_isotonic(np.array([0.0, np.nan]), 0.1)

    def test_monotone_derivative(self):
        x = np.array([1.0, 2.0, 4.0, 8.0])
        y = np.array([1.0, 3.0, 7.0, 19.0])
        np.testing.assert_allclose(
            compose._monotone_derivative(x, y),
            [2.0, 2.0, 54.0 / 23.0, 3.0],
            rtol=0.0,
            atol=2.0e-15,
        )
        with self.assertRaisesRegex(compose.ConversionError, "strictly increasing"):
            compose._monotone_derivative(x, np.array([1.0, 2.0, 2.0, 4.0]))

    def test_composition_projection(self):
        raw = np.array(
            [
                [0.70, -0.01, -0.50],
                [0.10, 0.40, 0.80],
                [0.10, 0.40, 0.80],
                [0.15, 0.30, -0.10],
            ],
            dtype=np.float64,
        )
        qheavy = np.array([0.25, 0.35, 0.25])
        ye = np.array([0.20, 0.45, 0.20])
        projected = compose._project_composition(raw, qheavy, ye)
        expected = np.array(
            [
                [0.66, 6551.0 / 41350.0, 31.0 / 60.0],
                [0.1228571428571429, 12869.0 / 82700.0, 0.0],
                [0.0914285714285714, 29939.0 / 82700.0, 19.0 / 60.0],
                [0.1257142857142857, 2679.0 / 8270.0, 1.0 / 6.0],
            ]
        )
        np.testing.assert_allclose(projected, expected, rtol=0.0, atol=2.0e-15)
        self.assertTrue(np.all(projected >= 0.0))
        np.testing.assert_allclose(projected.sum(axis=0), 1.0, rtol=0.0, atol=1.0e-13)
        np.testing.assert_allclose(
            projected[1] + 0.5 * projected[2] + qheavy * projected[3],
            ye,
            rtol=0.0,
            atol=1.0e-13,
        )
        with self.assertRaisesRegex(compose.ConversionError, "infeasible"):
            compose._project_composition(raw[:, :1], np.array([0.2]), np.array([1.1]))

    def test_conversion_writes_regularized_stellarcollapse_schema(self):
        report = compose.convert(self.input_dir, "test-analytic", self.output)
        self.assertEqual(report["profile"], "test-analytic")
        self.assertGreater(report["composition_max_change"], 0.0)
        for name in (
            "energy_shift_mev", "max_logenergy_change", "max_logpress_change",
            "max_entropy_change", "composition_max_change",
            "offgrid_charge_midpoint_max",
        ):
            self.assertTrue(np.isfinite(report[name]))
            self.assertGreaterEqual(report[name], 0.0)
        for name in ("bulk_clamped_low", "bulk_clamped_high"):
            self.assertIsInstance(report[name], int)
            self.assertGreaterEqual(report[name], 0)
        self.assertEqual(report["heavy_nucleus_absent_sentinels"], 1)
        self.assertEqual(report["heavy_charge_roundoff_clamps"], 1)
        self.assertGreater(report["heavy_charge_roundoff_max"], 0.0)

        nb, temperature, ye, thermo, _, _ = _fixture_arrays()
        rho = nb * 1.0e39 * M_REF_GRAM
        shape = (3, 4, 5)
        with h5py.File(self.output, "r") as h5:
            self.assertEqual(set(h5), set(OUTPUT_FIELDS) | {
                "energy_shift", "grhayl_compose", "have_rel_cs2", "logrho",
                "logtemp", "pointsrho", "pointstemp", "pointsye", "ye",
            })
            np.testing.assert_allclose(h5["logrho"], np.log10(rho), rtol=0.0, atol=2.0e-15)
            np.testing.assert_allclose(h5["logtemp"], np.log10(temperature), rtol=0.0, atol=2.0e-15)
            np.testing.assert_array_equal(h5["ye"], ye)
            self.assertEqual(h5["pointsrho"].dtype, np.dtype("int32"))
            self.assertEqual(h5["have_rel_cs2"].dtype, np.dtype("int32"))
            self.assertEqual(int(h5["have_rel_cs2"][0]), 1)
            for name, value in (("pointsrho", 5), ("pointstemp", 4), ("pointsye", 3)):
                self.assertEqual(h5[name].shape, (1,))
                self.assertEqual(h5[name].dtype, np.dtype("int32"))
                self.assertEqual(int(h5[name][0]), value)
            for name in ("logrho", "logtemp", "ye", "energy_shift"):
                self.assertEqual(h5[name].dtype, np.dtype("float64"))
                self.assertIsNone(h5[name].chunks)
            for name, expected_shape in (
                ("logrho", (5,)), ("logtemp", (4,)), ("ye", (3,)),
                ("energy_shift", (1,)), ("have_rel_cs2", (1,)),
            ):
                self.assertEqual(h5[name].shape, expected_shape)
                self.assertIsNone(h5[name].chunks)

            for name in OUTPUT_FIELDS:
                self.assertEqual(h5[name].shape, shape)
                self.assertEqual(h5[name].dtype, np.dtype("float64"))
                self.assertIsNone(h5[name].chunks)
                self.assertIsNone(h5[name].compression)
                self.assertTrue(np.all(np.isfinite(h5[name][...])))

            pressure = np.power(10.0, h5["logpress"][...])
            shifted_energy = np.power(10.0, h5["logenergy"][...])
            energy = shifted_energy - h5["energy_shift"][0]
            density = np.power(10.0, h5["logrho"][...])[None, None, :]
            enthalpy = C_SQUARED + energy + pressure / density
            self.assertTrue(np.all(pressure > 0.0))
            self.assertTrue(np.all(shifted_energy > 0.0))
            self.assertTrue(np.all(enthalpy > 0.0))
            self.assertTrue(np.all(h5["dedt"][...] > 0.0))
            self.assertTrue(np.all(density * h5["dpdrhoe"][...] >= 0.0))
            self.assertTrue(np.all(pressure / density * h5["dpderho"][...] >= 0.0))
            for values in (energy, pressure, h5["entropy"][...], enthalpy):
                self.assertTrue(np.all(np.diff(values, axis=1) > 0.0))

            eps_expected = (thermo[9] - M_REF_MEV) * MEV_TO_ERG / M_REF_GRAM
            mask = np.ones(shape, dtype=bool)
            mask[1, :, 2] = False
            np.testing.assert_allclose(energy[mask], eps_expected[mask], rtol=2.0e-9)
            np.testing.assert_allclose(
                pressure[mask],
                (thermo[0] * MEV_TO_ERG * 1.0e39)[mask],
                rtol=2.0e-9,
            )
            np.testing.assert_allclose(h5["mu_n"], thermo[2], rtol=0.0, atol=2.0e-13)
            np.testing.assert_allclose(h5["muhat"], -thermo[3], rtol=0.0, atol=2.0e-13)
            np.testing.assert_allclose(h5["mu_p"], thermo[2] + thermo[3], rtol=0.0, atol=2.0e-13)
            np.testing.assert_allclose(h5["munu"], thermo[4], rtol=0.0, atol=2.0e-13)
            np.testing.assert_allclose(h5["mu_e"], thermo[4] - thermo[3], rtol=0.0, atol=2.0e-13)

            mass = h5["Xn"][...] + h5["Xp"][...] + h5["Xa"][...] + h5["Xh"][...]
            charge = (
                h5["Xp"][...] + 0.5 * h5["Xa"][...]
                + (h5["Zbar"][...] / h5["Abar"][...]) * h5["Xh"][...]
            )
            np.testing.assert_allclose(mass, 1.0, rtol=0.0, atol=1.0e-12)
            np.testing.assert_allclose(
                charge,
                np.broadcast_to(ye[:, None, None], shape),
                rtol=0.0,
                atol=1.0e-12,
            )
            for name in ("Xn", "Xp", "Xa", "Xh"):
                self.assertTrue(np.all(h5[name][...] >= 0.0))
                self.assertTrue(np.all(h5[name][...] <= 1.0))
            self.assertEqual((h5["Abar"][2, 3, 4], h5["Zbar"][2, 3, 4], h5["Xh"][2, 3, 4]), (1.0, 0.0, 0.0))
            self.assertEqual(h5["Zbar"][2, 2, 3], h5["Abar"][2, 2, 3])
            np.testing.assert_allclose(
                [h5[name][2, 3, 4] for name in ("Xn", "Xp", "Xa", "Xh")],
                [0.5475, 0.3475, 0.105, 0.0],
                rtol=0.0,
                atol=2.0e-15,
            )

            H = density * (C_SQUARED + energy) + pressure
            bulk = density * h5["dpdrhoe"][...] + pressure / density * h5["dpderho"][...]
            np.testing.assert_allclose(h5["cs2"], C_SQUARED * bulk / H, rtol=2.0e-13)
            np.testing.assert_allclose(h5["gamma"], bulk / pressure, rtol=2.0e-13)
            self.assertTrue(np.all(h5["cs2"][...] >= 1.0e-14 * C_SQUARED))
            self.assertTrue(np.all(h5["cs2"][...] <= (1.0 - 1.0e-12) * C_SQUARED))

            # Independent field goldens on an untouched analytic ray.
            entropy_golden = 0.1 + 0.5 * np.log(temperature / temperature[0])
            dedt_golden = 0.5 * MEV_TO_ERG / M_REF_GRAM
            np.testing.assert_allclose(h5["entropy"][0, :, 0], entropy_golden, rtol=0.0, atol=5.0e-13)
            np.testing.assert_allclose(h5["dedt"][0, :, 0], dedt_golden, rtol=5.0e-13)
            np.testing.assert_allclose(h5["Abar"][0, :, 0], 80.0 + np.arange(4.0), rtol=0.0, atol=2.0e-13)
            np.testing.assert_allclose(h5["Zbar"][0, :, 0], 0.2 * (80.0 + np.arange(4.0)), rtol=0.0, atol=2.0e-13)
            np.testing.assert_allclose(h5["Xa"][0, :, 0], 4503.0 / 113500.0, rtol=0.0, atol=2.0e-13)

            composition_golden = np.array([140781.0 / 227000.0, 15909.0 / 113500.0, 4503.0 / 113500.0, 9079.0 / 45400.0])
            np.testing.assert_allclose(
                [h5[name][0, 0, 0] for name in ("Xn", "Xp", "Xa", "Xh")],
                composition_golden,
                rtol=0.0,
                atol=2.0e-15,
            )

            rho0 = 1.0e-5 * 1.0e39 * M_REF_GRAM
            pressure0 = 2.0e-5 * 1.01 * MEV_TO_ERG * 1.0e39
            energy0 = 2.005 * MEV_TO_ERG / M_REF_GRAM
            dpdrhoe0 = 1.01 * MEV_TO_ERG / M_REF_GRAM
            dpderho0 = 4.0 * rho0
            bulk0 = rho0 * dpdrhoe0 + pressure0 / rho0 * dpderho0
            H0 = rho0 * (C_SQUARED + energy0) + pressure0
            np.testing.assert_allclose(h5["dpdrhoe"][0, 0, 0], dpdrhoe0, rtol=3.0e-13)
            np.testing.assert_allclose(h5["dpderho"][0, 0, 0], dpderho0, rtol=3.0e-13)
            np.testing.assert_allclose(h5["gamma"][0, 0, 0], bulk0 / pressure0, rtol=3.0e-13)
            np.testing.assert_allclose(h5["cs2"][0, 0, 0], C_SQUARED * bulk0 / H0, rtol=3.0e-13)

            # Componentwise interpolation preserves mass, bounds, and used Xn/Xp.
            midpoint = {name: np.mean(h5[name][0:2, 0:2, 0:2]) for name in ("Abar", "Zbar", "Xn", "Xp", "Xa", "Xh")}
            self.assertGreater(midpoint["Abar"], 0.0)
            self.assertGreaterEqual(midpoint["Zbar"] / midpoint["Abar"], 0.0)
            self.assertLessEqual(midpoint["Zbar"] / midpoint["Abar"], 1.0)
            np.testing.assert_allclose(sum(midpoint[name] for name in ("Xn", "Xp", "Xa", "Xh")), 1.0, rtol=0.0, atol=1.0e-12)
            for name in ("Xn", "Xp", "Xa", "Xh"):
                self.assertGreaterEqual(midpoint[name], 0.0)
                self.assertLessEqual(midpoint[name], 1.0)
            midpoint_ye = np.mean(ye[0:2])
            charge_residual = midpoint["Xp"] + 0.5 * midpoint["Xa"] + midpoint["Zbar"] / midpoint["Abar"] * midpoint["Xh"] - midpoint_ye
            self.assertTrue(np.isfinite(charge_residual))
            def cell_midpoints(values):
                return (
                    values[:-1, :-1, :-1] + values[1:, :-1, :-1]
                    + values[:-1, 1:, :-1] + values[:-1, :-1, 1:]
                    + values[1:, 1:, :-1] + values[1:, :-1, 1:]
                    + values[:-1, 1:, 1:] + values[1:, 1:, 1:]
                ) / 8.0

            cell_values = {name: cell_midpoints(h5[name][...]) for name in ("Abar", "Zbar", "Xp", "Xa", "Xh")}
            cell_ye = 0.5 * (ye[:-1] + ye[1:])[:, None, None]
            all_charge_residuals = (
                cell_values["Xp"] + 0.5 * cell_values["Xa"]
                + cell_values["Zbar"] / cell_values["Abar"] * cell_values["Xh"]
                - cell_ye
            )
            expected_offgrid_max = float(np.max(np.abs(all_charge_residuals)))
            self.assertTrue(np.isfinite(expected_offgrid_max))
            self.assertEqual(report["offgrid_charge_midpoint_max"], expected_offgrid_max)

            # Exact strict-PAVA sentinel. The minimum step is profile-fixed.
            delta_e = compose._regularization_log_step(np.log10((np.array([2.0, 5.0, 3.0, 4.0]) + 20.0) * MEV_TO_ERG / M_REF_GRAM))
            raw_loge = np.log10((np.array([2.0, 5.0, 3.0, 4.0]) + 20.0) * MEV_TO_ERG / M_REF_GRAM)
            expected_delta_e = max(64.0 * np.finfo(float).eps * max(1.0, np.max(np.abs(raw_loge))), 8.0e-10 / np.log(10.0))
            self.assertEqual(delta_e, expected_delta_e)
            expected_loge = raw_loge.copy()
            expected_loge[1] = 0.5 * (raw_loge[1] + raw_loge[2]) - 0.5 * delta_e
            expected_loge[2] = 0.5 * (raw_loge[1] + raw_loge[2]) + 0.5 * delta_e
            np.testing.assert_allclose(h5["logenergy"][1, :, 2], expected_loge, rtol=0.0, atol=2.0e-14)

            raw_logp = np.log10(1.0e-2 * np.array([1.0, 4.0, 2.0, 3.0]) * MEV_TO_ERG * 1.0e39)
            delta_p = compose._pressure_regularization_log_step(raw_logp)
            loaded_logp = raw_logp * np.log(10.0) + np.log(CGS_TO_CODE_PRESSURE)
            expected_delta_p = max(
                64.0 * np.finfo(float).eps * max(1.0, np.max(np.abs(raw_logp))),
                8.0e-10 * max(1.0, np.max(np.abs(loaded_logp))) / np.log(10.0),
            )
            self.assertEqual(delta_p, expected_delta_p)
            expected_logp = raw_logp.copy()
            expected_logp[1] = 0.5 * (raw_logp[1] + raw_logp[2]) - 0.5 * delta_p
            expected_logp[2] = 0.5 * (raw_logp[1] + raw_logp[2]) + 0.5 * delta_p
            np.testing.assert_allclose(h5["logpress"][1, :, 2], expected_logp, rtol=0.0, atol=2.0e-14)

            self.assertAlmostEqual(h5["energy_shift"][0], 20.0 * MEV_TO_ERG / M_REF_GRAM, delta=2.0e-15 * h5["energy_shift"][0])

            manifest = json.loads(h5["grhayl_compose/manifest_json"].asstr()[()])
            self.assertEqual(h5["grhayl_compose/manifest_json"].shape, ())
            self.assertEqual(set(manifest), {
                "distortion", "m_ref_mu_mev", "m_ref_rho_mev", "manual_version",
                "nn", "pair_ids", "profile", "program_banner", "quadruple_ids",
                "regularization", "regularization_policy", "table_id",
                "thermo_ids", "ye_mapping",
            })
            self.assertEqual(manifest["profile"], "test-analytic")
            self.assertEqual(manifest["program_banner"], "2.17-test")
            self.assertEqual(manifest["manual_version"], "3.01-test")
            self.assertEqual(manifest["table_id"], 0)
            self.assertEqual((manifest["m_ref_rho_mev"], manifest["m_ref_mu_mev"]), (939.5654, 939.5654))
            self.assertEqual(manifest["ye_mapping"], "Y_e=Y_q; electron/positron-only local charge neutrality")
            self.assertEqual(manifest["regularization"], "strict-isotonic-derivative-projection-v1")
            self.assertEqual(manifest["thermo_ids"], list(THERMO_IDS))
            self.assertEqual(manifest["pair_ids"], list(PAIR_IDS))
            self.assertEqual(manifest["quadruple_ids"], [999])
            self.assertEqual(manifest["regularization_policy"], {
                "absent_heavy_nucleus": "A=1,Z=0,Xh=0 with three-species L2 projection",
                "bulk_fraction_bounds": [1.0e-14, 1.0 - 1.0e-12],
                "composition_projection": "minimum-L2 active-set",
                "density_derivative": "second-order nonuniform",
                "energy_shift_floor_mev": 20.0,
                "entropy_anchor": "max(raw_S(Tmin),0)",
                "entropy_integrator": "logarithmic-temperature mean",
                "isotonic_weights": "equal",
                "heavy_charge_roundoff": "clip Z to [0,A] within 1e-8*max(1,A)",
                "log_step": "pressure inverse margin uses 8e-10*max(1,abs(ln(P_code)))/ln(10)",
                "node_closure_tolerance": 1.0e-12,
                "offgrid_charge_closure": "diagnostic only; current consumers use Xn/Xp",
                "temperature_derivative": "PCHIP weighted harmonic",
            })

            raw_shifted_logenergy = np.log10(eps_expected + 20.0 * MEV_TO_ERG / M_REF_GRAM)
            raw_logpressure = np.log10(thermo[0] * MEV_TO_ERG * 1.0e39)
            raw_composition = np.stack(
                (_fixture_arrays()[4][0], _fixture_arrays()[4][1], 4.0 * _fixture_arrays()[4][2], _fixture_arrays()[5]["aav"][0] * _fixture_arrays()[5]["yav"][0]),
                axis=0,
            )
            cleaned_composition = np.stack((h5["Xn"][...], h5["Xp"][...], h5["Xa"][...], h5["Xh"][...]), axis=0)
            expected_distortion = {
                "bulk_clamped_high": 1,
                "bulk_clamped_low": 2,
                "composition_max_change": float(np.max(np.abs(cleaned_composition - raw_composition))),
                "energy_shift_mev": 20.0,
                "heavy_charge_roundoff_clamps": 1,
                "heavy_charge_roundoff_max": float(np.max(np.maximum(_fixture_arrays()[5]["zav"][0] - _fixture_arrays()[5]["aav"][0], 0.0))),
                "heavy_nucleus_absent_sentinels": 1,
                "max_entropy_change": float(np.max(np.abs(h5["entropy"][...] - thermo[1]))),
                "max_logenergy_change": float(np.max(np.abs(h5["logenergy"][...] - raw_shifted_logenergy))),
                "max_logpress_change": float(np.max(np.abs(h5["logpress"][...] - raw_logpressure))),
                "offgrid_charge_midpoint_max": expected_offgrid_max,
            }
            self.assertEqual(report, {"profile": "test-analytic", **expected_distortion})
            self.assertEqual(manifest["distortion"], expected_distortion)
            self.assertEqual(manifest["nn"], None)
            for forbidden in ("hash", "mtime", "date", "path"):
                self.assertNotIn(forbidden, json.dumps(manifest).lower())

    def test_separate_density_and_chemical_reference_masses(self):
        base = compose.PROFILES["test-analytic"]
        density_profile = base._replace(m_ref_rho_mev=base.m_ref_rho_mev + 1.0)
        chemical_profile = base._replace(m_ref_mu_mev=base.m_ref_mu_mev + 1.0)
        with mock.patch.dict(compose.PROFILES, {"rho-mutant": density_profile, "mu-mutant": chemical_profile}):
            base_output = self.root / "base.h5"
            rho_output = self.root / "rho.h5"
            mu_output = self.root / "mu.h5"
            compose.convert(self.input_dir, "test-analytic", base_output)
            compose.convert(self.input_dir, "rho-mutant", rho_output)
            compose.convert(self.input_dir, "mu-mutant", mu_output)
        with h5py.File(base_output, "r") as base_h5, h5py.File(rho_output, "r") as rho_h5, h5py.File(mu_output, "r") as mu_h5:
            self.assertFalse(np.array_equal(base_h5["logrho"], rho_h5["logrho"]))
            np.testing.assert_array_equal(base_h5["logrho"], mu_h5["logrho"])
            for field in ("mu_n", "mu_p", "mu_e", "muhat", "munu"):
                np.testing.assert_array_equal(base_h5[field], rho_h5[field])
            np.testing.assert_allclose(mu_h5["mu_n"], base_h5["mu_n"][...] - 1.0, rtol=0.0, atol=2.0e-13)
            np.testing.assert_allclose(mu_h5["mu_p"], base_h5["mu_p"][...] - 1.0, rtol=0.0, atol=2.0e-13)
            for field in ("mu_e", "muhat", "munu"):
                np.testing.assert_array_equal(base_h5[field], mu_h5[field])
            for field in (
                "logenergy", "energy_shift", "logpress", "entropy", "dedt",
                "dpderho", "dpdrhoe", "cs2", "gamma", "Abar", "Zbar",
                "Xa", "Xh", "Xn", "Xp",
            ):
                np.testing.assert_array_equal(base_h5[field], mu_h5[field])

    def test_energy_shift_and_bulk_projection(self):
        self.assertEqual(compose._energy_shift_from_extrema(1.0, 1.0), 20.0)
        for values in (np.array([-25.0, 1.0]), np.array([-1.0e12, 1.0]), np.array([-25.0, 1.0e12])):
            minimum = float(np.min(values))
            scale = max(1.0, float(np.max(np.abs(values))))
            expected = -minimum + 64.0 * np.finfo(float).eps * scale
            self.assertEqual(
                compose._energy_shift_from_extrema(minimum, scale), expected
            )

        u, v, low_count, high_count = compose._project_bulk(
            np.array([-2.0, 9.0]), np.array([1.0, 9.0]), np.array([10.0, 10.0])
        )
        np.testing.assert_allclose(u, [0.0, 5.0 * (1.0 - 1.0e-12)], rtol=0.0, atol=2.0e-15)
        np.testing.assert_allclose(v, [1.0e-13, 5.0 * (1.0 - 1.0e-12)], rtol=0.0, atol=2.0e-15)
        self.assertEqual((low_count, high_count), (1, 1))

    def test_profile_facts(self):
        profile = compose.PROFILES["sro-sly4-sna-141-regularized-v1"]
        self.assertEqual((profile.nrho, profile.nt, profile.nye), (391, 163, 66))
        self.assertEqual(profile.program_banner, "2.17")
        self.assertEqual(profile.manual_version, "3.01")
        self.assertEqual(profile.table_id, 141)
        self.assertEqual(profile.thermo_ids, tuple(THERMO_IDS))
        self.assertEqual(profile.pair_ids, tuple(PAIR_IDS))
        self.assertEqual(profile.quadruple_ids, (999,))
        self.assertEqual((profile.m_n_mev, profile.m_p_mev, profile.i_l), (939.5654, 938.2754, 1))
        self.assertEqual((profile.m_ref_rho_mev, profile.m_ref_mu_mev), (939.5654, 939.5654))
        self.assertEqual((profile.rho_axis, profile.temperature_axis, profile.ye_axis), ("log", "log", "linear"))
        self.assertFalse(profile.beta_equilibrium)
        self.assertFalse(profile.fixed_entropy)
        self.assertEqual(profile.nn, None)

    def test_rejects_bad_inputs_without_output_or_temp_leak(self):
        cases = []

        def duplicate_ids(h5):
            h5["Thermo_qty/index"][1] = h5["Thermo_qty/index"][0]
        cases.append((duplicate_ids, "identifier"))

        def reordered_ids(h5):
            values = h5["Thermo_qty/index"][...]
            values[[0, 1]] = values[[1, 0]]
            h5["Thermo_qty/index"][...] = values
        cases.append((reordered_ids, "identifier"))

        def wrong_dtype(h5):
            values = h5["Parameters/nb"][...].astype(np.float32)
            del h5["Parameters/nb"]
            h5["Parameters"].create_dataset("nb", data=values)
        cases.append((wrong_dtype, "dtype"))

        def extra_group(h5):
            h5.create_group("surprise")
        cases.append((extra_group, "root"))

        def missing_id(h5):
            values = h5["Thermo_qty/index"][...]
            values[-1] = 22
            h5["Thermo_qty/index"][...] = values
        cases.append((missing_id, "identifier"))

        def shortened_ids(h5):
            values = h5["Thermo_qty/index"][...][:-1]
            del h5["Thermo_qty/index"]
            h5["Thermo_qty"].create_dataset("index", data=values)
            h5["Thermo_qty"].attrs.modify("pointsqty", np.array([9], dtype=np.int32))
        cases.append((shortened_ids, "identifier"))

        def extended_pair_ids(h5):
            values = np.append(h5["Composition_pairs/index_yi"][...], np.int32(12))
            del h5["Composition_pairs/index_yi"]
            h5["Composition_pairs"].create_dataset("index_yi", data=values)
            h5["Composition_pairs"].attrs.modify("pointspairs", np.array([4], dtype=np.int32))
        cases.append((extended_pair_ids, "identifier"))

        def duplicate_pair(h5):
            h5["Composition_pairs/index_yi"][1] = 10
        cases.append((duplicate_pair, "identifier"))

        def wrong_rank(h5):
            values = h5["Composition_pairs/yi"][...]
            del h5["Composition_pairs/yi"]
            h5["Composition_pairs"].create_dataset("yi", data=values.reshape((3, -1)))
        cases.append((wrong_rank, "shape"))

        def chunked(h5):
            values = h5["Thermo_qty/thermo"][...]
            del h5["Thermo_qty/thermo"]
            h5["Thermo_qty"].create_dataset("thermo", data=values, chunks=True)
        cases.append((chunked, "layout"))

        def compressed(h5):
            values = h5["Composition_quadrupels/aav"][...]
            del h5["Composition_quadrupels/aav"]
            h5["Composition_quadrupels"].create_dataset("aav", data=values, compression="gzip")
        cases.append((compressed, "layout"))

        def wrong_order(h5):
            values = h5["Thermo_qty/thermo"][...].transpose(0, 1, 3, 2)
            del h5["Thermo_qty/thermo"]
            h5["Thermo_qty"].create_dataset("thermo", data=values)
        cases.append((wrong_order, "shape"))

        def malformed_attribute(h5):
            del h5["Parameters"].attrs["pointsnb"]
            h5["Parameters"].attrs.create("pointsnb", np.array(5, dtype=np.int32))
        cases.append((malformed_attribute, "attribute"))

        def extra_dataset(h5):
            h5["Parameters"].create_dataset("extra", data=np.array([1.0]))
        cases.append((extra_dataset, "dataset"))

        def inconsistent_neutron_number(h5):
            h5["Composition_quadrupels/nav"][0, 0, 0, 0] += 1.0
        cases.append((inconsistent_neutron_number, "N999"))

        def nonuniform_nb(h5):
            h5["Parameters/nb"][2] *= 1.1
        cases.append((nonuniform_nb, "uniform"))

        def nonfinite_temperature(h5):
            h5["Parameters/t"][1] = np.nan
        cases.append((nonfinite_temperature, "finite"))

        def nonuniform_ye(h5):
            h5["Parameters/yq"][1] += 0.01
        cases.append((nonuniform_ye, "uniform"))

        for index, (mutator, message) in enumerate(cases):
            with self.subTest(index=index):
                case_dir = self.root / f"bad-{index}"
                shutil.copytree(self.input_dir, case_dir)
                with h5py.File(case_dir / "eoscompose.h5", "r+") as h5:
                    mutator(h5)
                output = self.root / f"bad-{index}.h5"
                with self.assertRaisesRegex(compose.ConversionError, message):
                    compose.convert(case_dir, "test-analytic", output)
                self.assertFalse(output.exists())
                self.assertEqual(list(self.root.glob(f".{output.name}.*.tmp")), [])

        link_dir = self.root / "bad-link"
        shutil.copytree(self.input_dir, link_dir)
        with h5py.File(link_dir / "eoscompose.h5", "r+") as h5:
            del h5["Parameters/nb"]
            h5["Parameters/nb"] = h5py.SoftLink("/Parameters/t")
        with self.assertRaisesRegex(compose.ConversionError, "link"):
            compose.convert(link_dir, "test-analytic", self.root / "link.h5")
        self.assertFalse((self.root / "link.h5").exists())
        self.assertEqual(list(self.root.glob(".link.h5.*.tmp")), [])

        external_link_dir = self.root / "bad-external-link"
        shutil.copytree(self.input_dir, external_link_dir)
        with h5py.File(external_link_dir / "eoscompose.h5", "r+") as h5:
            del h5["Parameters/nb"]
            h5["Parameters/nb"] = h5py.ExternalLink("must-not-open.h5", "/nb")
        with self.assertRaisesRegex(compose.ConversionError, "must be a hard link"):
            compose.convert(external_link_dir, "test-analytic", self.root / "external-link.h5")
        self.assertFalse((self.root / "external-link.h5").exists())
        self.assertEqual(list(self.root.glob(".external-link.h5.*.tmp")), [])

        vds_dir = self.root / "bad-vds"
        shutil.copytree(self.input_dir, vds_dir)
        source_path = vds_dir / "vds-source.h5"
        with h5py.File(source_path, "w") as source:
            source.create_dataset("nb", data=np.logspace(-5.0, -1.0, 5))
        with h5py.File(vds_dir / "eoscompose.h5", "r+") as h5:
            del h5["Parameters/nb"]
            layout = h5py.VirtualLayout(shape=(5,), dtype=np.float64)
            layout[:] = h5py.VirtualSource(str(source_path), "nb", shape=(5,))
            h5["Parameters"].create_virtual_dataset("nb", layout)
        with self.assertRaisesRegex(compose.ConversionError, "virtual"):
            compose.convert(vds_dir, "test-analytic", self.root / "vds.h5")
        self.assertFalse((self.root / "vds.h5").exists())
        self.assertEqual(list(self.root.glob(".vds.h5.*.tmp")), [])

        external_dir = self.root / "bad-external-storage"
        shutil.copytree(self.input_dir, external_dir)
        with h5py.File(external_dir / "eoscompose.h5", "r+") as h5:
            values = h5["Parameters/nb"][...]
            del h5["Parameters/nb"]
            dataset = h5["Parameters"].create_dataset(
                "nb", shape=(5,), dtype=np.float64,
                external=[(str(external_dir / "nb-external.bin"), 0, h5py.h5f.UNLIMITED)],
            )
            dataset[...] = values
        with self.assertRaisesRegex(compose.ConversionError, "external storage"):
            compose.convert(external_dir, "test-analytic", self.root / "external-storage.h5")
        self.assertFalse((self.root / "external-storage.h5").exists())
        self.assertEqual(list(self.root.glob(".external-storage.h5.*.tmp")), [])

        control_dir = self.root / "bad-control"
        shutil.copytree(self.input_dir, control_dir)
        (control_dir / "eos.parameters").write_text("1 1 1\n", encoding="ascii")
        with self.assertRaisesRegex(compose.ConversionError, "control"):
            compose.convert(control_dir, "test-analytic", self.root / "control.h5")
        self.assertFalse((self.root / "control.h5").exists())
        self.assertEqual(list(self.root.glob(".control.h5.*.tmp")), [])

        for index, (filename, old, new) in enumerate((
            ("eos.parameters", "1 1 0\n", "1 1 1\n"),
            ("eos.quantities", "\n2\n", "\n1\n"),
            ("eos.thermo", "939.5654", "940.5654"),
        )):
            with self.subTest(control=filename):
                semantic_dir = self.root / f"bad-semantic-control-{index}"
                shutil.copytree(self.input_dir, semantic_dir)
                path = semantic_dir / filename
                path.write_text(path.read_text(encoding="ascii").replace(old, new), encoding="ascii")
                with self.assertRaisesRegex(compose.ConversionError, "control"):
                    compose.convert(semantic_dir, "test-analytic", self.root / f"semantic-{index}.h5")
                self.assertFalse((self.root / f"semantic-{index}.h5").exists())
                self.assertEqual(list(self.root.glob(f".semantic-{index}.h5.*.tmp")), [])

        os.unlink(self.input_dir / "eos.thermo")
        os.symlink("eos.parameters", self.input_dir / "eos.thermo")
        with self.assertRaisesRegex(compose.ConversionError, "regular file"):
            compose.convert(self.input_dir, "test-analytic", self.root / "symlink.h5")
        self.assertFalse((self.root / "symlink.h5").exists())
        self.assertEqual(list(self.root.glob(".symlink.h5.*.tmp")), [])

    def test_resource_limit_and_atomic_failures(self):
        self.assertEqual(compose._checked_cell_count(1, 1, 107_374_182), 107_374_182)
        with self.assertRaisesRegex(compose.ConversionError, "INT_MAX"):
            compose._checked_cell_count(1, 1, 107_374_183)
        with self.assertRaisesRegex(compose.ConversionError, "dimension"):
            compose._checked_cell_count(0, 2, 3)
        self.assertEqual(compose._checked_plane_cells(391, 163), 63_733)
        with self.assertRaisesRegex(compose.ConversionError, "plane"):
            compose._checked_plane_cells(400, 164)

        self.output.write_bytes(b"owner data")
        with self.assertRaisesRegex(compose.ConversionError, "exists"):
            compose.convert(self.input_dir, "test-analytic", self.output)
        self.assertEqual(self.output.read_bytes(), b"owner data")

        self.output.unlink()
        decoy = self.root / f".{self.output.name}.user.tmp"
        decoy.write_bytes(b"keep")
        with mock.patch.object(compose, "_write_plane", side_effect=OSError("injected write")):
            with self.assertRaisesRegex(compose.ConversionError, "injected write"):
                compose.convert(self.input_dir, "test-analytic", self.output)
        self.assertEqual([path for path in self.root.glob(f".{self.output.name}.*.tmp") if path != decoy], [])
        self.assertEqual(decoy.read_bytes(), b"keep")

        with mock.patch.object(compose, "_validate_output", side_effect=compose.ConversionError("reopen corruption")):
            with self.assertRaisesRegex(compose.ConversionError, "reopen corruption"):
                compose.convert(self.input_dir, "test-analytic", self.output)
        self.assertFalse(self.output.exists())
        self.assertEqual([path for path in self.root.glob(f".{self.output.name}.*.tmp") if path != decoy], [])

        def lose_race(source, destination):
            Path(destination).write_bytes(b"winner")
            raise FileExistsError("race")
        with mock.patch.object(compose.os, "link", side_effect=lose_race):
            with self.assertRaisesRegex(compose.ConversionError, "exists"):
                compose.convert(self.input_dir, "test-analytic", self.output)
        self.assertEqual(self.output.read_bytes(), b"winner")
        self.assertEqual([path for path in self.root.glob(f".{self.output.name}.*.tmp") if path != decoy], [])

        self.output.unlink()
        with mock.patch.object(compose.os, "link", side_effect=OSError("link failure")):
            with self.assertRaisesRegex(compose.ConversionError, "link failure"):
                compose.convert(self.input_dir, "test-analytic", self.output)
        self.assertFalse(self.output.exists())
        self.assertEqual([path for path in self.root.glob(f".{self.output.name}.*.tmp") if path != decoy], [])

    def test_reopen_validator_detects_corruption_classes(self):
        compose.convert(self.input_dir, "test-analytic", self.output)
        mutations = [
            ("scalar", lambda h5: h5["pointsrho"].__setitem__(0, 4)),
            ("axis", lambda h5: h5["logtemp"].__setitem__(1, h5["logtemp"][0])),
            ("manifest", lambda h5: h5["grhayl_compose/manifest_json"].__setitem__((), "{}")),
        ]
        mutations.extend(
            (field, lambda h5, field=field: h5[field].__setitem__((0, 0, 0), np.nan))
            for field in OUTPUT_FIELDS
        )
        for name, mutate in mutations:
            with self.subTest(name=name):
                corrupted = self.root / f"corrupt-{name}.h5"
                shutil.copy2(self.output, corrupted)
                with h5py.File(corrupted, "r+") as h5:
                    mutate(h5)
                with self.assertRaises(compose.ConversionError):
                    compose._validate_output(corrupted, compose.PROFILES["test-analytic"])

    def test_cli_exit_codes(self):
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(compose.main([]), 2)
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            self.assertEqual(compose.main(["--profile", "unknown", "--input-dir", str(self.input_dir), "--output", str(self.output)]), 2)
        self.assertIn("profile", stderr.getvalue())

        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            self.assertEqual(compose.main(["--profile", "test-analytic", "--input-dir", str(self.root / "missing"), "--output", str(self.output)]), 1)
        self.assertIn("input", stderr.getvalue())

        self.assertEqual(compose.main(["--profile", "test-analytic", "--input-dir", str(self.input_dir), "--output", str(self.output)]), 0)
        self.assertTrue(self.output.exists())


if __name__ == "__main__":
    unittest.main()
