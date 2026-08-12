import contextlib
import errno
import io
import json
import math
import os
from pathlib import Path
import runpy
import shutil
import tempfile
import unittest
from unittest import mock

import h5py
import numpy as np

from test_compose_to_grhayl import _fixture_arrays, _write_fixture, compose


class FailureCoverageTest(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.input_dir = self.root / "input"
        self.input_dir.mkdir()
        _write_fixture(self.input_dir)
        self.output = self.root / "output.h5"
        self.profile = compose.PROFILES["test-analytic"]

    def tearDown(self):
        self.temporary_directory.cleanup()

    def _mutated_input(self, name, mutation, message=None):
        directory = self.root / name
        shutil.copytree(self.input_dir, directory)
        with h5py.File(directory / "eoscompose.h5", "r+") as h5:
            mutation(h5)
        output = self.root / f"{name}.h5"
        context = self.assertRaisesRegex(compose.ConversionError, message) if message else self.assertRaises(compose.ConversionError)
        with context:
            compose.convert(directory, "test-analytic", output)
        self.assertFalse(output.exists())

    def _valid_output(self):
        compose.convert(self.input_dir, "test-analytic", self.output)

    def _corrupt_output(self, name, mutation, expected_distortion=None):
        corrupted = self.root / f"corrupt-{name}.h5"
        shutil.copy2(self.output, corrupted)
        with h5py.File(corrupted, "r+") as h5:
            mutation(h5)
        with self.assertRaises(compose.ConversionError):
            compose._validate_output(corrupted, self.profile, expected_distortion)

    def test_helper_rejections(self):
        with self.assertRaises(compose.ConversionError):
            compose._regularization_log_step([])
        with self.assertRaises(compose.ConversionError):
            compose._strict_isotonic([1.0, 2.0], 0.0)
        with mock.patch.object(compose.np, "diff", return_value=np.array([0.0])):
            with self.assertRaises(compose.ConversionError):
                compose._strict_isotonic([1.0, 2.0], 0.1)

        with self.assertRaises(compose.ConversionError):
            compose._monotone_derivative([1.0, 2.0], [1.0, 2.0])
        with self.assertRaises(compose.ConversionError):
            compose._monotone_derivative([1.0, np.nan, 3.0], [1.0, 2.0, 3.0])
        with self.assertRaises(compose.ConversionError):
            compose._energy_shift_mev([np.inf])
        with self.assertRaises(compose.ConversionError):
            compose._energy_shift_from_extrema(np.inf, 1.0)
        with np.errstate(over="ignore"):
            with self.assertRaises(compose.ConversionError):
                compose._energy_shift_from_extrema(-np.finfo(float).max, np.finfo(float).max)
        with self.assertRaises(compose.ConversionError):
            compose._inverse_temperature([1.0, 2.0, 3.0], [1.0, 2.0, 3.0], 4.0)
        with self.assertRaises(compose.ConversionError):
            compose._project_bulk([1.0], [1.0, 2.0], [1.0])
        with self.assertRaises(compose.ConversionError):
            compose._project_bulk([1.0], [1.0], [-1.0])
        with self.assertRaises(compose.ConversionError):
            compose._project_composition(np.ones((3, 1)), np.ones(1), np.ones(1))
        with self.assertRaises(compose.ConversionError):
            compose._project_composition(np.full((4, 1), np.nan), np.zeros(1), np.zeros(1))
        with mock.patch.object(compose.itertools, "combinations", return_value=[]):
            with self.assertRaises(compose.ConversionError):
                compose._project_composition(np.full((4, 1), 0.25), np.array([0.25]), np.array([0.25]))
        with mock.patch.object(compose.itertools, "combinations", return_value=[(0, 3)]):
            with self.assertRaises(compose.ConversionError):
                compose._project_composition(np.full((4, 1), 0.25), np.array([0.0]), np.array([0.25]))
        with mock.patch.object(compose.np, "max", return_value=1.0):
            with self.assertRaises(compose.ConversionError):
                compose._project_composition(np.full((4, 1), 0.25), np.array([0.25]), np.array([0.25]))

    def test_control_read_and_parse_rejections(self):
        with self.assertRaises(compose.ConversionError):
            compose._read_control_tokens(self.root / "missing-control")
        with self.assertRaises(compose.ConversionError):
            compose._require_regular_file(self.root / "missing-regular")

        oversized = self.root / "oversized-control"
        oversized.write_bytes(b"1" * (compose.MAX_CONTROL_BYTES + 1))
        with self.assertRaisesRegex(compose.ConversionError, "byte limit"):
            compose._read_control_tokens(oversized)

        too_many_tokens = self.root / "too-many-control-tokens"
        too_many_tokens.write_text(
            "1 " * (compose.MAX_CONTROL_TOKENS + 1), encoding="ascii"
        )
        with self.assertRaisesRegex(compose.ConversionError, "token limit"):
            compose._read_control_tokens(too_many_tokens)

        nonascii = self.root / "nonascii-control"
        nonascii.write_bytes(b"\xff")
        with self.assertRaisesRegex(compose.ConversionError, "could not be read"):
            compose._read_control_tokens(nonascii)

        cases = (
            ("parameters-nonnumeric", "eos.parameters", "bad\n"),
            ("quantities-nonintegral", "eos.quantities", "1.5\n"),
            ("thermo-overlong", "eos.thermo", "1" * 513),
            ("thermo-incomplete", "eos.thermo", "939.5654 938.2754\n"),
            ("thermo-malformed", "eos.thermo", "x 938.2754 1\n"),
        )
        for name, filename, contents in cases:
            with self.subTest(name=name):
                directory = self.root / name
                shutil.copytree(self.input_dir, directory)
                (directory / filename).write_text(contents, encoding="ascii")
                with self.assertRaises(compose.ConversionError):
                    compose.convert(directory, "test-analytic", self.root / f"{name}.h5")

        with mock.patch.object(Path, "open", side_effect=UnicodeError("bad encoding")):
            with self.assertRaises(compose.ConversionError):
                compose._validate_controls(self.input_dir, self.profile)

        parameter_tokens = compose._read_control_tokens(self.input_dir / "eos.parameters")
        quantity_tokens = compose._read_control_tokens(self.input_dir / "eos.quantities")

        def fixed_control_tokens(path):
            return parameter_tokens if path.name == "eos.parameters" else quantity_tokens

        with mock.patch.object(compose, "_read_control_tokens", side_effect=fixed_control_tokens), mock.patch.object(Path, "open", side_effect=UnicodeError("bad thermo encoding")):
            with self.assertRaises(compose.ConversionError):
                compose._validate_controls(self.input_dir, self.profile)

    def test_hdf_schema_rejection_families(self):
        invalid = self.root / "not-hdf.h5"
        invalid.write_text("not hdf", encoding="ascii")
        with self.assertRaises(compose.ConversionError):
            compose._open_validated_input(invalid, self.profile)

        self._mutated_input("root-not-group", lambda h5: (h5.__delitem__("metadata"), h5.create_dataset("metadata", data=np.array([1]))))
        self._mutated_input("thermo-extra", lambda h5: h5["Thermo_qty"].create_dataset("extra", data=[1]))
        self._mutated_input("pairs-extra", lambda h5: h5["Composition_pairs"].create_dataset("extra", data=[1]))
        self._mutated_input("quad-extra", lambda h5: h5["Composition_quadrupels"].create_dataset("extra", data=[1]))
        self._mutated_input("metadata-child", lambda h5: h5["metadata"].create_dataset("extra", data=[1]))
        self._mutated_input("metadata-bad-attr", lambda h5: h5["metadata"].attrs.__setitem__("date", np.int32(1)))
        self._mutated_input("parameter-extra-attr", lambda h5: h5["Parameters"].attrs.__setitem__("extra", np.array([1], dtype=np.int32)))
        self._mutated_input("quantity-extra-attr", lambda h5: h5["Thermo_qty"].attrs.__setitem__("extra", np.array([1], dtype=np.int32)))
        self._mutated_input("quad-id", lambda h5: h5["Composition_quadrupels/index_av"].__setitem__(0, 998), "identifier")
        self._mutated_input("missing-attribute", lambda h5: h5["Parameters"].attrs.__delitem__("pointsnb"), "attribute")
        self._mutated_input("tabulation-value", lambda h5: h5["Parameters"].attrs.__setitem__("tabulation_scheme", np.array([2], dtype=np.int32)), "attribute")

        def hard_link_alias(h5):
            del h5["Composition_quadrupels/nav"]
            h5["Composition_quadrupels/nav"] = h5["Composition_quadrupels/aav"]

        self._mutated_input("hard-link-alias", hard_link_alias, "aliases")
        self._mutated_input(
            "nb-lower",
            lambda h5: h5["Parameters/nb"].__setitem__(slice(None), np.logspace(math.log10(2.0e-5), -1.0, 5)),
            "lower",
        )
        self._mutated_input(
            "nb-upper",
            lambda h5: h5["Parameters/nb"].__setitem__(slice(None), np.logspace(-5.0, math.log10(2.0e-1), 5)),
            "upper",
        )

        def dataset_as_group(h5):
            del h5["Parameters/nb"]
            h5["Parameters"].create_group("nb")
        self._mutated_input("dataset-as-group", dataset_as_group, "dataset")

        def compact_dataset(h5):
            values = h5["Parameters/nb"][...]
            del h5["Parameters/nb"]
            space = h5py.h5s.create_simple(values.shape)
            properties = h5py.h5p.create(h5py.h5p.DATASET_CREATE)
            properties.set_layout(h5py.h5d.COMPACT)
            dataset = h5py.h5d.create(
                h5["Parameters"].id,
                b"nb",
                h5py.h5t.py_create(np.dtype("float64")),
                space,
                dcpl=properties,
            )
            dataset.write(h5py.h5s.ALL, h5py.h5s.ALL, values)
            dataset.close()

        self._mutated_input("compact", compact_dataset, "layout")

        class BrokenGroup:
            name = "/broken"

            def get(self, *args, **kwargs):
                raise RuntimeError("broken")

        with self.assertRaises(compose.ConversionError):
            compose._require_hard_link(BrokenGroup(), "x")
        with h5py.File(self.input_dir / "eoscompose.h5", "r") as h5:
            with self.assertRaises(compose.ConversionError):
                compose._validate_attribute(h5["Parameters"], "missing", 1)

    def test_energy_and_plane_rejections(self):
        with h5py.File(self.input_dir / "eoscompose.h5", "r") as h5:
            raw = compose._read_raw_plane(h5, 0, self.profile)
            rho = h5["Parameters/nb"][...] * 1.0e39 * self.profile.m_ref_rho_mev * compose.MEV_TO_GRAM
            temperature = h5["Parameters/t"][...]
        class NonfiniteThermo:
            def __getitem__(self, key):
                return np.full((4, 5), np.nan)

        with self.assertRaises(compose.ConversionError):
            compose._scan_energy_extrema(NonfiniteThermo(), self.profile)

        shift = 20.0 * compose.MEV_TO_ERG / (self.profile.m_ref_rho_mev * compose.MEV_TO_GRAM)

        def reject(mutator):
            changed = {name: np.array(values, copy=True) for name, values in raw.items()}
            mutator(changed)
            with self.assertRaises(compose.ConversionError):
                compose._regularize_plane(changed, rho, temperature, 0.2, shift, self.profile)

        reject(lambda values: values.__setitem__("pressure", values["pressure"][:-1]))
        reject(lambda values: values["entropy"].__setitem__((0, 0), np.nan))
        reject(lambda values: values["pressure"].__setitem__((0, 0), 0.0))
        reject(lambda values: values["energy"].__setitem__((0, 0), -shift))
        reject(lambda values: values["abar"].__setitem__((0, 0), 0.0))

        def invalid_present_mass(values):
            values["abar"][0, 0] = 0.0
            values["zbar"][0, 0] = 0.0
            values["neutron_number"][0, 0] = 0.0
            values["heavy_number_fraction"][0, 0] = 0.1

        reject(invalid_present_mass)

        def invalid_charge(values):
            values["zbar"][0, 0] = values["abar"][0, 0] + 1.0
            values["neutron_number"][0, 0] = -1.0

        reject(invalid_charge)

        with mock.patch.object(compose, "_project_bulk", return_value=(np.zeros((4, 5)), np.zeros((4, 5)), 0, 0)):
            changed = {name: np.array(values, copy=True) for name, values in raw.items()}
            changed["energy"][:] = -compose.C_SQUARED - changed["pressure"] / rho[None, :] - 1.0
            with self.assertRaises(compose.ConversionError):
                compose._regularize_plane(changed, rho, temperature, 0.2, shift + compose.C_SQUARED * 2.0, self.profile)

        with mock.patch.object(compose, "_project_composition", return_value=np.full((4, 4, 5), np.nan)):
            with self.assertRaises(compose.ConversionError):
                compose._regularize_plane(raw, rho, temperature, 0.2, shift, self.profile)

    def test_create_output_cleanup_and_exclusive_failures(self):
        rho = np.logspace(1.0, 5.0, 5)
        temperature = np.logspace(-2.0, 1.0, 4)
        ye = np.array([0.2, 0.3, 0.4])

        with mock.patch.object(compose.h5py, "File", side_effect=OSError("open failed")):
            with self.assertRaises(OSError):
                compose._create_output(self.root / "never.h5", self.profile, rho, temperature, ye, 1.0)

        class BadFile:
            def create_dataset(self, *args, **kwargs):
                raise RuntimeError("write failed")

            def close(self):
                pass

        owned = self.root / "owned.h5"

        def create_owned_file(*args, **kwargs):
            owned.write_bytes(b"converter partial")
            return BadFile()

        with mock.patch.object(compose.h5py, "File", side_effect=create_owned_file):
            with self.assertRaises(RuntimeError):
                compose._create_output(owned, self.profile, rho, temperature, ye, 1.0)
        self.assertFalse(owned.exists())

        with mock.patch.object(compose.h5py, "File", return_value=BadFile()):
            with self.assertRaises(RuntimeError):
                compose._create_output(self.root / "missing-partial.h5", self.profile, rho, temperature, ye, 1.0)

        with self.assertRaises(compose.ConversionError):
            compose._create_exclusive_output(self.root / "missing" / "out.h5", self.profile, rho, temperature, ye, 1.0)
        with mock.patch.object(compose, "_create_output", side_effect=OSError(errno.EIO, "disk")):
            with self.assertRaises(OSError):
                compose._create_exclusive_output(self.output, self.profile, rho, temperature, ye, 1.0)

        collision = self.root / ".output.h5.same.tmp"
        collision.write_bytes(b"racer")
        with mock.patch.object(compose.secrets, "token_hex", return_value="same"):
            with self.assertRaises(compose.ConversionError):
                compose._create_exclusive_output(self.output, self.profile, rho, temperature, ye, 1.0)
        self.assertEqual(collision.read_bytes(), b"racer")

    def test_reopen_rejection_families(self):
        self._valid_output()
        self._corrupt_output("root-extra", lambda h5: h5.create_dataset("extra", data=[1]))
        self._corrupt_output("scalar", lambda h5: h5["have_rel_cs2"].__setitem__(0, 0))
        self._corrupt_output("positive", lambda h5: h5["logpress"].__setitem__((0, 0, 0), -400.0))
        self._corrupt_output("inverse", lambda h5: h5["entropy"].__setitem__((0, 1, 0), h5["entropy"][0, 0, 0]))
        self._corrupt_output("cs-low", lambda h5: h5["cs2"].__setitem__((0, 0, 0), 0.0))
        self._corrupt_output("derivative-negative", lambda h5: h5["dpdrhoe"].__setitem__((0, 0, 0), -1.0))
        self._corrupt_output("identity", lambda h5: h5["gamma"].__setitem__((0, 0, 0), 0.0))
        self._corrupt_output("chemical-n", lambda h5: h5["muhat"].__setitem__((0, 0, 0), 0.0))
        self._corrupt_output("chemical-l", lambda h5: h5["munu"].__setitem__((0, 0, 0), 0.0))
        self._corrupt_output("fraction", lambda h5: h5["Xn"].__setitem__((0, 0, 0), -1.0))
        self._corrupt_output("closure", lambda h5: h5["Xp"].__setitem__((0, 0, 0), h5["Xp"][0, 0, 0] + 0.01))

        def activate_absent_sentinel(h5):
            h5["Xh"][2, 3, 4] = 0.1
            h5["Xn"][2, 3, 4] -= 0.1

        self._corrupt_output("active-sentinel", activate_absent_sentinel)
        self._corrupt_output("manifest-group", lambda h5: h5["grhayl_compose"].attrs.__setitem__("extra", 1))
        self._corrupt_output("manifest-json", lambda h5: h5["grhayl_compose/manifest_json"].__setitem__((), "not json"))

        def alter_manifest(h5, change):
            dataset = h5["grhayl_compose/manifest_json"]
            manifest = json.loads(dataset.asstr()[()])
            change(manifest)
            dataset[()] = json.dumps(manifest, separators=(",", ":"))

        changes = (
            ("facts", lambda value: value.__setitem__("profile", "wrong")),
            ("references", lambda value: value.__setitem__("m_ref_rho_mev", 1.0)),
            ("ids", lambda value: value.__setitem__("thermo_ids", [])),
            ("semantics", lambda value: value.__setitem__("ye_mapping", "wrong")),
            ("distortion-shape", lambda value: value.__setitem__("distortion", {})),
            ("distortion-value", lambda value: value["distortion"].__setitem__("max_entropy_change", -1.0)),
            ("distortion-count", lambda value: value["distortion"].__setitem__("bulk_clamped_low", 1.5)),
            ("heavy-count-type", lambda value: value["distortion"].__setitem__("heavy_charge_roundoff_clamps", 1.5)),
            ("sentinel-count", lambda value: value["distortion"].__setitem__("heavy_nucleus_absent_sentinels", 2)),
            ("distortion-shift", lambda value: value["distortion"].__setitem__("energy_shift_mev", 21.0)),
        )
        for name, change in changes:
            self._corrupt_output(name, lambda h5, change=change: alter_manifest(h5, change))

        with h5py.File(self.output, "r") as h5:
            distortion = json.loads(h5["grhayl_compose/manifest_json"].asstr()[()])["distortion"]
        wrong = dict(distortion)
        wrong["max_entropy_change"] += 1.0
        with self.assertRaises(compose.ConversionError):
            compose._validate_output(self.output, self.profile, wrong)

        with self.assertRaises(compose.ConversionError):
            compose._validate_output(self.root / "missing-output.h5", self.profile)

    def test_conversion_exception_and_entrypoint_paths(self):
        with self.assertRaises(compose.ConversionError):
            compose.convert(self.input_dir, "does-not-exist", self.output)
        with mock.patch.object(compose, "_open_validated_input", side_effect=RuntimeError("unexpected")):
            with self.assertRaisesRegex(compose.ConversionError, "unexpected"):
                compose.convert(self.input_dir, "test-analytic", self.output)

        with mock.patch.object(compose, "_write_plane", side_effect=compose.ConversionError("stop")), mock.patch.object(compose.os, "unlink", side_effect=FileNotFoundError):
            with self.assertRaisesRegex(compose.ConversionError, "stop"):
                compose.convert(self.input_dir, "test-analytic", self.output)

        stderr = io.StringIO()
        old_argv = list(os.sys.argv)
        try:
            os.sys.argv = [str(compose.MODULE_PATH)] if hasattr(compose, "MODULE_PATH") else ["compose_to_grhayl.py"]
            with contextlib.redirect_stderr(stderr), self.assertRaises(SystemExit) as exit_context:
                runpy.run_path(str(Path(compose.__file__)), run_name="__main__")
            self.assertEqual(exit_context.exception.code, 2)
        finally:
            os.sys.argv = old_argv


if __name__ == "__main__":
    unittest.main()
