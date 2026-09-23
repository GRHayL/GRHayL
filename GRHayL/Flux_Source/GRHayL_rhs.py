import argparse
import re
import sys
from pathlib import Path

import sympy as sp

script_dir = Path(__file__).resolve().parent
nrpy_dir = str(script_dir / "nrpy")
if str(script_dir) not in sys.path:
    sys.path.insert(0, str(script_dir))
if nrpy_dir not in sys.path:
    sys.path.insert(0, nrpy_dir)

import IGM_All_fluxes as fl
import IGM_Characteristic_Speeds as chsp
import IGM_All_Source_Terms as st


INCLUDES = ["ghl_flux_source.h"]
VARIANTS = ("hybrid", "hybrid_entropy", "tabulated", "tabulated_entropy")
SUPPORTED_SYMPY_VERSION = "1.11.1"


def manifest_sources(directory):
    text = (directory / "make.code.defn").read_text()
    return set(re.findall(r"\b[A-Za-z0-9_]+\.c\b", text))


def expected_outputs():
    expected = {Path(name) for name in manifest_sources(script_dir)}
    for variant in VARIANTS:
        expected.update(Path(variant) / name
                        for name in manifest_sources(script_dir / variant))
    return expected


def prepare_destination(raw_destination):
    destination = raw_destination.expanduser().resolve()
    if destination == script_dir or script_dir in destination.parents:
        raise ValueError("output directory must be outside GRHayL/Flux_Source")
    if destination.exists():
        if not destination.is_dir():
            raise ValueError("output path exists and is not a directory")
        if any(destination.iterdir()):
            raise ValueError("output directory must be empty")
    else:
        destination.mkdir(parents=True)
    return destination


def generate(destination):
    outcparams = "outCverbose=False,GoldenKernelsEnable=True"
    st.Cfunction__GRMHD_SourceTerms(
        destination, includes=INCLUDES, formalism="ADM",
        outCparams=outcparams)
    chsp.Cfunction__GRMHD_characteristic_speeds(
        destination, includes=INCLUDES, formalism="ADM",
        outCparams=outcparams)

    for variant in VARIANTS:
        variant_dir = destination / variant
        variant_dir.mkdir()
        fl.Cfunction__GRMHD_fluxes(
            variant_dir, variant, includes=INCLUDES, formalism="ADM",
            outCparams=outcparams,
            tabulated="tabulated" in variant,
            entropy="entropy" in variant)

    generated = {path.relative_to(destination)
                 for path in destination.rglob("*.c")}
    expected = expected_outputs()
    if generated != expected:
        missing = sorted(str(path) for path in expected - generated)
        extra = sorted(str(path) for path in generated - expected)
        raise RuntimeError(f"generated output does not match manifests; "
                           f"missing={missing}, extra={extra}")


def main():
    if sp.__version__ != SUPPORTED_SYMPY_VERSION:
        raise RuntimeError(
            f"Flux_Source generation requires SymPy {SUPPORTED_SYMPY_VERSION}; "
            f"found {sp.__version__}. Install GRHayL/Flux_Source/requirements.txt")
    parser = argparse.ArgumentParser(
        description="Generate Flux_Source C into a new or empty staging tree")
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    generate(prepare_destination(args.output_directory))


if __name__ == "__main__":
    main()
