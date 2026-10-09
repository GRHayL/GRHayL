# Fresh THC_M1 fixture production record

This page records the scoped fresh-THC_M1 production campaign behind the
Radiation M1 retained fixtures and its evidence limits. The removed
historical package had no verifiable producer-to-payload binding.
The accepted replacement was produced from fresh THC_M1 executions for all
five operation groups: pointwise, instantaneous source, covariant stress-energy,
Rusanov, four-point transport, and the standalone Jthick scalar. The candidate
is published in `GRHayL/TestData` as individual gzip-compressed `radiation/*.bin.gz`
members. The runner downloads those members at the full commit selected by
`.github/radiation-testdata-ref` and prepares them by default; the `run_m1`
action forwards its optional `fixture-dir` input to the runner as
`M1_FIXTURE_DIR` and performs no acquisition of its own. A nonempty
`fixture-dir` (or
`M1_FIXTURE_DIR` for the ordinary runner) instead replays a supplied raw or
gzip directory and skips downloads entirely.
The external producer-evidence archive's `fixture/README.md` maps every payload
to the captured producer records and exporter. Its manifest covers 14 payloads
and includes Jthick.
The published fixture families and binary layout are described by the maintained
[TestData provenance](https://github.com/GRHayL/TestData/blob/main/radiation/PROVENANCE.md).

## Scoped execution

The builds compiled Radiation and its required Core/Rusanov support, the
specific THC_M1 fragments, and the five fixture-consuming unit executables.
No Einstein Toolkit or whole-toolkit build was run. The pointwise, source, and
Jthick campaigns used the offline verifier executables captured in the
evidence archive's base snapshot. CL-04 stress-energy used the compiled
`assemble_rT` direct route. Rusanov used
the THC_M1 closure/face-flux fragments and recorded offline producer inputs.
Transport used the frozen `run_transport.py` captured in the evidence archive
because the current script adds a profile and changes the historical input
inventory. The fresh transport receipt records 63,504 constant-volume
agreements, 21,384 varying-volume agreements, and 42,120 varying-volume
differences; the exporter retained the documented agreeing operations.

The standalone producer evidence archive contains each campaign's completed
receipt, executable command, consumed input stream, raw THC_M1 output, source
snapshot, and exporter. The Rusanov receipts capture both the legacy and
current API runs; the Jthick command receipt captures a rerun of the THC_M1
helper. The package auditor reads each payload in full and verifies its
grammar, operation, policy, record count, and payload-only digest. Source files
are identified by their captured snapshots, not by stored source hashes.

## Replay and limits

The package auditor passed for 14 payloads and 123,702 records. The local
TestData preparation converts these 14 text payloads to versioned `.bin`
files without changing their recorded fields or binary64 values, and the
published package stores those members gzip-compressed. The GRHayL fixture
reader accepts the binary files while retaining text support for parser
checks.

The CI action runs the ordinary pinned download route by default: no
`fixture-dir` input is required. A nonempty `fixture-dir` instead selects the
supplied raw or gzip package and skips downloads. The runner supplies raw
members to the test executables on both routes, and external acquisition fails
rather than silently omitting replay if the pinned revision is unavailable.
The `.github/run_tests.sh m1` route refreshes all named `.bin.gz` members from
raw.githubusercontent.com at the revision in `.github/radiation-testdata-ref`
on every run and prepares them in a private `${TMPDIR:-/tmp}` directory
removed by its EXIT cleanup. Supplying `M1_FIXTURE_DIR` bypasses downloads and
replays the supplied raw or gzip package, preserving caller inputs.
Local preparation and replay do not establish a remote CI pass.
These fixtures establish the named discrete operations, not Cactus host
evolution or continuum validation.
