# THC_M1 transport fixtures

The published TestData package contains the three
`transport_four_point_d*.bin.gz` constant-volume shards, which expand to the
`transport_four_point_d*.bin` files the tests read. The values come from the
fresh transport THC_M1 campaign, whose completed receipt, source snapshot,
command, raw inputs, and producer outputs are in the producer evidence outside
this checkout. Unit tests do not invoke the producer. The
`.github/run_tests.sh m1` runner downloads the pinned members by default;
`M1_FIXTURE_DIR` supplies a local raw or gzip package instead.

- operation: `neutrino_four_point_transport_flux`
- retained records: 57384
- exact duplicate input pairs removed: 6120
- lossless physical-direction shards: {"0": 19128, "1": 19128, "2": 19128}
- comparison policy: `strict_relative_2e-12_propagated_response_v1`
- campaign receipt inventory: {"constant_volume": {"agreement": 63504, "difference": 0}, "varying_volume": {"agreement": 21384, "difference": 42120}}
- input admission: every constant-volume common-input pair, with only exact
  duplicate 50-field baseline/perturbed input tuples removed before outcome
  inspection; the complete tuple includes all operation controls and cell
  volumes
- exact branch inventory: {"direction": {"d0": 19128, "d1": 19128, "d2": 19128}, "mindiss": {"md0.0": 28692, "md1.0": 28692}, "opacity": {"packet": 19800, "transition": 18792, "zero": 18792}, "profile": {"constant": 16056, "monotone": 20664, "sawtooth": 20664}, "resolution": {"r128": 14346, "r256": 14346, "r32": 14346, "r64": 14346}, "theta": {"th0.0": 19128, "th1.0": 19128, "th2.0": 19128}}

The campaign receipt reports 21,384 varying-volume agreements and 42,120
differences. Those counts are the outcomes of the campaign's own audit; they do
not select records. Both constant and varying corpora retain every unique input
pair (6,120 and 1,512 exact duplicates are the only removals), and the
variable-volume fixture stores the exact THC_M1 outputs from the separately
transformed common-input stream. Current GRHayL outputs are not used to select
or regenerate them.

## Variable-volume prepared transport

The variable-volume corpus is a separate, versioned operation and file family:
`neutrino_four_point_prepared_transport_flux` in
`transport_four_point_varying_d*.bin`.  The three direction shards contain
only unique `changed_input` pairs so the prepared-transport consumer can
evaluate both current endpoints and compare the resulting response.  The unique
`input_invariant_control` pairs are retained separately in
`transport_four_point_varying_controls.bin`; they are fixed-consumed-input
controls, not perturbation-sensitivity evidence. Current GRHayL is evaluated at
both baseline and perturbed inputs; the stored THC output remains unchanged.

- unique variable pairs: 61992
- changed-input pairs: 38880
- invariant-control pairs: 23112
- exact duplicate variable input pairs removed: 1512
- changed-input direction shards: {"0": 12960, "1": 12960, "2": 12960}
- invariant-control file: `transport_four_point_varying_controls.bin`
- variable branch inventory: {"direction": {"d0": 12960, "d1": 12960, "d2": 12960}, "mindiss": {"md0.0": 19440, "md1.0": 19440}, "opacity": {"packet": 15120, "transition": 11880, "zero": 11880}, "profile": {"constant": 12960, "monotone": 12960, "sawtooth": 12960}, "resolution": {"r128": 9720, "r256": 9720, "r32": 9720, "r64": 9720}, "theta": {"th0.0": 12960, "th1.0": 12960, "th2.0": 12960}}
- record schema: versioned `M1THCMB1` binary records, 50 complete transformed
  transport inputs, and five exact THC_M1 outputs per role
- volume fields: face volume at input field 0; four cell volumes at fields
  6--9; states at fields 10--29; common transformed physical fluxes at
  fields 30--39
- reference values: only the common producer's `thcm1` outputs are exported;
  current GRHayL outputs are never used as fixture values

The exporter requires a completed campaign receipt, checks its
`source_inputs_match` result and recorded constant/varying inventory, binds the
final command receipt's stdout to the original raw output, and verifies the
common input stream's deterministic `fg -> ft` transform. The common command
receipt also binds its `argv`, transformed `stdin`, common `stdout`, and
successful return code. The exporter, receipts, and referenced raw streams are
in the producer evidence outside this checkout, and the regeneration package
replays them in an offline common producer built from the frozen
`current_official` transport source.

The producer template is `current_official/transport.cc`, with a THC_M1 face
fragment compiled from `THC_M1/src/thc_M1_calc_fluxes.cc`. Both sources and
their build records are held in the producer evidence. The campaign uses the
frozen `run_transport.py` captured there, because the current script adds a
profile and produces a different input inventory.

The accepted members are pinned by `.github/radiation-testdata-ref`. A
replacement TestData revision should add or replace only the affected members
and must not rewrite unrelated shards; it takes effect when the pin is updated.

The producer initializes with a valid `zeta_min` and then assigns the campaign
theta field directly to `minmod_theta`. The current unit-test adapter does the
same: it initializes valid parameters and assigns theta to `minmod_theta`
explicitly. The `th0.0`, `th1.0`, and `th2.0` labels consequently exercise
three distinct limiter controls in current replay against the retained outputs.

The campaign audit applies the exact strict-relative `2e-12`
endpoint/response relation,
`abs((A1-A0)-(T1-T0))/(den0+den1) <= 2e-12 + 4*DBL_EPSILON`, with
`denX=max(1e-300,abs(AX),abs(TX))`. The normal Unit_Tests replay evaluates the
current code at both baseline and perturbed inputs and applies the paired
endpoint/response comparator to current GRHayL and retained THC outputs.

Jthick is not part of this four-point transport fixture; see
[README.md](README.md) for the standalone Jthick fixture and the helper
campaign's recorded Jthick difference.
