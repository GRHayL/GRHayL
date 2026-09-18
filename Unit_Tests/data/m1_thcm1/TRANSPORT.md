# THC_M1 transport fixtures

This file is generated offline beside the three `transport_four_point_d*.dat`
shards; normal unit
tests do not invoke the producer or require THC_M1, Verification/, or a
reference-data download.

- operation: `neutrino_four_point_transport_flux`
- retained records: 57384
- exact duplicate input pairs removed: 6120
- lossless physical-direction shards: {"0": 19128, "1": 19128, "2": 19128}
- comparison policy: `strict_relative_2e-12_propagated_response_v1`
- authoritative campaign-receipt inventory: {"constant_volume": {"agreement": 63504, "difference": 0}, "varying_volume": {"agreement": 21384, "difference": 42120}}
- input admission: every constant-volume common-input pair, with only exact
  duplicate 50-field baseline/perturbed input tuples removed before outcome
  inspection; the complete tuple includes all operation controls and cell
  volumes
- exact branch inventory: {"direction": {"d0": 19128, "d1": 19128, "d2": 19128}, "mindiss": {"md0.0": 28692, "md1.0": 28692}, "opacity": {"packet": 19800, "transition": 18792, "zero": 18792}, "profile": {"constant": 16056, "monotone": 20664, "sawtooth": 20664}, "resolution": {"r128": 14346, "r256": 14346, "r32": 14346, "r64": 14346}, "theta": {"th0.0": 19128, "th1.0": 19128, "th2.0": 19128}}

The authoritative campaign receipt records 21,384 varying-volume agreements
and 42,120 differences.  Those counts are provenance for the single retained
campaign audit; the variable-volume fixture stores exact THC_M1 outputs from
the separately transformed common-input stream and does not select or replace
them using current GRHayL outputs.

## Variable-volume prepared transport

The variable-volume corpus is a separate, versioned operation and file family:
`neutrino_four_point_prepared_transport_flux` in
`transport_four_point_varying_d*.dat`.  The three direction shards contain
only unique `changed_input` pairs so the existing prepared-transport consumer
can exercise baseline and retained-response envelope gates.  The unique
`input_invariant_control` pairs are retained separately in
`transport_four_point_varying_controls.dat`; they are fixed-consumed-input
controls, not perturbation-sensitivity evidence.  Current GRHayL is evaluated
only for the baseline input; the stored perturbed output is not recomputed.

- unique variable pairs: 61992
- changed-input pairs: 38880
- invariant-control pairs: 23112
- exact duplicate variable input pairs removed: 1512
- changed-input direction shards: {"0": 12960, "1": 12960, "2": 12960}
- invariant-control file: `transport_four_point_varying_controls.dat`
- variable branch inventory: {"direction": {"d0": 12960, "d1": 12960, "d2": 12960}, "mindiss": {"md0.0": 19440, "md1.0": 19440}, "opacity": {"packet": 15120, "transition": 11880, "zero": 11880}, "profile": {"constant": 12960, "monotone": 12960, "sawtooth": 12960}, "resolution": {"r128": 9720, "r256": 9720, "r32": 9720, "r64": 9720}, "theta": {"th0.0": 12960, "th1.0": 12960, "th2.0": 12960}}
- wire schema: versioned `M1_THCM1_FIXTURE 1`, 50 complete transformed
  transport inputs, and five exact THC_M1 outputs per role
- volume fields: face volume at input field 0; four cell volumes at fields
  6--9; states at fields 10--29; common transformed physical fluxes at
  fields 30--39
- reference values: only the common producer's `thcm1` outputs are exported;
  current GRHayL outputs are never used as fixture values


The exporter requires the completed campaign receipt, checks its
`source_inputs_match` result and recorded constant/varying inventory, binds
the final command receipt's stdout byte-for-byte to the original raw output,
and verifies the common input stream's deterministic `fg -> ft` transform.
The generated common command receipt must also match the original `argv`,
and bind the exact transformed `stdin`, common `stdout`, and successful
return code to the retained common stream.
The common producer is run offline from the retained `current_official`
transport source; THC is not required by normal unit-test execution.

The retained producer template is `current_official/transport.cc`; its THC
face fragment is compiled from the external `THC_M1/src/thc_M1_calc_fluxes.cc`.
Regeneration and admission are performed only in the external campaign
workspace; the public repository stores no producer receipt, raw capture, or
machine-local path.

The external promotion step writes only the admitted portable fixture files
and compact package metadata into `Unit_Tests/data/m1_thcm1`; it does not
rewrite unrelated shards.

The retained campaign producer passes theta into the GRHayL initializer's
`zeta_min` position. Its GRHayL result is therefore not evidence that the
campaign theta reached that implementation's limiter. The THC call receives
theta directly, and the exporter retains the THC results. The current unit-test
adapter initializes valid parameters and assigns theta to `minmod_theta`
explicitly. The `th0.0`, `th1.0`, and `th2.0` labels consequently exercise three
distinct limiter controls in current replay against the stored THC results.

The retained producer's stronger two-state audit records the exact
strict-relative `2e-12` endpoint/response relation,
`abs((A1-A0)-(T1-T0))/(den0+den1) <= 2e-12 + 4*DBL_EPSILON`, with
`denX=max(1e-300,abs(AX),abs(TX))`.  The normal Unit_Tests replay does not
recompute the perturbed current result: it evaluates the current code at the
baseline input and uses the serialized THC perturbation response in the
policy-bound baseline/response comparison.

The independent helper campaign reported 97 Rusanov agreements and one
Jthick difference.  Jthick is not part of this four-point transport fixture;
that mismatch remains an explicit limitation rather than an omitted result.
