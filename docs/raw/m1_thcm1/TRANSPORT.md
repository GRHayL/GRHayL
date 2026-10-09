# THC_M1 transport fixtures

The retained sibling package contains the three
`transport_four_point_d*.dat` shards. Its manifest marks the package historical
and unadmitted; the receipt, source snapshot, command, raw inputs, and producer
outputs needed to verify the campaign binding were unavailable in the reviewed
workspace. The campaign counts and route details below are preserved from
historical family notes, not current admission evidence. Unit tests do not
invoke the producer. `.github/run_tests.sh` does not fetch M1 data; replay
requires an explicitly supplied `M1_FIXTURE_DIR`.

- operation: `neutrino_four_point_transport_flux`
- retained records: 57384
- exact duplicate input pairs removed: 6120
- lossless physical-direction shards: {"0": 19128, "1": 19128, "2": 19128}
- comparison policy: `strict_relative_2e-12_propagated_response_v1`
- historical note inventory: {"constant_volume": {"agreement": 63504, "difference": 0}, "varying_volume": {"agreement": 21384, "difference": 42120}}
- input admission: every constant-volume common-input pair, with only exact
  duplicate 50-field baseline/perturbed input tuples removed before outcome
  inspection; the complete tuple includes all operation controls and cell
  volumes
- exact branch inventory: {"direction": {"d0": 19128, "d1": 19128, "d2": 19128}, "mindiss": {"md0.0": 28692, "md1.0": 28692}, "opacity": {"packet": 19800, "transition": 18792, "zero": 18792}, "profile": {"constant": 16056, "monotone": 20664, "sawtooth": 20664}, "resolution": {"r128": 14346, "r256": 14346, "r32": 14346, "r64": 14346}, "theta": {"th0.0": 19128, "th1.0": 19128, "th2.0": 19128}}

The historical family notes report 21,384 varying-volume agreements and
42,120 differences. They attribute those counts to one campaign audit; the
receipt and raw streams were unavailable for independent review. The
variable-volume fixture stores the outputs present in the retained payload;
current GRHayL outputs are not used to regenerate them.

## Variable-volume prepared transport

The variable-volume corpus is a separate, versioned operation and file family:
`neutrino_four_point_prepared_transport_flux` in
`transport_four_point_varying_d*.dat`.  The three direction shards contain
only unique `changed_input` pairs so the prepared-transport consumer can
evaluate both current endpoints and compare the resulting response.  The unique
`input_invariant_control` pairs are retained separately in
`transport_four_point_varying_controls.dat`; they are fixed-consumed-input
controls, not perturbation-sensitivity evidence. Current GRHayL is evaluated at
both baseline and perturbed inputs; the stored THC output remains unchanged.

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


The historical exporter contract in the family notes requires a completed
campaign receipt, checks its `source_inputs_match` result and recorded
constant/varying inventory, binds the final command receipt's stdout to the
original raw output, and verifies the common input stream's deterministic
`fg -> ft` transform. The documented common command receipt also binds its
`argv`, transformed `stdin`, common `stdout`, and successful return code. The
exporter, receipt, and referenced raw streams are not staged with the retained
package, so these requirements have not been replayed here. The notes describe
an offline common producer using retained `current_official` transport source;
that source was not available for this review.

The historical notes name `current_official/transport.cc` as the producer
template and an external `THC_M1/src/thc_M1_calc_fluxes.cc` face fragment.
These sources and their build records are not available in the reviewed
workspace. Regeneration and admission remain pending recovery of those records.

No admitted package is currently pinned. A future TestData revision should
contain only admitted portable fixture files and compact package metadata; it
must not rewrite unrelated shards.

The historical notes report that the campaign producer passed theta into the
GRHayL initializer's `zeta_min` position. Its GRHayL result is therefore not
evidence that campaign theta reached that implementation's limiter. The notes
also report that the THC call received theta directly and that the exporter
retained THC results. The current unit-test adapter initializes valid
parameters and assigns theta to `minmod_theta` explicitly. The `th0.0`,
`th1.0`, and `th2.0` labels consequently exercise three distinct limiter
controls in current replay against the retained outputs.

The historical notes describe a stronger two-state audit with the exact
strict-relative `2e-12` endpoint/response relation,
`abs((A1-A0)-(T1-T0))/(den0+den1) <= 2e-12 + 4*DBL_EPSILON`, with
`denX=max(1e-300,abs(AX),abs(TX))`. The normal Unit_Tests replay evaluates the
current code at both baseline and perturbed inputs and applies the paired
endpoint/response comparator to current GRHayL and retained THC outputs.

Separate historical helper notes report 97 Rusanov agreements and one Jthick
difference. The underlying helper receipt and raw records were not available
in the reviewed workspace. Jthick is not part of this four-point transport
fixture; the reported mismatch remains a recorded limitation, not an admitted
campaign result.
