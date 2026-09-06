---
name: trialectic
description: Use for explicit GRHayL trialectic, three-seat independent review, "use tri", "run tri", or "engage tri", or when repository policy requires it. Mentioning, auditing, or editing this skill does not activate it. Three independent seats; delta-and-impact follow-ups; verified /work/ delivery.
---

# GRHayL Trialectic

Read [the shared protocol](../../review-protocol.md) with this skill; it owns modes,
independence, coverage, budgets, validation, and delivery. Distribute both and the
[GRHayL reference](../../grhayl-evidence.md); load only relevant reference sections,
not the other skill.

Use exactly three independent seats, excluding root. Root owns live writes and
delivery; seats never delegate or edit live targets. Default to `review`, including
documents and all KB work. Use `design` or `tri-build` only under protocol criteria.

## Roles

**Physics, numerics, and behavior.** Check relevant GRMHD/EOS/leakage/reconstruction/
induction assumptions, units, indices, signs, densitization, face/stagger ownership,
bounds, failure paths, stability, tolerances, and scientific oracles against source.

**API, build, and integration.** Check headers, structs/enums, definitions,
manifests/configuration/HDF5 guards, dispatch, installed exposure, compatibility,
consumers, and Core/gem boundaries. Prefer existing machinery and the simplest
sufficient change.

**Tests, docs, and downstream evidence.** Trace criteria to evidence; check fixture/
assertion strength, generated provenance, Doxygen/KB routes, CI selection, packaging,
and GRHayLib CCL/build/lifecycle impact. Catch coverage/delivery omissions and claims
that outrun evidence; do not audit unrelated repository areas.

For non-code work, use a first-principles analyst, a domain/implementation expert,
and an adversarial evidence editor. Roles add complementary scrutiny, not duplicate
suites.

## Completion

Follow-ups cover deltas, open findings, and affected code/docs, retaining valid
coverage. Root must [deliver and verify](../../review-protocol.md#delivery) ALL
intended authorized `/work/` destinations. A verdict, `DRAFT COMPLETE`, or exhausted
budget is not delivery. Install acceptable work; preserve blocked drafts separately.

If the protocol is missing, preserve drafts/companions in fresh
`/work/review-results/<task>/unapproved/` outside active source/KB/instruction discovery;
verify copies and report the missing file, never invent approval. Explicit outputs/
scoped write limits control; report inaccessible storage and actual retained paths.
