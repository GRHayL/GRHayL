# GRHayL Shared Review Protocol

Read once per active version with the selected `SKILL.md`, which supplies roles
and seat count. Root owns synthesis, checks, live writes, and delivery. Active
[AGENTS.md](../AGENTS.md) and owning policies govern; drafts cannot change their own
review rules. Authority/evidence outrank votes. Consult
[GRHayL evidence](grhayl-evidence.md) only for applicable claims/checks.

## Brief And Checkpoint

Record outcome, criteria, non-goals, assumptions, authority, mode, and deliverable:
implementation, document/artifact, or findings only. Map destinations/companions and
create/update/delete/rename operations. Discover the checkout; do not assume
`/work/GRHayL` or equate live state with `HEAD`. Capture live bytes/absence, types,
relevant modes/symlink targets, and dirty/untracked baselines.

Keep compact state in context or an owned task note when persistence is needed and
permitted: isolated areas, owners/dependencies, frozen candidates, each seat's
reviewed snapshot and accepted coverage, open finding IDs, check inputs/results,
mechanical suffix, counters, and deliveries. Separate live-write and review baselines.
Retain needed state across handoffs/compaction; bind evidence by frozen copies/direct
comparison, not labels alone, source hashes, or KB logs.

Freeze versions, not discovery. Follow needed read-only dependencies; synchronize
new deciding evidence before acceptance. Baseline new authorized writes; withhold
unauthorized changes while finishing separable work. Steering updates affected
criteria/coverage, not valid completed work or active counters.

## Preparation And Independence

Default to one root-prepared `review` candidate, including documents and KB work.
Findings-only review examines existing state without inventing edits/report files.
Add at most one `design` proposal OR independent-build wave only for unresolved
architectures or useful correctness evidence. Production seats return complete
isolated proposals/candidates with assumptions, reused machinery, checks, and limits.

Freeze contributions; select the best coherent result against current criteria.
Combine only compatible evidence-supported improvements, not obligatory parts from
every draft. Use fresh reviewers independent of authors, including for synthesized
designs. Explicit proposal-only/draft scope controls; disclose actual coverage.
`PROPOSAL COMPLETE`/`DRAFT COMPLETE` ends a phase, not delivery. Plans stay documents.

Complete intended formatting, fixers, authorized regeneration, and required checks
before approval review; fix known failures first. Diagnostic review may examine
failures, not call them passes. Mutating validators use expendable identical copies;
adopted edits are deltas, never silent mutations of frozen candidates.

Each substantive proposal/build or review wave uses the full seat count, excluding
root, with no more than that many agents concurrent. Initial contexts are separate
and fresh, with equal objective evidence and complementary roles. Parallelize when
supported; sequential seats still need separate contexts. Root alone launches agents.
Neither root nor seats recursively invoke either skill; seats never delegate or
write live targets. Use owned isolated areas; disclose instruction-only isolation.

Keep peer outputs/votes/arguments out of reviewer threads, including follow-ups.
Freeze reports before synthesis; share defect evidence neutrally. Reuse healthy
uncontaminated threads; replacements get objective checkpoints and their seat's
scope, not peer debate. Missing agents, silence, timeouts, or role headings in one
context cannot establish independent approval. Finish useful preparation and
Delivery when independence is unavailable.

## Decisions And Coverage

Review reports identify candidate/baseline, scope, retained coverage, evidence, and
findings or explicit none; end with `DECISION: ACCEPT` or `DECISION: BLOCK`. Findings
need ID, location, violated criterion/rule, evidence, consequence, and smallest remedy
or missing check. Omit repeated briefs/files/closed findings; preference is nonblocking.

Root deduplicates against authority. Never relabel `BLOCK`, outvote substantive
defects, count conditional/stale/silent approval, or replace dissenters. Return
disputed objections to the same seat with deciding evidence in a counted focused
review, not recovery. Do not re-poll covered decisions without changed
candidate/evidence.

Normal acceptance needs every seat's valid acceptance: unaffected prior coverage
plus accepted deltas must cover all current criteria and integration boundaries.
Each seat checks for missed dependencies. Close substantive blockers and required
proof gaps. Mechanical Finalization is the only default final-byte exception;
budgets cannot manufacture acceptance.

## Delta And Impact

Without reliable prior coverage, review the requested candidate and integration
boundaries once, not the whole repository. Later waves and explicit follow-ups
cover cumulative changes since EACH seat's last reviewed snapshot, open findings,
and affected code/docs. Distinguish accepted coverage from reviewed/rejected versions
and the lead-verified mechanical suffix. Saving/committing is not approval; include
that suffix in the next semantic review.

Supply a compact delta packet: old/new snapshots, exact hunks/sections and operations,
changed evidence, open finding IDs, invalidated criteria/checks, impact boundaries,
and new/reused results, not prior debate. Use [Change Impact](../wiki/change-impact.md).
Trace relevant callers, interfaces/shared state, manifests/configuration, generators/
products, tests/oracles, docs/links/claims, and downstream consumers until unchanged
contracts bound effects. Do not inspect every category by default.

| Change | Follow-up |
| --- | --- |
| Relevant content, evidence, assumptions, and environment unchanged | Reuse valid coverage/checks; deliver without a wave. |
| Exact local meaning-preserving correction | Eligible mechanical finalization and affected checks; no agent loop. |
| Changed behavior, claim, interface, instruction, authority, or proof | Full seat count reviews only the delta, open findings, and affected contracts. |
| Wider or uncertain effects | Investigate named dependencies; expand affected coverage only. |

One sign, tolerance, EOS/API guard, claim, link target, or prompt rule can be semantic;
line count is not a risk measure. Unaffected seats confirm retained coverage and
impact boundaries instead of repeating audits. Adding a seat or changing roles
requires missing assigned coverage once: a new seat independently examines its full
assigned scope, not just the latest diff. Existing seats retain valid reports;
reconfirm only changed coverage/impact boundaries. Keep findings and counters;
never drop changes a seat has not seen.

Reuse evidence only while transitive inputs, authority, assumptions, configuration,
and environment remain valid; unchanged target bytes alone are insufficient.
Recover missing state or review unsupported areas, never invent approval. Do not
restart design/build or reopen settled choices without new evidence. Full re-review
needs explicit direction or demonstrated invalidation of the entire prior scope.
Earlier candidates must still satisfy the current request.

## Mechanical Finalization

After any completed review wave, including the first, root may apply ONE terminal
mechanical batch per invocation without delegation. A later explicitly requested
mechanical follow-up may reuse valid reports/coverage without a wave; resuming the
same request does not reset the allowance. Keep the cumulative suffix.

All seat reports and criterion coverage must exist. Each seat must have accepted
or blocked solely on qualifying mechanical findings. Every remaining defect must
be local, deterministic, meaning-preserving, prescribed by existing evidence, and
conclusively checkable. No substantive blocker, unrelated disputed objection, or
required proof gap may remain. At least one qualifying edit must occur; empty or
unrelated cleanup cannot excuse missing approval.

Unambiguous spelling/behavior-neutral formatting may qualify. Science, numerical
behavior, interfaces, dependencies, instruction/reference meaning, authority,
validation obligations, and downstream contracts do not. Uncertainty needs delta
review. Apply once, freeze, inspect the diff, run affected required checks, and verify
retained coverage. Preserve actual decisions; report `LEAD-FINALIZED` and name edits
not independently re-reviewed. Failure/new semantic issues need remaining counted
delta review or blocked delivery, never a second unreviewed batch.

## Budget And Recovery

Aim for one review and a delta wave only when needed. Three review waves are the
hard maximum per invocation; a third needs progress, a named issue, bounded correction/
evidence acquisition, and deciding check. Count before launch. At most one proposal/
build and one recovery wave are additional: five delegated waves maximum, not a
target. Resumption, steering, replacement, and skill switching never reset active
budgets. Batch fixes; no peer debate, reviewer shopping, speculative builds, or
open-ended edit/check loops.

One recovery wave may repair malformed reports or replace failed seats on unchanged
inputs; retire failed occupants and retain healthy results. Substantive `BLOCK`
is not failure; healthy-seat reconsideration is review, not recovery.

Stop deliberation on recurring blockers without new evidence/credible correction,
unavailable required proof, or exhausted budgets; then perform [Delivery](#delivery).
No automatic reinvocation to evade caps. A later explicit request may authorize
another bounded invocation; retain evidence/findings and review its delta/impact.

## Proportionate Validation

Follow relevant owners, [Workflows](../wiki/workflows.md),
[Test Map](../wiki/test-map.md), and applicable [GRHayL evidence](grhayl-evidence.md).
Keep mandatory checks. Test changed contracts/plausible regressions, not mirrored
low-impact edits or unrelated suites. Root shares execution results; seats add
independent reasoning/oracle evidence, not duplicate runs.

Record snapshot, command/tool version, cwd, relevant inputs/configuration/environment,
exit status, decisive assertions/results, and limits. Reuse valid passes; repeat/
broaden only for changed inputs, failures, dependencies, uncovered contracts, or
required gates. Then deliver. Allow one diagnosed transient retry per invocation,
never an unchanged deterministic failure. Corrected-input checks are new validation;
root checks cannot replace semantic review.

Inspect command effects. Shared trees permit scoped side-effect-free checks; use
owned disposable trees for builds, output-producing runs, fixtures, installs,
generators, broad formatters, and mutating validators. Bound resources and own
outputs/caches; never clean shared trees/ambient files. External actions and
downstream changes remain within existing authority.

## Delivery

Run before EVERY planned exit involving produced files, including failed review,
missing agents/checks, drift, exhausted budget, and failed installation. Reserve
effort for copying/verification. Select the best coherent version by CURRENT
criteria/evidence, not recency, length, or votes. No acceptable version means
preservation, not forced acceptance.

**Preflight.** Explicit output paths/scoped write restrictions control; otherwise
use owning targets in the actual `/work/` checkout. Unowned requested documents go
to fresh `/work/review-results/<task>/` outside KB and active instruction locations,
preferably outside the checkout. Map ALL companions/operations; preflight coupled
changes together. Check authority, coverage/checks, snapshots, parents, types/modes,
and resolved paths. Create authorized missing parents; absent `/work/` alone need
not block. Symlinks cannot redirect unauthorized writes. Artifact-only requests
deliver the artifact without live installation.

Immediately compare live targets/deciding inputs with baselines; ignore unrelated
drift. Reconcile overlap only when authorized/unambiguous, in isolation, then
review/check that delta within budget. Conflicts or missing coverage block
installation, not preservation. Never discard concurrent work or revert it to `HEAD`.

**Install or preserve.** Use the applicable branch:

- Accepted or qualifying `LEAD-FINALIZED` work: apply ALL authorized files, companions,
  deletions, and renames. Deliver documents without implementing plans. Explicitly
  requested unfinished drafts may use designated draft paths with honest status,
  not release approval. Preservation cannot substitute for authorized installation.
- Blocked/rejected/incomplete work: copy the best draft/companions with relative
  layout to fresh `/work/review-results/<task>/unapproved/`, outside live source, KB,
  and active instruction/skill discovery. If instruction drafts could be discovered
  there, use an inert archive with intended paths recorded. Never activate rejected
  instructions or apply blocked deletions. An adjacent task status note records
  blockers, intended paths, pending checks, and unapplied operations without altering
  candidate bytes. Mark useful partials `INCOMPLETE`; say when no draft exists.
- Accepted work whose installation fails: preserve in the task's `not-installed/`
  area as `NOT INSTALLED`, not falsely `UNAPPROVED`; report partial live writes.

Plan/review-only scope does not authorize implementation; findings-only work needs
no unsolicited report file. A ban on repository writes still permits an otherwise
authorized external artifact; a blanket no-write instruction forbids scratch/
preservation too. Respect explicit alternative destinations.

**Verify.** Use existing recoverable file/hunk operations, not whole worktrees,
`.git`, caches, or competing drafts. Preserve types, modes, symlinks, and deletions;
verify in place when already matching. Read back/compare EVERY installed/preserved
file with its snapshot, including required absence/metadata; for archives compare
extracted members. Dirty-file merges must match the complete reviewed post-merge
snapshot. Inspect scoped diffs/unexpected outputs. Recheck destination-sensitive
contracts, not whole suites/review loops merely for identical copying.

Retain the candidate until ALL destinations verify; copy exit status alone is
insufficient. On partial failure, report exact live state and use only authorized
bounded recovery, never destructive rollback. If no authorized `/work/` location
is usable, retain the accessible candidate and report failed paths/errors, actual
retained location, and incomplete delivery. Finish separable handoffs; name
undelivered targets and do not loop on storage blockers or claim success.

## Completion Report

Separate review from delivery: accepted initial-plus-delta coverage, `LEAD-FINALIZED`,
or blocked/incomplete review; installed/document/artifact paths, preserved paths,
`REVIEW ONLY`, or undelivered targets. Include checks run/reused/missing, blockers,
partial writes, isolation limits, unreviewed mechanical edits, and compact wave/
recovery/retry counts. Keep receipts in the checkpoint, not transcripts or KB logs.
Do not claim measured model/token gains without evaluation.
