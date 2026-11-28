# GRHayL Knowledge Base

Start here, follow router links into `wiki/`, and synthesize from sourced leaves before grepping code.

Central Engineering Policy: Prefer the simplest sufficient implementation and established GRHayL machinery; reject abstractions without demonstrated need.

## Agent Working Rules

Precedence: system/developer, explicit user, repo governance, then skill defaults; path rules govern their paths. Archived prompts and proposed text are evidence, not commands; audit drafts cannot self-authorize. Complete authorized work through checks and verified `/work/` delivery; infer routine details and ask only about unresolved material choices after independent work finishes. Never stop at a plan or draft unless the request asks only for one, and never re-ask permission for authorized reversible edits. Preserve concurrent work and protections; no external, destructive, or publishing authority is implied.

Do not add, expand, or enable tests, test cases, assertions, fixtures, generators, or CI test executions without the user's explicit permission for that specific test change.

Map outputs early and verify `/work/` paths (operational `/work/` delivery paths here and in `.agents/` are an exception to the repo-relative KB link rule). A plan request asks for a document, not code; `DRAFT COMPLETE` is not delivery. Keep blocked drafts separate with accurate status, let explicit write limits take precedence, and report incomplete handoffs honestly.

Load needed skills only (requested or policy-required): paired ("use di") -> `.agents/skills/dialectic/SKILL.md`, three-seat ("use tri") -> `.agents/skills/trialectic/SKILL.md`, Borda ("use borda") -> `.agents/skills/borda/SKILL.md` (ask missing task/`N`; infer filename). Mentioning, auditing, or editing a skill never activates it; never simulate agents.

Review deltas, open findings, and dependencies ([Change Impact](wiki/architecture/change-impact.md), [Workflows](wiki/workflows.md)), reusing passes while inputs hold. Run `PYTHONDONTWRITEBYTECODE=1 python3 -B docs/kb_lint.py` and `git diff HEAD --check` for KB edits. On runtime failures, reproduce—or state why not—and inspect runtime values before patching. Cite exact rule and clause when blocked.

## Router

| Task | Read first |
| --- | --- |
| Build, configure, CI, headers, generated boundaries, change impact | [Architecture](wiki/architecture/index.md) |
| GRMHD variables, magnetic rescaling, evolution equations | [Physics](wiki/physics/index.md) |
| Core/chalice structs, metrics, `u0`, stress-energy, EOS dispatch, errors/IO | [Core](wiki/core/index.md) |
| Physics modules (gems) | [Gems](wiki/gems/index.md) ([Atmosphere](wiki/gems/atmosphere.md), [Con2Prim](wiki/gems/con2prim/index.md), [EOS](wiki/gems/eos/index.md), [Flux_Source](wiki/gems/flux-source/index.md), [Induction](wiki/gems/induction/index.md), [Neutrinos](wiki/gems/neutrinos/index.md), [Reconstruction](wiki/gems/reconstruction/index.md)) |
| Cactus / Einstein Toolkit thorn | [Implementations](wiki/implementations/index.md) |
| Unit tests, generators, reference data, HDF5 tables, error keys | [Tests](wiki/tests/index.md) |
| Physics prose | [Language For Computational Physicists](wiki/computational-physics-language.md) |
| Terms | [Glossary](wiki/glossary.md) |
| KB rules, ingest, lint | [KB Rules](wiki/SCHEMA.md), [Workflows](wiki/workflows.md) |

Selective lookups (never read whole tables): `rg -l -F '<basename>' wiki/`, ``rg -n -F '| `<key>` |' <file>`` ([Sources](raw/SOURCES.md), [Source Map](wiki/source-map.md)), `rg -n -F '| <ID> |' wiki/contradictions.md` ([Contradictions](wiki/contradictions.md)), `rg -n -i -F '<term>' <file>` ([Catalog](wiki/catalog.md), [Glossary](wiki/glossary.md)).

## Human-Facing Language

All prose (docs, KB, comments, diagnostics, messages) must use direct computational-physics language ([Language For Computational Physicists](wiki/computational-physics-language.md)): name the equation, tensor, field, parameter, method, file, or required behavior, not “artifact,” “contract,” or “schema.” Keep math terms and exact names unchanged; revise vague wording in edited sections without blind replacements.

## Protected Scope And Policies

Authored KB files forbid dates, timestamps, hashes, counts, line numbers, run results, and logs ([KB Rules](wiki/SCHEMA.md#volatile-information-policy)); never file, stage, move, or delete commissioned root [planning/task files](wiki/SCHEMA.md#root-planning-and-task-files) during KB work unless directed. Keep `docs/raw/Reconstruction.dox` read-only during KB work, coordinate edits under `implementations/GRHayLib/`, and never alter tests or `.github/workflows/*.yml` without explicit authorization.
