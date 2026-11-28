# Language For Computational Physicists

> Translate vague software-management wording into the direct language used by computational physicists. · Status: confirmed
> Up: [GRHayL Knowledge Base](../AGENTS.md)

## Summary

Write prose by naming the equation, physical quantity, numerical method, C struct, function, source file, or test result under discussion. This rule governs explanations read by people or agents across KB pages, documentation, comments, CLI help, diagnostics, and messages: state the mathematical problem, identify inputs and outputs, specify numerical assumptions, and name the compiled or tested C/C++/Python code directly.

## Detail

### Start From The Computation

Before writing a sentence, identify its concrete subject in GRHayL's technical vocabulary: primitive variable, conserved variable, ADM 3+1 spatial metric, lapse, shift, four-velocity, Lorentz factor, magnetic field, four-magnetic field, stress-energy tensor, hybrid or tabulated equation of state, chemical potential, fermion blocking factor, Riemann flux, characteristic speed, reconstruction stencil, vector-potential staggering, Lorenz gauge RHS, neutrino leakage rate, C struct, header declaration, recursive manifest, Cactus thorn parameter, unit-test driver, reference data file, and comparison tolerance.

State directly what is initialized, packed, unpacked, interpolated, inverted, reconstructed, evolved, rescaled, compiled, compared, or written, giving tensor indices, units, sign conventions, floor or ceiling bounds, inputs, outputs, and failure modes where relevant.

### Translation Table

Choose the specific physical, numerical, or file-level noun that matches the computation:

| Avoid in explanatory prose | Prefer |
| --- | --- |
| artifact | generated file, output file, compiled binary, library archive, reference table, test result, build output, or named source file |
| contract | required behavior, interface definition, mathematical assumption, numerical requirement, or expected result |
| schema | file format, data layout, field definitions, parameter list, table columns, or page layout |
| surface | public functions, public headers, command-line options, runtime parameters, or host interface |
| payload | values, input data, output data, struct fields, or array contents |
| plumbing | call sequence, data path, dispatch path, build steps, or input/output path |
| pipeline | ordered calculation stages, sequence of routines, or named numerical method |
| bridge | conversion function, host-transfer array, coupling term, or named interface |
| bundle | group of variables, component array, file set, or another named collection |
| gate | required check, acceptance condition, comparison, or test |
| harness | test program, test driver, comparison script, or test runner |
| fixture | test input, reference data, initial data, sample table, or expected output |
| consumer | calling function, downstream thorn, host code, reader, or downstream calculation |
| owner | defining module, defining leaf, implementing file, equation module, or responsible test |
| lifecycle | initialization, allocation, update, and cleanup steps, naming only the steps that apply |

Do not replace one vague noun with another. Rewrite "the test harness validates the EOS adapter contract" as a checkable statement such as "`unit_test_tabulated_eos_compose.c` checks that `NRPyEOS_read_table_set_EOS_params` loads a converted CompOSE HDF5 table and reproduces reference enthalpy and sound speed."

### Exact And Mathematical Names

Preserve exact code identifiers and established mathematical or physical terminology: "contract the upper index with the first lower index," apparent-horizon or ADM two-surfaces, normal or tangent bundles, named algorithms, and explicit C heap-memory or POSIX file ownership when the sentence names what is owned. Exact filenames, GitHub Action names (`actions/upload-artifact`), C macro or struct names, and direct quotations stay unchanged inside code spans or quotation marks. Outside exact names and quotations, reserve "artifact" for its archaeological meaning.

### Required Review And Conciseness

For every new or revised passage:

1. Replace left-column software metaphors (including inflections such as `owning`, `owned`, `owns`, `ownership`) with the concrete mathematical, numerical, or file-level subject.
2. Eliminate prose paragraphs that restate adjacent Markdown table cells; put per-row evidence directly in table cells with compact ``(Ev: `file.c::symbol`)`` tags when needed.
3. Omit repeated generic layered-availability disclaimers already covered by [KB Rules](SCHEMA.md#volatile-information-policy); state only module-specific unbuilt, unwired, or untested paths.
4. Compare recurring terminology with [Glossary](glossary.md) and review meaning rather than performing blind global word replacements.

Evidence (normative rule): `AGENTS.md::Human-Facing Language`.

## Sources

- [AGENTS.md](../AGENTS.md) - `## Human-Facing Language`
- [docs/raw/derivation.md](../docs/raw/derivation.md) - GRMHD variables, stress-energy notation, and evolution equations

## See Also

- Depends on: [Glossary](glossary.md)
- See also: [KB Rules](SCHEMA.md)
- Validated by: [Workflows](workflows.md#deterministic-and-review-checks)
