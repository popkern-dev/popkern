# popkern

BRep geometry kernel in Rust. AGPL-3.0 license.

This file is the working guide for anyone changing popkern, human or
coding agent. `CLAUDE.md` imports it.

## Project Structure

```
crates/popkern-math        Portable transcendentals: same bits natively and in wasm (leaf crate)
crates/popkern-core        Topology, geometry, validation (depends: math only)
crates/popkern-constraint  2D constraint solver (depends: math only)
crates/popkern-bool        Boolean operations (depends: core)
crates/popkern-draft       Draft/taper (depends: core)
crates/popkern-blend       Fillets and chamfers (depends: core)
crates/popkern-fastener    Parametric fasteners: threads, screws, nuts, holes (depends: core, sketch, sweep, bool, blend)
crates/popkern-hinge       Parametric hinges: butt, print-in-place, piano, living (depends: core, sketch, bool, fastener)
crates/popkern-io          STEP/IGES import/export (depends: core)
crates/popkern-mesh        Tessellation (depends: core)
crates/popkern-sketch      Sketch-to-BRep bridge (depends: core, constraint)
crates/popkern-sweep       Sweep and loft operations (depends: core, bool)
crates/popkern-shell       Shell and offset operations (depends: core, blend, bool)
crates/popkern-measure     Measurement queries (depends: core, mesh)
crates/popkern-viz         Debug visualization + headless render + surface analysis (depends: core, mesh, measure)
crates/popkern-part        Code-first part API (depends: core..shell, io, mesh)
crates/popkern-generative  Generative design as a part feature: density design -> faceted body (depends: part, mesh, analysis, analysis-brep, analysis-setup)
crates/popkern-history     Parametric feature tree (depends: most of the above)
crates/popkern-script      Scripting-surface recorder (depends: core..history, part)
crates/popkern-assembly    Assemblies: part definitions + placed instances, STEP product structure, .pkn v5 (depends: script, io)
crates/popkern-rhai        Embedded Rhai skin over the recorder (depends: script, assembly)
crates/popkern-project     Projects: folders of ordered migration scripts → one document (depends: rhai, assembly)
crates/popkern-server      `popkern serve`: popkern live self-hosted, projects on disk + native sessions (depends: live, project)
crates/popkern-cli         `popkern run` / `build` / `new` / `serve` commands (depends: rhai, project, server, script, viz)
crates/popkern-live        Browser engine for popkern live, wasm-bindgen (depends: rhai, script, mesh, measure)
crates/popkern-examples    Showcases + cross-crate integration tests (depends: everything)
crates/popkern-analysis           Thermal + structural analysis solver (depends: math only)
crates/popkern-analysis-brep      BRep -> AnalysisModel adapter, face selectors (depends: analysis, core, mesh, measure)
crates/popkern-analysis-setup     Analysis setup file format, TOML (depends: analysis)
crates/popkern-analysis-cli       `popkern-analysis run part.step setup.toml` (depends: the above, io, viz)
crates/popkern-analysis-fixtures  Frozen STEP fixture generator for the physics suite
```

## Dependency Rule

Dependencies flow forward only: Math → Core/Constraint →
Bool/Draft/Blend/IO/Mesh/Sketch → Sweep/Shell/Fastener/Measure/Viz →
Hinge → Part → Generative → History → Script → Assembly → Rhai →
Project → Live → Server → CLI/Examples. No reverse imports. Enforced
by `tools/structural-tests/tests/structural_dep_direction.rs`.

The analysis crates sit beside this chain: the solver at the Math tier,
the adapter at the Measure tier, and only generative (the Part-tier
feature), live, server, cli and examples may depend on them
(`structural_analysis_dep_direction.rs`).

## Quality Gates (run before every commit)

```
cargo fmt --check
cargo clippy -- -D warnings
cargo test
```

All three must pass. CI enforces the same checks.

**Two test suites:** the kernel suite (`tools/quality-check.sh --full`,
the default nextest profile) leaves out the analysis crates but keeps
the adapter's unit tests as a smoke set; the physics suite
(`tools/quality-check.sh --physics`, its own `target/physics`) runs
them. Run the physics suite when an analysis crate changes, when its
frozen fixtures (`tests/fixtures/analysis`) are regenerated, and before
a push touching popkern-io, -mesh, -measure or -core.
`-p popkern-analysis` on the default profile selects nothing: use
`--profile physics`.

`rust-toolchain.toml` pins the toolchain, so local clippy and rustfmt
match CI. Bump it deliberately, with the new version's clippy findings
fixed in the same commit.

## Working Style

- Research the codebase before editing. Never change code you haven't read.
- Aim for commercial-grade kernel behavior over quick wins. Commercial-grade
  means: complete behavior over narrow demos; robust handling of degeneracy,
  tolerances, and numeric edge cases; correct by construction when possible, by
  explicit repair otherwise; measured acceptance via validation or numeric
  gates (not visual inspection or a single favorable fixture); semantic
  ownership preserved (faces, edges, loops, generated features keep their roles
  and stable names); rollback/containment for risky operators; explicit
  fallback boundaries that record why the stronger path didn't apply;
  observability for hard geometry failures (structured errors, reproducible
  cases, debug artifacts).
- Quick fixes are acceptable only when they are also the commercial-grade
  solution, or are explicitly recorded as temporary containment with a
  follow-up issue.
- Do not defer hard kernel work merely because it is difficult. If uncertainty
  blocks implementation, plan explicit exploration and then revise the design.
- Use semantic domain names that describe the CAD/kernel concept: intersection
  graph, tolerance band, loop classifier, edge provenance, surface trim domain,
  repair action, rolling-ball contact, apex patch. Do not let planning labels
  ("Path A", "Plan B", "Phase 2", `temporary_builder`, `new_new_solver`) become
  code, API, or long-lived document names.
- When reviewing, lead with gaps that prevent commercial-grade quality:
  incomplete reasoning, hidden temporary assumptions, missing tests, missing
  observability, ambiguous naming. Summaries are secondary to finding flaws.

## Code Rules

- No .unwrap() or .expect() in library code — use ? or typed errors (thiserror)
- No unsafe anywhere in the workspace
- No source file > 500 lines — split into submodules
- All pub items need /// doc comments
- All errors use #[derive(thiserror::Error)]
- Geometric tolerance is always an explicit parameter, never global state
- Transcendentals go through `popkern_math::PortableMath` (`x.psin()`,
  `y.patan2(x)`, also re-exported as `popkern_core::geometry::PortableMath`),
  never the inherent `f64::sin` family: the platform libm differs in the
  last ulp natively vs wasm, and one script would build two different
  parts. Enforced by `structural_portable_math` (tests exempt).
- No print-trace as the artifact of kernel learning. Findings from probes must
  be promoted to structured diagnostics, classifier output, or repair reports
  that tests can assert. Local text traces are fine for transient exploration
  but must be removed or promoted before commit.

## Design Invariants

1. Topology and geometry are separate layers (topology references geometry, not reverse)
2. Data model is a DAG — no ownership cycles
3. Boolean operations are transactional (succeed completely or fail)
4. Euler operators are the only topology modification path
5. Geometric tolerance is explicit (no global state)

**Every accepted topology/geometry change must state or encode:**
- Invariants it preserves
- Validation and postconditions it must pass
- Rollback or fallback behavior on failure
- Diagnostics produced when it declines or fails
- Why the solution is not fixture-specific
- How it preserves semantic ownership and stable handles/names
- Which tolerance is used and why that tolerance is local and explicit

## Testing

- 100% line coverage target for all crates
- Reference models in tests/fixtures/ with known-correct topology
- Geometric validation: Euler characteristic, manifold edges, vertex valence
- Structural tests (`tools/structural-tests/tests/structural_*.rs`) enforce
  architecture rules

## Rendering & Visualization

**External viewer (f3d):** validates STEP exports against OpenCASCADE:
```bash
f3d --output out.png --resolution 800,600 --edges --up +Z file.step
# Requires a display server or OSMesa.
```

**Built-in CPU renderer** (no GPU required):

```bash
# Render a showcase example
cargo run --bin boolean_showcase    # outputs examples/output/boolean_showcase.{step,png}

# Headless render any BRep (from code)
popkern_viz::headless::render_to_png(&brep, &camera, 800, 600, path, true)
# render_to_png auto-falls-back from GPU to CPU software renderer
```

If you change geometry output, render it and look at the result.

## Parameter Range Conventions

`SurfaceGeometry::parameter_range()` returns `(±∞, ±∞)` for planes/cylinders. Never pass these raw to the intersection solver. Three bounded alternatives exist in `popkern_core::geometry::uv_projection`:

| Function | Padding | Use case |
|---|---|---|
| `bounded_parameter_range` | 50%, min 1.0 (scaled down for faces under 0.5) | Boolean intersection solvers (need headroom) |
| `tight_bounded_parameter_range` | 5%, min 0.05 (scaled down for faces under 0.5) | Tessellation, AABB sampling |
| `exact_bounded_parameter_range` | 0% | Self-intersection validation |

All three take face vertex positions from `BRep::face_outer_loop_positions(face_handle)`.

## Known Gotchas

- **Seam edges** — edges where `face_1 == face_2` (e.g., cylinder seam). Can't be filleted. The fillet analysis rejects them, but other operations should also check.
- **`validate()` vs `validate_all()`** — `validate()` is the fast path (structural + geometric, no self-intersection). `validate_all()` adds the expensive O(n^2) self-intersection check. Don't call `validate_all()` in hot paths.
- **`find_starting_points` queue cap** — capped at 50k pairs. Near-coincident NURBS surfaces may produce incomplete results (missing branches). The coincident fast path in `intersect_surfaces` handles the common case.
- **`no_file_exceeds_500_lines`** — structural test enforced. When adding code, check `wc -l` before committing. Extract tests to `*_tests.rs` files with `#[path]` attribute if needed.

## Local Checks

- `tools/quality-check.sh [--fast|--full|--physics]` — shared quality gate
- `tools/check-semantic-names.sh` — flags planning labels in long-lived code
- `tools/health-check.sh [--pre-commit]` — health metrics and line budgets

## Task Classification (Adaptive Rigor)

Scale process to the blast radius of a change.

| Signal | Rigor | Process |
|--------|-------|---------|
| Bug fix, single file | Light | Fix, test, commit |
| Refactor within one crate | Medium | Check cross-crate impact, run full test suite |
| New public API surface | Heavy | Design review first, write a usage example, then implement |
| New crate or cross-crate change | Full | Design doc, implementation plan, staged commits |
| Any change to `popkern-core` public types | Full | Most other crates depend on it: grep for cross-crate impact |

**Escalation triggers** (bump rigor up one level):
- Change touches more than 3 files across 2+ crates
- Change modifies a public trait or error type
- Change affects serialization, tolerance handling, or Euler operators
- You're unsure which level applies

Every new public API should be designed with an eventual GUI in mind:
selectable inputs, a feasibility precheck, structured failure reasons.

**Design docs** (Heavy/Full) should cover: problem and user-visible goal;
relevant existing code; invariants and validity model; data model and API
shape; algorithms and fallback/repair paths; degenerate and adversarial
cases; error/diagnostic behavior; performance expectations; tests,
fixtures and acceptance criteria; migration or integration plan.
Remaining uncertainty is labeled as exploration with its own acceptance
criteria.

## Work Classification

Before substantial edits, state which of three a change is:

1. **Implementation path** — intended to be accepted code because it satisfies
   the commercial-grade standard.
2. **Diagnostic probe** — temporary exploration to learn facts. Must be gated
   or isolated, measured, and have its findings promoted into structured
   diagnostics, classifications, or repair facts. Removed unless promoted
   through the implementation standard.
3. **Temporary containment** — bounded mitigation that preserves correctness or
   fallback behavior while a stronger solution is pending. Must be recorded
   with a follow-up issue naming the commercial-grade replacement.

Do not pursue shortcuts framed as implementation.

## Hard Kernel Work

For hard kernel gates, work in stages. Record the baseline, current
metric, target gate, and next invariant. Make one layer true by
construction before advancing: topology ownership, G0 identity, G1 field
agreement, G2 rows, then determinism if results vary. Instrument before
broadening construction logic. A failed iteration only counts as progress
when it records the numeric delta, worst sample or handle, rejected
hypothesis, and next mechanism. When local repairs plateau, stop tuning
and revise the design before building the next structural substrate.
