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
crates/popkern-io          STEP import/export (depends: core)
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
- **Non-manifold boolean results decline** — operands touching only
  along an edge, a line or at a point make a result with zero-thickness
  contacts. `boolean()` returns `BoolError::NonManifoldResult` with the
  contact locations instead of serving coincident duplicate edges or
  vertices, or separate shells that touch. No single body holds lumps
  that touch without sharing, and nothing is recorded (PK-1115):
  `NonManifoldPolicy::SplitToManifold` serves separate manifold bodies,
  each valid alone, from the `_bodies` entry points (`boolean_bodies`,
  `nary_boolean_bodies`, `nary_boolean_where_bodies`,
  `Partition::select_bodies`, `symmetric_difference_bodies`; Part's
  `boolean_bodies_with`), the single-body ones declining with
  `BoolError::SplitIntoBodies` (a solid touching itself has no separate
  bodies and declines; in the feature tree each piece is its own body,
  see N-body operations). `NonManifoldPolicy::General` serves
  one general body: the split pieces glued where they touch, found
  from their geometry alone (`BRep::glue_touching`, PK-1117: every
  piece given a vertex at each contact point, its own, its edge split
  there or a vertex loop of its face; each contact curve, an edge or
  tangent ruling lying across a face, imprinted there by
  `BRep::imprint_edge`; then `glue_coincident`: one vertex, one ring in
  radial order; a glue that declines is
  `BoolError::GeneralBodyDeclined`, a `GeneralDeclined` shadow line),
  and general operands enter as their manifold pieces
  (`BRep::separate_at_contacts`, an intermediate; `split_to_manifold`
  returns them as separate bodies), in two-body and N-ary booleans
  alike (Rhai `non_manifold: "general"`). An old document's stored
  contact list loads and is ignored. The cell complex
  (`popkern-bool/src/complex/`) decides membership by region and is
  the primary classifier, falling back to the cellular classifier on a
  decline (a `PrimaryFallback` shadow line). `POPKERN_COMPLEX_PRIMARY=0`
  makes the cellular classifier primary with the complex as its rescue,
  `POPKERN_COMPLEX=0` disables the complex altogether,
  `POPKERN_COMPLEX_SHADOW=1` (with `POPKERN_COMPLEX_SHADOW_LOG=<file>`,
  under `POPKERN_COMPLEX_PRIMARY=0`) records its agreement with the
  cellular classifier on every boolean,
  and `POPKERN_COMPLEX_ARTIFACT=<dir>` dumps a declined complex as JSON.
- **Vertex table** — `popkern-bool/src/tolerance_zone/` resolves every
  intersection's candidate vertices together and decides the cellular
  graph's record-end vertices and paves. `POPKERN_VERTEX_TABLE=0`
  restores the graph's snap-then-adopt identity; a declined resolution
  falls back to it for that boolean (a `VertexTableFallback` shadow
  line, and its candidates under `POPKERN_COMPLEX_ARTIFACT`).
- **N-body operations** — `popkern_bool::partition` builds one cell
  complex over N solids (`popkern-bool/src/partition/`); `select` /
  `nary_boolean` keep the regions a `SetExpr` accepts, `region_solids`
  and `split` keep each region as its own solid with the faces
  neighbours share (`Interface`), `imprint` classifies a target's faces
  (`FaceClass`). Operations declare accepted `BodyKind`s
  (`BoolError::UnsupportedBodyKind`). Coplanar faces between bodies
  are held (operand-qualified `FacePair`s); a region bounded by several
  shells (a body floating in another, an operand's cavity) keeps them
  together (`partition/region_nesting.rs`); sheets enter as laminae
  that cut regions without bounding volume (PK-978); record crossings of
  different pairs adopt a node within 20 tolerances
  (`partition/crossings.rs`), so crowded triple points share a vertex.
  When a stage still declines, `nary_boolean` falls back to
  two-operand booleans, logged as a `PartitionFallback` shadow line.
  In the feature tree they are `Feature::Split { tools, keep }`,
  `NaryBoolean { tools, expr, non_manifold }`, `Imprint`, and
  `KeepPieces` / `SeparatePieces { source, rule }` (whole solids taken
  out by `popkern_bool::pieces`, names kept; a body made by
  `SeparatePieces` is chosen over the body its names point to). A
  split serves each kept piece as its own body, and a split-policy
  boolean each touching piece (PK-1120: no body holds lumps that
  touch): the node records its bodies (`FeatureNode::output_bodies`,
  its own body first), each keyed by a face name only its piece
  carried, names with lineage before positional generated ones; a
  replay serves the piece carrying the key as that body, a body whose
  key is gone takes the next unclaimed piece in order (a `Naming`
  fallback record), a gone piece leaves its body absent and a new
  piece gets a new body (`replay_outputs`); `FeatureTree::served_bodies`,
  script `Body::split` (returns the bodies) and `Body::served_bodies`,
  Rhai `b.split(..)` (an array) and `b.served_bodies()`; an
  imprint's face classes are recorded by stable name in its report row
  (`FeatureTree::imprint_class`, Rhai `m.imprint_class(name)`). Named by
  `popkern_bool::name_from_operands` (the target's names carry, the
  rest by the operand entities they are pieces of). The shortcut oracle
  `POPKERN_SHORTCUT_ORACLE=1` runs the general path on every touching /
  containment shortcut and logs agreement (`ShortcutOracle`).
- **Body kinds are a contract** — `BRep::body_kind` reads the kind off
  regions, rings and vertex fans (`Empty`, `Acorn`, `Wire`, `Sheet`,
  `Solid`, `General`; a solid opened mid-edit is `General` until
  `settle_open_shells`). Operations gate at entry with
  `BRep::require_kind(operation, accepted)` (`CoreError::UnsupportedBodyKind`,
  carried by each crate's `Core`/`Topology` error) or, in booleans and
  the operations built on them, `popkern_bool::require_operand_kind`
  (`BoolError::UnsupportedBodyKind` with the operand index): chamfers,
  fillet edits, drafts, sweeps into a body, holes, hinges and booleans
  take solids; edge fillets (open caps at a sheet's boundary), corner
  blends, face blends, sheet booleans, partition and imprint solids or
  sheets; shell, offset and plane splits solids (or empty, answered
  with `NoFaces` or no pieces), thicken and the maximum safe offset
  solids or sheets (or empty);
  centroid and inertia solids (or empty); volume solids, general
  bodies (the sum of their solid regions, `region_volumes`) or empty;
  tessellation and area any kind.
- **Non-manifold Euler operators** — `euler/nonmanifold/`: acorns
  (`make_acorn`/`kill_acorn`), wire edges (`make_wire_edge`, joining
  the shells its ends lie in; `kill_wire_edge` splits a shell it
  disconnects and leaves a lone end an acorn), a lamina spanning a wire
  cycle (`make_face_from_wire_cycle`, geometry the caller's), a region
  behind a closed surface (`make_face_region(shell, side, kind)`), a vertex on a face's interior
  held as a vertex loop (`make_vertex_loop`/`kill_vertex_loop`: a loop
  of no fins naming its vertex, none of the face's boundary loops,
  `BRep::face_vertex_loops`, a `VertexUse::Loop`; STEP `VERTEX_LOOP`
  bounds; PK-1116), an edge imprinted across a face between two of its
  vertices (`imprint_edge`/`kill_imprinted_edge`: the face split in two
  between points of its outer loop, a hole joined to it, a spur into a
  vertex loop, a slit between two; PK-1117), and
  the glues that make one general body of touching pieces
  (`glue_vertices`, `glue_edges` with the ring in radial order and
  shells of one region joined across new wedges, `glue_faces` merging a
  facing pair, its edges and the gap between in one step), each with
  its inverse. Each checks every precondition before it changes the
  body (`EulerError`, carried by `CoreError::Euler`); counts move by the
  documented delta over `NonManifoldCounts`. Radial order is core's
  (`topology/radial_order/`, generic over `RadialFins`;
  `BRep::radial_order`, `order_ring`); the cell complex adapts its edges
  to it, and `validate_all` checks a ring of more than two fins is in
  radial order (`RingOrder`, `RingOrderUndecided`).
- **General bodies downstream** — `popkern_mesh::tessellate_general`
  meshes any kind: faces tagged with the regions behind and in front
  (`FaceRegions`: `two_sided` for sheets and laminae, `internal` between
  two solids), wire edges as polylines, acorns as points,
  `region_triangles` a solid region's boundary faced out of it;
  `popkern_measure::region_volumes` (and `body_volume` of a general
  body, their sum) sums each face's divergence term into the region
  behind and out of the region in front. STEP: `StepImportOptions::
  non_manifold = StepNonManifold::General` imports a self-touching
  `CLOSED_SHELL` as one general body (its pieces glued), and a general
  body exports as items in one representation (`export_general`):
  touching solids as OCC writes them (one `CLOSED_SHELL`, shared edges
  reused, a face between two solids listed once; cavities, touching ones
  one shell, in a `BREP_WITH_VOIDS`), sheets and laminae in a
  `SHELL_BASED_SURFACE_MODEL` (a closed sheet around a void a
  `CLOSED_SHELL` there), wires and acorns as a `GEOMETRIC_CURVE_SET`
  (trimmed curves, points), which import as wire bodies; under the
  general option the items touching or held in a solid merge back into
  one general body (`import_general::merge_touching_items`: glued where
  they coincide, a sheet glued on a solid's material side adopted as its
  lamina by `glue_edges`, free sheets, wires and acorns placed in the
  cell holding them by `fill_enclosed_cells`; PK-1114). Sewing with `allow_nonmanifold_edges` glues
  each pair into a ring in radial order (a third sheet onto a sewn edge,
  `GlueDeclined` when it cannot). `BRep::contacts()` lists where a body
  touches itself, from its topology (rings of more than two fins,
  vertices of several uses, a contact edge's ends).
- **Sheet booleans** — `boolean_with_sheets` decides sheet ∩ solid,
  sheet − solid, solid − sheet (pieces) and the union's appendage on
  the cell complex (`sheet/complex_route.rs`, sheets as laminae:
  `partition/lamina_part.rs`), falling back to the planar split or
  imprint-and-trim path with a `PartitionFallback sheet` shadow line;
  `boolean()` still declines sheets. A sheet's free edge is canonical:
  a lone fin (its radial ring links null, `BRep::is_free_edge`), the
  survivor the edge's canonical half-edge with the curve run its way.
  The lamina view the complex reads mirrors each free edge with a
  twin on the same face, which on a periodic face reads as a seam: the
  view's open edges and lamina operands are passed explicitly
  (`pipeline_view`, `ComplexInput::laminae`), never re-read off it.
  The constructors' two-face sheets are converted at entry. OCC
  oracle: `tools/occ-oracle/occ_sheet.py`.
- **Wire booleans and contacts** — `boolean_with_wire` cuts a wire
  body where its curves cross the solid's faces and keeps pieces by
  midpoint (∩ keeps inside and on-face pieces, − the outside, solid −
  wire is the solid, ∪ a mixed body); on-face pieces are tested by
  face distance first, because `point_in_solid_final` calls a point on
  a face outside. A partition `Interface` carries `patch`, the shared
  surface as a one-face sheet: region extraction merges co-surfacial
  faces, so `faces` may reach past the patch, and contact or bonded
  areas are measured on `patch`. Assembly interference reports a
  touching pair's extent and area when asked (`contact_area`, off by
  default: a partition per touching pair).
- **STEP non-manifold bodies** — a `CLOSED_SHELL` that touches itself
  (a contact edge with four face uses, or a vertex two closed sets of
  faces share) imports as separate manifold bodies, one per piece, each
  with its own copies of the contact edges and vertices
  (`popkern-io/src/step/mapping/topology/shell_pieces.rs`,
  `import_general::separate_touching_pieces`; healing welds vertices
  within a piece only).
  The pieces are the face sets joined through two-use edges; a piece
  using an edge four times and odd uses (an internal face) decline
  under the default option, and under `StepNonManifold::General` import
  face by face, glued into one general body whose enclosed cells
  become its solids (`BRep::fill_enclosed_cells`: a closed piece of
  positive signed volume a solid when it holds a back side, else a
  bounded void; a negative piece an inner shell of the smallest cell
  containing it by winding number, `BRep::winding_about`; PK-1070,
  PK-1071). Pieces touching without sharing entities (cavities touching
  along an edge, PK-1072) glue where they coincide
  (`BRep::glue_coincident`) under the same option. Fixtures
  are OCC-written (`tools/occ-oracle/occ_write_nonmanifold.py`).
- **Tangency and contacts (PK-1123..1130)** — where two surfaces
  touch without crossing, `contact_classify` reads the relative
  curvature `D = II₁ − II₂` in one tangent frame (a touch, a crossing
  with its branch directions, a tangent curve, or `Undecided`, an
  eigenvalue under `2ε/L²` decided by the gap on a ring), and
  `intersect_surfaces_typed` reports `contacts` beside the unchanged
  curves (`contact_discovery`: critical points of the signed gap). An
  arm meeting a touch emits no curve, as every arm does: contacts are
  typed beside. A boolean whose operands meet in no curve classifies
  only by probes that agree off the other body's boundary and declines
  `BoolError::RelationUndecided` otherwise. Records meeting inside one
  another on a face are split there (`tolerance_zone/meetings.rs`),
  the one rule for records of one face pair passing within the
  tolerance (the vertex table's old branch rule is retired). A key
  runs its span's way by its ends when open and its leaving tangent
  when closed
  (`assemble_span_unify::runs_same_way`: a closed span's ends are one
  node). Touches of an edge with a face land on the root of `C′ · N`,
  not a minimum of distance values (`surface_curve_touch.rs`). Exact
  arms: sphere × sphere, sphere × coaxial cone, cones sharing an apex,
  and a NURBS surface of revolution (either parameter the angle,
  clamped or periodic form) × a coaxial plane, cylinder, cone or
  sphere (`revolved_nurbs.rs`). Manifoldness has one authority,
  `BRep::vertex_neighbourhood(s)` / `edge_neighbourhood` (a face
  meeting a vertex twice round one fan is a pinch use, still a
  manifold; the boolean's pinched-vertex test is its client;
  `imprint_edge` from a pinch vertex splices into the wedge the curve
  enters, read off chords since a tangency pinch's wedges are cusps,
  and declines when no single wedge holds it); fillets,
  chamfers at a pinch vertex, shells, offsets and thickens decline
  `CoreError::PinchUse` before building. Separate shells touching are
  found by alternating projection between face surfaces, as critical
  points of the signed gap where projection stalls at a tangent touch
  (`find_surface_contacts` over each face's window, a face closed round
  a seam spanning its period whole), and by each edge walked against
  the other face (`postprocess/shell_contact*.rs`: an edge resting on a
  face is a contact). Naming scheme 6 (PK-1160): an edgeless contact (a
  vertex loop) is `RoleCategory::Contact`, named by the faces it lies
  on; a general two-operand boolean is named after its glue; a face
  lying on several operand faces is confirmed on each at a second
  point (one point cannot tell a coincidence from a touch). A shell whose result
  touches itself (a wall closed exactly) declines
  `ShellError::SelfContact` (`popkern-shell/src/self_contact.rs`,
  tangent contacts of faces sharing no vertex); under
  `ShellOptions::with_non_manifold(General)` it serves one general
  body glued along the contact (`glue_touching` also glues a piece
  touching itself). A draft whose drafted faces meet across a strip
  face serves their ridge (`popkern-draft/src/ridge.rs`: the strip's
  short edges contracted, the strip merged into a side, the draft run
  again; `DraftWarning::FaceConsumedToRidge`); another consumed face
  declines `DraftError::FaceConsumed`. The linear tolerance is the
  resolution (PK-1131): within it is inclusive with a rounding
  allowance (`Tolerance::within_linear`, `linear_reach`), everywhere
  two things are decided one place. An operand or a result thinner
  than it everywhere has no inside (the other operand or empty), a
  union closes a gap at or under it between facing faces (an N-ary
  selection too, wherever it keeps both sides), and a merge flattens a
  step at or under it between faces facing one way; each is recorded
  as a `popkern_bool::sub_tolerance::SubToleranceFeature` (history and
  Part report it as a warning). Vertex identity reads the same
  tolerance (PK-1165): an operand vertex or an edge-face hit has half
  the linear reach as its radius (two are one place within the
  tolerance), a record end its faces' tolerances plus its fit
  deviation, a candidate joins the nearest compatible vertex, and the
  graph's snap and the assembler's coincident-span band are the reach
  (a record's band plus its fit deviation). Tolerance tests on edges
  are lengths, never parameter fractions, and a trim test's boundary
  band (ten tolerances) is not "on the face": contacts confirm it in
  space (`point_near_face_boundary`).
  `boolean_with_diagnostics` reports where the operands touch
  (`BoolDiagnostics::contacts`, each `OperandContact` classified);
  a declined vertex resolution's candidates, record meetings among
  them, dump under `POPKERN_COMPLEX_ARTIFACT`. Corpora:
  `popkern-examples/tests/tangency_corpus.rs` (placements, scales,
  self-turns, per policy, closed forms) and
  `popkern-rhai/tests/tangency_bench_scripts.rs`; both ratchet (a gate
  must pass, a target that starts passing must be promoted);
  `pinch_use_consumers.rs` puts a pinch use through every consumer.
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
