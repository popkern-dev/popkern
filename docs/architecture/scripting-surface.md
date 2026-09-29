# Scripting surface — normative specification

**Status:** accepted (PK-381, 2026-08-30) · **Owner:** scripted-modeling
campaign (`docs/architecture/scripted-modeling-design.md`, decisions
D1–D7) · **Skins bound by this document:** `popkern-script` (Rust
recorder, phase 2), `popkern-rhai` (phase 3), `popkern-py` (phase 5),
MCP tool schema (phase 6), and any future GUI command vocabulary.

This is the one language-agnostic contract behind every scripting
front-end (D7: *one surface, many skins*). A skin may adapt spelling
to its host language's conventions (method chaining, keyword
arguments, exceptions) but may not add operations, parameters, or
semantics that are not here, and every operation here has exactly one
canonical spelling that a migration script or a GUI can emit
(**the GUI-migration invariant**: every tree operation is readable,
writable code). Changes to this document go through the change log in
§13; a `preview` → `stable` promotion is a deliberate entry there,
never a side effect of shipping.

The kernel-side types this surface projects are cited inline as
`crate::path` so the spec can be diffed against the code; the
authority for wire forms is `crates/popkern-history` (`.pkn`,
format-compatibility ADR `docs/decisions/2026-08-30-pkn-format-compatibility-policy.md`).

---

## 1. Conventions

### 1.1 Stability marks

| Mark | Meaning |
|---|---|
| `stable` | Frozen for the launch train: name, parameter shape, and semantics change only through a change-log entry with a deprecation period. The launch `stable` set is deliberately small: profiles, extrude/revolve, booleans, fillet/chamfer, shell, parameters, labels, save/open/export, replay + gates. |
| `preview` | Shipped and usable, subject to change; docs and the playground label it as such. |
| `reserved` | Named here so skins do not squat the spelling; not implemented in this phase. |

### 1.2 Naming

- Operations are `snake_case` verbs or nouns: `extrude`, `fillet`,
  `linear_pattern`, `to_step`. Skins may camelCase for host
  conventions (`toStep`) but must document the mapping.
- Feature operation identifiers (the `operation` string on report
  rows and the auto-label prefix) are exactly
  `Feature::operation_label()`: `sketch`, `raw_body`, `extrude`,
  `fillet`, `chamfer`, `boolean`, `shell`, `transform`, `thicken`,
  `offset_solid`, `draft`, `sweep`, `loft`, `pipe`, `helical_sweep`,
  `revolve`, `reference`, `linear_pattern`, `circular_pattern`,
  `curve_pattern`, `mirror`.
- Auto-labels are `{operation}-{n}` with `n` document-unique and never
  reused (`popkern_history::labels`). A user-supplied label that is
  already taken is replaced by an auto-label at record time
  (recording never fails on a name); `rename` is the strict path.

### 1.3 Scalars, parameters, expressions

Every numeric argument of a feature op is a **scalar binding**
(`popkern_history::feature_sketch::ScalarBinding` generalised):

| Form | Spelling | Records as |
|---|---|---|
| literal | `0.08` | a feature-scoped parameter `{feature}.{role}` holding the literal (so it is editable later without changing the script) |
| named parameter | `p("wall")` / a `Param` handle | `ParameterId` reference |
| expression | `expr("2*wall + clearance")` | a feature-scoped parameter with `Expression` (§4.3) |

Units: parameters carry a `Unit` (`Meters`, `Millimeters`,
`Centimeters`, `Inches`, `Radians`, `Degrees`, dimensionless). Length
arguments default to the document's length unit; angle arguments to
degrees at the surface (stored in radians). Dimension mismatch is a
record-time error (`HistoryError::DimensionMismatch`).

### 1.4 Geometry literals

- Points and vectors: `[x, y, z]` arrays in document units.
- Planes: `Plane::XY | XZ | YZ`, `Plane::xy_at(z)` (and `xz_at`,
  `yz_at`), `Plane::on_face(sel)`, `Plane::on_face_offset(sel, d)`,
  `Plane::custom(origin, normal)` — projecting
  `feature::SketchPlaneDesc`, which is the authority. (The
  `popkern-part` facade currently exposes only `XY`/`xy_at`; phase 2b
  compiles the full set through `popkern-sketch`.)
- Directions: `Dir::X | Y | Z`, a literal `[x, y, z]`,
  `face_normal(sel)`, `edge_direction(sel)`, and `-dir` (negation) —
  projecting `feature::Vec3Desc`.
- Axes: `axis(edge_sel)` or `axis_at(origin, direction)` —
  `feature::AxisDesc`.

### 1.5 Tolerance

The document carries one build tolerance
(`FeatureTree::build_tolerance`; `effective_build_tolerance()` is the
kernel default when unset). `model(...)` accepts `tolerance:` once;
no op takes a tolerance argument. Skins never expose global state.

### 1.6 Transactionality

Every recording call either appends exactly one feature (and, for
`Boss`/`Cut` modes, the auto-merge boolean it implies) and replays it
clean, or appends nothing and raises the structured error (§10).
Migration verbs are transactional per call in the same sense. This is
boolean invariant 3 extended to the whole surface.

---

## 2. Model / session

| Op | Mark | Signature (canonical) | Semantics |
|---|---|---|---|
| `model` | stable | `model(name?, tolerance?, unit?) -> Model` | Begin a new document with an empty timeline. |
| `open` | stable | `open(path) -> Model` | Load a `.pkn` (versions `MIN_SUPPORTED..=CURRENT`; else `UnsupportedVersion { file_version, min_supported, max_supported }`). Labels are made unique on load. No replay is performed until requested. |
| `save` | stable | `model.save(path)` | Write `.pkn` at `CURRENT_VERSION` with the naming-scheme stamp. Never bumps the version for additive fields (compat ADR). |
| `replay` | stable | `model.replay() -> ReplayReport` | Full replay under the document tolerance. Returns the report (§11); raises `ReplayError` (§10.3) on the first failed feature (the tree's replay stops there today). |
| `replay_lenient` | preview | `model.replay_lenient() -> ReplayReport` | Failed features are reported (`outcome: failed`) and replay continues with the body up to the failure; requires the tree to gain continue-on-failure replay. |
| `to_brep` | stable | `model.to_brep() -> BRep` | Ephemeral exit: replay (incremental, from the cache) and return the current body of the active `BodyId`. |
| `to_part` | stable | `model.to_part() -> Part` | As above wrapped in `popkern_part::Part` (validated at the serving boundary). |
| `to_step` / `to_stl` / `to_obj` | stable / preview / preview | `model.to_step(path)` … | Export the current body. |
| `assert_watertight` | stable | `model.assert_watertight()` | Gate: 0 open mesh edges; raises `NotWatertight { open_edges }`. |
| `assert_step_roundtrip` | stable | `model.assert_step_roundtrip(path)` | Gate: export, re-import, face count identical. |
| `bodies` | preview | `model.bodies() -> [BodyId]` | Bodies alive at the end of the timeline. |
| `body` | preview | `model.body(id_or_label) -> BodyId` | Resolve a body by id or label (§9). |
| `active_body` / `set_active_body` | preview | | The body new features attach to; default = the body of the last feature. |

`Model` is the recorder session: it owns the `FeatureTree`, the
`BRepCache`, and the label table, and every feature op below is a
method on it (or on a `Body` view of it).

---

## 3. Sketching

Two dialects (D2), never silently converted.

### 3.1 Analytic profiles (`profile`) — stable

```
model.profile(plane)            -> ProfileBuilder
  .rect(w, h) | .circle(r) | .polygon([[x,y], ...])
  | .path([x0,y0]).line_to([x,y]).arc_to(via=[x,y], to=[x,y]).close()
  .at(cx, cy)                    // placement on the plane, default [0,0]
  .sketch(label?)                -> Feature   // records Feature::Sketch { Analytic }
```

Records `Feature::Sketch { plane, sketch_data: Analytic(AnalyticProfile<ScalarBinding>) }`;
every scalar is a scalar binding (§1.3), so a profile's dimensions are
parameters. Compilation is `popkern_sketch::profile_compile`, shared
with `popkern-part` — the recorder and the direct API produce
identical faces (parity fixtures, design §13).

### 3.2 Constraint sketches (`sketch`) — preview

```
model.sketch(plane)             -> SketchBuilder
  .point(name, [x,y]) .line(name, a, b) .circle(name, c, r) .arc(name, c, a, b)
  .coincident(a, b) .horizontal(l) .vertical(l) .parallel(a, b) .perpendicular(a, b)
  .tangent(a, b) .equal(a, b) .fixed(p)
  .distance(a, b, binding) .angle(a, b, binding) .radius(c, binding)   // driving dims
  .sketch(label?)                -> Feature   // records Feature::Sketch { Constraint(json) }
```

Driving dimensions bind to parameters exactly like feature scalars.
Payload is `popkern_constraint::Sketch` JSON, re-solved on replay.

### 3.3 Promotion — reserved

`model.feature(sketch).promote_to_constraint()` (design §5.4): one-way
analytic → constraint; synthesises a fully constrained sketch (0 DOF)
and records the representation change in the timeline. Reserved
until §11.E3 (constraint arc exactness) resolves.

---

## 4. Parameters and configurations

### 4.1 Parameters — stable

| Op | Signature | Records |
|---|---|---|
| `param` | `model.param(name, value, unit?, bounds?) -> Param` | `ParameterTable::add_global` (global scope; display-unit variant when `unit` given) |
| `param_expr` | `model.param_expr(name, "text", unit?) -> Param` | global parameter with `Expression` |
| `set_param` | `model.set_param(name_or_param, value)` | `set_value` (bounds clamp reported as `EvaluationWarning::BoundsClamped`) |
| `set_expr` / `clear_expr` | | `set_expression` / `clear_expression` |
| `set_bounds` / `clear_bounds` | | |
| `rename_param` | `model.rename_param(old, new)` | `ParameterTable::rename` (qualified names re-keyed) |
| `params` | `model.params() -> [ParamInfo]` | name, qualified name, scope, value, unit, expression text, bounds, read-only, last evaluation warning |

Feature-scoped parameters are named `{feature_qualified_prefix}.{role}`
(`fillet-3.radius`, `extrude-1.distance`); `qualified_prefix` is the
one normalisation rule (lower-case, spaces → `_`). Skins expose them
through `model.feature(label).param(role)`.

### 4.2 Expression language

Text form parsed to `popkern_history::expression::Expression`:
literals, parameter references by qualified name, `+ - * / ^`,
unary functions (`sin cos tan sqrt abs floor ceil min max …`),
comparisons and `if(cond, a, b)`, and measurement terms
`measure(<query>)` (§9). Evaluation warnings (`MathError`,
`BoundsClamped`) ride the replay result, never abort replay; cycles
and unknown references are record-time errors.

### 4.3 Configurations — preview

Mirror of `config_api`: `add_config(name, description)`,
`add_derived_config(name, parent)`, `duplicate_config`,
`delete_config`, `rename_config`, `set_active_config(name|None)`,
`set_override(config, param, value)`, `clear_override`,
`set_suppression(config, feature, bool)`, `clear_suppression`,
`configs()`, `config(name)`. Replay runs under the active
configuration.

---

## 5. Feature operations

Each row: canonical spelling; every argument is either a selection
(§6), a scalar binding (§1.3), a geometry literal (§1.4), or an enum
spelled exactly as listed. `mode` is `new_body | boss | cut`
(`feature::OperationMode`; default `new_body`).

**Receivers and return values.** Feature ops are methods on a `Body`
view (`model.active()` or the `Body` returned by an earlier op); the
few that need no body (`profile`, `sketch`, `raw_body`, `reference`,
and `loft`, whose inputs are section sketch bodies) are on `Model`. Every op returns the `Body` it served so calls chain
in timeline order; the recorded feature is `body.last()` / the
`Feature` handle (`FeatureId` + label) for §7, and `model.feature(label)`
always works. `new_body`-mode ops start a new body and return it;
`boss`/`cut` return the receiving body.

| Op | Mark | Canonical signature | Records | Notes |
|---|---|---|---|---|
| `extrude` | stable | `extrude(face, end, direction?, mode?, draft_angle?, second?)` | `Feature::Extrude` | `end`: `blind(d)`, `through_all()`, `up_to_surface(face)`, `offset_from_surface(face, d)`, `up_to_vertex(v)`, `mid_plane(half_d)`. `second`: `{direction, end, draft_angle?}` (`ExtrudeSecondDir`). `face` defaults to the last sketch's face. |
| `revolve` | stable | `revolve(face, axis, end, mode?)` | `Feature::Revolve` | `end`: `blind(angle)`, `full_circle()`, `mid_plane(half_angle)`, `up_to_surface(face)`. |
| `sweep` | stable | `sweep(face, path_body, path_edge, twist?, mode?)` | `Feature::Sweep` | `path_edge` is a single edge of `path_body` (the profile sketch body has no path edges; the path body is referenced, not consumed). The profile contributes its shape; placement comes from the path. `twist` in degrees, default 0. |
| `loft` | stable | `m.loft(sections, closed?, mode?)` | `Feature::Loft` | `sections` is an ordered list of sketch bodies, each contributing its most recent sketch face; the loft records on the first section's body and consumes the others. Mixed open/closed sections are rejected. |
| `pipe` | stable | `pipe(path_edge, radius, mode?)` | `Feature::Pipe` | `path_edge` is a single edge of the receiving body. |
| `helical_sweep` | stable | `helical_sweep(face, axis, radius, pitch, turns, taper_to?, start_angle?, left_hand?, mode?)` | `Feature::HelicalSweep` | Sweeps the profile face along a helix about `axis` (springs, coils, auger flights, thread forms). `radius`/`pitch` are lengths; `turns` is a count (fractional allowed); `taper_to` is the end radius of a conical helix; `start_angle` (degrees, default 0) phases the helix about the axis — what places multi-start threads and multiple blades apart; `left_hand` (default false) mirrors the winding. The profile contributes its shape; placement comes from the helix. |
| `bolt` | preview | `m.bolt(spec, length, head?, thread_length?, left_hand?, starts?)` | `Feature::Fastener` | Preview until the PK-460 walk-class campaign closes (heads/nuts/holes exercise merge classes still being rooted; `threaded_rod` is the stable representative). A hex / socket-cap / countersunk screw with a real modeled thread, built in the insertion frame (axis +Z, origin at the bearing plane, threads toward +Z) and placed with the transform verbs. `spec`: a designation string (`"M6"`, `"M6x0.75"`, `"1/4-20"`, `"#8-32"` — millimetre dimensions regardless of the model unit) or `#{ major, pitch }` in model units; `head`: `"hex"` (default) `\| "socket" \| "countersunk"`; `thread_length` limits the threaded stretch (default full). |
| `threaded_rod` | stable | `m.threaded_rod(spec, length, left_hand?, starts?)` | `Feature::Fastener` | Full-length external thread, both ends free-running, tip lead-in chamfer. |
| `nut` | preview | `m.nut(spec, clearance?, thickness?, left_hand?, starts?)` | `Feature::Fastener` | Hex nut with a modeled internal thread; `clearance` is the printed-fit grow of the internal form (0.15–0.3 typical for FDM); `thickness` overrides the ISO 4032 nominal. |
| `tapped_hole` | preview | `tapped_hole(axis, spec, depth, clearance?, left_hand?, starts?, entry_chamfer?)` | `Feature::TappedHole` | Bore at the internal minor diameter plus a real modeled internal thread, entered at the axis origin, cutting along the axis direction (INTO the material). `depth`: `blind(d)` or `through_all()`. A fastener built on the same axis mates by construction. |
| `clearance_hole` | preview | `clearance_hole(axis, spec, depth, fit?, diameter?, countersink?, counterbore?, entry_chamfer?)` | `Feature::ClearanceHole` | Plain hole for a screw of `spec`: `fit`: `"close" \| "normal"` (default) `\| "loose"` (ISO 273) or `diameter` explicit; `countersink`/`counterbore` sink the matching head standard plus the given extra room. |
| `union` / `subtract` / `intersect` | stable | `body.union(tool_body)` … | `Feature::Boolean { operation, tool_body }` | `tool_body` is a `BodyId`/label of the same tree (§9); result replaces the receiving body. |
| `fillet` | stable | `fillet(edges, radius, continuity?, cap?, profile?, interference?, law?, corner?, corner_cap?, corner_cap_strict?)` | `Feature::Fillet { edges, radius, recipe }` | `continuity`: `g0 | g1 | g2`; `cap`: `open | flat | rounded | tangent_match`; `profile`: `circular | conic(rho) | curvature_continuous`; `interference`: `decline | serve_non_interfering`; `law`: `constant | linear(start, end) | by_stations([(t, r)…])`; `corner`: `smooth | mitre`. Preview (PK-603): `corner_cap`: `auto | setback | fitted | classic` — which construction closes a corner where three or more rounds meet (`recipe.corner_cap`); `corner_cap_strict`: `true` serves that construction or fails the feature, never falling back (needs `setback` or `fitted`). `corner_cap` with `corner: "mitre"` is refused. The report's `blend.corners` names what served each corner. |
| `fillet_each` | preview | `fillet_each(groups, continuity?, cap?, profile?, interference?, corner?, corner_cap?, corner_cap_strict?)` | `Feature::Fillet { edges, radius, recipe }` with `recipe.edge_radii` | One fillet feature with several radius groups (PK-602): `groups` is `[[edge_query, radius], …]`; each group's edges take its radius and the corner where they meet is one vertex blend. Groups never share an edge (refused, naming each overlapping pair); a group that selects nothing is dropped when its query allows empty; radii are constant (no `law`). Group 0's radius is the feature's `radius`, group `k`'s is the parameter `radius_k`. |
| `chamfer` | stable | `chamfer(edges, distance, kind?)` | `Feature::Chamfer` | `kind`: `equal` (default), `two_distance(face, d2)`, `distance_angle(face, angle)`. |
| `shell` | stable | `shell(faces_to_remove, thickness, direction?, join?, trim?)` | `Feature::Shell` | `direction`: `inward` (default) `| outward | midplane`; `join`: `arc | intersection`; `trim`: `trim | decline`. Empty `faces_to_remove` = closed hollow. |
| `thicken` | preview | `thicken(thickness, join?, trim?)` | `Feature::Thicken` | symmetric ±t/2 (PK-323). |
| `offset_solid` | preview | `offset_solid(distance, join?, trim?)` | `Feature::OffsetSolid` | |
| `draft` | preview | `draft(faces, angle, neutral)` | `Feature::Draft` | `neutral`: `face(sel)` or `plane(origin, normal)`. |
| `translate` / `rotate` / `transform` | stable / stable / preview | `translate(dx, dy, dz)`, `rotate(axis, angle)`, `transform(matrix4x4)` | `Feature::Transform { motion }` | Rigid motions of the active body. |
| `linear_pattern` | preview | `linear_pattern(seeds, primary, secondary?, geometry_pattern?, overrides?)` | `Feature::LinearPattern` | `primary`/`secondary`: `{direction, count, total_distance}`. Seeds: contiguous features of one body starting with a sketch and containing a boss/cut (`HistoryError` seed rules). |
| `circular_pattern` | preview | `circular_pattern(seeds, axis, count, total_angle, geometry_pattern?, overrides?)` | `Feature::CircularPattern` | |
| `curve_pattern` | preview | `curve_pattern(seeds, curve, count, equal_spacing?, align_to_curve?, geometry_pattern?, overrides?)` | `Feature::CurvePattern` | |
| `mirror` | preview | `mirror(seeds, plane, merge_result?)` | `Feature::Mirror` | `plane`: `xy | xz | yz | on_face(sel) | custom(origin, normal)`. |
| `raw_body` | preview | `raw_body(brep, origin?)` | `Feature::RawBody { snapshot, origin }` | `origin`: `import(source) | polyhedron | smooth_mesh | engrave | other(summary)`. Document-as-truth: the snapshot is embedded; file-referencing is a non-goal. `import_step(path)`, `polyhedron(faces)`, `smooth_mesh(cage)`, `engrave(...)` are conveniences that produce a `raw_body`. |
| `reference` | preview | `reference(name, query)` | `Feature::Reference` | Publishes a measurement (§9) as a read-only parameter `name`. |

`instance_overrides` for patterns: `{index: {suppressed?, params: {qualified_name: value}}}`.

---

## 6. Selection

### 6.1 Grammar (data, not closures)

A selection is a **query** evaluated against the live body at
record time; the document stores the result, never the query as
code. The grammar is data so every skin (and the MCP schema) can
carry it verbatim:

```
edge_query  := edges() ( "." predicate )* ( ".or(" edge_query ")" )* ( ".allow_empty()" )?
face_query  := faces() ( "." predicate )* ( ".allow_empty()" )?

edge predicate :=
    circular() | straight() | conic() | nurbs()
  | between(SurfaceKind, SurfaceKind)
  | at_height(z) | at_vertex([x,y,z]) | direction(Dir | [x,y,z])
  | chord_length_between(min, max) | not_seam()
  | labeled(label) | named(entity_name)

face predicate :=
    kind(SurfaceKind) | facing(Dir) | facing_neg(Dir) | at_height(z)
  | labeled(label) | named(entity_name)

SurfaceKind := plane | cylinder | cone | sphere | torus | nurbs
```

Predicates compose by conjunction; `.or(...)` disjoins whole
conjunctions (`EdgeQuery::or`). The vocabulary is
`popkern_part::select::{EdgeQuery, FaceQuery}` today
(`circular`, `straight`, `between`, `at_height`, `direction`,
`chord_length_between`, `not_seam`; `kind`, `facing`, `facing_neg`,
`at_height`) — all `stable` — plus four predicates the recorder adds
over the same `EdgeInfo`/`FaceInfo` data and that `popkern-part` gains
in phase 2c: `conic()`, `nurbs()` (the remaining `CurveKind`s),
`labeled(label)`, `named(entity_name)` — `preview`; and `at_vertex([x,y,z])`
(PK-602: edges with an end vertex at the point, within 1e-6 — the edges
meeting at one corner) — `preview`. `filter(closure)`
exists in the Rust skin only and is **not** part of the surface (a
closure has no canonical spelling — the recorder records its resolved
names and `describe()` text like any other query, flagged `opaque`).

Single-entity arguments (`face` of `extrude`, `path_edge` of `sweep`,
`vertex` of `up_to_vertex`) take a query that must resolve to exactly
one entity (`SelectionCardinality { expected: 1, got }` otherwise).

### 6.2 Recording semantics (D3)

- Resolution happens once, at the call. Seam edges are filtered out
  and the **post-filter** membership is stored.
- Stored form per selection: `{ names: [EntityName], describe: text,
  labels_used: [label] }`. The `describe()` text is the query's
  canonical spelling and is diagnostics only; it is never re-evaluated
  at replay.
- Empty resolution is an error (`EmptySelection { selection, entity,
  available }`) unless `.allow_empty()`, in which case the feature is
  recorded as a no-op with a `SkipRecord`.
- Replay failures cite label and text:
  `'mounting_boss' (fillet-3, authored as: circular ∧ between(cylinder, plane)) — 2 of 9 edges no longer resolve`.

### 6.3 Labels

| Op | Mark | Semantics |
|---|---|---|
| `label(name, selection)` | stable | Bind an entity label to the resolved (single) entity; `DuplicateLabel` if taken. Entity and feature labels are separate namespaces. |
| `unlabel(name)` | stable | |
| `rename_label(old, new)` | stable | Edits the label table only; zero references change. |
| `labels()` | stable | `[(label, EntityName)]` |
| feature labels | stable | `feature.label`, `model.rename(feature, new_label)` (strict: `DuplicateLabel`, `EmptyLabel`). |

---

## 7. Timeline / migration verbs

All verbs operate on `Model` and are the canonical code form of the
GUI's history-tree edits.

| Verb | Mark | Signature | Records / effect |
|---|---|---|---|
| `feature` | stable | `model.feature(id_or_label) -> Feature` | `feature_by_label` / id lookup; `LabelNotFound`. |
| `features` | stable | `model.features() -> [FeatureInfo]` | timeline order: id, label, operation, suppressed, body, params. |
| `insert_before` / `insert_after` | preview | `model.insert_before(anchor, |m| m.fillet(...))` | `insert_feature_at` with dependency validation (`DependencyAfterInsertion`, `DependencyCycle`); the closure records exactly one feature at the rollback position. |
| `modify` | stable | `feature.modify(|f| f.radius(0.1).continuity(g2))` | `modify_feature` with the same op's parameter shape; unknown fields for the op are errors. Selection edits use `.edges(query)` and re-resolve at call time. |
| `reorder` | preview | `model.reorder(feature, index)` | `reorder` (validated). |
| `suppress` / `unsuppress` | stable | | `suppress` / `unsuppress`. |
| `delete` | preview | `model.delete(feature)` | removes the feature and its parameters; `ParameterReferenced` if others depend on its params. |
| `rollback_to` | preview | `model.rollback_to(feature)` / `rollback_to_end()` | `rollback_to` (feature index); subsequent recording inserts at the rollback position. |
| `describe` | stable | `feature.describe() -> text` | One plain-language line (`"fillet 0.08 on 9 edges: circular ∧ between(cylinder, plane)"`) — the GUI history-tree caption and the `describe()` text of §6.2; derived from the recorded feature, never stored. |
| `dependents` / `dependencies` | preview | `feature.dependents() -> [Feature]` | parameter and entity-provenance dependency queries. |
| `preview` | reserved | `model.preview(|m| m.fillet(...)) -> PreviewResult` | Record into a scratch position, replay, return body + report, discard (UI preview mode). Reserved until 2b lands rollback-position recording. |
| `undo` / `redo` | preview | | `TreeCommand` stack (rename, label, tolerance, timeline edits are all undoable). |

---

## 8. Interrogation

| Op | Mark | Semantics |
|---|---|---|
| `measure(query)` | preview | Evaluate a `MeasureQuery` against the current body: `edge_length(e)`, `edge_tangent(e, x|y|z)`, `face_area(f)`, `face_perimeter(f)`, `face_normal(f, c)`, `face_centroid(f, c)`, `face_radius(f)`, `face_min_curvature(f)`, `face_max_curvature(f)`, `vertex_position(v, c)`, `min_distance(a, b)`, `dihedral_angle(a, b)`, `edge_angle(a, b)`, `body_volume()`, `body_surface_area()`, `body_centroid(c)`, `body_bound(min_x…size_z)`, `body_moment_of_inertia(axis)`. |
| `max_fillet_radius(edges)` | preview | Slider bound: the largest radius the precheck admits for the selection (UI rule: every numeric argument has a queryable range). |
| `precheck_fillet(edges, radius)` | preview | `FilletFeasibility` (issues, `max_radius_estimate`, convexity, seam). |
| `param_range(feature, role)` | reserved | Generic slider bounds for any scalar role (`{min, max, default, continuous}`); today only fillet radius has a measured bound (`max_fillet_radius`), the rest fall back to declared parameter bounds. |
| `report()` | stable | The last `ReplayReport` (§11). |
| `feature.evidence()` | stable | That feature's `FeatureReport`. |

---

## 9. Bodies

- A `BodyId` is minted when a `new_body`-mode feature (or `raw_body`)
  serves; `boss`/`cut` modes fold the new body into the receiving
  body through the implied boolean (`auto_merge`).
- Tool models for explicit booleans are recorded **into the same
  tree** under their own `BodyId`: `let tool = model.tool(|m| { m.profile(...).circle(r).extrude(...) })`
  returns the body handle (`model.body(id_or_label)` is the lookup);
  `plate.subtract(tool)` records
  `Feature::Boolean { operation: Subtraction, tool_body }` and the tool
  body is consumed. Bodies can be labelled (`label_body(name, body)`,
  preview) and are stable across save/open (labels are document
  data).
- Assemblies/product structure are out of scope (design §1).

---

## 10. Error payload schema

One taxonomy at every skin (design §8). Each error is a structured
payload with the fields below; skins present it natively (Rust enum,
Python exception attributes, Rhai error map, MCP `tool-error` JSON)
and `message` is always the Rust `Display` string.

### 10.1 Envelope

```
{ kind: <class>, message, feature?: {id, label, operation}, selection?: {describe, labels_used}, ...payload }
```

### 10.2 Classes (record-time; projections of `popkern_part::PartError` and `popkern_history::HistoryError`)

| Class | Payload fields | Source |
|---|---|---|
| `BlendDeclined` | `selection`, `blend_kind` (BlendError variant), `remediation {action, parameter, value, rationale}`, plus per-variant: `fillet_a/fillet_b`; `mechanism/detail`; `parameter/value/reason`; `detail`; `failed/total/failures[]` | `PartError::Blend` |
| ↳ `failures[]` | `{edge, reason_kind, requested?, max_radius?, reason?, vertex?, spine_parameter?, location [x,y,z]?, remediation}` | `BlendFailure` + `remediation_hint()` |
| `EmptySelection` | `selection`, `entity` (`edges|faces`), `available` | `PartError::EmptySelection` |
| `SelectionCardinality` | `selection`, `expected`, `got` | surface |
| `BooleanFailed` | `operation`, `detail`, `stage?`, `remediation?` | `PartError::Boolean` (`BoolError`) |
| `SweepFailed` | `feature`, `detail` | `PartError::Sweep` |
| `ValidationFailed` | `feature`, `detail`, `issues[]` | `PartError::Validation` |
| `NotWatertight` | `open_edges` | `PartError::NotWatertight` |
| `IoFailed` | `path`, `detail` | `PartError::Io` |
| `NotApplicable` | `detail` | `PartError::NotApplicable` |
| `Cancelled` | — | `PartError::Cancelled` |
| `InvalidParameter` | `parameter` (**surface name**, e.g. `radius`, with `kernel_path`, e.g. `specs[0].radius`, kept beside it), `value`, `reason` | `BlendError::InvalidParameter` and equivalents — the skin translates kernel spec paths to surface names (PK-370 finding) |
| `LabelError` | `label`, `reason` (`duplicate|empty|not_found`) | `HistoryError::{DuplicateLabel, EmptyLabel, LabelNotFound}` |
| `ParameterError` | `parameter`, `reason` (`cycle|unknown_reference|referenced|read_only|duplicate|dimension_mismatch {expected, got}`) | `HistoryError` parameter variants |
| `TimelineError` | `reason` (`feature_not_found|invalid_position|dependency_cycle|dependency_after_insertion|body_not_found|seed_*`) | `HistoryError` timeline/pattern variants |
| `DocumentError` | `reason` (`unsupported_version {file_version, min, max}|unsupported_naming_scheme {file_scheme, supported}|serialization`) , `path` | `HistoryError` document variants |
| `ArgumentError` | `argument`, `reason` | skin-level (bad enum spelling, wrong type) — raised before the kernel runs |

### 10.3 Replay-time (`popkern_history::ReplayError`)

| Class | Payload |
|---|---|
| `FeatureFailed` | `feature {id, label, operation}`, `error` (a `FeatureFailure` — §10.2 class carried as data: `kind` = `BlendDeclined` (`blend_kind`, `message`, `remediation`, `failed/total/failures[]` with `edge`, `edge_name`, `reason_kind`, `remediation`) \| `BooleanFailed` (`operation`, `bool_kind`, `stage?`, `detail`) \| `SweepFailed`/`ShellFailed` (`feature`, `detail`) \| `DraftFailed` \| `ValidationFailed` \| `InvalidParameter` \| `Cancelled` \| `Other` (`detail`)), `report` (the partial `ReplayReport` up to the failure; the recorder's `Model::report()` after the call). Landed PK-388: `popkern_history::FeatureFailure`, carried by `ReplayError::FeatureFailed`, `FeatureResult::Failed`, `ReplayOutcome::Failed`, `ScriptError::FeatureFailed`; `Display` is the kernel's message. |
| `ResolutionFailed` | `feature`, `name` (EntityName), `selection` (describe/labels), `resolved`/`expected` counts |
| `CoreError` | `detail` |

### 10.4 Rules

- Every class above is a distinct type in every skin (Python: one
  exception class each, all subclasses of `PopkernError`; Rhai: the
  `kind` string; MCP: `error.kind`). Verified by the E2 probe
  (`probes/pyo3-diagnostics`).
- Fields are never flattened into the message; the message is
  redundant with the fields.
- Floats cross bit-exactly (serde_json `float_roundtrip` where JSON
  is the carrier).
- One projection: `popkern_script::ScriptError::payload()` is the
  envelope every skin presents (PK-403); a skin adapts the carrier
  (map, exception attributes, JSON), never the field set. Parity is
  gated per class in `crates/popkern-rhai/tests/error_parity.rs`.

---

## 11. Report schema

`ReplayReport` (`popkern_history::replay_report`) is the observability
spine. Surface form:

```
ReplayReport { features: [FeatureReport], is_clean, warnings[], fallbacks[], skips[], downgraded_blends[], failures[] }
FeatureReport { feature: {id, label, operation}, outcome: served|failed{error}|suppressed{reason}, elapsed_ms?,
                warnings: [{source, message}], blend?: BlendEvidenceRecord, fallbacks: [{source, mechanism, reason}],
                skipped: [{source, what, reason}] }
BlendEvidenceRecord { edges_requested, corner_patches, requested_continuity, continuity?: {measured…}, downgraded }
```

`ReplayResult` additionally carries `tolerance`, `elapsed`, and
`eval_warnings: [(parameter, warning)]`. A fillet always carries
measured continuity evidence; a downgrade is a first-class fact
(`downgraded == true`), never a green checkmark.

---

## 12. Canonical spellings (worked examples)

The boss-plate showcase (`crates/popkern-examples/src/flagship_boss_plate_api.rs`)
as recorded — this is the phase-2 acceptance fixture:

```
let m = model("boss_plate", unit: mm);
let plate = m.profile(Plane::XY).rect(2.4, 1.6).sketch("plate_sketch")
             .extrude(blind(0.4));                              // extrude-1, body A
let boss  = m.tool(|m| m.profile(Plane::xy_at(0.2)).circle(0.35).sketch()
             .extrude(blind(0.8)));                             // extrude-2, body B
plate.union(boss);                                              // boolean-1
plate.fillet(
    edges().between(cylinder, plane).at_height(0.4)
        .or(edges().straight().at_height(0.4))
        .or(edges().straight().direction(Dir::Z)),
    0.08, continuity: g1, label: "rounds");                     // fillet-1 → label "rounds"
m.assert_watertight(); m.to_step("boss_plate.step"); m.save("boss_plate.pkn");
```

Migration, as a GUI would emit it:

```
let m = open("boss_plate.pkn");
m.feature("rounds").modify(|f| f.radius(0.1).continuity(g2));
m.set_param("extrude-1.distance", 0.5);
let report = m.replay();
assert!(report.is_clean());
m.save("boss_plate.pkn");
```

Rule for skins: a script re-generated from a document (the GUI's
"show as code") must consist only of spellings from §§2–9, in
timeline order, one statement per feature, with every scalar as its
current parameter binding.

Rust skin note (`popkern-script`, phase 2b): `Body` is a copyable
handle and every op takes the `Model` explicitly
(`plate.fillet(&mut m, edges()…, 0.08)`), `profile(...).sketch(label)`
returns the `Body` the sketch belongs to, and `End::Blind(0.4.into())`
spells `blind(0.4)`. The canonical spellings above are the
host-language-neutral form; the Rust skin's mapping is one-to-one.

---

## 13. Change log

| Date | Change | Marks affected |
|---|---|---|
| 2026-09-29 | Corner-cap preference (PK-603): `fillet` / `fillet_each` take `corner_cap` (`auto | setback | fitted | classic`) and `corner_cap_strict`, recorded as `recipe.corner_cap` (additive, `serde(default)`, omitted on disk when default). A strict preference serves its construction or fails the feature (`BlendFailure::CornerCapUnavailable`). `corner_cap` with a mitre, or strict without `setback`/`fitted`, is refused. The blend report gains `corners` (construction, sides, achieved, patches, centroid per corner). | `corner_cap`, `corner_cap_strict` added preview |
| 2026-09-29 | Radius groups (PK-602): `fillet_each(groups, …)` records one `Feature::Fillet` whose `recipe.edge_radii` gives each group's edges its own radius (additive, `serde(default)`; single-radius trees are unchanged on disk), so a corner's edges can carry different radii and be built as one vertex blend. Overlapping groups are refused; group radii are parameters `radius_k`; `modify(edges)` on a grouped fillet is refused (re-record with `fillet_each`); emission spells `fillet_each`. New edge predicate `at_vertex([x,y,z])` (`EdgeInfo` gains `endpoints`). | `fillet_each` added preview; `at_vertex` added preview |
| 2026-09-04 | Sweep family landed (PK-440): `sweep`, `loft`, `pipe` promoted preview → stable and bound in both skins; new stable op `helical_sweep` records the new `Feature::HelicalSweep` (additive variant, no format bump). `sweep` takes its path edge from an explicit `path_body` (referenced, not consumed); `loft` moves to `Model` (sections are sketch bodies, consumed like boolean tools); profiles sweep/loft their EXACT analytic wires (the PK-408 rule). Kernel: `HelixSpec` gains `start_angle` + `left_handed` (exact phased arcs); arc-path sweeps no longer promote curved walls to planes (collinear-rail gate); multi-turn sweeps no longer trip the self-intersection scan (knot-span v-density); mixed open/closed loft sections are rejected. Emission: `loft`, `pipe`, `helical_sweep` are spelled; `sweep` still declines by name (its cross-body path query has no verified re-spelling yet). | `sweep`, `loft`, `pipe` preview → stable; `helical_sweep` added stable |
| 2026-09-02 | Conformance checklist (PK-405): `docs/reference/scripting-surface-checklist.json` is generated from this document's op tables (§§2, 4.1, 5, 6.3, 7, 8), the §3 builder code blocks, the §4.3 prose, the §6.1 grammar, and the §10 class tables (`popkern_rhai::checklist_from_spec`); `crates/popkern-rhai/tests/conformance.rs` fails when the checked-in file is stale (regenerate with `POPKERN_UPDATE_CHECKLIST=1`), audits the Rhai engine (registered names + binding tables: positional / option-key / documented-default coverage of every spec parameter), the Rust skin (`pub fn` per stable op), and the shared payload projection (every §10 class and field), pins a negative fixture, and runs every stable op in one script. Preview ops the Rhai skin does not bind are allow-listed with a reason. Verbs landed to close the audit: `set_expr/clear_expr/set_bounds/clear_bounds/rename_param`, `add_derived_config/duplicate_config/delete_config/rename_config/clear_suppression/config`, `.close()`; `FeatureFailed.report` and `ResolutionFailed.resolved/expected` are now carried. | none |
| 2026-09-02 | Canonical emission + migration verbs in Rhai (PK-404): `popkern_rhai::to_rhai(&mut Model)` / script-side `m.to_rhai()` shows a document as §§2–9 spellings (timeline order, one statement per feature, scalars as `p("name")` / `expr(..)` / literal); a recorded query is re-emitted only when it still resolves to the recorded name snapshot on the body the feature saw, else the snapshot is spelled `edges().named("label/Role/i").or(..)` (D3). New spellings: entity names as text `label/Category/index[@ancestors]` (`Model::spell_entity`/`parse_entity`, `named(text)` on both selectors, `label_entity(name, text)`); `label` option on `union/subtract/intersect/translate/rotate` (`Body::boolean_with/translate_with/rotate_with`); `revolve`/`rotate` axes by query (`AxisSpec`, recorded as the `axis` selection); single-entity arguments now record their query (`plane`, `end_face`, `face`, `axis` roles). The selection grammar has a text reader (`popkern_script::parse_edges/parse_faces`). Rhai verbs: `insert_before/insert_after(anchor, \|m, b\| ..)`, `reorder`, `rollback_to/rollback_to_end`, `undo/redo`, `add_config/set_active_config/active_config/set_override/clear_override/set_suppression/configs`; report maps carry `is_clean`, `failures`, `downgraded_blends`. Not spelled (emission declines by name): raw_body, constraint sketches, sweep/loft/pipe/thicken/offset_solid/draft/patterns/mirror/reference, `up_to_vertex`, extrude direction/draft/second direction, revolve boss/cut, `suppress_if`, parameter bounds, derived configurations, expression overrides, stale multi-face queries. | none |
| 2026-09-02 | Error payload projection is one function for every skin (PK-403): `ScriptError::payload()` in `popkern-script` yields the §10.1 envelope as JSON; `InvalidParameter` records carry the **surface** `parameter` with the kernel spec path beside it as `kernel_path` (`specs[0].radius` → `radius`, `options.profile.rho` → `profile`, `dist_b` → `distance_2`, `radius_end`/`….end` → `law`, `specs[i].edge` → `edges`; `popkern_script::surface_parameter`). Rhai skin (PK-402/403): the runtime error's value IS the payload map, positioned at the raising call (a method call is positioned at its method name), `try`/`catch` receives the same map; `ScriptRunError::render()` prints the class, `script:line:column`, the message, and one remediation line per hint. | none |
| 2026-08-30 | Full-circle profiles extrude to one periodic wall on every route (PK-387): constraint-dialect circles, every end condition, second direction; draft records a fallback. Exact-route names follow the extrude contract (`Modified 0`/`Cap 0`/`Lateral i`). Every sketch is its own body; Boss/Cut merge into the host named by the on-face reference or the up-to/offset target; `Body::extrude_with` returns the served host. | none |
| 2026-08-30 | Structured failure payload landed (PK-388): `FeatureFailed.error` is a `FeatureFailure` (§10.3); the PyO3 probe's JSON is serde of the same value (`kind`/`blend_kind`/`reason_kind`, `edge_name`). | none |
| 2026-08-30 | Bodies landed (PK-385): `tool`, `raw_body`, `label_body`/`body`/`body_label`, `alive_bodies`, `set_active`, `Body::is_alive`; booleans consume the tool and reject dead bodies (`BodyNotAlive`). | none |
| 2026-08-30 | Migration verbs landed (PK-384): `modify` via `Edit` (scalar roles, `continuity`, `edges`/`faces`), `insert_before/after`, `reorder`, `suppress/unsuppress`, `delete`, `rollback_to/rollback_to_end`, `dependents`, `undo/redo`, `describe`, `params`, `param_expr`, configurations. Documents now persist `SelectionRecord`s per feature (§6.2 describe text). | none |
| 2026-08-30 | Rust skin landed (PK-382): `model/save/open/replay/to_brep/to_part/to_step`, gates, `profile` (rect/circle/polygon/path, base + offset + on-face planes), `extrude` (all ends), `revolve`, `fillet` (recipe), `chamfer` (equal), `shell`, booleans, `translate/rotate`, `param/set_param`, feature/entity labels, `feature/features/evidence/report`. `FeatureFailed.error` was a string until PK-388. | none |
| 2026-08-30 | Initial normative surface (PK-381). Stable set: `model/open/save/replay/to_brep/to_part/to_step`, gates, `profile`, `extrude/revolve`, booleans, `fillet/chamfer`, `shell`, `translate/rotate`, parameters (§4.1), labels, `feature/features/modify/suppress/unsuppress`, `report/evidence`. Everything else `preview`; `promote_to_constraint` reserved. | — |

Promotion rule: `preview → stable` requires (1) a parity fixture
exercising the op through record→replay, (2) an error-payload row in
§10 exercised by a test, (3) an entry in this table.

---

## 14. Design review record (2026-08-30, `/design-review`)

1. **Consumer perspective.** §12 is the README example. First draft
   was awkward in two places and was changed: feature ops now return
   the served `Body` (chaining in timeline order) instead of a bare
   `Feature`, and the tool-model closure is `model.tool(|m| …)` so it
   no longer collides with the `model.body(id)` lookup.
2. **Error UX.** Every class in §10 carries the values a dialog needs
   (`BlendDeclined.failures[].max_radius` + `remediation.value`,
   `EmptySelection.available`, `SelectionCardinality.expected/got`,
   `ResolutionFailed` citing label + authored text + counts).
   `InvalidParameter.parameter` is the **surface** name with the kernel
   path beside it — the E2 probe measured `specs[0].radius` leaking.
   Messages are redundant with fields, never the only carrier.
3. **Composability.** Patterns/mirror take `[Feature]` seeds; every
   selection argument is the same query grammar; bodies are first-class
   for booleans; migration verbs accept labels everywhere ids are
   accepted.
4. **Naming.** Operation identifiers are exactly
   `Feature::operation_label()`; predicates are exactly
   `popkern_part::select`; enums are the `*Desc`/`*Kind` variants in
   `snake_case`. Gap fixed: `profile`/`sketch` were both called "sketch"
   in prose — the analytic builder is `profile`, the constraint builder
   is `sketch`, both record `Feature::Sketch`.
5. **UI implications.** Preview mode → `preview` (reserved, §7);
   slider-friendly → `max_fillet_radius` + declared bounds now,
   `param_range` reserved (§8); undoable → every verb is a
   `TreeCommand`; nameable → `describe` added (§7, stable);
   progressive disclosure → every op has a short form plus keyword
   options; selection-aware → all ops take queries; measurable →
   `measure` (§8). Error messages: §10.
6. **Commercial calibration.** Commercial kernels' and OpenCASCADE's scripting layers all
   converge on: persistent naming for selections across edits, an
   explicit feasibility precheck per op, transactional feature
   application, and structured (not string) failure reasons. Matched
   here by D3 name snapshots, `precheck_*`/`*_range`, §1.6, and §10.
   Where this surface is deliberately narrower: one document = one
   part (no assemblies), analytic profiles are not silently promoted
   to constraint sketches (D2).
7. **Verdict.** API shape: `Model`, `Body`, `Feature`, `Param` handles;
   the op table of §5; the data selector grammar of §6; verbs of §7;
   payload classes of §10. Key choices: literals become feature
   parameters (every number is editable later); selectors are data,
   closures are opaque and Rust-only; tool bodies live in the same
   tree; stable set kept small (§13). UI readiness: satisfied — undo,
   naming, selection, measurement, errors; deferred — `preview`,
   `param_range`, `promote_to_constraint` (reserved with named
   unblockers). Recorded in `docs/journal/decisions.md`.
