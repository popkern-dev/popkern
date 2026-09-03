# popkern scripting guide

One page that gets a person — or a language model given this page as
context — productive in popkern's Rhai dialect. Everything here is
executable: every code block on this page runs in the test suite and
must serve (or raise exactly the error class its comment names), and
every stable verb, selector predicate and error class of the
[scripting surface](../architecture/scripting-surface.md) is named
here by its exact spelling (`crates/popkern-rhai/tests/guide.rs`
checks both). Try any block in the playground (`web/playground/`) or
with `popkern run file.rhai --png out.png --step out.step`.

## 1. The mental model in five sentences

1. A script **records a feature tree**; the tree is the truth, the
   script is one way of writing it down. `save("part.pkn")` stores the
   tree; `open("part.pkn")` gets it back; "show as code" prints the
   canonical script for any tree.
2. Every operation appends **one feature** and replays it immediately.
   A feature the kernel declines is not recorded: the model is exactly
   as it was before the call, and the error tells you why.
3. Selections are **queries evaluated once**, at the call. The
   document stores the entities' stable names, never the query as
   code — an edited script re-resolves, a saved document does not.
4. Every number you pass becomes a **parameter** (`extrude-1.distance`,
   `rounds.radius`) that sliders and `set_param` can edit later
   without changing the script.
5. Booleans, fillets and shells are **transactional and honest**:
   the replay report says what was measured (a G2 request that
   achieved G1 is reported as downgraded, never as a green tick).

## 2. Rhai in ninety seconds (for Python and JavaScript hands)

```rhai
let m = model();                      // `let`, not `var`; no `const` needed
let width = 2.4;                      // numbers are f64; 2 and 2.0 differ (2 is an int)
let plate = m.profile(xy()).rect(width, 1.6).sketch("plate").extrude(blind(0.4));
print(plate.label);                   // print writes to the run's output, not stdout
let rounds = plate.fillet(edges().straight().direction("z"), 0.1, #{ label: "rounds" });
m
```

- **Options ride in a trailing `#{ key: value }` map.** Unknown keys
  are an `ArgumentError`, so a typo cannot silently do nothing.
- **Enums are strings**: `"g2"`, `"cylinder"`, `"z"`, `"inward"`,
  `"cut"`. Recipes are constructor calls: `conic(0.7)`, `linear(0.2)`,
  `two_distance(faces(), 0.05)`.
- **Closures are `|args| expr`**, and they receive the model when a
  verb needs a nested recording: `m.tool(|m| ...)`,
  `m.insert_before(anchor, |m, b| ...)`.
- **The script's last expression is its value.** End with `m` so the
  host receives the model (the playground and `popkern run` use it).
- **Errors are maps.** `try { ... } catch (e) { print(e.kind) }` gives
  you the same structured payload the host prints; every error is
  positioned at the raising call.
- Not available on purpose: `eval`, `sleep`, `exit`, timestamps, file
  I/O beyond `open`/`save`/`to_step`. Scripts are pure functions of
  their text.
- Integers vs floats matter: `rect(2, 1)` is an `ArgumentError`
  (`w` expects a number); write `rect(2.0, 1.0)`.

## 3. The boss plate, line by line

```rhai
let m = model(#{ unit: "m" });
// A plate: a rectangle on the XY plane, extruded 0.4 up.
let plate = m.profile(xy()).rect(2.4, 1.6).sketch("plate_sketch")
             .extrude(blind(0.4));
// A boss: a circle 0.2 above the base, extruded 0.8, recorded as a
// tool body so the union can consume it.
let boss = m.tool(|m| m.profile(xy_at(0.2)).circle(0.35).sketch("boss_sketch")
             .extrude(blind(0.8)));
plate = plate.union(boss);
// Rounds: the boss-to-plate junction, the plate's top edges, the
// plate's vertical corners — all by geometry, none by index.
plate.fillet(
    edges().between("cylinder", "plane").at_height(0.4)
        .or(edges().straight().at_height(0.4))
        .or(edges().straight().direction("z")),
    0.08, #{ continuity: "g1", label: "rounds" });
m.assert_watertight();
print(m.feature("rounds").describe());
m
```

`plate = plate.union(boss)` reads oddly to Rust eyes: every verb
returns the body it served, and the union's result replaces the
receiving body, so rebinding keeps the variable pointing at the live
body.

## 4. Planes, profiles, sketches

A profile is a closed outline on a plane; `sketch(label)` records it
as a feature and starts a new body.

| Spelling | Meaning |
|---|---|
| `xy()`, `xz()`, `yz()` | the base planes |
| `xy_at(z)`, `xz_at(y)`, `yz_at(x)` | a base plane offset along its normal. In-plane coordinates `(u, v)` are `(x, y)` on XY, `(x, z)` on XZ and `(y, z)` on YZ; `plane_at` uses the frame of its normal (u along the least-aligned world axis) |
| `on_face(faces()...)`, `on_face_offset(faces()..., d)` | the plane of a planar face (the query must pick exactly one face) |
| `plane_at([ox, oy, oz], [nx, ny, nz])` | a custom plane by origin and normal |
| `.rect(w, h)`, `.circle(r)`, `.polygon([[x, y], ...])` | closed outlines centred on the plane origin |
| `.path(x0, y0).line_to(x, y).arc_to([vx, vy], [x, y]).close()` | a path of lines and three-point arcs, closed back to its start |
| `.at(cx, cy)` | place the outline's centre |

```rhai
let m = model();
let knob = m.profile(xz())
    .path(0.0, 0.0).line_to(0.45, 0.0).line_to(0.45, 0.25).line_to(0.30, 0.25)
    .line_to(0.30, 0.55).line_to(0.0, 0.55).close()
    .sketch("knob_section")
    .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
m
```

Profiles work on every plane, not only XY — an L-bracket sketched on
`yz()` extrudes along +x.

## 5. Features

Every feature verb is a method on a body; every one returns the body
it served; every one takes an optional trailing options map with at
least `label`.

| Verb | Positional arguments | Options | Notes |
|---|---|---|---|
| `extrude(end)` | `end`: `blind(d)`, `through_all()`, `up_to_surface(faces()...)`, `offset_from_surface(faces()..., d)`, `mid_plane(half_d)` | `label`, `mode` (`"new_body"` default, `"boss"`, `"cut"`) | `boss`/`cut` fold the new body into the body the sketch sits on (`on_face` sketches); a base-plane sketch is its own body — combine it with `union`/`subtract`. |
| `revolve(axis, end)` | `axis`: `axis_at([o], [d])` or `axis(edges()...)`; `end`: `blind(angle_deg)`, `full_circle()`, `mid_plane(half_angle)` | `label` | Revolved spans become exact planes, cylinders, cones and tori. |
| `union(tool)`, `subtract(tool)`, `intersect(tool)` | a body from `m.tool(...)` | `label` | The tool body is consumed; rebinding (`plate = plate.union(boss)`) keeps the live body. |
| `fillet(edges, radius)` | `edges()...` query; `radius` scalar | `label`, `continuity` (`"g0"`, `"g1"`, `"g2"`), `cap` (`"open"`, `"flat"`, `"rounded"`, `"tangent_match"`), `profile` (`circular()`, `conic(rho)`, `curvature_continuous()`), `interference` (`"decline"`, `"serve_non_interfering"`), `law` (`constant()`, `linear(end_radius)`, `by_stations([[t, r], ...])`) | A radius the edge cannot carry is a `FeatureFailed` with the maximum radius and a `ReduceRadius` remediation. |
| `chamfer(edges, distance)` | as fillet | `label`, `kind` (`equal()` default, `two_distance(faces()..., d2)`, `distance_angle(faces()..., angle_deg)`) | The face query of a kind must resolve on the chamfered body. |
| `shell(faces_to_remove, thickness)` | `faces()...` (may be `.allow_empty()` for a closed hollow) | `label`, `direction` (`"inward"` default, `"outward"`, `"midplane"`), `join` (`"arc"`, `"intersection"`), `trim` (`"trim"`, `"decline"`) | |
| `translate(dx, dy, dz)`, `rotate(axis, angle_deg)` | rigid motions; the distances and the angle are scalars (`p()`/`expr()` apply) | `label` | Usually applied to tool bodies before a boolean. Selector arguments (`at_height(z)`, `chord_length_between`) are plain numbers, not parameters. |

```rhai
let m = model();
let block = m.profile(xy()).rect(1.2, 0.8).sketch("block_sketch").extrude(blind(0.5));
let bore = m.tool(|m| m.profile(xy_at(-0.1)).circle(0.18).sketch("bore_sketch").extrude(blind(0.7)));
block = block.subtract(bore);
m.label_edge("bore_rim", edges().circular().at_height(0.5));
block = block.fillet(edges().labeled("bore_rim"), 0.05, #{ label: "rim_round" });
block = block.chamfer(edges().straight().at_height(0.0), 0.04,
    #{ kind: distance_angle(faces().facing_neg("z"), 30.0), label: "lead_in" });
m.assert_watertight();
m
```

A shell hollows a body from a removed face; round the corners first,
and keep the wall thinner than any round it must follow:

```rhai
let m = model();
let box = m.profile(xy()).rect(1.6, 1.0).sketch("outline").extrude(blind(0.6));
box = box.fillet(edges().straight().direction("z"), 0.15, #{ label: "corners" });
box = box.shell(faces().facing("z").at_height(0.6), 0.06, #{ direction: "inward", label: "hollow" });
m.assert_watertight();
m
```

Scalars: any numeric argument accepts a literal, `p("name")` for a
declared parameter, or `expr("text")` for an expression over
parameters.

## 6. Selecting edges and faces

Queries are data: `edges()` or `faces()` followed by predicates, all
of which must hold; `.or(...)` joins whole queries; `.allow_empty()`
turns "matched nothing" from an error into a recorded no-op. A query
that must name exactly one entity (a `label_edge`, an `on_face`
plane, a chamfer's reference face) raises `SelectionCardinality` when
it matches more.

| Edge predicate | Meaning |
|---|---|
| `circular()`, `straight()`, `conic()`, `nurbs()` | the edge's curve kind |
| `between("cylinder", "plane")` | the two adjacent surface kinds (`"plane"`, `"cylinder"`, `"cone"`, `"sphere"`, `"torus"`, `"nurbs"`) |
| `at_height(z)` | the edge lies at that z |
| `direction("z")` / `direction([x, y, z])` | a straight edge along that direction |
| `chord_length_between(min, max)` | endpoint distance in that range (a closed circle's chord is 0; use it to tell straight edges apart by length, and `at_height` plus the feature order to tell concentric rims apart) |
| `not_seam()` | exclude periodic seams (they are excluded from blends anyway) |
| `labeled("name")` | an edge labelled with `label_edge` |
| `named("rounds/Generated/3")` | the entity's stable name (what "show as code" emits when a query went stale) |

| Face predicate | Meaning |
|---|---|
| `kind("plane")` | surface kind |
| `facing("z")`, `facing_neg("z")` | outward normal along ± an axis (or `[x, y, z]`) |
| `at_height(z)` | a planar face at that z |
| `labeled("name")`, `named("plate_sketch/Generated/0")` | as for edges |

```rhai
let m = model();
// A triangular prism: its top face has exactly one edge along x.
let prism = m.profile(xy()).polygon([[0.0, 0.0], [1.0, 0.0], [0.5, 0.8]]).sketch("tri").extrude(blind(1.0));
m.label_face("top", faces().facing("z").at_height(1.0));
m.label_edge("front_top", edges().straight().at_height(1.0).direction("x").chord_length_between(0.9, 1.1));
prism.fillet(edges().straight().direction("z").or(edges().labeled("front_top")), 0.08, #{ label: "rounds" });
prism.chamfer(edges().circular().at_height(9.0).allow_empty(), 0.02, #{ label: "no_rims" });   // a no-op, recorded as such
print(m.labels().len());
m
```

Labels: `label_edge(name, query)`, `label_face(name, query)`,
`unlabel(name)`, `rename_label(old, new)`, `labels()`. Entity labels
and feature labels are separate namespaces; a duplicate is a
`LabelError`.

## 7. Parameters, expressions, configurations

```rhai
let m = model();
m.param("wall", 0.3, "m");
m.param("span", 2.0, "m");
m.param_expr("lip", "wall * 0.5", "m");
let plate = m.profile(xy()).rect(p("span"), 1.2).sketch("plate_sketch").extrude(blind(p("wall")));
plate = plate.fillet(edges().straight().direction("z"), 0.1, #{ label: "corners" });
let rib = m.tool(|m| m.profile(xy_at(p("wall"))).rect(expr("span * 0.8"), 0.12).sketch("rib_sketch").extrude(blind(p("lip"))));
plate = plate.union(rib);
m.set_bounds("wall", 0.1, 1.0);
m.set_expr("lip", "wall * 0.4");
m.clear_expr("lip");
m.clear_bounds("wall");
m.rename_param("lip", "rib_height");
m.set_param("wall", 0.35);
m.add_config("heavy", "double wall");
m.set_override("heavy", "wall", 0.6);
m.set_active_config("heavy");
m.replay();
m.set_active_config(());
print(m.params().len());
print(m.configs().len());
m
```

Units: `"m"`, `"mm"`, `"cm"`, `"in"`, `"deg"`, `"rad"`, `"ratio"`.
Every literal you pass to a feature becomes `feature.role`
(`extrude-1.distance`, `corners.radius`) — `params()` lists them and
`set_param` edits them. Configurations (`add_config`,
`add_derived_config(name, parent)`, `duplicate_config`,
`delete_config`, `rename_config`, `set_override`, `clear_override`,
`set_suppression(config, feature, true)`, `clear_suppression`,
`set_active_config`, `configs()`, `config(name)`) are preview verbs:
they work, and their spellings may still move.

## 8. Reading the model back

```rhai
let m = model();
let cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));
cube.fillet(edges().straight().direction("z"), 0.1, #{ label: "rounds" });
let rounds = m.feature("rounds");
print(rounds.describe());               // "rounds — fillet 0.1 on 4 edges: ..."
print(rounds.evidence().outcome);       // the feature's replay row
print(m.features().len());              // timeline order
print(m.report().failures.len());       // spec §11 report: features, warnings, blend evidence
print(m.bodies().len());
print(m.active().label);
m.label_face("top", faces().facing("z").at_height(1.0));
print(m.entity_spelling("top"));        // "extrude-1/Cap/0": the label's stable entity name
print(m.to_rhai());                     // the document as canonical Rhai ("show as code")
m
```

`report()` carries, per feature, the outcome, the wall time, kernel
warnings, and for blends the measured continuity beside the requested
one. `assert_watertight()` and `assert_step_roundtrip(path)` are
gates: they raise `NotWatertight` / `ValidationFailed` instead of
returning.

## 9. Editing a saved document (migration verbs)

The same verbs a GUI's history tree would use. Stable:
`feature(label)`, `features()`, `describe()`, `modify(#{...})`,
`suppress()`, `unsuppress()`, `evidence()`. Preview: `insert_before`,
`insert_after`, `reorder`, `delete`, `rollback_to`,
`rollback_to_end`, `undo`, `redo`, `dependents`, `to_rhai`.

```rhai
let m = model();
let cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));
cube = cube.fillet(edges().straight().direction("z"), 0.1, #{ label: "rounds" });
let rounds = m.feature("rounds");
rounds.modify(#{ radius: 0.12 });
rounds.suppress();
rounds.unsuppress();
m.insert_before(rounds, |m, b| b.chamfer(edges().straight().at_height(0.0), 0.05, #{ label: "foot_break" }));
m.reorder(m.feature("foot_break"), 2);
m.rollback_to(m.feature("foot_break"));
m.rollback_to_end();
m.rename(rounds, "corner_rounds");
print(m.feature("corner_rounds").dependents().len());
m.replay();
m
```

`insert_before(anchor, |m, b| ...)` hands the closure the model and
the body as it was **before** the anchor; the closure records exactly
one feature there. Names are the truth: a feature inserted earlier
keeps every later selection valid because later features reference
entities by stable name, not by index.

## 10. When it goes wrong: the error classes

Every failure is a map with a `kind` and the fields below; the host
prints it as `Kind at file:line:column: message` plus one remediation
line per hint. In a script, `try`/`catch` receives the same map.

| Class | When | Fields |
|---|---|---|
| `ArgumentError` | wrong type, unknown option key, unknown parameter name, bad enum spelling | `argument`, `reason` |
| `EmptySelection` | a query matched nothing (add `.allow_empty()` if a no-op is intended) | `selection`, `entity` |
| `SelectionCardinality` | a single-entity argument matched several | `selection`, `expected`, `got` |
| `FeatureFailed` | the kernel declined the feature (its `error` carries the class below) | `feature`, `error`, `report` |
| `BlendDeclined` | inside `FeatureFailed`: a fillet/chamfer the geometry cannot carry — each failure names the edge, the `max_radius` and a remediation (`ReduceRadius`) | `selection`, `blend_kind`, `detail`, `failures` |
| `InvalidParameter` | inside `FeatureFailed`: a negative radius, a distance out of range | `parameter`, `kernel_path`, `value`, `reason` |
| `BooleanFailed`, `SweepFailed`, `ShellFailed`, `DraftFailed` | a boolean, extrude/revolve, shell or draft declined | `operation`/`feature`, `detail` |
| `ResolutionFailed` | on replay, a stored name no longer resolves (after a migration) | `feature`, `name`, `selection`, `resolved`, `expected` |
| `LabelError` | unknown or duplicate label | `label`, `reason` |
| `ParameterError` | bad bounds, dimension mismatch, dangling expression | `parameter`, `reason` |
| `TimelineError` | a consumed tool body, an invalid reorder, an insertion that breaks a dependency | `reason` |
| `DocumentError`, `IoFailed`, `NotApplicable` | a document that will not open, a path that will not write, an op the body cannot take | `reason` / `path`, `detail` |
| `NotWatertight`, `ValidationFailed` | a gate failed | `open_edges` / `feature`, `detail`, `issues` |
| `CoreError` | a kernel invariant tripped (report it) | `detail` |
| `Cancelled`, `BudgetExceeded` | the host stopped the run (playground and CLI budgets: `budget` is `operations` or `wall_clock`) | `budget` |
| `SyntaxError` | Rhai could not parse the script | position |

```rhai
// raises: FeatureFailed
let m = model();
let cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));
cube.fillet(edges(), 0.9);
```

```rhai
// raises: EmptySelection
let m = model();
let cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));
cube.fillet(edges().circular(), 0.1);
```

```rhai
// raises: ArgumentError
let m = model();
let cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));
cube.fillet(edges(), 0.1, #{ continuity: "g1", colour: "red" });
```

```rhai
let m = model();
let cube = m.profile(xy()).rect(1.0, 1.0).sketch("cube").extrude(blind(1.0));
try {
    cube.fillet(edges(), 0.9);
} catch (e) {
    print(e.kind);                      // "FeatureFailed"
    print(e.error.kind);                // "BlendDeclined"
    print(e.error.failures[0].max_radius);
}
m
```

## 11. Budgets, size, speed

- The playground runs a script under an operation budget and a wall
  clock budget (about 20 s); `popkern run --max-ops N --timeout S`
  does the same. A script that overruns is stopped with
  `BudgetExceeded`; the last good model stays on screen.
- G2 fillets with corners are the slow path today (seconds each,
  tens of seconds for a fully rounded cube); G1 rounds and chamfers
  are fast. Prefer G1 unless curvature continuity is the point.
- The wasm engine runs about three times slower than native; a
  twenty-feature part builds in well under a second natively and in a
  few seconds in the browser. A 200-feature document (a hundred
  bodies) took ten seconds in the browser and ten megabytes of wasm
  memory: tens of features stay interactive, hundreds take seconds.
- Cut-mode extrudes (`#{ mode: "cut" }`) need a sketch on a face of
  the target body (`on_face(...)`); from a base plane, record the
  cutter with `m.tool(...)` and `subtract` it.
- A single expression may nest to 64 levels at the top level and 32
  inside a `fn` body (the same limit in every build — CLI, tests, and
  the wasm playground). Computed-geometry arithmetic rarely reaches
  it, but a deep chain like `a + b*c/d - f(g(x))` repeated in a point
  literal can; if you hit `Expression exceeds maximum complexity`,
  bind sub-results to `let` locals (`let ang = rot - half; [r*cos(ang),
  r*sin(ang)]`), which also reads better.

## 12. Cheat-sheet

```text
model(#{ unit: "mm" })            open(path)   save(path)   to_step(path)   to_rhai()
m.profile(PLANE).rect|circle|polygon|path…close().at(cx, cy).sketch(label)
PLANE: xy() xz() yz() xy_at(z) xz_at(y) yz_at(x) on_face(F) on_face_offset(F, d) plane_at(o, n)
END:   blind(d) through_all() up_to_surface(F) offset_from_surface(F, d) mid_plane(h) full_circle()
AXIS:  axis_at([ox,oy,oz], [dx,dy,dz])   axis(E)
body.extrude(END, #{label, mode})        body.revolve(AXIS, END, #{label})
body.union(tool) subtract(tool) intersect(tool)      tool = m.tool(|m| …)
body.fillet(E, r, #{label, continuity, cap, profile, interference, law})
body.chamfer(E, d, #{label, kind})       body.shell(F, t, #{label, direction, join, trim})
body.translate(dx, dy, dz)               body.rotate(AXIS, deg)
E: edges().circular().straight().conic().nurbs().between(k, k).at_height(z).direction(d)
   .chord_length_between(a, b).not_seam().labeled(l).named(n).or(E).allow_empty()
F: faces().kind(k).facing(d).facing_neg(d).at_height(z).labeled(l).named(n).allow_empty()
m.param(name, v, unit) param_expr(name, text, unit) set_param set_expr clear_expr
  set_bounds clear_bounds rename_param params()      p("name")  expr("a * 2")
m.label_edge label_face unlabel rename_label labels()
m.feature(label) features() report() bodies() active() assert_watertight() assert_step_roundtrip(path)
feature.describe() evidence() modify(#{…}) suppress() unsuppress() dependents()
m.insert_before(f, |m, b| …) insert_after reorder(f, i) delete(f) rollback_to(f) rollback_to_end() undo() redo()
m.add_config add_derived_config duplicate_config delete_config rename_config set_active_config
  set_override clear_override set_suppression clear_suppression configs() config(name)
```

Gallery scripts that exercise all of this: `crates/popkern-rhai/examples/`
(`boss_plate`, `revolved_knob`, `shelled_enclosure`,
`chamfered_bracket`, `bored_block`, `tapered_fillet`,
`parametric_plate`, `keyed_shaft`).
