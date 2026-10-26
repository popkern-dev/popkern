# popkern scripting guide

One page that gets a person — or a language model given this page as
context — productive in popkern's Rhai dialect. Everything here is
executable: every code block on this page runs in the test suite and
must serve (or raise exactly the error class its comment names), and
every stable verb, selector predicate and error class of the
[scripting surface](../architecture/scripting-surface.md) is named
here by its exact spelling (`crates/popkern-rhai/tests/guide.rs`
checks both). Try any block in popkern live (`web/live/`) or
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
  host receives the model (popkern live and `popkern run` use it).
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
| `fillet(edges, radius)` | `edges()...` query; `radius` scalar | `label`, `continuity` (`"g0"`, `"g1"`, `"g2"`), `cap` (`"open"`, `"flat"`, `"rounded"`, `"tangent_match"`), `profile` (`circular()`, `conic(rho)`, `curvature_continuous()`), `interference` (`"decline"`, `"serve_non_interfering"`), `law` (`constant()`, `linear(end_radius)`, `by_stations([[t, r], ...])`), `corner` (`"smooth"`, `"mitre"`: trim every strip pair at a corner to its mutual intersection), `corner_cap` (`"auto"`, `"setback"`, `"fitted"`, `"classic"`, `"rolling_ball"`: which construction closes a corner where three or more rounds meet; `"rolling_ball"` gives the exact rolling-ball sphere where three rounds of one radius meet on flat faces), `corner_cap_strict` (true: that construction or a failure, never a fallback; needs `"setback"`, `"fitted"` or `"rolling_ball"`), `setback` (a positive distance: each round ends that far from a corner where three or more rounds meet), `setback_ends` (a list of [edges, at, distance] groups: the end of each edge of the first selection at the corner it shares with an edge of the second, overriding `setback`) | A radius the edge cannot carry is a `FeatureFailed` with the maximum radius and a `ReduceRadius` remediation. A mitre with no point shared by every strip is a `FeatureFailed` naming the mitre's reason. A strict `corner_cap` that cannot serve a corner is a `FeatureFailed` naming the corner and every decline. A `setback` the corner cannot take is a `FeatureFailed` naming the admissible range. |
| `fillet_each(groups)` | `[[edges()..., radius], ...]`: one radius per group of edges | `label`, `continuity`, `cap`, `profile`, `interference`, `corner`, `corner_cap`, `corner_cap_strict`, `setback`, `setback_ends` (as for fillet; radius laws are not taken) | One fillet whose groups carry different radii: the edges meeting at a corner blend together as one corner even when their radii differ. An edge selected by two groups is refused (the message names every overlapping pair). The first group's radius is the feature's `radius` parameter, later groups' are `radius_1`, `radius_2`, … |
| `chamfer(edges, distance)` | as fillet | `label`, `kind` (`equal()` default, `two_distance(faces()..., d2)`, `distance_angle(faces()..., angle_deg)`) | The face query of a kind must resolve on the chamfered body. |
| `shell(faces_to_remove, thickness)` | `faces()...` (may be `.allow_empty()` for a closed hollow) | `label`, `direction` (`"inward"` default, `"outward"`, `"midplane"`), `join` (`"arc"`, `"intersection"`), `trim` (`"trim"`, `"decline"`) | |
| `translate(dx, dy, dz)`, `rotate(axis, angle_deg)` | rigid motions; the distances and the angle are scalars (`p()`/`expr()` apply) | `label` | Usually applied to tool bodies before a boolean. Selector arguments (`at_height(z)`, `chord_length_between`) are plain numbers, not parameters. |
| `sweep(path_body, path_edge)` | the body carrying the path edge, then an edge query on it (exactly one edge) | `label`, `mode`, `twist` | The receiver's last sketch is the profile; it contributes its shape — placement comes from the path. `twist` in degrees over the full path. The path body is referenced, not consumed. |
| `loft(sections)` | on `Model`: `m.loft([s0, s1, …])` — ordered sketch bodies, each contributing its most recent sketch face | `label`, `mode`, `closed` | Records on the first section's body and consumes the others. Sections must be all closed or all open; matching point counts keep the surface smooth (circle sections stay exact circles). |
| `pipe(path_edge, radius)` | an edge query on the receiving body (exactly one edge), then the tube radius | `label`, `mode` | The tube is a new body; the path body keeps its edge. |
| `bolt(spec, length)` | on `Model`: `spec` is a designation string (`"M6"`, `"M6x0.75"`, `"1/4-20"`, `"#8-32"` — millimetre-magnitude numbers bound in model units) or `#{ major: .., pitch: .. }`; `length` is under-head (hex/socket) or overall (countersunk) | `label`, `head` (`"hex"` default, `"socket"`, `"countersunk"`), `thread_length`, `left_hand`, `starts`, `tip_chamfer` | Real modeled threads (PK-460). Bodies build in the insertion frame — axis +Z, origin at the bearing plane, threads toward +Z — and place with `translate`/`rotate`. Author fastener models in a default-unit session with millimetre-magnitude numbers. |
| `threaded_rod(spec, length)` | as `bolt` | `label`, `left_hand`, `starts`, `tip_chamfer` | Full-length external thread, both ends free-running. |
| `nut(spec)` | as `bolt` (no length) | `label`, `clearance`, `thickness`, `left_hand`, `starts`, `corner_chamfer`, `mouth_chamfer` | `clearance` grows the internal form for printed fits (0.15–0.3 typical); `thickness` overrides the ISO 4032 nominal; `corner_chamfer: false` / `mouth_chamfer: false` drop the 30° hex corner and 45° thread-mouth chamfers (both on by default). |
| `tapped_hole(axis, spec, depth)` | `axis_at(origin, direction)` with the direction pointing INTO the material; `spec` as above; `depth`: `blind(d)` or `through_all()` | `label`, `clearance`, `left_hand`, `starts`, `entry_chamfer` | Bore at the internal minor diameter plus a modeled internal thread; a fastener on the same axis mates by construction. Blind taps thread to ¾ lead above the bottom (the tap-relief shoulder). |
| `clearance_hole(axis, spec, depth)` | as `tapped_hole` | `label`, `fit` (`"close"`, `"normal"` default, `"loose"`), `diameter`, `countersink`, `counterbore`, `entry_chamfer` | ISO 273 hole series, or an explicit `diameter`; `countersink`/`counterbore` sink the matching head standard plus the given extra room. |
| `hinge(kind, length)` | on `Model`: `kind` is `"butt"` (two leaves and a pin; many knuckles make a piano hinge) or `"print_in_place"` (the pin is part of leaf A, gaps printed open); `length` along the pin axis | `label`, `leaf_width`, `thickness`, `knuckles`, `pin`, `knuckle_diameter`, `clearance`, `gap`, `angle`, `holes`, `screw`, `countersink`, `only` | Returns `#{ leaf_a, leaf_b, pin }` built in the hinge frame — pin axis +Z from z = 0, leaf A toward −X, leaf B at `angle` (180 flat, 0 folded shut). Unset sizes follow the length; `screw: "M3"` adds countersunk holes (`holes:` per leaf, default 2); `only: "leaf_b"` records one leaf (an assembly part). One parameter set per hinge. |
| `hinge_component(owner, component)` | `owner` is the label of a hinge's first component; `component` `"leaf_a"`, `"leaf_b"` or `"pin"` | `label` | Another component of that hinge, sharing its parameters. |
| `living_hinge(style, span, length)` | `style` `"straight"` (a thin web with rounded roots) or `"slotted"` (serpentine kerf slots, for rigid materials); `span` is the flex zone's width, `length` along the bend axis (+Z) | `label`, `panel_width`, `thickness`, `web`, `root_radius`, `slot_pitch`, `slot_width`, `bridge` | One body: two panels and the flexure, mounting face y = 0. The style's own sizes only (`web`/`root_radius` straight, `slot_*`/`bridge` slotted). |
| `helical_sweep(axis, radius, pitch, turns)` | `axis_at([o], [d])` or `axis(edges()...)`; then helix radius, pitch (advance per turn), turns (fractional allowed) | `label`, `mode`, `taper_to`, `start_angle`, `left_hand` | Springs, coils, auger flights, thread forms. `taper_to` is the end radius of a conical helix; `start_angle` (degrees) phases multiple helical bodies apart — multi-start threads, blades; `left_hand: true` mirrors the winding. The profile contributes its shape; placement comes from the helix. |

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

`fillet_each` gives each group of edges its own radius in one
feature, and `at_vertex` picks the edges meeting at one corner. Here
a block's top corner meets three different radii; with `continuity:
"g2"` the corner is closed by a curvature-continuous cap:

```rhai
let m = model();
let block = m.profile(xy()).rect(1.0, 1.0).sketch("block").extrude(blind(1.0));
let corner = [0.5, 0.5, 1.0];
block.fillet_each([
    [edges().at_vertex(corner).direction("x"), 0.05],
    [edges().at_vertex(corner).direction("y"), 0.2],
    [edges().at_vertex(corner).direction("z"), 0.3],
], #{ continuity: "g2", label: "mixed_corner" });
m.assert_watertight();
m
```

A corner where three or more rounds meet is closed by one of several
constructions. By default the kernel tries the setback cap first, then
the fitted cap (a variational fit on the same corner), then the classic
corner patches. `corner_cap` states a preference: `"fitted"` tries the
fitted cap first, `"setback"` the setback cap, and `"classic"` skips
both. `"rolling_ball"` closes a corner where three rounds of one radius
meet on flat faces with the exact sphere of the ball that rolls them (a
G1 corner, smooth but not curvature-continuous); other corners keep the
setback cap. With `corner_cap_strict: true` the fillet uses that construction
or fails, never falling back. The report's `blend.corners` names what
closed each corner:

```rhai
let m = model();
let block = m.profile(xy()).rect(1.0, 1.0).sketch("block").extrude(blind(1.0));
block.fillet(edges().at_vertex([0.5, 0.5, 1.0]), 0.2, #{ continuity: "g2", corner_cap: "fitted", label: "fitted_corner" });
m.assert_watertight();
m
```

`setback` sets how far each round stops short of such a corner before
the cap takes over. The default is the shortest setback the corner
allows. A longer one spreads the corner over more of the part. The
corner is built at exactly that setback, or the fillet fails with the
range the corner admits. `setback_ends` sets one round's end at one
corner, naming the corner by the edges that meet there, for example
`setback_ends: [[edges().at_vertex(c).direction("z"), edges().at_vertex(c), 0.4]]`:

```rhai
let m = model();
let block = m.profile(xy()).rect(1.0, 1.0).sketch("block").extrude(blind(1.0));
block.fillet(edges().at_vertex([0.5, 0.5, 1.0]), 0.1, #{ setback: 0.3, label: "set_back_corner" });
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

The sweep family: `helical_sweep` coils a profile about an axis, and
`loft` skins ordered section sketches (each sketch is its own body;
the loft consumes them):

```rhai
let m = model();
let spring = m.profile(xy()).circle(0.04).sketch("wire_section")
    .helical_sweep(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), 0.3, 0.12, 6.5,
        #{ label: "coil" });
let s0 = m.profile(xy_at(2.0)).circle(0.3).sketch("horn_base");
let s1 = m.profile(xy_at(2.35)).circle(0.18).sketch("horn_waist");
let s2 = m.profile(xy_at(2.7)).circle(0.26).sketch("horn_lip");
let horn = m.loft([s0, s1, s2], #{ label: "horn" });
m.assert_watertight();
m
```

Hinges, like fasteners, are one call. A butt hinge with countersunk
screw holes, and a slotted living hinge beside it:

```rhai
let m = model();
let h = m.hinge("butt", 60.0, #{ knuckles: 5, screw: "M3", label: "door_hinge" });
let flex = m.living_hinge("slotted", 12.0, 40.0, #{ label: "flex" });
flex = flex.translate(0.0, 30.0, 0.0);
print(h.leaf_a.is_alive() && h.leaf_b.is_alive() && h.pin.is_alive());
m.assert_watertight();
m
```

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
| `at_vertex([x, y, z])` | an edge with an end at that point: the edges meeting at one corner |
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

### Projects: a part built by a series of scripts

A real part often grows in steps. A **project** keeps those steps as
ordered scripts, the way a database keeps migrations:

```text
my-bracket/
  popkern.toml            [project] name = "Bracket"
  migrations/
    001_plate.rhai        let m = model(); … m
    002_boss.rhai         let m = previous(); … m
    003_rounds.rhai       let m = previous(); … m
```

**The first migration** builds the document as any script does.

**Every later migration** starts from `previous()`, the document the
earlier migrations built, and edits it:
- it can add features;
- it can change parameters;
- it can insert before, reorder or suppress earlier features.

It must end with that same document (return it, or nothing). Starting
over with `model()` is refused (`DocumentError`, reason
`document_replaced`). Other models a migration creates are tools.

**Building:**
- `popkern build my-bracket` runs the migrations in order and reports
  one line per migration: the features it recorded and the earlier
  ones it edited.
- `popkern build my-bracket --pkn bracket.pkn` saves the document with
  a record of the migrations it applied. Building against that file
  later runs only the new migrations.
- If an applied migration changed since, the document is rebuilt from
  the folder (`--strict` refuses instead).
- `popkern new my-bracket` creates a project to start from.

**When an upstream edit breaks a later migration** (a label it selects
no longer exists, say), the failure names that migration's line and
lists the upstream migrations that changed since it last built.

**Paths.** File paths in a project (`open`, `save`, `to_step`) stay
inside the project folder.

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
| `Cancelled`, `BudgetExceeded` | the host stopped the run (popkern live and CLI budgets: `budget` is `operations` or `wall_clock`) | `budget` |
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

## 11. Inspecting surface quality

- The **Shade** buttons in the top bar switch between CAD (the default: grey faces with black edges, or a part's own colour and finish; the face or edge under the mouse is highlighted, and a click names it in the viewport's corner), zebra, and Gauss K / Mean H (gaussian and mean curvature). On CAD shading the **Color** buttons paint one colour per part or per face instead. On a phone the same modes are the chip row under the top bar. Edges are the faces' true boundaries, so tangent fillet seams are drawn too. Zebra draws reflection-line stripes computed per pixel from the surface normals: a tangency (G1) break kinks the stripes across a seam; a curvature (G2) break kinks how fast they bend while staying connected — orbit the model and watch the stripes flow across fillet seams.
- The curvature views paint analytic per-vertex curvature through a diverging palette: blue negative, white zero, red positive, ranged at the median nonzero magnitude so one spike cannot wash the map out. A constant-radius fillet reads as a uniform band in mean curvature; a value step across a seam is a curvature discontinuity.

## 12. Budgets, size, speed

- popkern live stops a build at its **time limit**: 24 s by default,
  set in the Build button's menu from 10 s up to no limit and
  remembered in this browser. The status line counts the seconds
  while a build runs;
  **Stop** always cancels. A script that overruns is stopped with
  `BudgetExceeded` and the last good model stays on screen. The
  interpreter's operation budget, a guard against runaway loops,
  scales with the limit and is lifted with it.
  `popkern run --max-ops N --timeout S` sets the same budgets from
  the command line.
- The default document and the gallery examples are **pre-built**,
  and the page always opens on one (an edited script from your last
  visit is kept aside: a notice above the editor restores it).
  Opening one shows the model, its feature tree (continuity already
  certified) and its parameters at once, without waiting for the
  kernel; the status line says served (pre-built) and how long it
  takes to build natively. Editing the script, moving a slider,
  exporting, show-as-code and the zebra view build the model live
  first; a heavy example may then need a longer time limit.
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
  the browser build). Computed-geometry arithmetic rarely reaches
  it, but a deep chain like `a + b*c/d - f(g(x))` repeated in a point
  literal can; if you hit `Expression exceeds maximum complexity`,
  bind sub-results to `let` locals (`let ang = rot - half; [r*cos(ang),
  r*sin(ang)]`), which also reads better.

## 13. Assemblies: parts built once, placed many times

A model with repeated parts (a bearing's balls, a gearbox's planets)
is an **assembly**: part definitions, each recorded by a closure on a
fresh model and built **once**, and named instances that place them.
An instance moves with `translate` / `rotate` (the later call applies
after the earlier, as on bodies) or `frame(origin, z_dir, x_dir)`;
`seat(axis_at(origin, direction))` places a fastener on a hole cut
along that axis, in the hole's own insertion frame, so its thread phase
matches the hole's. Placements are rigid (no scale or mirror). A script
may return an assembly instead of a model.

```rhai
let a = assembly("pins");
let pin = a.part("pin", |m| {
    m.param("radius", 0.1, "m");
    m.profile(xy()).circle(p("radius")).sketch("pin").extrude(blind(1.0));
});
for i in 0..4 {
    a.place(pin, "pin_" + i).translate(i * 0.5, 0.0, 0.0);
}
let pair = a.subassembly("pair", |s| {
    s.place(pin, "low");
    s.place(pin, "high").translate(0.0, 0.0, 1.5);
});
a.place(pair).rotate(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), 90.0).translate(0.0, 2.0, 0.0);
a.assert_watertight();
a
```

A **part family** is one parametric part: declare its default
parameters, and ask for other sizes with `variant` — each distinct
parameter set is built once, however often it is placed, and its name
spells the changed values (`block[w=3]`):

```rhai
let a = assembly("blocks");
let block = a.part("block", #{ w: 1.0, h: 0.5 }, |m, p| {
    m.profile(xy()).rect(p.w, p.h).sketch("block").extrude(blind(0.2));
});
a.place(block);
a.place(block.variant(#{ w: 3.0 })).translate(0.0, 2.0, 0.0);
a
```

An **assembly parameter** drives a parameter in several parts at
once — one slider for a value the parts share:

```rhai
let a = assembly("plates");
let thin = a.part("thin", |m| {
    m.param("t", 0.1, "m");
    m.profile(xy()).rect(1.0, 1.0).sketch("thin").extrude(blind(p("t")));
});
let wide = a.part("wide", |m| {
    m.param("t", 0.1, "m");
    m.profile(xy()).rect(2.0, 1.0).sketch("wide").extrude(blind(p("t")));
});
a.place(thin);
a.place(wide).translate(0.0, 2.0, 0.0);
a.param("t", 0.25, "m");
a.drive("t", [thin, wide]);       // both parts' own "t" now follow it
a
```

A hinge is an assembly of its leaves and pin, mated about the pin
axis, with its open angle as an assembly parameter (a slider in popkern
live); `a.swing` sweeps it and reports where the motion first collides:

```rhai
let a = assembly("lid_box");
let h = a.hinge("lid", "print_in_place", 60.0, #{ angle: 110.0 });
a.limits("lid_angle", 0.0, 270.0);
let swing = a.swing("lid_angle", #{ from: 180.0, to: 270.0, step: 15.0 });
print(swing.first_contact);
a
```

Any angle or distance mate can follow an assembly parameter the same
way: `a.angle(x, y, p("turn"), #{ about: base.datum("z") })` (with
`scale:` / `offset:` to map the parameter). A box and its lid placed as
separate parts hinge along their seam with
`a.hinge_bodies("lid", "print_in_place", box, lid, axis_at([x, y, z], [1.0, 0.0, 0.0]), #{ length: 50.0 })`:
the lid then opens by `lid_angle`. Mates can also take a line or plane
given in a part's own frame: `box.axis([o], [d])`, `box.plane([o], [n])`.

- STEP export writes the product structure: a product per part, an
  occurrence per instance with its placement and name, each part's
  solids once.
- `.pkn` saves the whole tree.
- Show-as-code emits every part's own recorded script and the
  placements.
- An assembly's parameters are qualified by their part: `pin:radius`
  (a CLI `--param pin:radius=0.2`, or the part's slider in popkern
  live). Changing one rebuilds that part once and moves every
  instance of it.
- In popkern live the inspector lists the instance tree:
  - hide and show instances;
  - select one;
  - **open** a part to see its features and parameters.
- Picks name the occurrence: `pin_2:pin/…`.
- The `part colours` view tints each part.

Instance and part names cannot contain `/` or `:` (the instance path
and entity separators). Placing a sub-assembly inside itself, a
duplicate sibling name, or a part whose closure built no solid each
raise an `AssemblyError` naming the instance, the part and the
reason. An error inside a part's closure is that operation's own
error, with the closure's line and `in part '…'`.

## 14. Films: presentation videos of a part or an assembly

A **film** is a script that animates the built document into a short
presentation video: captioned steps, the camera turning and moving,
parts sliding and turning into place, joints sweeping, a section
sweeping through. It reads the document and never changes it, so a
film reruns in a moment and never rebuilds the geometry. A project
keeps any number of them as `films/<name>.rhai` beside `migrations/`
(a single script `stem.rhai` takes `stem.<name>.film.rhai` beside it);
in popkern live, **+ film** adds one and its tab plays it, and Export
writes MP4 (Slides, Discord, social), vertical MP4, WebM, GIF or a PNG.
`popkern film <project | script> --out guide.mp4` renders one from the
command line.

The film below plays on the pins assembly above: `document()` is the
built document, `storyboard()` the film, and each step's closure adds
what moves during it.

```rhai-film
let doc = document();
let film = storyboard(#{ look: "toon", camera: #{ azimuth: -40.0, elevation: 25.0 } });
film.title("Pins", "four pins and a pair", 2.0, |s| { s.orbit(30.0); });
film.step("The pair drops in", "from 2 units up", 3.0, |s| {
    s.move("pair_1", #{ from: [0.0, 0.0, 2.0], duration: 2.0 });
    s.highlight("pair_1/high", #{ at: 2.0, duration: 1.0 });
});
film.step("One pin turns", "", 2.0, |s| {
    s.turn("pin_3", axis_at([1.5, 0.0, 0.0], [0.0, 0.0, 1.0]), 0.0, 180.0);
    s.camera(#{ zoom: 1.4, target: "pin_3" });
});
film.step("Inside", "a section across the row", 2.5, |s| {
    s.section(#{ from: "y", from_at: 1.0, to_at: 0.5 });
    s.fade("pair_1", 1.0, 0.3);
});
film
```

Every action takes `#{ at, duration, ease }`: `at` seconds into its
step (default 0), `duration` (default: to the step's end), `ease`
(`"inout"` default, `"in"`, `"out"`, `"linear"`). Before an action
starts it holds its first value; after it ends, its last, so moves
that come `from` an offset end on the real assembly.

| Call | What it does |
|---|---|
| `storyboard(#{ size, fps, look, background, camera })` | the film; `size: [1920, 1080]`, `fps: 30`, `look: "toon"` (flat colours, outlines, ground shadow), `"studio"` or `"cad"`; `background: "#ffffff"`; `camera: #{ azimuth, elevation, zoom, target }` |
| `film.step(title, subtitle, seconds, body)` | a captioned step; `body` is the closure the step's actions go in (`film.step(title, seconds, body)` without a subtitle) |
| `film.title(title, subtitle, seconds, body)` | a step with a large caption |
| `film.hold(seconds)` | a pause with no caption |
| `s.orbit(degrees)` | turn the camera about its target |
| `s.camera(#{ azimuth, elevation, zoom, target })` | move the camera; unset fields keep their value; `target` an occurrence or `"all"` |
| `s.move(target, #{ from: [x, y, z] })` / `#{ to: … }` | slide from an offset into place, or out to an offset (document units) |
| `s.turn(target, axis_at(p, d), from, to)` | rotate about an axis, degrees |
| `s.joint(name, from, to)` | sweep an assembly parameter that only moves occurrences (a revolute's angle, a mate's distance) |
| `s.explode(from, to)` | the explode factor |
| `s.section(#{ from, to, from_at, to_at, combine })` | a section plane: `"x"`, `"y"`, `"z"` at 0..1 across the model, `#{ normal, offset }`, or `"off"` |
| `s.fade(target, from, to)` / `s.show(target, visible)` / `s.highlight(target)` | opacity, visibility, emphasis |
| `s.evolve(label, #{ from, to })` | a generative design converging: the design of the generative feature `label` from update `from` (default 0) to `to` (default `"end"`), over its design space drawn faint; the time spreads over the designs the solve recorded, which are denser early. Its layer is the target `"design label"`; the intermediate designs are display only, never bodies |
| `doc.instances()`, `doc.joints()`, `doc.value(name)`, `doc.bodies`, `doc.kind` | what the document holds |

A target is an occurrence path (`"pair_1"` moves everything under it,
`"pair_1/high"` one pin), `"body N"` in a part document, `"design label"` for a generative
design's layer, or `"all"`.
A name the document does not have is an `ArgumentError` at its line
that lists the names it does have.

## 15. Cheat-sheet

```text
model(#{ unit: "mm" })            open(path)   save(path)   to_step(path)   to_rhai()
m.profile(PLANE).rect|circle|polygon|path…close().at(cx, cy).sketch(label)
PLANE: xy() xz() yz() xy_at(z) xz_at(y) yz_at(x) on_face(F) on_face_offset(F, d) plane_at(o, n)
END:   blind(d) through_all() up_to_surface(F) offset_from_surface(F, d) mid_plane(h) full_circle()
AXIS:  axis_at([ox,oy,oz], [dx,dy,dz])   axis(E)
body.extrude(END, #{label, mode})        body.revolve(AXIS, END, #{label})
body.union(tool) subtract(tool) intersect(tool)      tool = m.tool(|m| …)
body.fillet(E, r, #{label, continuity, cap, profile, interference, law, corner, corner_cap, corner_cap_strict, setback, setback_ends})
body.fillet_each([[E, r], …], #{…as fillet, no law})
body.chamfer(E, d, #{label, kind})       body.shell(F, t, #{label, direction, join, trim})
body.translate(dx, dy, dz)               body.rotate(AXIS, deg)
body.sweep(path_body, E, #{label, mode, twist})     body.pipe(E, r, #{label, mode})
body.helical_sweep(AXIS, r, pitch, turns, #{label, mode, taper_to, start_angle, left_hand})
m.loft([s0, s1, …], #{label, mode, closed})
m.bolt(SPEC, len, #{label, head, thread_length, left_hand, starts})   m.threaded_rod(SPEC, len, #{…})
m.hinge("butt"|"print_in_place", len, #{knuckles, pin, clearance, gap, angle, screw, holes, only})   m.living_hinge(STYLE, span, len)
m.nut(SPEC, #{label, clearance, thickness, left_hand, starts, corner_chamfer, mouth_chamfer})
body.tapped_hole(AXIS, SPEC, DEPTH, #{label, clearance, left_hand, starts, entry_chamfer})
body.clearance_hole(AXIS, SPEC, DEPTH, #{label, fit, diameter, countersink, counterbore, entry_chamfer})
SPEC: "M6" "M6x0.75" "1/4-20" "#8-32" (mm) or #{ major, pitch } (model units)   DEPTH: blind(d) through_all()
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
assembly(name)   a.part(name, |m| …)   a.part(name, #{defaults}, |m, p| …)   def.variant(#{…})
a.subassembly(name, |s| …)   a.place(def, name?).translate(dx,dy,dz).rotate(AXIS, deg).frame(o, z, x).seat(AXIS).hide()
a.param(name, v, unit)   a.drive(name, [defs])   a.drive(name, def, part_param)   a.set_param(name, v)
a.hinge(name, KIND, len, #{angle, …})   a.limits(name, min, max)   a.swing(name, #{from, to, step})
a.assert_watertight()   a.save(path)   a.to_step(path)
```

Gallery scripts that exercise all of this: `crates/popkern-rhai/examples/`
(`boss_plate`, `revolved_knob`, `shelled_enclosure`,
`chamfered_bracket`, `bored_block`, `tapered_fillet`,
`parametric_plate`, `keyed_shaft`); the assemblies are
`ball_bearing`, `planetary_gearbox` and `fastener_kit`.
