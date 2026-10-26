# popkern live worker protocol — version 1

The page never calls the wasm module directly. It posts messages to a
Web Worker (`worker.js`) that owns one `LiveSession`, and the
worker answers with the engine's envelope plus the body buffers as
transferable typed arrays. The same shapes are the seed of the daemon's
WebSocket schema (design §6.4, phase 6): a GUI that speaks this to a
worker will speak it to a socket.

Every request carries a client-chosen `id`; every reply echoes it. The
worker answers requests one at a time in order; a request that arrives
while the engine is busy waits (the wasm call is synchronous). A page
that wants to abandon a run terminates the worker (`WorkerHost.cancel`)
and starts a fresh one — the engine has no cooperative cancel beyond
its budgets, and a hostile script must not be able to wedge the page.

## Requests

| type | fields | reply |
|---|---|---|
| `init` | `options: { max_operations: number \| null, timeout_ms: number \| null, continuity_certification?: "Inline" \| "Deferred" }` | `{ ok, protocol_version, clock_available }` |
| `set_budget` | `options: { max_operations: number \| null, timeout_ms: number \| null }` | `{ ok: true, budget: { max_operations, timeout_ms } }` or `{ ok: false, outcome }` — replaces the budgets for later runs (the page's time limit); the certification mode and the model are kept |
| `run` | `source: string`, `name: string`, `generative?: "solve"` | envelope reply — with `generative: "solve"` this one build may solve generative designs (see below) and runs free of the budgets |
| `seed_generative` | `designs: string` (a bake's `generative`, JSON text) | `{ ok: true, added }` or `{ ok: false, outcome }` — adds solved generative designs to the worker, so builds serve them without solving (PK-956) |
| `open_step` | `bytes: ArrayBuffer`, `name: string` | envelope reply — a STEP file's product structure as a view-only assembly (assemblies.md §8.3): imported parts with no recipe, re-exportable and savable; `printed` carries the import summary and anything left out |
| `set_param` | `name: string`, `value: number` | envelope reply |
| `set_joints` | `values: { [name]: number }` | envelope reply — several placement-only assembly parameters (joints) set with one mate solve and one placement envelope (a film frame's pose); when any name is not a joint, each is set in turn as by `set_param` and the reply is the last envelope or the first failure |
| `certify` | `slice_ms: number` | `{ ok: true, step: { certified, remaining, report: ReplayReport \| null, failure: outcome \| null } }` |
| `to_rhai` | — | `{ ok: true, code: string }` or `{ ok: false, outcome }` |
| `export_step` | `name: string` | `{ ok: true, bytes: Uint8Array }` or `{ ok: false, outcome }` |
| `export_pkn` | — | `{ ok: true, bytes: Uint8Array }` or `{ ok: false, outcome }` |
| `section` | `normal: [x, y, z]`, `offset: number` (a clipping plane, world SI: the side `normal · x > offset` removed) | `{ ok: true, section: { svg, dxf, bodies: [{ body, area, centroid }], frame } }` or `{ ok: false, outcome }` — the built document's EXACT section by the plane (PK-856): each body's loops cut by the kernel, drawn as SVG (hatched, holes even-odd) and DXF (closed polylines, a layer per body), measured in document units in the plane's 2-D frame `frame = [x_dir, y_dir]` |
| `interference` | `clearance: number` | `{ ok: true, report: { pairs: [{ a, b, status: "clear" \| "touching" \| "interfering", distance?, lower_bound?, refined?, volume?, decline?, closest?, faces? }], occurrences, candidate_pairs, tolerance, clearance } }` or `{ ok: false, outcome }` — every pair of assembly part occurrences within `clearance` (assemblies.md §8.1): certified clearances, contact, overlap volumes (null with `decline` when the boolean declined) |
| `set_visible` | `path: string`, `visible: boolean` | `{ ok: true }` or `{ ok: false, outcome }` — show or hide an assembly occurrence (`stage/planet_1`) in the document (saved with `.pkn`); no geometry changes, nothing is re-meshed |
| `analysis_mesh` | — | `{ ok: true, stats: { facet_deg, triangles, mesh_ms }, buffers }` or `{ ok: false, outcome }` — the live bodies re-meshed at the zebra analysis bound (6°), one `{ positions, normals, indices, face_ids }` per body in envelope order; cached until the next `run`/`set_param` |
| `simulate` | `setup: string` (an analysis setup file: structural, or `analysis = "thermal"`; docs/architecture/analysis-integration.md) | `{ ok: true, simulation: { status: "solved" \| "declined", analysis, mm_per_unit, model?, result?, solve?, error?: { kind, message } }, meshes: [{ body, field_names, positions, indices, displacement, fields: { [name]: Float32Array } }] }` or `{ ok: false, outcome }` — solves the setup on the part document (PK-948): faces are referenced by their stable names (`face_labels`), so a setup survives a parametric edit; one field mesh per solid, positions in model units, fields per vertex in mm, MPa or °C (`von_mises_mpa`, `max_principal_mpa`, `displacement_mm`, `zz_error_mpa`, or `temperature_c`); a setup or solver decline is `status: "declined"` with the typed kind |
| `materials` | — | `{ ok: true, materials: [{ name, youngs_modulus_mpa, printed }] }` — the presets a `simulate` setup may name |
| `run_project` | `files: { [path]: string }`, `base64?: { [path]: string }` | envelope reply with `project: { name, versions: [{ name, status: "built" \| "cached" \| "failed" \| "stale", elapsed_ms, features, edited, printed }], error }` — builds a project folder (live-projects-and-history.md §5.4) held in the page: text files by path, binary files (a STEP base) base64; checkpoints persist across calls, so only migrations from the first changed one rerun |
| `film` | `source: string`, `name: string` (`films/guide.rhai`), `context?` (the envelope's `film_context`, used when the worker holds no live document: a pre-built example) | `{ ok: true, film, printed, context }` or `{ ok: false, outcome }` — runs a film script against the built document without rebuilding it (docs/architecture/presentation-films.md): `film` is the storyboard (`{ name, size, fps, look, background, duration, camera, steps, tracks }`, lengths SI), `context` what it could name (`{ kind, name, unit_scale, paths, bodies, joints, params }`); an unknown name is an `ArgumentError` at the film's line |
| `timeline` | — | `{ ok: true, parts: [{ part, steps: [{ index, feature, label, operation, migration, suppressed, failed, bodies_changed, source: { script, line, column } \| null }] }], assembly: { steps: [{ index, kind: "place" \| "move" \| "mate" \| "solve", label, instances, migration, source }], recorded } \| null }` or `{ ok: false, outcome }` — the document's history, one row per timeline feature of every part (`part` is null for a part document); for an assembly, `assembly` also lists its construction steps (`recorded` is false for one opened from a file: its instances at their stored placements, no mate steps) |
| `assembly_step` | `index: number` | `{ ok: true, instances: [...] }` or `{ ok: false, outcome }` — the occurrences after the first `index` construction steps, in the run's `assembly.instances` shape (an occurrence not yet placed is absent); placements only, nothing is meshed |
| `step_meshes` | `part: string` (`""` for a part document), `index: number`, `known?: string[]` | `{ ok: true, rows: [{ key, body, created_faces, buffers, mesh, positions?, normals?, indices?, face_ids? }] }` or `{ ok: false, outcome }` — the bodies alive after the first `index` features; `key` (`body:feature`, the step that last changed the body) names identical meshes, so rows whose key is in `known` carry no buffers; `created_faces` indexes `mesh.face_labels` for the faces this step created |
| `stats` | — | `{ ok: true, wasm_memory_bytes: number \| null, body_count }` — the module's linear memory size (grows, never shrinks) |

`init` must be the first message; it loads the wasm module, registers
the host clock (`performance.now`) and creates the session under the
budgets. `clock_available` is `false` only where the worker has no
`performance` object — then the wall-clock budget cannot fire and the
envelope's `wall_budget` says `unavailable`.

## Generative designs (PK-956)

A generative design takes minutes to solve in the browser, so the
worker never solves one on its own: a `run` or `set_param` that meets a
design it has not solved declines it (`kernel_decline`, `FeatureFailed`
with `payload.error.kind = "GenerativeDeclined"` and
`payload.error.generative_kind = "NeedsSolve"`), and the page keeps the
model it shows. The page's re-solve is a `run` with `generative:
"solve"`. While it solves, the worker posts messages that are not
replies, one per optimisation update:

    { type: "progress", id, progress: { count, max, compliance_nmm, volume_fraction, change } }

`id` is the request's. A design once solved is kept by the worker, so
later builds whose design inputs are unchanged serve it without
solving. A gallery bake carries the designs its native run solved
(`generative` in its header); the page hands them to the worker with
`seed_generative`, so building the baked example live (a slider on
another feature, Simulate, an export) solves nothing. A native session
(`popkern serve`) solves as it builds.

A solve also keeps its design's history (PK-959): the density field at
the first dozen updates, then about every 15%, and the final design. A
film showing the design converge (`s.evolve`) reads it by the design's
id (its feature report's `generative.design`, also in the film
context's `designs`):

| type | fields | reply |
| --- | --- | --- |
| `seed_design_history` | `bytes: ArrayBuffer` (a bake's `<stem>.designs.bin`, listed as `designs` in the bake index) | `{ ok: true, added }` or `{ ok: false, outcome }` |
| `design_history` | `design: string` | `{ ok: true, history: { design, updates, cell, frames: [{ update, compliance_nmm, volume_fraction }] } }`, or `{ ok: false, outcome }` when the worker holds no history for it |
| `design_mesh` | `design: string`, `at: number` (a position among the records: 0 the first, `frames.length - 1` the final design, between them two records blended), `space: boolean` (the design space instead) | `{ ok: true, mesh: { positions, normals, indices } }` in the design body's units, or `{ ok: false, outcome }` |

A gallery bake writes its designs' histories beside it; the page
fetches the file only when a film needs it, and hands it to every
worker it starts later.

## Deferred continuity certification

Certifying a fillet's seams (interval bounds over each whole seam) can
cost far more than building the fillet. With `continuity_certification:
"Deferred"` in `init`, runs and parameter edits serve on the sampled
continuity evidence: each blend row's `blend.continuity` carries
`deferred_seams > 0` and `certified_achieved: null` (no certified claim
yet). The page then sends `certify` requests; each certifies seams for
about `slice_ms` (at least one) and answers how many are `remaining`.
When `remaining` reaches 0 the step carries the upgraded `report` (same
shape as the envelope's); certification never changes geometry, so no
buffers are sent. A newer `run`/`set_param` replaces the report and its
deferred seams, so a page abandons a certify loop simply by not sending
the next slice. `failure` is set (and `remaining` is 0) when
certification could not continue. Without the option the default is
`"Inline"`: every run certifies before it answers, as before.

## Envelope reply

```
{
  id,
  ok: true,
  envelope: {
    outcome: "served" | "script_error" | "kernel_decline" | "gate_failure",
    kind?, rendered?, payload?,        // script_error / kernel_decline
    step?, detail?,                    // gate_failure
    printed: string[],
    elapsed_ms: number | null,
    clock_available: boolean,
    wall_budget: "none" | "armed" | "unavailable",
    report: ReplayReport | null,       // spec §11 + is_clean/failures/downgraded_blends
    features: number,
    params: [{ name, value, unit, expression, bounds: [lo, hi] | null, read_only, explicit }],   // explicit: declared by name, not a feature argument
    bodies: [{ id, label, vertex_count, triangle_count, face_count,
               face_labels: string[], bounds: [[x,y,z],[x,y,z]], open_edges }],
    bodies_unchanged: boolean,         // true: no bodies sent, keep the last ones (below)
    assembly: null | {                 // an assembly document (assemblies.md §6.5)
      name,
      definitions: [{ name, kind: "part" | "assembly", bodies: number[], features }],
      instances: [{ path, name, depth, definition, world: number[16] /* column-major */,
                    visible, is_part }]
    }
  },
  // For an assembly `bodies` are the part definitions' meshes, each
  // drawn once per occurrence under its `world` matrix; an assembly's
  // `params` are qualified "part:name", and `set_param` takes the same.
  buffers: [{ positions: Float32Array, normals: Float32Array,
              indices: Uint32Array, face_ids: Uint32Array,
              gaussian: Float32Array, mean: Float32Array }]    // one per body, transferred
}
```

`buffers[i]` belongs to `envelope.bodies[i]`. `face_ids[t]` indexes
`face_labels` for triangle `t`. A failed run (`ok: true`, outcome not
`served`) carries no bodies: the page keeps the last good model on
screen and shows the outcome beside it. A `set_param` of an assembly
parameter that drives no part parameter (a joint) only moves
occurrences: its envelope has `bodies_unchanged: true`, empty `bodies`
and `buffers`, `report: null`, only the assembly's own `params` rows,
and the last served definitions; the page keeps its bodies, report and
part rows and re-places the meshes it holds at the new `instances`
matrices (the analysis and step meshes stay valid). The classes and payload fields
are the scripting-surface spec §10; the partition
(`script_error`/`kernel_decline`/`gate_failure`) is the one `popkern
run` maps to exit codes 2/3/4.

## Errors

A malformed request, a message before `init`, or an exception inside
the worker answers `{ id, ok: false, outcome: { outcome:
"gate_failure", step: "worker", detail } }` — the same shape as an
engine gate failure, so the page has one error path.

## Socket transport (popkern serve, PK-665)

Under `popkern serve` with native compute, `host.json` says
`{ "compute": "native" }` and the page speaks this protocol to
`./api/session` on its own origin, through `socket-worker.js`, which is
shaped like a Worker.
- Requests are the same JSON with two differences: `open_step` sends
  `bytes_b64` (base64) instead of an ArrayBuffer, and there is one more
  request, `{ "type": "cancel" }`. Cancel has no reply. It stops the
  running request's script at its next step: that request answers with
  the `Cancelled` class, and the session keeps its document.
- Each reply is a text frame (the reply JSON plus `attachments: [{ path,
  kind: "f32" | "u32" | "u8" }]`), followed by one binary frame per
  attachment, in order. The page puts each typed array at its dotted
  path (`buffers.0.positions`, `rows.2.indices`, `bytes`).

## Versioning

`protocol_version` is the engine's `PROTOCOL_VERSION`. The page refuses
to speak to a worker whose version it does not know. Bumped when a
request or reply shape changes; additive fields and new requests
(`run_project`, `timeline`, `step_meshes`, PK-660/661; `assembly_step`
and the `timeline` reply's `assembly`, PK-666) do not bump.
