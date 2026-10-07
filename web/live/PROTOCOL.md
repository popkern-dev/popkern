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
| `run` | `source: string`, `name: string` | envelope reply |
| `open_step` | `bytes: ArrayBuffer`, `name: string` | envelope reply — a STEP file's product structure as a view-only assembly (assemblies.md §8.3): imported parts with no recipe, re-exportable and savable; `printed` carries the import summary and anything left out |
| `set_param` | `name: string`, `value: number` | envelope reply |
| `certify` | `slice_ms: number` | `{ ok: true, step: { certified, remaining, report: ReplayReport \| null, failure: outcome \| null } }` |
| `to_rhai` | — | `{ ok: true, code: string }` or `{ ok: false, outcome }` |
| `export_step` | `name: string` | `{ ok: true, bytes: Uint8Array }` or `{ ok: false, outcome }` |
| `export_pkn` | — | `{ ok: true, bytes: Uint8Array }` or `{ ok: false, outcome }` |
| `interference` | `clearance: number` | `{ ok: true, report: { pairs: [{ a, b, status: "clear" \| "touching" \| "interfering", distance?, lower_bound?, refined?, volume?, decline?, closest?, faces? }], occurrences, candidate_pairs, tolerance, clearance } }` or `{ ok: false, outcome }` — every pair of assembly part occurrences within `clearance` (assemblies.md §8.1): certified clearances, contact, overlap volumes (null with `decline` when the boolean declined) |
| `set_visible` | `path: string`, `visible: boolean` | `{ ok: true }` or `{ ok: false, outcome }` — show or hide an assembly occurrence (`stage/planet_1`) in the document (saved with `.pkn`); no geometry changes, nothing is re-meshed |
| `analysis_mesh` | — | `{ ok: true, stats: { facet_deg, triangles, mesh_ms }, buffers }` or `{ ok: false, outcome }` — the live bodies re-meshed at the zebra analysis bound (6°), one `{ positions, normals, indices, face_ids }` per body in envelope order; cached until the next `run`/`set_param` |
| `run_project` | `files: { [path]: string }`, `base64?: { [path]: string }` | envelope reply with `project: { name, versions: [{ name, status: "built" \| "cached" \| "failed" \| "stale", elapsed_ms, features, edited, printed }], error }` — builds a project folder (live-projects-and-history.md §5.4) held in the page: text files by path, binary files (a STEP base) base64; checkpoints persist across calls, so only migrations from the first changed one rerun |
| `timeline` | — | `{ ok: true, parts: [{ part, steps: [{ index, feature, label, operation, migration, suppressed, failed, bodies_changed, source: { script, line, column } \| null }] }], assembly: { steps: [{ index, kind: "place" \| "move" \| "mate" \| "solve", label, instances, migration, source }], recorded } \| null }` or `{ ok: false, outcome }` — the document's history, one row per timeline feature of every part (`part` is null for a part document); for an assembly, `assembly` also lists its construction steps (`recorded` is false for one opened from a file: its instances at their stored placements, no mate steps) |
| `assembly_step` | `index: number` | `{ ok: true, instances: [...] }` or `{ ok: false, outcome }` — the occurrences after the first `index` construction steps, in the run's `assembly.instances` shape (an occurrence not yet placed is absent); placements only, nothing is meshed |
| `step_meshes` | `part: string` (`""` for a part document), `index: number`, `known?: string[]` | `{ ok: true, rows: [{ key, body, created_faces, buffers, mesh, positions?, normals?, indices?, face_ids? }] }` or `{ ok: false, outcome }` — the bodies alive after the first `index` features; `key` (`body:feature`, the step that last changed the body) names identical meshes, so rows whose key is in `known` carry no buffers; `created_faces` indexes `mesh.face_labels` for the faces this step created |
| `stats` | — | `{ ok: true, wasm_memory_bytes: number \| null, body_count }` — the module's linear memory size (grows, never shrinks) |

`init` must be the first message; it loads the wasm module, registers
the host clock (`performance.now`) and creates the session under the
budgets. `clock_available` is `false` only where the worker has no
`performance` object — then the wall-clock budget cannot fire and the
envelope's `wall_budget` says `unavailable`.

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
screen and shows the outcome beside it. The classes and payload fields
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
