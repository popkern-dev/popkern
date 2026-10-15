// The worker's request handler, independent of the message transport
// so the node harness can drive it exactly as the browser worker does.
// See PROTOCOL.md for the shapes.

export const KNOWN_PROTOCOLS = ["1"];

function gateFailure(id, step, detail) {
  return { id, ok: false, outcome: { outcome: "gate_failure", step, detail } };
}

/**
 * Build a handler over a loaded wasm package.
 * `pkg` is the wasm-bindgen module namespace (`init`, `LiveSession`,
 * `protocol_version`); `loadModule` performs the one-time wasm
 * instantiation (browser: `pkg.default()`, node: `initSync`).
 */
export function createHandler(pkg, loadModule) {
  let session = null;
  let loaded = false;
  let wasmMemory = null;

  async function init(message) {
    if (!loaded) {
      wasmMemory = (await loadModule()) ?? null;
      loaded = true;
    }
    const version = pkg.protocol_version();
    if (!KNOWN_PROTOCOLS.includes(version)) {
      return gateFailure(message.id, "worker", `unknown protocol version ${version}`);
    }
    const clock = pkg.init();
    const options = message.options ?? {};
    session = new pkg.LiveSession(JSON.stringify(options));
    return { id: message.id, ok: true, protocol_version: version, clock_available: clock };
  }

  function envelopeReply(id, json) {
    const envelope = JSON.parse(json);
    const buffers = [];
    const transfer = [];
    // A joint moved placements only: the page keeps the meshes it has.
    const count = envelope.bodies_unchanged ? 0 : session.body_count();
    for (let i = 0; i < count; i++) {
      const body = {
        positions: session.body_positions(i),
        normals: session.body_normals(i),
        indices: session.body_indices(i),
        face_ids: session.body_face_ids(i),
        gaussian: session.body_curvature(i, "gaussian"),
        mean: session.body_curvature(i, "mean"),
      };
      buffers.push(body);
      transfer.push(
        body.positions.buffer,
        body.normals.buffer,
        body.indices.buffer,
        body.face_ids.buffer,
        body.gaussian.buffer,
        body.mean.buffer,
      );
    }
    return { reply: { id, ok: true, envelope, buffers }, transfer };
  }

  function bytesReply(id, produce) {
    try {
      const bytes = produce();
      return { reply: { id, ok: true, bytes }, transfer: [bytes.buffer] };
    } catch (raised) {
      return { reply: { id, ok: false, outcome: parseOutcome(raised) }, transfer: [] };
    }
  }

  function parseOutcome(raised) {
    if (typeof raised === "string") {
      try {
        return JSON.parse(raised);
      } catch {
        return { outcome: "gate_failure", step: "worker", detail: raised };
      }
    }
    return { outcome: "gate_failure", step: "worker", detail: String(raised?.message ?? raised) };
  }

  /** Handle one request. Resolves to `{ reply, transfer }`. */
  return async function handle(message) {
    const id = message?.id ?? null;
    try {
      if (message?.type === "init") {
        return { reply: await init(message), transfer: [] };
      }
      if (!session) {
        return { reply: gateFailure(id, "worker", "send init before any other request"), transfer: [] };
      }
      switch (message.type) {
        case "run":
          return envelopeReply(id, session.run(String(message.source ?? ""), String(message.name ?? "script.rhai")));
        case "set_budget":
          try {
            const budget = JSON.parse(session.set_budget(JSON.stringify(message.options ?? {})));
            return { reply: { id, ok: true, budget }, transfer: [] };
          } catch (raised) {
            return { reply: { id, ok: false, outcome: parseOutcome(raised) }, transfer: [] };
          }
        case "open_step":
          return envelopeReply(id, session.open_step(new Uint8Array(message.bytes ?? []), String(message.name ?? "model")));
        case "set_param":
          return envelopeReply(id, session.set_param(String(message.name), Number(message.value)));
        case "certify":
          return { reply: { id, ok: true, step: JSON.parse(session.certify(Number(message.slice_ms ?? 0))) }, transfer: [] };
        case "to_rhai":
          try {
            return { reply: { id, ok: true, code: session.to_rhai() }, transfer: [] };
          } catch (raised) {
            return { reply: { id, ok: false, outcome: parseOutcome(raised) }, transfer: [] };
          }
        case "export_step":
          return bytesReply(id, () => session.export_step(String(message.name ?? "model")));
        case "export_pkn":
          return bytesReply(id, () => session.export_pkn());
        case "interference":
          try {
            return { reply: { id, ok: true, report: JSON.parse(session.interference(Number(message.clearance ?? 0))) }, transfer: [] };
          } catch (raised) {
            return { reply: { id, ok: false, outcome: parseOutcome(raised) }, transfer: [] };
          }
        case "set_visible":
          try {
            session.set_visible(String(message.path ?? ""), Boolean(message.visible));
            return { reply: { id, ok: true }, transfer: [] };
          } catch (raised) {
            return { reply: { id, ok: false, outcome: parseOutcome(raised) }, transfer: [] };
          }
        case "analysis_mesh": {
          // Zebra analysis mesh (PK-577): the live bodies at the 6°
          // facet bound, same order as the last envelope's bodies.
          const stats = JSON.parse(session.analysis_mesh());
          if (!stats.ok) {
            return { reply: { id, ok: false, outcome: stats.outcome }, transfer: [] };
          }
          const buffers = [];
          const transfer = [];
          for (let i = 0; i < session.body_count(); i++) {
            const body = {
              positions: session.analysis_f32(i, "positions"),
              normals: session.analysis_f32(i, "normals"),
              indices: session.analysis_u32(i, "indices"),
              face_ids: session.analysis_u32(i, "face_ids"),
            };
            buffers.push(body);
            transfer.push(body.positions.buffer, body.normals.buffer, body.indices.buffer, body.face_ids.buffer);
          }
          delete stats.ok;
          return { reply: { id, ok: true, stats, buffers }, transfer };
        }
        case "run_project":
          // A project (live-projects-and-history.md): text files by
          // path in `files`, binary ones base64 in `base64`.
          return envelopeReply(id, session.run_project(JSON.stringify({ files: message.files ?? {}, base64: message.base64 ?? {} })));
        case "film": {
          // A film script against the built document (no rebuild).
          const reply = JSON.parse(session.film(String(message.source ?? ""), String(message.name ?? "film.rhai"), message.context ? JSON.stringify(message.context) : ""));
          if (reply.outcome) return { reply: { id, ok: false, outcome: reply.outcome }, transfer: [] };
          return { reply: { id, ok: true, ...reply }, transfer: [] };
        }
        case "timeline": {
          const parts = JSON.parse(session.timeline());
          if (!Array.isArray(parts)) return { reply: { id, ok: false, outcome: parts.outcome }, transfer: [] };
          const assembly = JSON.parse(session.assembly_timeline());
          return { reply: { id, ok: true, parts, assembly }, transfer: [] };
        }
        case "assembly_step": {
          // The occurrences after `index` construction steps (placements
          // only: the page holds every definition's meshes).
          const instances = JSON.parse(session.assembly_step(Number(message.index ?? 0)));
          if (!Array.isArray(instances)) return { reply: { id, ok: false, outcome: instances.outcome }, transfer: [] };
          return { reply: { id, ok: true, instances }, transfer: [] };
        }
        case "step_meshes": {
          // The bodies after `index` features of `part`; rows whose key
          // the page lists in `known` carry no buffers.
          const rows = JSON.parse(session.step_meshes(String(message.part ?? ""), Number(message.index ?? 0), JSON.stringify(message.known ?? [])));
          if (!Array.isArray(rows)) return { reply: { id, ok: false, outcome: rows.outcome }, transfer: [] };
          const transfer = [];
          rows.forEach((row, i) => {
            if (!row.buffers) return;
            row.positions = session.step_f32(i, "positions");
            row.normals = session.step_f32(i, "normals");
            row.indices = session.step_u32(i, "indices");
            row.face_ids = session.step_u32(i, "face_ids");
            transfer.push(row.positions.buffer, row.normals.buffer, row.indices.buffer, row.face_ids.buffer);
          });
          return { reply: { id, ok: true, rows }, transfer };
        }
        case "stats":
          return {
            reply: {
              id,
              ok: true,
              wasm_memory_bytes: wasmMemory?.buffer?.byteLength ?? null,
              body_count: session.body_count(),
            },
            transfer: [],
          };
        default:
          return { reply: gateFailure(id, "worker", `unknown request type ${String(message?.type)}`), transfer: [] };
      }
    } catch (raised) {
      return { reply: gateFailure(id, "worker", String(raised?.message ?? raised)), transfer: [] };
    }
  };
}
