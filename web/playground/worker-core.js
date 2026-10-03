// The worker's request handler, independent of the message transport
// so the node harness can drive it exactly as the browser worker does.
// See PROTOCOL.md for the shapes.

export const KNOWN_PROTOCOLS = ["1"];

function gateFailure(id, step, detail) {
  return { id, ok: false, outcome: { outcome: "gate_failure", step, detail } };
}

/**
 * Build a handler over a loaded wasm package.
 * `pkg` is the wasm-bindgen module namespace (`init`, `Playground`,
 * `protocol_version`); `loadModule` performs the one-time wasm
 * instantiation (browser: `pkg.default()`, node: `initSync`).
 */
export function createHandler(pkg, loadModule) {
  let playground = null;
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
    playground = new pkg.Playground(JSON.stringify(options));
    return { id: message.id, ok: true, protocol_version: version, clock_available: clock };
  }

  function envelopeReply(id, json) {
    const envelope = JSON.parse(json);
    const buffers = [];
    const transfer = [];
    for (let i = 0; i < playground.body_count(); i++) {
      const body = {
        positions: playground.body_positions(i),
        normals: playground.body_normals(i),
        indices: playground.body_indices(i),
        face_ids: playground.body_face_ids(i),
        gaussian: playground.body_curvature(i, "gaussian"),
        mean: playground.body_curvature(i, "mean"),
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
      if (!playground) {
        return { reply: gateFailure(id, "worker", "send init before any other request"), transfer: [] };
      }
      switch (message.type) {
        case "run":
          return envelopeReply(id, playground.run(String(message.source ?? ""), String(message.name ?? "script.rhai")));
        case "set_budget":
          try {
            const budget = JSON.parse(playground.set_budget(JSON.stringify(message.options ?? {})));
            return { reply: { id, ok: true, budget }, transfer: [] };
          } catch (raised) {
            return { reply: { id, ok: false, outcome: parseOutcome(raised) }, transfer: [] };
          }
        case "set_param":
          return envelopeReply(id, playground.set_param(String(message.name), Number(message.value)));
        case "certify":
          return { reply: { id, ok: true, step: JSON.parse(playground.certify(Number(message.slice_ms ?? 0))) }, transfer: [] };
        case "to_rhai":
          try {
            return { reply: { id, ok: true, code: playground.to_rhai() }, transfer: [] };
          } catch (raised) {
            return { reply: { id, ok: false, outcome: parseOutcome(raised) }, transfer: [] };
          }
        case "export_step":
          return bytesReply(id, () => playground.export_step(String(message.name ?? "model")));
        case "export_pkn":
          return bytesReply(id, () => playground.export_pkn());
        case "analysis_mesh": {
          // Zebra analysis mesh (PK-577): the live bodies at the 6°
          // facet bound, same order as the last envelope's bodies.
          const stats = JSON.parse(playground.analysis_mesh());
          if (!stats.ok) {
            return { reply: { id, ok: false, outcome: stats.outcome }, transfer: [] };
          }
          const buffers = [];
          const transfer = [];
          for (let i = 0; i < playground.body_count(); i++) {
            const body = {
              positions: playground.analysis_f32(i, "positions"),
              normals: playground.analysis_f32(i, "normals"),
              indices: playground.analysis_u32(i, "indices"),
              face_ids: playground.analysis_u32(i, "face_ids"),
            };
            buffers.push(body);
            transfer.push(body.positions.buffer, body.normals.buffer, body.indices.buffer, body.face_ids.buffer);
          }
          delete stats.ok;
          return { reply: { id, ok: true, stats, buffers }, transfer };
        }
        case "stats":
          return {
            reply: {
              id,
              ok: true,
              wasm_memory_bytes: wasmMemory?.buffer?.byteLength ?? null,
              body_count: playground.body_count(),
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
