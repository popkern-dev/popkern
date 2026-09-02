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

  async function init(message) {
    if (!loaded) {
      await loadModule();
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
      };
      buffers.push(body);
      transfer.push(body.positions.buffer, body.normals.buffer, body.indices.buffer, body.face_ids.buffer);
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
        case "set_param":
          return envelopeReply(id, playground.set_param(String(message.name), Number(message.value)));
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
        default:
          return { reply: gateFailure(id, "worker", `unknown request type ${String(message?.type)}`), transfer: [] };
      }
    } catch (raised) {
      return { reply: gateFailure(id, "worker", String(raised?.message ?? raised)), transfer: [] };
    }
  };
}
