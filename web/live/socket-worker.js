// The worker protocol over a socket (PK-665): `popkern serve` with
// native compute answers the page's requests from a native session at
// `./api/session` on the page's own origin.

/** The worker protocol over `./api/session` (a native session behind
 *  `popkern serve`), shaped like a Worker for WorkerHost: postMessage,
 *  message/error events, terminate — plus `cancelRun()`, which stops
 *  the running script and keeps the session (a Worker can only be
 *  terminated). Typed arrays travel as binary frames after each reply. */
export class SocketWorker {
  constructor() {
    const url = new URL("./api/session", location.href);
    url.protocol = url.protocol === "https:" ? "wss:" : "ws:";
    this.socket = new WebSocket(url);
    this.socket.binaryType = "arraybuffer";
    this.listeners = { message: [], error: [], messageerror: [] };
    this.queue = [];
    this.pending = null;
    this.socket.addEventListener("open", () => {
      for (const text of this.queue) this.socket.send(text);
      this.queue = [];
    });
    this.socket.addEventListener("message", (event) => this.frame(event.data));
    this.socket.addEventListener("error", () => this.emit("error", { message: "the session socket failed" }));
    this.socket.addEventListener("close", () => {
      if (!this.closing) this.emit("error", { message: "the session socket closed" });
    });
  }

  addEventListener(type, listener) {
    this.listeners[type]?.push(listener);
  }

  emit(type, event) {
    for (const listener of this.listeners[type] ?? []) listener(event);
  }

  postMessage(message) {
    const out = { ...message };
    if (out.bytes) {
      const bytes = new Uint8Array(out.bytes);
      let text = "";
      for (let i = 0; i < bytes.length; i += 0x8000) text += String.fromCharCode(...bytes.subarray(i, i + 0x8000));
      out.bytes_b64 = btoa(text);
      delete out.bytes;
    }
    const text = JSON.stringify(out);
    if (this.socket.readyState === WebSocket.OPEN) this.socket.send(text);
    else this.queue.push(text);
  }

  /** Stop the running script; the session and its document stay. */
  cancelRun() {
    this.postMessage({ type: "cancel" });
  }

  frame(data) {
    if (typeof data === "string") {
      const reply = JSON.parse(data);
      const attachments = reply.attachments ?? [];
      delete reply.attachments;
      if (attachments.length === 0) return this.emit("message", { data: reply });
      this.pending = { reply, attachments, index: 0 };
      return;
    }
    const p = this.pending;
    if (!p) return;
    const { path, kind } = p.attachments[p.index++];
    const array = kind === "f32" ? new Float32Array(data) : kind === "u32" ? new Uint32Array(data) : new Uint8Array(data);
    const keys = path.split(".");
    let at = p.reply;
    for (const key of keys.slice(0, -1)) at = at[key] ??= {};
    at[keys.at(-1)] = array;
    if (p.index === p.attachments.length) {
      this.pending = null;
      this.emit("message", { data: p.reply });
    }
  }

  terminate() {
    this.closing = true;
    this.socket.close();
  }
}
