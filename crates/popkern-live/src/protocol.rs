//! The worker protocol (`web/live/PROTOCOL.md`) answered natively: the
//! same requests and reply shapes the browser worker
//! (`web/live/worker-core.js` over the wasm bindings) produces, for a
//! host that runs the kernel itself (`popkern serve`, PK-665). Bulk
//! buffers travel beside the JSON as [`Attachment`]s named by their path
//! in the reply (`buffers.0.positions`), so a transport can send them as
//! binary frames and the page rebuilds the typed arrays in place.

use serde_json::{Value, json};

use popkern_core::time::monotonic_clock_available;

use crate::envelope::{BodyMesh, Envelope};
use crate::session::{RunOptions, Session};
use crate::session_project::{base64_decode, decode_project_files};

/// One bulk buffer of a reply.
#[derive(Debug, Clone, PartialEq)]
pub enum Buffer {
    /// A `Float32Array`.
    F32(Vec<f32>),
    /// A `Uint32Array`.
    U32(Vec<u32>),
    /// A `Uint8Array`.
    U8(Vec<u8>),
}

impl Buffer {
    /// The buffer's bytes (little-endian, as typed arrays hold them).
    pub fn bytes(&self) -> Vec<u8> {
        match self {
            Buffer::F32(v) => v.iter().flat_map(|x| x.to_le_bytes()).collect(),
            Buffer::U32(v) => v.iter().flat_map(|x| x.to_le_bytes()).collect(),
            Buffer::U8(v) => v.clone(),
        }
    }

    /// `"f32"`, `"u32"` or `"u8"`.
    pub fn kind(&self) -> &'static str {
        match self {
            Buffer::F32(_) => "f32",
            Buffer::U32(_) => "u32",
            Buffer::U8(_) => "u8",
        }
    }
}

/// A buffer and where it belongs in the reply (`buffers.0.positions`).
#[derive(Debug, Clone, PartialEq)]
pub struct Attachment {
    /// Dotted path into the reply JSON.
    pub path: String,
    /// The buffer.
    pub buffer: Buffer,
}

/// One answer: the JSON reply and its buffers.
#[derive(Debug, Clone)]
pub struct Reply {
    /// The reply (`id`, `ok`, …) without the buffers.
    pub json: Value,
    /// The buffers, by path.
    pub attachments: Vec<Attachment>,
}

impl Reply {
    fn plain(json: Value) -> Reply {
        Reply {
            json,
            attachments: Vec::new(),
        }
    }
}

fn gate_failure(id: &Value, step: &str, detail: &str) -> Reply {
    Reply::plain(
        json!({ "id": id, "ok": false, "outcome": { "outcome": "gate_failure", "step": step, "detail": detail } }),
    )
}

fn body_attachments(prefix: &str, bodies: &[BodyMesh], curvature: bool) -> Vec<Attachment> {
    let mut out = Vec::new();
    for (i, b) in bodies.iter().enumerate() {
        let at = |field: &str, buffer| Attachment {
            path: format!("{prefix}.{i}.{field}"),
            buffer,
        };
        out.push(at("positions", Buffer::F32(b.positions.clone())));
        out.push(at("normals", Buffer::F32(b.normals.clone())));
        out.push(at("indices", Buffer::U32(b.indices.clone())));
        out.push(at("face_ids", Buffer::U32(b.face_ids.clone())));
        if curvature {
            out.push(at("gaussian", Buffer::F32(b.gaussian.clone())));
            out.push(at("mean", Buffer::F32(b.mean.clone())));
        }
    }
    out
}

/// A session answering protocol requests.
pub struct Protocol {
    session: Option<Session>,
    /// Bodies of the last envelope (the `stats` request's count).
    last_bodies: usize,
}

impl Default for Protocol {
    fn default() -> Self {
        Protocol::new()
    }
}

impl Protocol {
    /// A protocol endpoint awaiting `init`.
    pub fn new() -> Protocol {
        Protocol {
            session: None,
            last_bodies: 0,
        }
    }

    /// The session, once `init` created it.
    pub fn session(&mut self) -> Option<&mut Session> {
        self.session.as_mut()
    }

    fn envelope(&mut self, id: &Value, envelope: Envelope) -> Reply {
        // A placement-only envelope keeps the bodies the page holds.
        if !envelope.bodies_unchanged {
            self.last_bodies = envelope.bodies.len();
        }
        let attachments = body_attachments("buffers", &envelope.bodies, true);
        let buffers: Vec<Value> = envelope.bodies.iter().map(|_| json!({})).collect();
        Reply {
            json: json!({ "id": id, "ok": true, "envelope": envelope.to_json(), "buffers": buffers }),
            attachments,
        }
    }

    /// Answer one request.
    pub fn handle(&mut self, message: &Value) -> Reply {
        let id = message.get("id").cloned().unwrap_or(Value::Null);
        let kind = message.get("type").and_then(Value::as_str).unwrap_or("");
        if kind == "init" {
            let options: RunOptions = message
                .get("options")
                .cloned()
                .and_then(|o| serde_json::from_value(o).ok())
                .unwrap_or_default();
            self.session = Some(Session::new(options.host()));
            return Reply::plain(json!({
                "id": id, "ok": true,
                "protocol_version": crate::PROTOCOL_VERSION,
                "clock_available": monotonic_clock_available(),
            }));
        }
        let Some(session) = self.session.as_mut() else {
            return gate_failure(&id, "worker", "send init before any other request");
        };
        let text = |key: &str, default: &str| {
            message
                .get(key)
                .and_then(Value::as_str)
                .unwrap_or(default)
                .to_string()
        };
        let number = |key: &str| message.get(key).and_then(Value::as_f64).unwrap_or(0.0);
        let failed = |outcome: popkern_rhai::Outcome| {
            Reply::plain(json!({ "id": id, "ok": false, "outcome": outcome }))
        };
        match kind {
            "run" => {
                let e = session.run(&text("source", ""), &text("name", "script.rhai"));
                self.envelope(&id, e)
            }
            "open_step" => match base64_decode(&text("bytes_b64", "")) {
                Ok(bytes) => {
                    let e = session.open_step(&bytes, &text("name", "model"));
                    self.envelope(&id, e)
                }
                Err(e) => gate_failure(&id, "worker", &format!("bytes: {e}")),
            },
            "set_param" => {
                let e = session.set_param(&text("name", ""), number("value"));
                self.envelope(&id, e)
            }
            "run_project" => {
                let transfer = json!({
                    "files": message.get("files").cloned().unwrap_or(json!({})),
                    "base64": message.get("base64").cloned().unwrap_or(json!({})),
                });
                match decode_project_files(&transfer.to_string()) {
                    Ok(files) => {
                        let e = session.run_project(files);
                        self.envelope(&id, e)
                    }
                    Err(detail) => gate_failure(&id, "project", &detail),
                }
            }
            "set_budget" => {
                let options: RunOptions = message
                    .get("options")
                    .cloned()
                    .and_then(|o| serde_json::from_value(o).ok())
                    .unwrap_or_default();
                let host = options.host();
                let current = session.options();
                session.set_options(popkern_rhai::HostOptions {
                    max_operations: host.max_operations,
                    timeout: host.timeout,
                    ..current
                });
                Reply::plain(json!({ "id": id, "ok": true, "budget": {
                    "max_operations": options.max_operations, "timeout_ms": options.timeout_ms } }))
            }
            "certify" => {
                let ms = number("slice_ms");
                let slice = if ms.is_finite() && ms > 0.0 {
                    std::time::Duration::from_secs_f64(ms / 1000.0)
                } else {
                    std::time::Duration::ZERO
                };
                Reply::plain(json!({ "id": id, "ok": true, "step": session.certify(slice) }))
            }
            "to_rhai" => match session.to_rhai() {
                Ok(code) => Reply::plain(json!({ "id": id, "ok": true, "code": code })),
                Err(outcome) => failed(outcome),
            },
            "export_step" | "export_pkn" => {
                let bytes = if kind == "export_step" {
                    session.export_step(&text("name", "model"))
                } else {
                    session.export_pkn()
                };
                match bytes {
                    Ok(bytes) => Reply {
                        json: json!({ "id": id, "ok": true }),
                        attachments: vec![Attachment {
                            path: "bytes".into(),
                            buffer: Buffer::U8(bytes),
                        }],
                    },
                    Err(outcome) => failed(outcome),
                }
            }
            "interference" => match session.interference(number("clearance")) {
                Ok(report) => Reply::plain(json!({ "id": id, "ok": true, "report": report })),
                Err(outcome) => failed(outcome),
            },
            "set_visible" => {
                let visible = message
                    .get("visible")
                    .and_then(Value::as_bool)
                    .unwrap_or(true);
                match session.set_visible(&text("path", ""), visible) {
                    Ok(()) => Reply::plain(json!({ "id": id, "ok": true })),
                    Err(outcome) => failed(outcome),
                }
            }
            "analysis_mesh" => {
                let started = popkern_core::time::Stopwatch::start();
                match session.analysis_meshes() {
                    Ok(bodies) => {
                        let triangles: usize = bodies.iter().map(|b| b.indices.len() / 3).sum();
                        let mesh_ms = started.elapsed().map(|d| d.as_secs_f64() * 1000.0);
                        Reply {
                            json: json!({ "id": id, "ok": true,
                                "stats": { "facet_deg": crate::mesh_buffers::ANALYSIS_MAX_FACET_ANGLE_DEG,
                                           "triangles": triangles, "mesh_ms": mesh_ms },
                                "buffers": bodies.iter().map(|_| json!({})).collect::<Vec<_>>() }),
                            attachments: body_attachments("buffers", &bodies, false),
                        }
                    }
                    Err(outcome) => failed(outcome),
                }
            }
            "stats" => Reply::plain(json!({ "id": id, "ok": true, "wasm_memory_bytes": null,
                "body_count": self.last_bodies })),
            "timeline" => match session
                .timeline()
                .and_then(|parts| Ok((parts, session.assembly_timeline()?)))
            {
                Ok((parts, assembly)) => Reply::plain(
                    json!({ "id": id, "ok": true, "parts": parts, "assembly": assembly }),
                ),
                Err(outcome) => failed(outcome),
            },
            "assembly_step" => {
                let index = message.get("index").and_then(Value::as_u64).unwrap_or(0) as usize;
                match session.assembly_step(index) {
                    Ok(instances) => {
                        Reply::plain(json!({ "id": id, "ok": true, "instances": instances }))
                    }
                    Err(outcome) => failed(outcome),
                }
            }
            "step_meshes" => {
                let part = text("part", "");
                let known: Vec<String> = message
                    .get("known")
                    .cloned()
                    .and_then(|k| serde_json::from_value(k).ok())
                    .unwrap_or_default();
                let index = message.get("index").and_then(Value::as_u64).unwrap_or(0) as usize;
                let part = (!part.is_empty()).then_some(part.as_str());
                match session.step_bodies(part, index, &known) {
                    Ok(rows) => {
                        let mut attachments = Vec::new();
                        for (i, row) in rows.iter().enumerate() {
                            if let Some(mesh) = row.mesh.as_ref().filter(|_| row.buffers) {
                                let mut b =
                                    body_attachments("rows", std::slice::from_ref(mesh), false);
                                for a in &mut b {
                                    a.path = a.path.replacen("rows.0.", &format!("rows.{i}."), 1);
                                }
                                attachments.extend(b);
                            }
                        }
                        Reply {
                            json: json!({ "id": id, "ok": true, "rows": rows }),
                            attachments,
                        }
                    }
                    Err(outcome) => failed(outcome),
                }
            }
            other => gate_failure(&id, "worker", &format!("unknown request type {other}")),
        }
    }
}

#[cfg(test)]
#[path = "protocol_tests.rs"]
mod tests;
