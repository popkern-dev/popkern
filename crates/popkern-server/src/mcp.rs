//! The Model Context Protocol endpoint (`/mcp`): an agent drives the
//! same projects the page edits, and builds them on the native kernel.
//!
//! Transport: MCP "Streamable HTTP" without the optional server stream.
//! Each `POST /mcp` carries one JSON-RPC message (or a batch) and is
//! answered with `application/json`; notifications and client responses
//! are acknowledged with `202 Accepted`. `GET /mcp` answers 405 (the
//! server never pushes). There are no sessions: every call reads the
//! projects folder afresh, so the page and an agent see one truth, and
//! the page reloads a project an agent changed (projects.js polls the
//! open project's revision).
//!
//! Like the rest of `popkern serve` this is loopback only with no
//! authentication; a request whose `Origin` is not a loopback origin is
//! refused (403) so a web page elsewhere cannot drive the endpoint
//! through DNS rebinding.

use serde_json::{Value, json};

use crate::Context;
use crate::mcp_tools;

/// Protocol revisions this endpoint speaks, newest first.
const PROTOCOL_VERSIONS: [&str; 3] = ["2025-06-18", "2025-03-26", "2024-11-05"];

/// What the server tells a connecting agent about itself.
const INSTRUCTIONS: &str = "popkern is a B-rep CAD kernel driven by Rhai model scripts. \
A project is a folder holding popkern.toml and ordered migrations (migrations/NNN_name.rhai): \
the first builds the document (`let m = model(); ...; m`), each later one starts with \
`let m = previous();`, edits it and ends with `m`. Use list_projects and read_project to see \
what exists, write_file to change a migration, build_project to run it on the kernel (the \
outcome, failures, parameters and bodies come back), and export_step for a STEP file. \
Edits appear in the popkern live page that serves the same folder.";

/// The answer to one `POST /mcp`: an HTTP status and an optional body.
pub(crate) struct McpReply {
    /// 200 with a body, 202 without.
    pub(crate) status: u16,
    /// The JSON-RPC response(s).
    pub(crate) body: Option<Value>,
}

/// Whether `origin` (an `Origin` header) names this machine.
pub(crate) fn loopback_origin(origin: &str) -> bool {
    let Some(rest) = origin
        .strip_prefix("http://")
        .or_else(|| origin.strip_prefix("https://"))
    else {
        return false;
    };
    let host = if let Some(v6) = rest.strip_prefix('[') {
        v6.split(']').next().unwrap_or("")
    } else {
        rest.split([':', '/']).next().unwrap_or("")
    };
    host == "localhost"
        || host.ends_with(".localhost")
        || host
            .parse::<std::net::IpAddr>()
            .is_ok_and(|ip| ip.is_loopback())
}

fn error(id: Value, code: i64, message: &str) -> Value {
    json!({ "jsonrpc": "2.0", "id": id, "error": { "code": code, "message": message } })
}

fn result(id: Value, value: Value) -> Value {
    json!({ "jsonrpc": "2.0", "id": id, "result": value })
}

/// Answer the body of one `POST /mcp`.
pub(crate) fn post(body: &[u8], context: &Context) -> McpReply {
    let Ok(message) = serde_json::from_slice::<Value>(body) else {
        return McpReply {
            status: 200,
            body: Some(error(Value::Null, -32700, "parse error: not JSON")),
        };
    };
    let reply = match message {
        Value::Array(batch) => {
            let answers: Vec<Value> = batch.iter().filter_map(|m| answer(m, context)).collect();
            (!answers.is_empty()).then_some(Value::Array(answers))
        }
        single => answer(&single, context),
    };
    McpReply {
        status: if reply.is_some() { 200 } else { 202 },
        body: reply,
    }
}

/// One message's response; `None` for a notification or a client
/// response (nothing to answer).
fn answer(message: &Value, context: &Context) -> Option<Value> {
    let method = message.get("method").and_then(Value::as_str);
    let id = message.get("id").cloned();
    let (Some(method), Some(id)) = (method, id) else {
        // A notification (`notifications/initialized`, a cancel) or a
        // response to a server request this server never sends.
        return None;
    };
    let params = message.get("params").cloned().unwrap_or(json!({}));
    Some(match method {
        "initialize" => {
            let asked = params
                .get("protocolVersion")
                .and_then(Value::as_str)
                .unwrap_or("");
            let version = PROTOCOL_VERSIONS
                .iter()
                .find(|v| **v == asked)
                .unwrap_or(&PROTOCOL_VERSIONS[0]);
            result(
                id,
                json!({
                    "protocolVersion": version,
                    "capabilities": { "tools": { "listChanged": false } },
                    "serverInfo": { "name": "popkern", "version": env!("CARGO_PKG_VERSION") },
                    "instructions": INSTRUCTIONS,
                }),
            )
        }
        "ping" => result(id, json!({})),
        "tools/list" => result(id, json!({ "tools": mcp_tools::definitions() })),
        "tools/call" => {
            let Some(name) = params.get("name").and_then(Value::as_str) else {
                return Some(error(id, -32602, "tools/call needs a tool name"));
            };
            let arguments = params.get("arguments").cloned().unwrap_or(json!({}));
            match mcp_tools::call(name, &arguments, context) {
                Some(outcome) => result(id, outcome.to_json()),
                None => error(id, -32602, &format!("unknown tool {name}")),
            }
        }
        // No resources or prompts: say so politely rather than erroring
        // clients that probe every list.
        "resources/list" => result(id, json!({ "resources": [] })),
        "prompts/list" => result(id, json!({ "prompts": [] })),
        _ => error(id, -32601, &format!("method {method} is not supported")),
    })
}

#[cfg(test)]
mod tests {
    use super::loopback_origin;

    #[test]
    fn only_loopback_origins_are_this_machine() {
        for ok in [
            "http://localhost:7690",
            "http://127.0.0.1",
            "http://[::1]:7690",
            "https://app.localhost",
        ] {
            assert!(loopback_origin(ok), "{ok}");
        }
        for refused in [
            "http://evil.example",
            "http://localhost.evil.example",
            "http://192.168.1.4:7690",
            "null",
            "file://",
        ] {
            assert!(!loopback_origin(refused), "{refused}");
        }
    }
}
