//! One HTTP connection: the page's static files, `host.json`, the
//! projects API, and the session WebSocket upgrade. Minimal HTTP/1.1:
//! one request per connection (`Connection: close`), bodies by
//! `Content-Length`.

use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::net::TcpStream;
use std::path::{Component, Path};

use serde_json::{Value, json};

use crate::store::{ProjectsRoot, StoreError};
use crate::{Compute, Context};

/// Requests larger than this are refused (a project is small text).
const MAX_BODY: usize = 32 * 1024 * 1024;

struct Request {
    method: String,
    path: String,
    headers: Vec<(String, String)>,
}

impl Request {
    fn header(&self, name: &str) -> Option<&str> {
        self.headers
            .iter()
            .find(|(k, _)| k.eq_ignore_ascii_case(name))
            .map(|(_, v)| v.as_str())
    }
}

/// Read the request head without consuming it (a WebSocket upgrade
/// hands the untouched stream to the handshake).
fn peek_head(stream: &TcpStream) -> Option<(Request, usize)> {
    let mut buf = vec![0u8; 16 * 1024];
    for _ in 0..200 {
        let n = stream.peek(&mut buf).ok()?;
        let mut headers = [httparse::EMPTY_HEADER; 64];
        let mut req = httparse::Request::new(&mut headers);
        match req.parse(&buf[..n]) {
            Ok(httparse::Status::Complete(len)) => {
                let request = Request {
                    method: req.method?.to_string(),
                    path: req.path?.to_string(),
                    headers: req
                        .headers
                        .iter()
                        .map(|h| {
                            (
                                h.name.to_string(),
                                String::from_utf8_lossy(h.value).into_owned(),
                            )
                        })
                        .collect(),
                };
                return Some((request, len));
            }
            Ok(httparse::Status::Partial) if n < buf.len() => {
                std::thread::sleep(std::time::Duration::from_millis(5));
            }
            _ => return None,
        }
    }
    None
}

fn respond(stream: &mut TcpStream, status: u16, content_type: &str, body: &[u8]) {
    let reason = match status {
        200 => "OK",
        201 => "Created",
        204 => "No Content",
        400 => "Bad Request",
        404 => "Not Found",
        409 => "Conflict",
        413 => "Payload Too Large",
        _ => "Error",
    };
    let head = format!(
        "HTTP/1.1 {status} {reason}\r\nContent-Type: {content_type}\r\nContent-Length: {}\r\nCache-Control: no-cache\r\nX-Content-Type-Options: nosniff\r\nConnection: close\r\n\r\n",
        body.len()
    );
    let _ = stream.write_all(head.as_bytes());
    let _ = stream.write_all(body);
}

fn json_reply(stream: &mut TcpStream, status: u16, value: &Value) {
    respond(
        stream,
        status,
        "application/json",
        value.to_string().as_bytes(),
    );
}

fn content_type(path: &Path) -> &'static str {
    match path.extension().and_then(|e| e.to_str()).unwrap_or("") {
        "html" => "text/html; charset=utf-8",
        "js" | "mjs" => "text/javascript; charset=utf-8",
        "css" => "text/css; charset=utf-8",
        "json" => "application/json",
        "wasm" => "application/wasm",
        "svg" => "image/svg+xml",
        "png" => "image/png",
        "rhai" | "toml" | "md" | "txt" => "text/plain; charset=utf-8",
        _ => "application/octet-stream",
    }
}

/// A page file for URL path `path`, never outside the page folder.
fn static_file(page: &Path, path: &str) -> Option<std::path::PathBuf> {
    let path = path.split(['?', '#']).next().unwrap_or("/");
    let relative = path.trim_start_matches('/');
    let relative = if relative.is_empty() {
        "index.html"
    } else {
        relative
    };
    let candidate = Path::new(relative);
    if candidate
        .components()
        .any(|c| !matches!(c, Component::Normal(_)))
    {
        return None;
    }
    let full = page.join(candidate);
    full.is_file().then_some(full)
}

/// Serve one connection.
pub(crate) fn connection(mut stream: TcpStream, context: &Context) {
    let Some((request, head_len)) = peek_head(&stream) else {
        respond(&mut stream, 400, "text/plain", b"bad request");
        return;
    };
    let path = request.path.split('?').next().unwrap_or("/").to_string();
    if path == "/api/session" {
        if context.options.compute == Compute::Native
            && request
                .header("upgrade")
                .is_some_and(|u| u.eq_ignore_ascii_case("websocket"))
        {
            crate::socket::session(stream);
        } else {
            let mut head = vec![0u8; head_len];
            let _ = stream.read_exact(&mut head);
            respond(
                &mut stream,
                404,
                "text/plain",
                b"this host runs the kernel in the page",
            );
        }
        return;
    }
    // Consume the head; read the body by Content-Length.
    let mut head = vec![0u8; head_len];
    if stream.read_exact(&mut head).is_err() {
        return;
    }
    let length: usize = request
        .header("content-length")
        .and_then(|v| v.trim().parse().ok())
        .unwrap_or(0);
    if length > MAX_BODY {
        respond(&mut stream, 413, "text/plain", b"too large");
        return;
    }
    let mut body = vec![0u8; length];
    if stream.read_exact(&mut body).is_err() {
        return;
    }
    if path == "/host.json" {
        let value = json!({
            "storage": "server",
            "compute": context.options.compute.as_str(),
            "version": env!("CARGO_PKG_VERSION"),
        });
        json_reply(&mut stream, 200, &value);
        return;
    }
    if let Some(rest) = path.strip_prefix("/api/projects") {
        projects(&mut stream, &request.method, rest, &body, context);
        return;
    }
    if request.method != "GET" && request.method != "HEAD" {
        respond(&mut stream, 404, "text/plain", b"not found");
        return;
    }
    match static_file(&context.options.page, &path) {
        Some(file) => match std::fs::read(&file) {
            Ok(bytes) => respond(&mut stream, 200, content_type(&file), &bytes),
            Err(_) => respond(&mut stream, 404, "text/plain", b"not found"),
        },
        None => respond(&mut stream, 404, "text/plain", b"not found"),
    }
}

fn store_error(stream: &mut TcpStream, error: StoreError) {
    match error {
        StoreError::NotFound(id) => {
            json_reply(stream, 404, &json!({ "error": format!("no project {id}") }))
        }
        StoreError::Conflict(current) => json_reply(
            stream,
            409,
            &json!({ "error": "conflict", "project": current }),
        ),
        StoreError::BadPath(p) => json_reply(
            stream,
            400,
            &json!({ "error": format!("{p} is not a project file") }),
        ),
        StoreError::Io(detail) => json_reply(stream, 500, &json!({ "error": detail })),
    }
}

fn files_of(body: &Value) -> BTreeMap<String, String> {
    body.get("files")
        .and_then(Value::as_object)
        .map(|m| {
            m.iter()
                .filter_map(|(k, v)| v.as_str().map(|s| (k.clone(), s.to_string())))
                .collect()
        })
        .unwrap_or_default()
}

/// `/api/projects[/<id>]`.
fn projects(stream: &mut TcpStream, method: &str, rest: &str, body: &[u8], context: &Context) {
    let store = ProjectsRoot::new(&context.options.root);
    let id = rest.trim_start_matches('/');
    let body: Value = serde_json::from_slice(body).unwrap_or(Value::Null);
    let result = match (method, id.is_empty()) {
        ("GET", true) => store.list().map(|list| (200, json!(list))),
        ("POST", true) => {
            let name = body
                .get("name")
                .and_then(Value::as_str)
                .unwrap_or("Project");
            store
                .create(name, &files_of(&body))
                .map(|p| (201, json!(p)))
        }
        ("GET", false) => store.get(id).map(|p| (200, json!(p))),
        ("PUT", false) => {
            let revision = body.get("revision").and_then(Value::as_str).unwrap_or("");
            store
                .save(id, &files_of(&body), revision)
                .map(|p| (200, json!(p)))
        }
        ("DELETE", false) => store.remove(id).map(|()| (204, Value::Null)),
        _ => {
            json_reply(stream, 400, &json!({ "error": "unsupported" }));
            return;
        }
    };
    match result {
        Ok((204, _)) => respond(stream, 204, "text/plain", b""),
        Ok((status, value)) => json_reply(stream, status, &value),
        Err(error) => store_error(stream, error),
    }
}
