//! PK-665: `popkern serve` end to end — loopback-only binding, the page
//! and host.json, the projects API with revision-checked saves, and the
//! session socket (protocol replies with binary attachments; a cancel
//! stops a runaway script without losing the session's document).

use std::io::{Read, Write};
use std::net::{SocketAddr, TcpStream};
use std::path::{Path, PathBuf};
use std::time::Duration;

use popkern_server::{Compute, ServeError, ServeOptions, bind};
use serde_json::{Value, json};
use tungstenite::Message;

fn page() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../web/live")
}

fn root(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("serve-{name}"));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("root");
    dir
}

fn start(name: &str, compute: Compute) -> SocketAddr {
    let (listener, handle) = bind(ServeOptions {
        root: root(name),
        page: page(),
        bind: "127.0.0.1:0".parse().expect("addr"),
        compute,
    })
    .expect("bind");
    let addr = listener.local_addr().expect("addr");
    std::thread::spawn(move || handle.run(listener));
    addr
}

/// One HTTP request; (status, body).
fn http(addr: SocketAddr, method: &str, path: &str, body: Option<&Value>) -> (u16, Vec<u8>) {
    let mut stream = TcpStream::connect(addr).expect("connect");
    let payload = body.map(Value::to_string).unwrap_or_default();
    write!(
        stream,
        "{method} {path} HTTP/1.1\r\nHost: localhost\r\nContent-Length: {}\r\nContent-Type: application/json\r\n\r\n{payload}",
        payload.len()
    )
    .expect("write");
    let mut response = Vec::new();
    stream.read_to_end(&mut response).expect("read");
    let split = response
        .windows(4)
        .position(|w| w == b"\r\n\r\n")
        .expect("head")
        + 4;
    let status = String::from_utf8_lossy(&response[9..12])
        .parse()
        .expect("status");
    (status, response[split..].to_vec())
}

fn json_body(bytes: &[u8]) -> Value {
    serde_json::from_slice(bytes).expect("json")
}

#[test]
fn only_loopback_addresses_are_served() {
    let refused = bind(ServeOptions {
        root: root("refused"),
        page: page(),
        bind: "0.0.0.0:0".parse().expect("addr"),
        compute: Compute::Native,
    });
    assert!(matches!(refused, Err(ServeError::NonLoopback { .. })));
}

#[test]
fn the_page_and_host_json() {
    let addr = start("page", Compute::Native);
    let (status, body) = http(addr, "GET", "/", None);
    assert_eq!(status, 200);
    assert!(String::from_utf8_lossy(&body).contains("popkern live"));
    let (status, body) = http(addr, "GET", "/host.json", None);
    assert_eq!(status, 200);
    assert_eq!(
        json_body(&body),
        json!({ "storage": "server", "compute": "native", "version": env!("CARGO_PKG_VERSION") })
    );
    assert_eq!(http(addr, "GET", "/../Cargo.toml", None).0, 404);
    assert_eq!(http(addr, "GET", "/%2e%2e/Cargo.toml", None).0, 404);
}

#[test]
fn projects_save_with_revisions() {
    let addr = start("projects", Compute::Native);
    let files = json!({ "popkern.toml": "[project]\nname = \"Plate\"\n", "migrations/001_base.rhai": "let m = model();\nm\n" });
    let (status, body) = http(
        addr,
        "POST",
        "/api/projects",
        Some(&json!({ "name": "Plate", "files": files })),
    );
    assert_eq!(status, 201);
    let created = json_body(&body);
    assert_eq!(created["id"], "plate");
    let rev = created["revision"].as_str().expect("rev").to_string();
    let mut edited = files.clone();
    edited["migrations/001_base.rhai"] = json!("// edited\nlet m = model();\nm\n");
    let (status, body) = http(
        addr,
        "PUT",
        "/api/projects/plate",
        Some(&json!({ "files": edited, "revision": rev })),
    );
    assert_eq!(status, 200);
    assert_ne!(json_body(&body)["revision"], json!(rev));
    // The same (now stale) revision again: 409 with the stored project.
    let (status, body) = http(
        addr,
        "PUT",
        "/api/projects/plate",
        Some(&json!({ "files": files, "revision": rev })),
    );
    assert_eq!(status, 409);
    assert!(
        json_body(&body)["project"]["files"]["migrations/001_base.rhai"]
            .as_str()
            .unwrap_or("")
            .starts_with("// edited")
    );
    let (_, list) = http(addr, "GET", "/api/projects", None);
    assert_eq!(json_body(&list).as_array().map(Vec::len), Some(1));
    assert_eq!(http(addr, "DELETE", "/api/projects/plate", None).0, 204);
    assert_eq!(http(addr, "GET", "/api/projects/plate", None).0, 404);
}

/// The next reply on `ws`: its JSON and binary frames.
fn reply(ws: &mut tungstenite::WebSocket<TcpStream>) -> (Value, Vec<Vec<u8>>) {
    loop {
        match ws.read().expect("read") {
            Message::Text(text) => {
                let head: Value = serde_json::from_str(&text).expect("json");
                let n = head["attachments"].as_array().map(Vec::len).unwrap_or(0);
                let mut frames = Vec::new();
                while frames.len() < n {
                    if let Message::Binary(b) = ws.read().expect("binary") {
                        frames.push(b);
                    }
                }
                return (head, frames);
            }
            _ => continue,
        }
    }
}

#[test]
fn the_session_socket_runs_and_cancels_natively() {
    let addr = start("session", Compute::Native);
    let stream = TcpStream::connect(addr).expect("connect");
    let (mut ws, _) = tungstenite::client(format!("ws://{addr}/api/session"), stream).expect("ws");
    let send = |ws: &mut tungstenite::WebSocket<TcpStream>, v: Value| {
        ws.send(Message::Text(v.to_string())).expect("send")
    };
    send(&mut ws, json!({ "id": 1, "type": "init", "options": {} }));
    let (init, _) = reply(&mut ws);
    assert_eq!(init["ok"], true);
    assert_eq!(init["clock_available"], true);
    let plate = "let m = model();\nm.profile(xy()).rect(1.0, 1.0).sketch(\"plate\").extrude(blind(0.2));\nm\n";
    send(
        &mut ws,
        json!({ "id": 2, "type": "run", "source": plate, "name": "plate.rhai" }),
    );
    let (run, frames) = reply(&mut ws);
    assert_eq!(run["envelope"]["outcome"], "served");
    assert_eq!(frames.len(), 6, "one body: six buffers");
    assert_eq!(
        run["attachments"][0],
        json!({ "path": "buffers.0.positions", "kind": "f32" })
    );
    let vertices = run["envelope"]["bodies"][0]["vertex_count"]
        .as_u64()
        .expect("vertices") as usize;
    assert_eq!(
        frames[0].len(),
        vertices * 3 * 4,
        "positions as little-endian f32"
    );

    // A runaway script: cancel stops it, and the session keeps the plate.
    send(
        &mut ws,
        json!({ "id": 3, "type": "run", "source": "let m = model();\nloop {}\n", "name": "runaway.rhai" }),
    );
    std::thread::sleep(Duration::from_millis(300));
    let started = std::time::Instant::now();
    send(&mut ws, json!({ "type": "cancel" }));
    let (cancelled, _) = reply(&mut ws);
    assert!(
        started.elapsed() < Duration::from_secs(5),
        "cancel answers promptly"
    );
    assert_eq!(cancelled["id"], 3);
    assert_eq!(cancelled["envelope"]["kind"], "Cancelled", "{cancelled}");
    send(&mut ws, json!({ "id": 4, "type": "export_pkn" }));
    let (pkn, frames) = reply(&mut ws);
    assert_eq!(pkn["ok"], true, "the session kept its document");
    let doc: Value = serde_json::from_slice(&frames[0]).expect("pkn");
    assert!(
        doc["tree"]["timeline"]
            .as_array()
            .is_some_and(|t| t.len() == 2),
        "the plate's two features"
    );
}

#[test]
fn browser_compute_has_no_session_socket() {
    let addr = start("browser", Compute::Browser);
    let (_, body) = http(addr, "GET", "/host.json", None);
    assert_eq!(json_body(&body)["compute"], "browser");
    let stream = TcpStream::connect(addr).expect("connect");
    assert!(tungstenite::client(format!("ws://{addr}/api/session"), stream).is_err());
}
