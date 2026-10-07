//! `popkern serve`'s MCP endpoint end to end: the JSON-RPC handshake,
//! notifications acknowledged without a body, the tool list, a project
//! created, edited and built natively (served, parameters applied, a
//! script error located), a STEP export on disk, the explicit/implicit
//! parameter split, and the transport's refusals (foreign origin, GET).

use std::io::{Read, Write};
use std::net::{SocketAddr, TcpStream};
use std::path::{Path, PathBuf};

use popkern_server::{Compute, ServeError, ServeOptions, bind};
use serde_json::{Value, json};

fn root(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("mcp-{name}"));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("root");
    dir
}

fn start(root: PathBuf) -> SocketAddr {
    let (listener, handle) = bind(ServeOptions {
        root,
        page: Path::new(env!("CARGO_MANIFEST_DIR")).join("../../web/live"),
        bind: "127.0.0.1:0".parse().expect("addr"),
        compute: Compute::Browser,
    })
    .expect("bind");
    let addr = listener.local_addr().expect("addr");
    std::thread::spawn(move || handle.run(listener));
    addr
}

/// One request to `/mcp`; (status, body).
fn request(addr: SocketAddr, method: &str, origin: Option<&str>, body: &str) -> (u16, String) {
    let mut stream = TcpStream::connect(addr).expect("connect");
    let origin = origin
        .map(|o| format!("Origin: {o}\r\n"))
        .unwrap_or_default();
    write!(
        stream,
        "{method} /mcp HTTP/1.1\r\nHost: localhost\r\n{origin}Content-Type: application/json\r\nAccept: application/json, text/event-stream\r\nContent-Length: {}\r\n\r\n{body}",
        body.len()
    )
    .expect("write");
    let mut response = String::new();
    stream.read_to_string(&mut response).expect("read");
    let status = response[9..12].parse().expect("status");
    let body = response
        .split_once("\r\n\r\n")
        .map(|(_, b)| b.to_string())
        .unwrap_or_default();
    (status, body)
}

fn rpc(addr: SocketAddr, id: u64, method: &str, params: Value) -> Value {
    let message = json!({ "jsonrpc": "2.0", "id": id, "method": method, "params": params });
    let (status, body) = request(addr, "POST", None, &message.to_string());
    assert_eq!(status, 200, "{body}");
    let reply: Value = serde_json::from_str(&body).expect("json");
    assert_eq!(reply["id"], id);
    reply
}

/// A tool call's structured result and its error flag.
fn call(addr: SocketAddr, tool: &str, arguments: Value) -> (Value, bool) {
    let reply = rpc(
        addr,
        7,
        "tools/call",
        json!({ "name": tool, "arguments": arguments }),
    );
    let result = &reply["result"];
    let text = result["content"][0]["text"].as_str().expect("text content");
    let structured: Value = serde_json::from_str(text).expect("text is the JSON");
    assert_eq!(structured, result["structuredContent"]);
    (structured, result["isError"].as_bool().expect("isError"))
}

#[test]
fn the_handshake_and_the_tool_list() {
    let addr = start(root("handshake"));
    let reply = rpc(
        addr,
        1,
        "initialize",
        json!({ "protocolVersion": "2025-03-26", "capabilities": {}, "clientInfo": { "name": "test", "version": "0" } }),
    );
    assert_eq!(reply["result"]["protocolVersion"], "2025-03-26");
    assert_eq!(reply["result"]["serverInfo"]["name"], "popkern");
    assert!(reply["result"]["capabilities"]["tools"].is_object());
    // An unknown revision is answered with the newest this server speaks.
    let reply = rpc(
        addr,
        2,
        "initialize",
        json!({ "protocolVersion": "1999-01-01" }),
    );
    assert_eq!(reply["result"]["protocolVersion"], "2025-06-18");
    // A notification has no answer: 202 and no body.
    let (status, body) = request(
        addr,
        "POST",
        None,
        r#"{"jsonrpc":"2.0","method":"notifications/initialized"}"#,
    );
    assert_eq!((status, body.as_str()), (202, ""));
    let tools = rpc(addr, 3, "tools/list", json!({}));
    let names: Vec<&str> = tools["result"]["tools"]
        .as_array()
        .expect("tools")
        .iter()
        .map(|t| t["name"].as_str().expect("name"))
        .collect();
    for expected in [
        "list_projects",
        "read_project",
        "create_project",
        "write_file",
        "delete_file",
        "build_project",
        "run_script",
        "export_step",
    ] {
        assert!(names.contains(&expected), "{expected} in {names:?}");
    }
    assert_eq!(rpc(addr, 4, "ping", json!({}))["result"], json!({}));
    assert_eq!(rpc(addr, 5, "no/such", json!({}))["error"]["code"], -32601);
    // A batch answers each request in it.
    let (status, body) = request(
        addr,
        "POST",
        None,
        r#"[{"jsonrpc":"2.0","id":8,"method":"ping"},{"jsonrpc":"2.0","method":"notifications/initialized"}]"#,
    );
    assert_eq!(status, 200);
    let batch: Value = serde_json::from_str(&body).expect("json");
    assert_eq!(batch.as_array().map(Vec::len), Some(1));
    // Not JSON at all: a parse error, not a dropped connection.
    let (_, body) = request(addr, "POST", None, "{nope");
    assert_eq!(
        serde_json::from_str::<Value>(&body).expect("json")["error"]["code"],
        -32700
    );
}

#[test]
fn foreign_origins_and_get_are_refused() {
    let addr = start(root("refusals"));
    let ping = r#"{"jsonrpc":"2.0","id":1,"method":"ping"}"#;
    assert_eq!(
        request(addr, "POST", Some("http://evil.example"), ping).0,
        403
    );
    assert_eq!(
        request(addr, "POST", Some("http://127.0.0.1:7690"), ping).0,
        200
    );
    assert_eq!(request(addr, "GET", None, "").0, 405);
}

#[test]
fn a_project_is_created_edited_built_and_exported() {
    let dir = root("build");
    let addr = start(dir.clone());
    let (created, failed) = call(
        addr,
        "create_project",
        json!({ "name": "Bracket", "source": "let m = model();\nm.param(\"wall\", 0.01, \"m\");\nm.profile(xy()).rect(0.08, 0.05).sketch(\"plate\").extrude(blind(p(\"wall\")));\nm\n" }),
    );
    assert!(!failed, "{created}");
    assert_eq!(created["id"], "bracket");
    assert!(dir.join("bracket/migrations/001_base.rhai").is_file());

    let (list, _) = call(addr, "list_projects", json!({}));
    assert_eq!(list["projects"][0]["id"], "bracket");
    assert_eq!(
        list["projects"][0]["path"],
        json!(
            dir.canonicalize()
                .expect("root")
                .join("bracket")
                .display()
                .to_string()
        )
    );

    // A second migration (bare file name) rounds the vertical edges.
    let (written, failed) = call(
        addr,
        "write_file",
        json!({ "project": "bracket", "path": "002_round.rhai", "text": "let m = previous();\nm.bodies()[0].fillet(edges().straight().direction(\"z\"), 0.005);\nm\n" }),
    );
    assert!(!failed, "{written}");
    assert_eq!(written["written"], "migrations/002_round.rhai");

    let (built, failed) = call(addr, "build_project", json!({ "project": "bracket" }));
    assert!(!failed, "{built:#}");
    assert_eq!(built["outcome"], "served");
    assert_eq!(built["bodies"][0]["watertight"], true);
    assert_eq!(
        built["project"]["versions"].as_array().map(Vec::len),
        Some(2)
    );
    // The declared wall is explicit; the fillet radius is implicit.
    let params = built["params"].as_array().expect("params");
    let wall = params.iter().find(|p| p["name"] == "wall").expect("wall");
    assert_eq!(wall["explicit"], true);
    assert!(
        params.iter().any(|p| p["explicit"] == false),
        "an implicit feature parameter in {params:?}"
    );
    let height = |b: &Value| {
        b["bodies"][0]["bounds"][1][2].as_f64().expect("z")
            - b["bodies"][0]["bounds"][0][2].as_f64().expect("z")
    };
    assert!((height(&built) - 0.01).abs() < 1e-6);

    // Parameters apply after the build.
    let (thick, failed) = call(
        addr,
        "build_project",
        json!({ "project": "bracket", "params": { "wall": 0.02 } }),
    );
    assert!(!failed, "{thick:#}");
    assert!((height(&thick) - 0.02).abs() < 1e-6);

    // A broken migration is reported with where it broke.
    call(
        addr,
        "write_file",
        json!({ "project": "bracket", "path": "migrations/003_broken.rhai", "text": "let m = previous();\nm.no_such_call();\nm\n" }),
    );
    let (broken, failed) = call(addr, "build_project", json!({ "project": "bracket" }));
    assert!(failed);
    assert_eq!(broken["outcome"], "script_error", "{broken:#}");
    assert_eq!(broken["project"]["error"]["migration"], "003_broken.rhai");
    let (_, failed) = call(
        addr,
        "delete_file",
        json!({ "project": "bracket", "path": "migrations/003_broken.rhai" }),
    );
    assert!(!failed);

    let (exported, failed) = call(addr, "export_step", json!({ "project": "bracket" }));
    assert!(!failed, "{exported}");
    let step = std::fs::read_to_string(exported["path"].as_str().expect("path")).expect("step");
    assert!(step.starts_with("ISO-10303-21;"));
    assert!(
        dir.join("bracket/.popkern/exports/bracket.step").is_file(),
        "the export lands in the reserved build folder"
    );

    // The page's API sees the agent's edits (one folder, one truth).
    let (read, _) = call(addr, "read_project", json!({ "project": "bracket" }));
    assert_eq!(
        read["files"]
            .as_object()
            .map(|f| f.keys().cloned().collect::<Vec<_>>()),
        Some(vec![
            "migrations/001_base.rhai".to_string(),
            "migrations/002_round.rhai".to_string(),
            "popkern.toml".to_string(),
        ])
    );
}

#[test]
fn a_loose_script_runs_without_saving() {
    let dir = root("script");
    let addr = start(dir.clone());
    let (ran, failed) = call(
        addr,
        "run_script",
        json!({ "source": "let m = model();\nm.profile(xy()).circle(0.01).sketch(\"s\").extrude(blind(0.02));\nm\n", "face_names": true }),
    );
    assert!(!failed, "{ran:#}");
    assert_eq!(ran["bodies"][0]["faces"], 3);
    assert_eq!(
        ran["bodies"][0]["face_names"].as_array().map(Vec::len),
        Some(3)
    );
    assert_eq!(std::fs::read_dir(&dir).expect("dir").count(), 0);
    let (unknown, failed) = call(addr, "read_project", json!({ "project": "nope" }));
    assert!(failed);
    assert!(unknown["error"].as_str().expect("error").contains("nope"));
}

#[test]
fn a_project_folder_is_not_a_projects_root() {
    let dir = root("single");
    std::fs::write(dir.join("popkern.toml"), "[project]\nname = \"x\"\n").expect("manifest");
    let refused = bind(ServeOptions {
        root: dir,
        page: Path::new(env!("CARGO_MANIFEST_DIR")).join("../../web/live"),
        bind: "127.0.0.1:0".parse().expect("addr"),
        compute: Compute::Browser,
    });
    assert!(matches!(refused, Err(ServeError::ProjectFolder { .. })));
}
