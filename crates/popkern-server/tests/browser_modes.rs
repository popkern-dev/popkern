//! PK-665: the popkern live browser harness (web/live/test/browser.mjs)
//! against `popkern serve` in both compute modes — projects on the
//! server, the kernel native or in the page. Skips like the page test
//! when the web bundle, node or a headless browser is missing.

use std::path::{Path, PathBuf};
use std::process::Command;

use popkern_server::{Compute, ServeOptions, bind};
use serde_json::Value;

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// The harness output for a fresh server in `compute` mode, or None
/// when the environment cannot run it.
fn harness(compute: Compute) -> Option<Value> {
    let page = root().join("web/live");
    if !page.join("pkg/popkern_live.js").exists()
        || Command::new("node").arg("--version").output().is_err()
    {
        eprintln!("[pk665] SKIP: web bundle or node absent");
        return None;
    }
    let projects =
        Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("browser-{}", compute.as_str()));
    let _ = std::fs::remove_dir_all(&projects);
    std::fs::create_dir_all(&projects).expect("projects");
    let (listener, handle) = bind(ServeOptions {
        root: projects.clone(),
        page,
        bind: "127.0.0.1:0".parse().expect("addr"),
        compute,
    })
    .expect("bind");
    let addr = listener.local_addr().expect("addr");
    std::thread::spawn(move || handle.run(listener));
    let shot =
        Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("serve-{}.png", compute.as_str()));
    let output = Command::new("node")
        .arg(root().join("web/live/test/browser.mjs"))
        .arg(&shot)
        .env("POPKERN_PAGE_ORIGIN", format!("http://{addr}"))
        .output()
        .expect("node");
    if output.status.code() == Some(3) {
        eprintln!("[pk665] SKIP: no headless browser");
        return None;
    }
    let out: Value = serde_json::from_slice(&output.stdout).unwrap_or_else(|e| {
        panic!(
            "harness output: {e}\n{}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        )
    });
    // A harness that stopped early says why before anything else.
    assert!(out["error"].is_null(), "harness error: {}", out["error"]);
    // Projects the page made live in the server's folder.
    assert!(
        projects.join("scratch/popkern.toml").is_file(),
        "Scratch is a folder on the server"
    );
    Some(out)
}

fn assert_common(out: &Value, host: &str) {
    assert!(out["error"].is_null(), "{out}");
    // Edit-to-refresh (PK-827, the PK-425 gate): a first-feature
    // parameter edit on the boss plate replays the tree, crosses the
    // socket and draws within a second (fastest of six: suite load
    // only adds to a sample).
    let e = &out["edit_latency"];
    eprintln!(
        "[pk827] {host}: edit-to-refresh {} ms (kernel {} ms, transfer + draw {} ms), samples {}",
        e["fastest_ms"], e["kernel_ms"], e["rest_ms"], e["samples"]
    );
    assert!(
        e["outcomes"]
            .as_array()
            .is_some_and(|o| o.iter().all(|x| x == "served")),
        "every edit serves: {e}"
    );
    assert!(
        e["fastest_ms"].as_f64().unwrap_or(f64::MAX) <= 1000.0,
        "edit-to-refresh within a second: {e}"
    );
    assert_eq!(out["exceptions"], serde_json::json!([]), "{out}");
    assert_eq!(out["host"], host, "{out}");
    assert_eq!(out["after_cancel"], "served", "{out}");
    assert_eq!(out["scrubber"]["ticks"], 6, "{}", out["scrubber"]);
    assert_eq!(out["scrubber"]["camera_kept"], true, "{}", out["scrubber"]);
    let m = &out["migrations"];
    assert_eq!(m["failed_badge"], "failed", "{m}");
    assert_eq!(m["overlay_names_upstream"], true, "{m}");
    assert!(
        out["foreign_requests"]
            .as_array()
            .is_some_and(Vec::is_empty),
        "{out}"
    );
}

#[test]
fn the_page_works_against_popkern_serve_with_native_compute() {
    let Some(out) = harness(Compute::Native) else {
        return;
    };
    assert_common(&out, "server-native");
    // Cancel stops the script and keeps the session (no restart).
    assert_eq!(out["cancel"]["restarts"], 0, "{}", out["cancel"]);
    assert!(
        out["cancel"]["status"]
            .as_str()
            .is_some_and(|s| s.contains("kept its model")),
        "{}",
        out["cancel"]
    );
    eprintln!(
        "[pk665] native: first render {} ms, history revisit {} ms",
        out["first_render_ms"], out["scrubber"]["revisit_ms"]
    );
}

#[test]
fn the_page_works_against_popkern_serve_with_browser_compute() {
    let Some(out) = harness(Compute::Browser) else {
        return;
    };
    assert_common(&out, "server-browser");
    // The kernel is in the page: cancel replaces the worker.
    assert_eq!(out["cancel"]["restarts"], 1, "{}", out["cancel"]);
}
