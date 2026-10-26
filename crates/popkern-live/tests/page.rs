//! PK-411: the static page's contracts, checked without a browser.
//!
//! - The worker protocol (`web/live/PROTOCOL.md`) is driven end
//!   to end under node's worker_threads by `web/live/test/
//!   host.test.mjs`: init, run, set_param, exports, show-as-code, a
//!   malformed request, and cancel of a script that ignores its
//!   budget (the worker is terminated and replaced; the next run
//!   answers). The answers are compared with the native engine.
//! - Privacy posture: the page's own JavaScript reaches the network
//!   only to load same-origin assets (`fetch` with a relative literal
//!   path); `XMLHttpRequest`, `WebSocket`, `sendBeacon`, `EventSource`
//!   and dynamic `fetch` targets are absent. The CSP meta pins every
//!   source to `'self'`.
//! - Supply chain: every vendored file's SHA-256 matches
//!   `web/live/vendor/LICENSES.md`.
//! - PROTOCOL.md names every request type the worker core handles.
//!
//! The harness needs the web bundle from `tools/build-live.sh
//! --web`; without it the protocol test records why it skipped.

use std::path::{Path, PathBuf};
use std::process::Command;

#[path = "support/page_rows.rs"]
mod page_rows;
#[path = "support/sim_rows.rs"]
mod sim_rows;
mod support;

use popkern_live::mesh_buffers::fnv1a_u32;
use popkern_live::{HostOptions, Session};
use serde_json::Value;
use support::{base64_decode, fnv1a_bytes, sha256_hex};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("workspace root")
}

fn page() -> PathBuf {
    root().join("web/live")
}

fn read(path: impl AsRef<Path>) -> String {
    let path = path.as_ref();
    std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

#[test]
fn worker_protocol_round_trips_under_node_and_cancel_replaces_a_wedged_worker() {
    let bundle = page().join("pkg/popkern_live.js");
    if !bundle.exists() {
        eprintln!(
            "[pk411] SKIP: no web bundle at {} (run tools/build-live.sh --web)",
            bundle.display()
        );
        return;
    }
    if Command::new("node").arg("--version").output().is_err() {
        eprintln!("[pk411] SKIP: node is not installed");
        return;
    }
    let script = root().join("crates/popkern-rhai/examples/boss_plate.rhai");
    let output = Command::new("node")
        .arg(page().join("test/host.test.mjs"))
        .arg(&script)
        .output()
        .expect("node");
    assert!(
        output.status.success(),
        "host.test.mjs failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let out: Value = serde_json::from_slice(&output.stdout).expect("harness JSON");
    assert_eq!(out["init"]["protocol_version"], "1");
    assert_eq!(out["init"]["clock_available"], true);
    assert_eq!(out["run"]["ok"], true, "{out}");
    assert_eq!(out["run"]["outcome"], "served");
    assert_eq!(out["run"]["features"], 6);
    page_rows::assert_history_rows(&out);
    assert_eq!(out["run"]["bodies"], 1);
    assert_eq!(
        out["run"]["transferred"][0], true,
        "typed arrays transferred"
    );
    // Same answers as the native engine.
    let mut session = Session::new(HostOptions::default());
    let native = session.run(&read(&script), "boss_plate.rhai");
    assert!(native.outcome.is_served());
    assert_eq!(
        out["run"]["indices_fnv"][0].as_str(),
        Some(fnv1a_u32(&native.bodies[0].indices).to_string().as_str())
    );
    assert_eq!(
        out["run"]["face_ids_fnv"][0].as_str(),
        Some(fnv1a_u32(&native.bodies[0].face_ids).to_string().as_str())
    );
    assert_eq!(
        out["run"]["face_labels"][0],
        serde_json::to_value(&native.bodies[0].face_labels).expect("json")
    );
    assert_eq!(out["run"]["params"], native.params.len());
    // The zebra analysis mesh (PK-577): 6°, finer than the display
    // mesh, the same buffers as the native engine's.
    let analysis = &out["analysis"];
    assert_eq!(analysis["ok"], true, "{analysis}");
    assert_eq!(analysis["facet_deg"], 6.0);
    assert_eq!(analysis["triangles"], analysis["buffer_triangles"]);
    assert!(
        analysis["triangles"].as_f64() > analysis["display_triangles"].as_f64(),
        "{analysis}"
    );
    let native_analysis = session.analysis_meshes().expect("native analysis mesh");
    assert_eq!(
        analysis["indices_fnv"][0].as_str(),
        Some(fnv1a_u32(&native_analysis[0].indices).to_string().as_str())
    );
    // Parameter edit replays; a bad name is the ArgumentError class.
    assert_eq!(out["set_param"]["outcome"], "served");
    assert_eq!(out["set_param"]["value"], 0.5);
    assert_eq!(out["set_param"]["bodies"], 1);
    assert_eq!(out["bad_param"]["outcome"], "script_error");
    assert_eq!(out["bad_param"]["kind"], "ArgumentError");
    // Exports and show-as-code agree with the native engine after the same edit.
    let edited = session.set_param("extrude-1.distance", 0.5);
    assert!(edited.outcome.is_served());
    let pkn = session.export_pkn().expect("pkn");
    let rhai = session.to_rhai().expect("rhai");
    assert_eq!(out["exports"]["step_ok"], true);
    assert_eq!(
        out["exports"]["pkn_fnv"].as_str(),
        Some(fnv1a_bytes(&pkn).to_string().as_str()),
        ".pkn bytes"
    );
    assert_eq!(
        out["exports"]["rhai_fnv"].as_str(),
        Some(fnv1a_bytes(rhai.as_bytes()).to_string().as_str()),
        "canonical Rhai"
    );
    // One error shape for a malformed request.
    let budget = &out["set_budget"];
    assert_eq!(budget["ok"], true, "{out}");
    assert_eq!(budget["echoed"]["max_operations"], 1000);
    assert!(budget["echoed"]["timeout_ms"].is_null());
    assert_eq!(budget["stopped_kind"], "BudgetExceeded", "{out}");
    assert_eq!(budget["finished_outcome"], "served", "{out}");
    assert_eq!(out["bogus"]["outcome"], "gate_failure");
    assert_eq!(out["bogus"]["step"], "worker");
    // Cancel: the pending request resolves as a gate failure within
    // the budget, the worker is replaced, and the next run serves.
    let cancel = &out["cancel"];
    let budget = cancel["budget_ms"].as_f64().expect("budget");
    let cancel_ms = cancel["cancel_ms"].as_f64().expect("cancel_ms");
    let restart_ms = cancel["restart_ms"].as_f64().expect("restart_ms");
    assert!(
        cancel_ms <= 2.0 * budget,
        "cancel took {cancel_ms} ms against a {budget} ms budget"
    );
    assert_eq!(cancel["runaway_outcome"], "gate_failure");
    assert_eq!(cancel["runaway_step"], "cancelled");
    assert_eq!(cancel["restarted"], true);
    assert_eq!(cancel["recovered"], true);
    eprintln!(
        "[pk411] cancel {cancel_ms:.1} ms, worker restart {restart_ms:.0} ms (budget {budget} ms)"
    );
}

#[test]
fn page_serves_the_boss_plate_in_headless_chrome_when_a_browser_is_present() {
    if !page().join("pkg/popkern_live.js").exists()
        || Command::new("node").arg("--version").output().is_err()
    {
        eprintln!("[pk411] SKIP: web bundle or node absent");
        return;
    }
    let shot = root().join("target/live/boss-plate-headless.png");
    let _ = std::fs::create_dir_all(shot.parent().expect("dir"));
    let output = Command::new("node")
        .arg(page().join("test/browser.mjs"))
        .arg(&shot)
        .output()
        .expect("node");
    if output.status.code() == Some(3) {
        eprintln!("[pk411] SKIP: no headless browser on this machine");
        return;
    }
    let out: Value = serde_json::from_slice(&output.stdout).unwrap_or_else(|e| {
        panic!(
            "browser harness output: {e}\n{}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        )
    });
    assert!(out["error"].is_null(), "{out}");
    let status = out["status_after_first_run"].as_str().unwrap_or("");
    // A baked bundle shows the default script pre-built: a refresh
    // never waits on the kernel.
    let expected = if out["bakes_present"] == true {
        "served (pre-built): 6 feature(s), 1 body(ies)"
    } else {
        "served: 6 feature(s), 1 body(ies)"
    };
    assert!(status.starts_with(expected), "{status}");
    if out["first_run_live"] == false {
        assert_eq!(out["live_run"]["outcome"], "served", "{out}");
    }
    assert!(status.contains("watertight"), "{status}");
    assert_eq!(out["webgl"], true, "WebGL context created");
    assert_eq!(out["bodies_on_screen"], 1);
    assert_eq!(
        out["dialog_hidden"], true,
        "the code dialog is closed on load"
    );
    assert_eq!(
        out["overlay_hidden"], true,
        "no error overlay after a served run"
    );
    assert_eq!(out["report_rows"], 6);
    assert!(out["param_controls"].as_u64().unwrap_or(0) >= 8, "{out}");
    assert!(
        out["pick_text"]
            .as_str()
            .is_some_and(|t| t.contains("rounds/")),
        "face pick names a stable entity: {}",
        out["pick_text"]
    );
    assert_eq!(out["set_param"]["outcome"], "served");
    assert_eq!(out["set_param"]["value"], 0.5);
    page_rows::assert_scrubber(&out);
    page_rows::assert_migrations(&out);
    page_rows::assert_assembly_history(&out);
    page_rows::assert_insert_here(&out);
    page_rows::assert_folder(&out);
    page_rows::assert_part_insert(&out);
    page_rows::assert_time_limit_and_bakes(&out);
    page_rows::assert_assembly(&out);
    page_rows::assert_appearance(&out);
    page_rows::assert_arm_rows(&out);
    page_rows::assert_open_step(&out);
    page_rows::assert_view_rows(&out);
    sim_rows::assert_simulate(&out);
    // Broken script: overlay with class + position, line marked, last good model kept.
    assert_eq!(out["broken"]["outcome"], "script_error");
    assert_eq!(out["broken"]["kind"], "SyntaxError");
    assert_eq!(out["broken"]["overlay_visible"], true);
    assert_eq!(out["broken"]["overlay_class"], "script error");
    assert_eq!(out["broken"]["marked_line"], 2);
    assert_eq!(out["broken"]["bodies_still_on_screen"], 1);
    assert!(
        out["broken"]["status"]
            .as_str()
            .is_some_and(|s| s.contains("last good model")),
        "{}",
        out["broken"]["status"]
    );
    // Kernel declines: class, position and remediation in the overlay;
    // a generative design's re-solve (PK-956).
    page_rows::assert_overlay_rows(&out);
    // Blend rows state measured continuity and how it was measured
    // (sampled while certification is pending, certified, or measured
    // where a certificate could not decide); a downgrade is flagged.
    assert!(
        out["report_rows_text"]
            .as_array()
            .is_some_and(|rows| rows.iter().any(|r| r.as_str().is_some_and(|t| [
                "G1 sampled",
                "G1 certified",
                "G1 measured"
            ]
            .iter()
            .any(|how| t.contains(how))))),
        "{}",
        out["report_rows_text"]
    );
    // PK-561: the served fillet certifies after render and its row
    // upgrades to the certified claim.
    assert!(
        out["status_after_certify"]
            .as_str()
            .is_some_and(|s| s.contains("continuity certified")),
        "{}",
        out["status_after_certify"]
    );
    assert!(
        out["certified_rows_text"]
            .as_array()
            .is_some_and(|rows| rows
                .iter()
                .any(|r| r.as_str().is_some_and(|t| t.contains("G1 certified")))),
        "{}",
        out["certified_rows_text"]
    );
    assert_eq!(out["downgraded_row"]["css"], "cell-warn");
    assert!(
        out["downgraded_row"]["text"]
            .as_str()
            .is_some_and(|t| t.contains("downgraded") && t.contains("G0") && t.contains("G2")),
        "{}",
        out["downgraded_row"]
    );
    // Exports and show-as-code answer through the page's host: the
    // STEP passes the OCC oracle, the .pkn reopens and replays clean,
    // and the emitted script re-runs to the same mesh.
    assert!(out["export_step_bytes"].as_i64().unwrap_or(0) > 1000);
    assert!(out["export_pkn_bytes"].as_i64().unwrap_or(0) > 100);
    assert!(out["show_code_lines"].as_i64().unwrap_or(0) >= 5);
    let step = base64_decode(out["export_step_b64"].as_str().expect("step b64"));
    let pkn = base64_decode(out["export_pkn_b64"].as_str().expect("pkn b64"));
    let step_path = root().join("target/live/page-download.step");
    std::fs::write(&step_path, &step).expect("write step");
    let oracle = Path::new(&std::env::var("HOME").unwrap_or_default()).join("venvs/ocp/bin/python");
    let check = root().join("tools/occ-oracle/occ_check.py");
    if oracle.exists() && check.exists() {
        let output = Command::new(&oracle)
            .arg(&check)
            .arg(&step_path)
            .output()
            .expect("oracle");
        let report = String::from_utf8_lossy(&output.stdout);
        assert!(
            report.contains("\"passed\": true"),
            "OCC oracle on the page's STEP: {report}"
        );
    } else {
        eprintln!("[pk411] OCC oracle absent; page STEP checked by size only");
    }
    let mut reopened = popkern_script::Model::from_pkn_bytes(&pkn).expect("reopen the page's .pkn");
    assert_eq!(reopened.replay().expect("replay").failures().count(), 0);
    assert_eq!(reopened.features().len(), 6);
    let code = out["show_code"].as_str().expect("code");
    let mut session = Session::new(HostOptions::default());
    let boss = session.run(
        &read(root().join("crates/popkern-rhai/examples/boss_plate.rhai")),
        "boss.rhai",
    );
    let mut again = Session::new(HostOptions::default());
    let rerun = again.run(code, "emitted.rhai");
    assert!(rerun.outcome.is_served(), "{:?}", rerun.outcome);
    assert_eq!(
        rerun.bodies[0].triangle_count, boss.bodies[0].triangle_count,
        "show-as-code re-runs to the same mesh"
    );
    assert_eq!(
        fnv1a_u32(&rerun.bodies[0].indices),
        fnv1a_u32(&boss.bodies[0].indices)
    );
    // Cancel stops a runaway and the page recovers.
    assert_eq!(out["cancel"]["restarts"], 1);
    assert_eq!(out["after_cancel"], "served");
    // Nothing left the origin; nothing threw.
    assert_eq!(out["foreign_requests"], serde_json::json!([]));
    assert_eq!(out["exceptions"], serde_json::json!([]));
    eprintln!(
        "[pk411] headless chrome: kernel ready {} ms, first render {} ms, cancel {} ms, history revisit {} ms; screenshot {}",
        out["kernel_ready_ms"],
        out["first_render_ms"],
        out["cancel"]["elapsed_ms"],
        out["scrubber"]["revisit_ms"],
        shot.display()
    );
}

#[test]
fn page_scripts_reach_the_network_only_for_same_origin_assets() {
    // Every page module (each top-level `.js` file; `vendor/` is
    // pinned by hash below), so a new module is vetted without being
    // listed here.
    let mut own: Vec<String> = std::fs::read_dir(page())
        .expect("page dir")
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.file_name().to_string_lossy().into_owned())
        .filter(|name| name.ends_with(".js"))
        .collect();
    own.sort();
    for core in ["app.js", "local-store.js", "worker.js", "worker-host.js"] {
        assert!(own.iter().any(|f| f == core), "{core} missing: {own:?}");
    }
    for file in &own {
        let source = read(page().join(file));
        let mut violations = page_rows::network_violations(&source);
        if file == "socket-worker.js" {
            // The self-hosted session socket (PK-665): same origin only,
            // and only when `popkern serve`'s host.json asks for it.
            assert!(
                source.contains(r#"new URL("./api/session", location.href)"#),
                "socket-worker.js must build its URL from the page's own origin"
            );
            assert!(
                !source.contains("ws://") && !source.contains("wss://"),
                "no absolute socket URLs"
            );
            violations.retain(|v| !v.ends_with(": WebSocket"));
        }
        assert!(violations.is_empty(), "{file}: {violations:?}");
    }
    let html = read(page().join("index.html"));
    let csp = html
        .lines()
        .find(|l| l.contains("Content-Security-Policy"))
        .expect("CSP meta present");
    for directive in [
        "default-src 'self'",
        "script-src 'self' 'wasm-unsafe-eval'",
        "worker-src 'self'",
        "connect-src 'self'",
        "style-src 'self'",
        "object-src 'none'",
        "base-uri 'self'",
        "form-action 'none'",
    ] {
        assert!(csp.contains(directive), "CSP lacks `{directive}`: {csp}");
    }
    assert!(
        !html.contains("https://"),
        "the page references no remote resource"
    );
    assert!(html.contains("<meta name=\"referrer\" content=\"no-referrer\">"));
    // Every script tag is a same-origin module.
    for line in html.lines().filter(|l| l.contains("<script")) {
        assert!(
            line.contains("src=\"./"),
            "script tag not same-origin: {line}"
        );
    }
    // Drafts persist only under try/catch and only locally.
    page_rows::assert_storage_owned(&own, |file| read(page().join(file)));
}

#[test]
fn vendored_files_match_their_recorded_hashes() {
    let licenses = read(page().join("vendor/LICENSES.md"));
    let mut checked = 0;
    for line in licenses.lines().filter(|l| l.starts_with("| `")) {
        let cells: Vec<&str> = line.split('|').map(str::trim).collect();
        let file = cells[1].trim_matches('`');
        let hash = cells[5].trim_matches('`');
        let bytes = std::fs::read(page().join("vendor").join(file))
            .unwrap_or_else(|e| panic!("vendor/{file}: {e}"));
        assert_eq!(sha256_hex(&bytes), hash, "vendor/{file} hash");
        checked += 1;
    }
    assert!(checked >= 1, "LICENSES.md lists the vendored files");
    // Every file in vendor/ is listed.
    for entry in std::fs::read_dir(page().join("vendor")).expect("vendor dir") {
        let name = entry
            .expect("entry")
            .file_name()
            .to_string_lossy()
            .into_owned();
        if name == "LICENSES.md" {
            continue;
        }
        assert!(
            licenses.contains(&format!("`{name}`")),
            "vendor/{name} is not licensed"
        );
    }
}

#[test]
fn protocol_document_names_every_request_the_worker_handles() {
    let protocol = read(page().join("PROTOCOL.md"));
    let core = read(page().join("worker-core.js"));
    for verb in [
        "init",
        "run",
        "set_param",
        "set_joints",
        "to_rhai",
        "export_step",
        "export_pkn",
    ] {
        assert!(
            protocol.contains(&format!("| `{verb}` |")),
            "PROTOCOL.md lacks {verb}"
        );
        assert!(
            core.contains(&format!("\"{verb}\"")),
            "worker-core.js does not handle {verb}"
        );
    }
    assert!(protocol.contains("version 1"));
    assert!(core.contains("KNOWN_PROTOCOLS = [\"1\"]"));
}
