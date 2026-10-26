//! Bake the gallery: run the page's default script
//! (`web/live/default.rhai`) and every example in
//! `web/live/examples.json` natively with no budgets, certify every deferred seam, and write
//! `web/live/bakes/<stem>.bin` plus `bakes/index.json` (see
//! `popkern_live::bake`). An existing bake is reused when the
//! index records the same kernel identity and source key.
//!
//! `cargo run --profile live -p popkern-live --example bake_gallery -- <kernel-id> [jobs]`
//!
//! `<kernel-id>` is the wasm bundle's SHA-256 (tools/bake-live.sh
//! passes it); bakes made under another identity are discarded.

use std::path::Path;
use std::sync::Mutex;
use std::time::Duration;

use popkern_core::time::Stopwatch;
use popkern_live::bake::{self, BakeTimings, FORMAT};
use popkern_live::{ContinuityCertification, HostOptions, Session};
use serde_json::{Value, json};

fn ms(watch: &Stopwatch) -> f64 {
    watch.elapsed().map_or(f64::NAN, |d| d.as_secs_f64() * 1e3)
}

/// Run, certify and encode one example; the index record (with
/// `bake` set on success, `declined` otherwise).
fn bake_one(out: &Path, file: &str, source: &str) -> Value {
    let stem = file.trim_end_matches(".rhai");
    let key = bake::source_key(source);
    let mut session = Session::new(HostOptions {
        max_operations: None,
        timeout: None,
        continuity_certification: ContinuityCertification::Deferred,
        track_sources: false,
    });
    let run = Stopwatch::start();
    // The designs this run solves go into its bake (PK-956).
    let _ = popkern_generative::export_touched();
    let envelope = session.run(source, file);
    let generative = popkern_generative::export_touched();
    let run_ms = ms(&run);
    if !envelope.outcome.is_served() {
        let outcome = serde_json::to_value(&envelope.outcome).unwrap_or(Value::Null);
        return json!({ "file": file, "source_key": key, "declined": outcome });
    }
    let certify = Stopwatch::start();
    let report = loop {
        let step = session.certify(Duration::from_secs(3600));
        if let Some(failure) = step.failure {
            break Err(serde_json::to_value(failure).unwrap_or(Value::Null));
        }
        if step.remaining == 0 {
            break Ok(step
                .report
                .or_else(|| envelope.report.clone())
                .unwrap_or(Value::Null));
        }
    };
    let certify_ms = ms(&certify);
    let report = match report {
        Ok(report) => report,
        Err(failure) => {
            return json!({ "file": file, "source_key": key, "declined": { "certify": failure } });
        }
    };
    let timings = BakeTimings { run_ms, certify_ms };
    let bytes = match bake::encode(source, &envelope, report, timings, generative.as_deref()) {
        Ok(bytes) => bytes,
        Err(error) => {
            return json!({ "file": file, "source_key": key, "declined": error.to_string() });
        }
    };
    let name = format!("{stem}.bin");
    if let Err(error) = std::fs::write(out.join(&name), &bytes) {
        return json!({ "file": file, "source_key": key, "declined": error.to_string() });
    }
    json!({
        "file": file,
        "source_key": key,
        "bake": name,
        "bytes": bytes.len(),
        "run_ms": run_ms,
        "certify_ms": certify_ms,
    })
}

fn main() {
    let mut args = std::env::args().skip(1);
    let Some(kernel) = args.next() else {
        eprintln!("usage: bake_gallery <kernel-id> [jobs]");
        std::process::exit(2);
    };
    let jobs: usize = args.next().and_then(|j| j.parse().ok()).unwrap_or_else(|| {
        std::thread::available_parallelism()
            .map_or(4, |n| n.get())
            .min(8)
    });
    let page = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../web/live");
    let dir = page.join("examples");
    let out = page.join("bakes");
    if let Err(error) = std::fs::create_dir_all(&out) {
        eprintln!("cannot create {}: {error}", out.display());
        std::process::exit(1);
    }
    let manifest: Value = std::fs::read_to_string(page.join("examples.json"))
        .ok()
        .and_then(|t| serde_json::from_str(&t).ok())
        .unwrap_or(Value::Null);
    // The page's default script, then the gallery: (index name, path).
    let mut files: Vec<(String, std::path::PathBuf)> =
        vec![("default.rhai".to_string(), page.join("default.rhai"))];
    if let Some(entries) = manifest["examples"].as_array() {
        files.extend(
            entries
                .iter()
                .filter_map(|e| e["file"].as_str())
                .map(|f| (f.to_string(), dir.join(f))),
        );
    }
    // Reuse bakes from the same kernel whose source is unchanged.
    let previous: Value = std::fs::read_to_string(out.join("index.json"))
        .ok()
        .and_then(|t| serde_json::from_str(&t).ok())
        .unwrap_or(Value::Null);
    let same_kernel = previous["kernel"].as_str() == Some(kernel.as_str())
        && previous["format"].as_u64() == Some(u64::from(FORMAT));
    let reusable = |file: &str, key: &str| -> Option<Value> {
        if !same_kernel {
            return None;
        }
        let record = previous["entries"]
            .as_array()?
            .iter()
            .find(|r| r["file"] == file && r["source_key"] == key)?;
        match record["bake"].as_str() {
            Some(name) => out.join(name).exists().then(|| record.clone()),
            None => Some(record.clone()),
        }
    };
    let queue = Mutex::new(files.iter().enumerate().collect::<Vec<_>>());
    let results = Mutex::new(Vec::<(usize, Value)>::new());
    let total = Stopwatch::start();
    std::thread::scope(|scope| {
        for _ in 0..jobs.max(1) {
            scope.spawn(|| {
                loop {
                    let next = queue.lock().ok().and_then(|mut q| q.pop());
                    let Some((position, (file, path))) = next else {
                        break;
                    };
                    let Ok(source) = std::fs::read_to_string(path) else {
                        eprintln!("{file}: unreadable");
                        continue;
                    };
                    let key = bake::source_key(&source);
                    let record = match reusable(file, &key) {
                        Some(record) => {
                            eprintln!("{file}: reused");
                            record
                        }
                        None => {
                            let record = bake_one(&out, file, &source);
                            match record.get("bake") {
                                Some(_) => eprintln!(
                                    "{file}: baked in {:.1} s + certified in {:.1} s, {:.1} MiB",
                                    record["run_ms"].as_f64().unwrap_or(0.0) / 1e3,
                                    record["certify_ms"].as_f64().unwrap_or(0.0) / 1e3,
                                    record["bytes"].as_f64().unwrap_or(0.0) / 1048576.0
                                ),
                                None => eprintln!("{file}: NOT baked: {}", record["declined"]),
                            }
                            record
                        }
                    };
                    if let Ok(mut results) = results.lock() {
                        results.push((position, record));
                    }
                }
            });
        }
    });
    let mut results = results.into_inner().unwrap_or_default();
    results.sort_by_key(|(position, _)| *position);
    // Bakes no entry refers to any more (a removed example, a kernel
    // change) are deleted, so the directory holds exactly the index.
    let keep: Vec<String> = results
        .iter()
        .filter_map(|(_, record)| record["bake"].as_str().map(str::to_string))
        .collect();
    if let Ok(entries) = std::fs::read_dir(&out) {
        for path in entries.filter_map(|e| e.ok().map(|e| e.path())) {
            let name = path.file_name().and_then(|n| n.to_str()).unwrap_or("");
            if name.ends_with(".bin") && !keep.iter().any(|k| k == name) {
                let _ = std::fs::remove_file(&path);
            }
        }
    }
    let baked = keep.len();
    let index = json!({
        "format": FORMAT,
        "kernel": kernel,
        "entries": results.into_iter().map(|(_, record)| record).collect::<Vec<_>>(),
    });
    let text = serde_json::to_string_pretty(&index).unwrap_or_default();
    if let Err(error) = std::fs::write(out.join("index.json"), text) {
        eprintln!("cannot write the index: {error}");
        std::process::exit(1);
    }
    eprintln!(
        "{baked}/{} examples baked in {:.0} s",
        files.len(),
        ms(&total) / 1e3
    );
}
