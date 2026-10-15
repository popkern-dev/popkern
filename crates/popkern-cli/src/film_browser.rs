//! `popkern film --out`: render a film with the page's own export
//! (web/live/film-export.js) in headless Chrome. The project is served
//! in-process by `popkern serve` (native kernel) on a free loopback
//! port; Chrome is driven over the DevTools protocol (the browser
//! endpoint from its `DevToolsActivePort` file, one attached page
//! target, `Runtime.evaluate`), and the bytes come back as base64.

use std::io::ErrorKind;
use std::net::TcpStream;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::time::Duration;

use popkern_core::time::Stopwatch;

use serde_json::{Value, json};
use tungstenite::stream::MaybeTlsStream;
use tungstenite::{Message, WebSocket};

use crate::args::Cli;
use crate::film::FilmSource;

/// How long Chrome may take to open its DevTools port.
const STARTUP: Duration = Duration::from_secs(30);
/// How long one render may take (the arm's film at 1080p is minutes).
const RENDER: Duration = Duration::from_secs(3600);

/// The Chrome to use: `--chrome`, `POPKERN_CHROME`, or the first of the
/// usual names on PATH.
fn chrome(cli: &Cli) -> Result<PathBuf, String> {
    if let Some(path) = &cli.chrome {
        return Ok(path.clone());
    }
    if let Some(path) = std::env::var_os("POPKERN_CHROME") {
        return Ok(PathBuf::from(path));
    }
    let path = std::env::var_os("PATH").unwrap_or_default();
    for name in [
        "google-chrome",
        "chromium",
        "chromium-browser",
        "google-chrome-stable",
        "chrome-headless-shell",
    ] {
        for dir in std::env::split_paths(&path) {
            let candidate = dir.join(name);
            if candidate.is_file() {
                return Ok(candidate);
            }
        }
    }
    Err("no headless Chrome: pass --chrome <path> or set POPKERN_CHROME".into())
}

/// A scratch folder removed on drop.
struct Scratch(PathBuf);

impl Scratch {
    fn new(label: &str) -> Result<Scratch, String> {
        let stamp = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_nanos())
            .unwrap_or(0);
        let dir = std::env::temp_dir().join(format!(
            "popkern-film-{label}-{}-{stamp}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        Ok(Scratch(dir))
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// Chrome, killed on drop.
struct Browser(Child);

impl Drop for Browser {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}

/// The projects root to serve: always a scratch copy of the project's
/// files. Serving the project's own parent folder let the page write
/// there (its Scratch project appeared beside the user's projects).
fn projects_root(source: &FilmSource, scratch: &Scratch) -> Result<PathBuf, String> {
    let project = scratch.0.join("projects").join(&source.id);
    for (path, bytes) in &source.files {
        let target = project.join(path);
        if let Some(parent) = target.parent() {
            std::fs::create_dir_all(parent).map_err(|e| e.to_string())?;
        }
        std::fs::write(&target, bytes).map_err(|e| e.to_string())?;
    }
    Ok(scratch.0.join("projects"))
}

/// The browser-level DevTools URL Chrome writes into its profile.
fn devtools_url(profile: &Path, child: &mut Child) -> Result<String, String> {
    let started = Deadline::new(STARTUP);
    let file = profile.join("DevToolsActivePort");
    while !started.passed() {
        if let Ok(text) = std::fs::read_to_string(&file) {
            let mut lines = text.lines();
            if let (Some(port), Some(path)) = (lines.next(), lines.next()) {
                return Ok(format!("ws://127.0.0.1:{}{}", port.trim(), path.trim()));
            }
        }
        if let Ok(Some(status)) = child.try_wait() {
            return Err(format!("Chrome exited at startup ({status})"));
        }
        std::thread::sleep(Duration::from_millis(100));
    }
    Err("Chrome did not open its DevTools port".into())
}

/// A wall-clock limit (a `Stopwatch`: `std::time::Instant` is barred
/// from library code, PK-674's wasm rule).
struct Deadline {
    watch: Stopwatch,
    limit: Duration,
}

impl Deadline {
    fn new(limit: Duration) -> Deadline {
        Deadline {
            watch: Stopwatch::start(),
            limit,
        }
    }

    fn passed(&self) -> bool {
        self.watch.elapsed().is_some_and(|e| e > self.limit)
    }
}

/// One DevTools session over the browser socket.
struct Devtools {
    socket: WebSocket<MaybeTlsStream<TcpStream>>,
    next: u64,
}

impl Devtools {
    fn call(
        &mut self,
        session: Option<&str>,
        method: &str,
        params: Value,
        deadline: &Deadline,
    ) -> Result<Value, String> {
        self.next += 1;
        let id = self.next;
        let mut message = json!({ "id": id, "method": method, "params": params });
        if let Some(session) = session {
            message["sessionId"] = json!(session);
        }
        self.socket
            .send(Message::Text(message.to_string()))
            .map_err(|e| format!("DevTools send: {e}"))?;
        loop {
            if deadline.passed() {
                return Err(format!("DevTools {method}: no answer in time"));
            }
            let reply = match self.socket.read() {
                Ok(Message::Text(text)) => text,
                Ok(_) => continue,
                Err(tungstenite::Error::Io(e))
                    if matches!(e.kind(), ErrorKind::WouldBlock | ErrorKind::TimedOut) =>
                {
                    continue;
                }
                Err(e) => return Err(format!("DevTools read: {e}")),
            };
            let value: Value = serde_json::from_str(&reply).map_err(|e| e.to_string())?;
            if value["id"].as_u64() != Some(id) {
                continue;
            }
            if let Some(error) = value.get("error") {
                return Err(format!("DevTools {method}: {error}"));
            }
            return Ok(value["result"].clone());
        }
    }
}

/// The in-page export: open the project, play the film, export it,
/// return the bytes as base64 (JSON text).
fn export_script(id: &str, path: &str, preset: &str, at: f64) -> String {
    format!(
        r#"(async () => {{
  const wait = async (ready, ms, what) => {{
    const t0 = performance.now();
    while (!(await ready())) {{
      if (performance.now() - t0 > ms) throw new Error("timed out waiting for " + what);
      await new Promise((r) => setTimeout(r, 100));
    }}
  }};
  await wait(() => window.popkernLive && window.popkernLive.kernelReadyMs, 120000, "the kernel");
  const p = window.popkernLive;
  p.setTimeLimit(null);
  await p.projects.openStored({id});
  await wait(() => p.runDone && !p.running, 3600000, "the build");
  await p.runDone;
  if (!p.lastGood) throw new Error("the document did not build in the page");
  p.projects.select({path});
  const reply = await p.run();
  if (!reply || !reply.ok || !p.film.active) throw new Error("the film did not run: " + JSON.stringify(reply && reply.outcome));
  await p.film.seek({at});
  const {{ exportFilm }} = await import("./film-export.js");
  const out = await exportFilm(p, {preset});
  let text = "";
  for (let i = 0; i < out.bytes.length; i += 0x8000) text += String.fromCharCode.apply(null, out.bytes.subarray(i, i + 0x8000));
  return JSON.stringify({{ ok: true, base64: btoa(text), frames: out.frames }});
}})().catch((e) => JSON.stringify({{ ok: false, error: String((e && e.stack) || e) }}))"#,
        id = json!(id),
        path = json!(path),
        preset = json!(preset),
        at = json!(at),
    )
}

/// Render film `path` of `source` with the CLI's preset: the file bytes.
pub fn render(cli: &Cli, source: &FilmSource, path: &str) -> Result<Vec<u8>, String> {
    let chrome = chrome(cli)?;
    let scratch = Scratch::new("render")?;
    let root = projects_root(source, &scratch)?;
    let page = cli.page.clone().unwrap_or_else(crate::default_page);
    let options = popkern_server::ServeOptions {
        root,
        page,
        bind: "127.0.0.1:0".parse().map_err(|_| "loopback address")?,
        compute: popkern_server::Compute::Native,
    };
    let (listener, handle) =
        popkern_server::bind(options).map_err(|e| format!("serving the project: {e}"))?;
    let addr = listener.local_addr().map_err(|e| e.to_string())?;
    std::thread::spawn(move || handle.run(listener));
    let profile = scratch.0.join("chrome");
    std::fs::create_dir_all(&profile).map_err(|e| e.to_string())?;
    let child = Command::new(&chrome)
        .args([
            "--headless=new",
            "--no-sandbox",
            "--use-gl=angle",
            "--use-angle=swiftshader",
            "--enable-unsafe-swiftshader",
            "--ignore-gpu-blocklist",
            "--disable-dev-shm-usage",
            "--remote-debugging-port=0",
        ])
        .arg(format!("--user-data-dir={}", profile.display()))
        .arg("about:blank")
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .map_err(|e| format!("{}: {e}", chrome.display()))?;
    let mut browser = Browser(child);
    let url = devtools_url(&profile, &mut browser.0)?;
    let (socket, _) =
        tungstenite::connect(url.as_str()).map_err(|e| format!("DevTools connect: {e}"))?;
    if let MaybeTlsStream::Plain(stream) = socket.get_ref() {
        stream
            .set_read_timeout(Some(Duration::from_secs(5)))
            .map_err(|e| e.to_string())?;
    }
    let mut devtools = Devtools { socket, next: 0 };
    let deadline = &Deadline::new(RENDER);
    let target = devtools.call(
        None,
        "Target.createTarget",
        json!({ "url": format!("http://{addr}/") }),
        deadline,
    )?;
    let target_id = target["targetId"]
        .as_str()
        .ok_or("no page target")?
        .to_string();
    let attached = devtools.call(
        None,
        "Target.attachToTarget",
        json!({ "targetId": target_id, "flatten": true }),
        deadline,
    )?;
    let session = attached["sessionId"]
        .as_str()
        .ok_or("no page session")?
        .to_string();
    // The target starts on about:blank and navigates: wait until the
    // served page has loaded (an evaluate during the navigation finds
    // its context destroyed).
    let loaded = Deadline::new(STARTUP);
    loop {
        let probe = devtools.call(
            Some(&session),
            "Runtime.evaluate",
            json!({ "expression": "location.protocol + ' ' + document.readyState", "returnByValue": true }),
            deadline,
        );
        if probe.is_ok_and(|r| r["result"]["value"] == json!("http: complete")) {
            break;
        }
        if loaded.passed() {
            return Err(format!("the page at http://{addr}/ did not load"));
        }
        std::thread::sleep(Duration::from_millis(200));
    }
    let result = devtools.call(
        Some(&session),
        "Runtime.evaluate",
        json!({ "expression": export_script(&source.id, path, &cli.preset, cli.at), "awaitPromise": true, "returnByValue": true }),
        deadline,
    )?;
    let text = result["result"]["value"]
        .as_str()
        .ok_or_else(|| format!("the page answered {result}"))?;
    let reply: Value = serde_json::from_str(text).map_err(|e| e.to_string())?;
    if reply["ok"] != json!(true) {
        return Err(format!(
            "the page's export failed: {}",
            reply["error"].as_str().unwrap_or("unknown")
        ));
    }
    popkern_live::session_project::base64_decode(reply["base64"].as_str().unwrap_or_default())
}
