//! PK-662: a project exported from popkern live as a .zip unzips with
//! standard tools (Python's zipfile, and `unzip` when installed) and
//! re-imports losslessly; archives of the folder itself are stripped to
//! the project root; a compressed archive is refused with a reason.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;

use serde_json::Value;

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn have(tool: &str) -> bool {
    Command::new(tool).arg("--version").output().is_ok()
}

fn node(args: &[&str]) -> Value {
    let out = Command::new("node")
        .arg(root().join("crates/popkern-live/tests/node/zip_roundtrip.mjs"))
        .args(args)
        .output()
        .expect("node");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    if out.stdout.is_empty() {
        Value::Null
    } else {
        serde_json::from_slice(&out.stdout).expect("json")
    }
}

fn folder(dir: &Path, from: &Path, out: &mut BTreeMap<String, Vec<u8>>) {
    for entry in std::fs::read_dir(dir).expect("dir") {
        let path = entry.expect("entry").path();
        if path.is_dir() {
            folder(&path, from, out);
        } else {
            let key = path
                .strip_prefix(from)
                .expect("inside")
                .to_string_lossy()
                .replace('\\', "/");
            out.insert(key, std::fs::read(&path).expect("read"));
        }
    }
}

fn decoded(v: &Value) -> BTreeMap<String, Vec<u8>> {
    v["files"]
        .as_object()
        .expect("files")
        .iter()
        .map(|(k, b64)| {
            let bytes = popkern_live::session_project::base64_decode(b64.as_str().expect("b64"))
                .expect("b64");
            (k.clone(), bytes)
        })
        .collect()
}

fn python(script: &str) -> std::process::Output {
    Command::new("python3")
        .arg("-c")
        .arg(script)
        .output()
        .expect("python3")
}

#[test]
fn exported_projects_open_in_standard_tools_and_reimport_losslessly() {
    if !have("node") || !have("python3") {
        eprintln!("[pk662] SKIP: node or python3 missing");
        return;
    }
    let project = root().join("crates/popkern-project/gallery/boss_plate_in_steps");
    let zip = Path::new(env!("CARGO_TARGET_TMPDIR")).join("project.zip");
    node(&[
        "write",
        project.to_str().expect("utf8"),
        zip.to_str().expect("utf8"),
    ]);

    // Python's zipfile: CRCs check, and the names and bytes match.
    let check = python(&format!(
        "import zipfile,sys,json\nz=zipfile.ZipFile({:?})\nassert z.testzip() is None\nprint(json.dumps(sorted(z.namelist())))",
        zip.display().to_string()
    ));
    assert!(
        check.status.success(),
        "{}",
        String::from_utf8_lossy(&check.stderr)
    );
    let mut original = BTreeMap::new();
    folder(&project, &project, &mut original);
    let names: Vec<String> = serde_json::from_slice(&check.stdout).expect("names");
    assert_eq!(names, original.keys().cloned().collect::<Vec<_>>());
    if have("unzip") {
        let t = Command::new("unzip")
            .arg("-t")
            .arg(&zip)
            .output()
            .expect("unzip");
        assert!(
            t.status.success(),
            "unzip -t: {}",
            String::from_utf8_lossy(&t.stdout)
        );
    }

    // Back through the page's reader: identical files.
    let back = node(&["read", zip.to_str().expect("utf8")]);
    assert_eq!(back["ok"], true, "{back}");
    assert_eq!(decoded(&back), original);

    // An OS "compress folder" archive (STORED, one top folder) imports
    // with the folder stripped; a deflated one is refused.
    let top = Path::new(env!("CARGO_TARGET_TMPDIR")).join("folder.zip");
    let deflated = Path::new(env!("CARGO_TARGET_TMPDIR")).join("deflated.zip");
    let make = python(&format!(
        "import zipfile,os\nsrc={:?}\nfor out,mode in (({:?},zipfile.ZIP_STORED),({:?},zipfile.ZIP_DEFLATED)):\n  z=zipfile.ZipFile(out,'w',mode)\n  for d,_,fs in os.walk(src):\n    for f in fs:\n      p=os.path.join(d,f); z.write(p, os.path.join('my-bracket', os.path.relpath(p, src)))\n  z.close()",
        project.display().to_string(),
        top.display().to_string(),
        deflated.display().to_string()
    ));
    assert!(
        make.status.success(),
        "{}",
        String::from_utf8_lossy(&make.stderr)
    );
    let stripped = node(&["read", top.to_str().expect("utf8")]);
    assert_eq!(decoded(&stripped), original, "top folder stripped");
    let refused = node(&["read", deflated.to_str().expect("utf8")]);
    assert_eq!(refused["ok"], false);
    assert!(
        refused["error"]
            .as_str()
            .unwrap_or("")
            .contains("compressed"),
        "{refused}"
    );
}
