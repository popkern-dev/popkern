//! Shared oracle for multi-body unions with no closed-form volume
//! (PK-758, PK-753): popkern's served body against OpenCASCADE's OWN
//! fuse of the same operands (`tools/occ-oracle/occ_fuse.py`), the two
//! kernels sharing no boolean code.

use popkern_rhai::ScriptHost;

/// Run `body` after `prelude` in a fresh model and return the model.
pub fn run(prelude: &str, body: &str, name: &str) -> Result<popkern_script::Model, String> {
    let source = format!("{prelude}\nlet m = model();\n{body}\nm\n");
    let mut host = ScriptHost::new();
    let outcome = host
        .run(&source, name)
        .map_err(|e| format!("{name}: {}", e.render()))?;
    outcome.model.ok_or_else(|| "no model".into())
}

/// OCC's verdict on `served` (one valid solid) and its own fuse of the
/// `operands` (each a body expression over `m`, after `prelude`):
/// `(served volume, fused volume)`, or `None` without the OCP venv.
pub fn occ_volumes(
    prelude: &str,
    operands: &[&str],
    served: &mut popkern_script::Model,
    name: &str,
) -> Option<(f64, f64)> {
    let home = std::path::PathBuf::from(std::env::var("HOME").unwrap_or_default());
    let python = home.join("venvs/ocp/bin/python");
    let tools = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occ-oracle");
    if !python.exists() || !tools.exists() {
        eprintln!("[union-oracle] OCC oracle absent; volumes not cross-checked");
        return None;
    }
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("union_oracle");
    std::fs::create_dir_all(&dir).expect("dir");
    let step = dir.join(format!("{name}.step"));
    served.to_step(&step).expect("STEP export");
    let oracle = |script: &str, args: &[&std::path::Path], extra: &[&str]| -> serde_json::Value {
        let out = std::process::Command::new(&python)
            .arg(tools.join(script))
            .args(args)
            .args(extra)
            .output()
            .expect("oracle");
        serde_json::from_slice(&out.stdout).expect("oracle json")
    };
    let check = oracle("occ_check.py", &[&step], &[]);
    assert_eq!(
        check["passed"], true,
        "{name}: OCC rejects the STEP: {check}"
    );
    let solids = oracle("occ_solids.py", &[&step], &["--adaptive"]);
    let rows = solids["solids"].as_array().expect("solids");
    assert_eq!(rows.len(), 1, "{name}: one solid {solids}");
    let served_volume = rows[0]["volume"].as_f64().expect("volume");
    let mut parts = Vec::new();
    for (k, body) in operands.iter().enumerate() {
        let mut model = run(
            prelude,
            &format!("let d = {body};"),
            &format!("{name}_{k}.rhai"),
        )
        .expect("operand");
        let path = dir.join(format!("{name}_operand_{k}.step"));
        model.to_step(&path).expect("operand STEP");
        parts.push(path);
    }
    let refs: Vec<&std::path::Path> = parts.iter().map(|p| p.as_path()).collect();
    let fused = oracle("occ_fuse.py", &refs, &[]);
    assert_eq!(fused["done"], true, "{name}: OCC fuse {fused}");
    Some((
        served_volume,
        fused["volume"].as_f64().expect("fused volume"),
    ))
}
