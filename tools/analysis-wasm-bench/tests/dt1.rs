//! DT1, native half: the three bench cases' result hashes are pinned
//! here, and `run.sh` checks that the wasm build reproduces them in
//! node. `export_cases` (ignored) writes the case files for the
//! runners.

use std::path::{Path, PathBuf};

use popkern_analysis::element::vector_boundary::PieceSupport;
use popkern_analysis::model::{AnalysisModel, Surface};
use popkern_analysis::structural::assemble::FaceMechanics;
use popkern_analysis_wasm_bench::{encode, run};

fn repo(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
}

fn plane(m: &AnalysisModel, n: [f64; 3], d: f64) -> Vec<usize> {
    m.faces
        .iter()
        .enumerate()
        .filter(|(_, f)| match f.surface {
            Surface::Plane { normal, offset } => {
                (0..3).all(|k| (normal[k] - n[k]).abs() < 1e-9) && (offset - d).abs() < 1e-9
            }
            _ => false,
        })
        .map(|(i, _)| i)
        .collect()
}

/// A plane `n · x = d` as (n, d).
type Plane = ([f64; 3], f64);

/// The three cases: (name, bytes). Clamped on one plane face, loaded
/// by a traction on another.
fn cases() -> Vec<(&'static str, Vec<u8>)> {
    let setups: [(&str, Plane, Plane, [f64; 3]); 3] = [
        (
            "block_two_bores",
            ([0.0, 0.0, -1.0], 0.0),
            ([0.0, 0.0, 1.0], 10.0),
            [1.0, 0.0, 0.0],
        ),
        (
            "block_dimple",
            ([0.0, 0.0, -1.0], 0.0),
            ([1.0, 0.0, 0.0], 10.0),
            [0.0, 0.0, -1.0],
        ),
        (
            "l_bracket",
            ([1.0, 0.0, 0.0], 40.0),
            ([0.0, 1.0, 0.0], 30.0),
            [1.0, 0.0, 0.0],
        ),
    ];
    setups
        .iter()
        .map(|(name, clamp, load, t)| {
            let m = popkern_analysis_cli::step_model::model_from_step_file(&repo(&format!(
                "tests/fixtures/analysis/{name}.step"
            )))
            .unwrap_or_else(|e| panic!("{name}: {e}"));
            let mut faces = vec![FaceMechanics::Free; m.faces.len()];
            let clamped = plane(&m, clamp.0, clamp.1);
            let loaded = plane(&m, load.0, load.1);
            assert!(!clamped.is_empty() && !loaded.is_empty(), "{name}");
            for f in clamped {
                faces[f] = FaceMechanics::Support(PieceSupport::Fixed {
                    u0: [0.0; 3],
                    gradient: [[0.0; 3]; 3],
                });
            }
            for f in loaded {
                faces[f] = FaceMechanics::Traction(*t);
            }
            (*name, encode(&m, &faces, 2000.0, 0.35))
        })
        .collect()
}

/// DT1 budget: small enough to run in every gate pass.
const DT1_CELLS: usize = 4000;

/// Native hashes at `DT1_CELLS`, which the wasm build must reproduce.
/// They move with the kernel's STEP import and tessellation, which feed
/// the model the solve runs on: re-pinned 2026-10-25 at the move into
/// popkern (wasm reproduced all three; block_two_bores kept its 3,165
/// cells, 11,808 dofs and 11 iterations, block_dimple's sphere sense was
/// corrected, docs/architecture/analysis-integration.md fixture log).
const PINNED: [(&str, u64); 3] = [
    ("block_two_bores", 0xe74b_1c78_84ae_bba3),
    ("block_dimple", 0x7785_3ae9_0b50_5a5e),
    ("l_bracket", 0xb952_afcf_f3d8_f6a5),
];

#[test]
fn native_hashes_are_pinned() {
    for ((name, bytes), (pinned_name, pinned)) in cases().iter().zip(PINNED) {
        assert_eq!(*name, pinned_name);
        let o = run(bytes, DT1_CELLS).unwrap_or_else(|e| panic!("{name}: {e}"));
        // A second native run gives the same bits.
        assert_eq!(run(bytes, DT1_CELLS), Ok(o), "{name}");
        assert_eq!(o.hash, pinned, "{name}: {o:?}");
    }
}

#[test]
#[ignore]
fn export_cases() {
    let dir = repo("target/wasm-bench");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let mut native = String::new();
    for (name, bytes) in cases() {
        std::fs::write(dir.join(format!("{name}.bin")), &bytes).unwrap_or_else(|e| panic!("{e}"));
        let o = run(&bytes, DT1_CELLS).unwrap_or_else(|e| panic!("{name}: {e}"));
        native += &format!(
            "{name} {:016x} {} {} {}\n",
            o.hash, o.cells, o.dof, o.iterations
        );
    }
    std::fs::write(dir.join("native.txt"), native).unwrap_or_else(|e| panic!("{e}"));
}

#[test]
#[ignore]
fn native_e3_timing() {
    let (name, bytes) = cases().swap_remove(0);
    // The browser's numbers come from run.sh --e3; these are the same
    // solves natively, for the wasm/native ratio.
    for cells in [15_000, 20_000, 30_000] {
        let t = std::time::Instant::now();
        let o = run(&bytes, cells).unwrap_or_else(|e| panic!("{name}: {e}"));
        eprintln!(
            "native {name} budget {cells}: {:.2} s {o:?}",
            t.elapsed().as_secs_f64()
        );
    }
}
