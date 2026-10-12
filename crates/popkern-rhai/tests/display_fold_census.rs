//! Display-mesh fold ratchet over the gallery scripts that carry
//! in-face folds (PK-687..PK-690).
//!
//! A fold is the angle between two adjacent facets of ONE face (welded
//! at 1e-9), so it measures the triangulation, not the model's own
//! creases. The display preset bounds vertex-normal spread at 12°;
//! acceptance for the fold tasks is no in-face fold over 12.5° away
//! from chart singularities. Each row pins the current count of folds
//! over 12.5° and the worst fold as ceilings: a change that adds folds
//! fails here, and a task that removes them tightens its row.

use popkern_core::geometry::{Point3, Vec3};
use std::collections::HashMap;

/// `(script, body, max folds over 12.5°, max worst fold in degrees)`.
const CEILINGS: &[(&str, usize, usize, f64)] = &[
    // PK-690: the screw thread's end-of-knot seed made an 82° needle.
    ("angled_fitting", 1, 2, 19.0),
    ("countersunk_joint", 1, 3, 15.0),
    ("hex_bolt", 0, 3, 19.0),
    ("threaded_stud", 0, 1, 18.5),
    // PK-687 (crease flips took the coil from 12228/2736/12879 to 266/116/242).
    ("coil_spring", 0, 266, 49.5),
    ("coil_spring", 1, 116, 25.5),
    ("coil_spring", 2, 242, 23.0),
    // PK-688 (crease flips: 3771 -> 597, worst 84.9° -> 30.8°).
    ("propeller", 0, 597, 31.0),
    // PK-689 (crease flips: 85 -> 1 and 144 -> 15).
    ("mixed_radius_corner", 0, 1, 21.0),
    ("pyramid_apex", 0, 15, 29.0),
];

fn facet_normal(p: [Point3; 3]) -> Vec3 {
    (p[1] - p[0]).cross(&(p[2] - p[0]))
}

/// `(folds over 12.5°, worst fold)` over every face of `mesh`.
fn in_face_folds(mesh: &popkern_mesh::TriMesh) -> (usize, f64) {
    let (mut over, mut worst) = (0usize, 0.0f64);
    for fm in &mesh.face_meshes {
        let mut weld: HashMap<(i64, i64, i64), usize> = HashMap::new();
        let ids: Vec<usize> = fm
            .vertices
            .iter()
            .map(|v| {
                let p = v.position;
                let k = (
                    (p.x * 1e9).round() as i64,
                    (p.y * 1e9).round() as i64,
                    (p.z * 1e9).round() as i64,
                );
                let n = weld.len();
                *weld.entry(k).or_insert(n)
            })
            .collect();
        let mut edges: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for (ti, t) in fm.triangles.iter().enumerate() {
            let w = [ids[t[0]], ids[t[1]], ids[t[2]]];
            if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
                continue;
            }
            for (x, y) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
                edges.entry((x.min(y), x.max(y))).or_default().push(ti);
            }
        }
        let normal = |ti: usize| {
            let t = fm.triangles[ti];
            facet_normal([
                fm.vertices[t[0]].position,
                fm.vertices[t[1]].position,
                fm.vertices[t[2]].position,
            ])
        };
        for tris in edges.values() {
            let [a, b] = tris.as_slice() else { continue };
            let (n0, n1) = (normal(*a), normal(*b));
            if n0.length() < 1e-300 || n1.length() < 1e-300 {
                continue;
            }
            let fold = (n0.dot(&n1) / (n0.length() * n1.length()))
                .clamp(-1.0, 1.0)
                .acos()
                .to_degrees();
            if fold > 12.5 {
                over += 1;
            }
            worst = worst.max(fold);
        }
    }
    (over, worst)
}

#[test]
fn gallery_display_folds_stay_under_their_ceilings() {
    let root = format!("{}/examples", env!("CARGO_MANIFEST_DIR"));
    let mut failures = Vec::new();
    let mut scripts: Vec<&str> = CEILINGS.iter().map(|c| c.0).collect();
    scripts.dedup();
    for stem in scripts {
        let source = std::fs::read_to_string(format!("{root}/{stem}.rhai")).expect("script");
        let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
        let run = host.run(&source, stem).expect("script runs");
        let mut model = run.model.expect("script returns a model");
        let parts = model.to_parts_all().expect("parts");
        for &(_, body, max_over, max_worst) in CEILINGS.iter().filter(|c| c.0 == stem) {
            let mesh = popkern_mesh::tessellate(
                parts[body].brep(),
                &popkern_mesh::TessellationOptions::display(),
            )
            .expect("display mesh");
            let (over, worst) = in_face_folds(&mesh);
            eprintln!("{stem}#{body}: folds over 12.5° {over} worst {worst:.1}°");
            if over > max_over || worst > max_worst {
                failures.push(format!(
                    "{stem}#{body}: {over} folds (ceiling {max_over}), worst {worst:.1}° (ceiling {max_worst})"
                ));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
