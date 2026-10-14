//! Display-mesh facet-tilt ratchet over the gallery scripts that carry
//! leaning facets (PK-687..PK-690, PK-809).
//!
//! The contract (docs/architecture/display-tessellation-in-face-folds.md):
//! every display facet's normal lies within 12° of the surface normal at
//! each of its vertices, chart-singular vertices excepted. Each row pins
//! the current count of facets tilted past 12° and the worst tilt as
//! ceilings: a change that adds leaning facets fails here, and a task
//! that removes them tightens its row. In-face folds (the dihedral
//! between adjacent facets of one face, welded at 1e-9) are printed as a
//! quality metric: neighbours that each honour the tilt bound can meet
//! at up to about twice it.

use popkern_core::geometry::{Point3, Vec3};
use std::collections::HashMap;

/// `(script, body, max facets tilted over 12°, max worst tilt in degrees)`.
const CEILINGS: &[(&str, usize, usize, f64)] = &[
    // E1 (cross-face boundary caps) then E4 steps (3D Delaunay flips,
    // shear-aware chart Delaunay, obtuse-facet bisection): from the
    // 2026-10-13 baseline 0/1/2/1/188/48/124/502/10/19.
    ("angled_fitting", 1, 0, 11.5),
    ("countersunk_joint", 1, 0, 12.0),
    ("hex_bolt", 0, 0, 12.0),
    ("threaded_stud", 0, 0, 12.0),
    ("coil_spring", 0, 6, 14.4),
    ("coil_spring", 1, 0, 12.0),
    ("coil_spring", 2, 0, 12.0),
    ("propeller", 0, 0, 12.0),
    ("mixed_radius_corner", 0, 0, 12.0),
    ("pyramid_apex", 0, 0, 12.0),
];

/// Edges of the globally welded display mesh (bit-identical positions)
/// used by one triangle (open) or more than two (non-manifold).
fn welded_defects(mesh: &popkern_mesh::TriMesh) -> (usize, usize) {
    let mut weld: HashMap<[u64; 3], usize> = HashMap::new();
    let ids: Vec<usize> = mesh
        .vertices
        .iter()
        .map(|v| {
            let n = weld.len();
            let p = v.position;
            *weld
                .entry([p.x.to_bits(), p.y.to_bits(), p.z.to_bits()])
                .or_insert(n)
        })
        .collect();
    let mut uses: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        for k in 0..3 {
            let (a, b) = (ids[t[k]], ids[t[(k + 1) % 3]]);
            if a != b {
                *uses.entry((a.min(b), a.max(b))).or_default() += 1;
            }
        }
    }
    (
        uses.values().filter(|&&c| c == 1).count(),
        uses.values().filter(|&&c| c > 2).count(),
    )
}

/// Facet tilt bound of the display contract, in degrees.
const TILT_BOUND: f64 = 12.0;

/// `(facets tilted over the bound, worst tilt)`: each facet's largest
/// angle to the surface normal at one of its vertices, skipping facets
/// with a chart-singular vertex (no unique normal) or zero area.
fn facet_tilts(brep: &popkern_core::topology::BRep, mesh: &popkern_mesh::TriMesh) -> (usize, f64) {
    use popkern_core::geometry::SurfaceEval;
    let (mut over, mut worst) = (0usize, 0.0f64);
    for fm in &mesh.face_meshes {
        let Some(surface) = brep.faces[fm.face].surface.as_ref().map(|s| &s.geometry) else {
            continue;
        };
        for t in &fm.triangles {
            let n = facet_normal(t.map(|i| fm.vertices[i].position));
            if n.length() < 1e-300 {
                continue;
            }
            let mut tilt = 0.0f64;
            let mut singular = false;
            for &i in t {
                let v = &fm.vertices[i];
                let (su, sv) = (
                    surface.partial_u(v.uv.0, v.uv.1),
                    surface.partial_v(v.uv.0, v.uv.1),
                );
                if su.cross(&sv).length() < 1e-9 * (su.length() * sv.length()).max(1e-300) {
                    singular = true;
                    break;
                }
                let c = n.dot(&v.normal) / (n.length() * v.normal.length());
                tilt = tilt.max(c.clamp(-1.0, 1.0).acos().to_degrees());
            }
            if singular {
                continue;
            }
            if tilt > TILT_BOUND {
                over += 1;
            }
            worst = worst.max(tilt);
        }
    }
    (over, worst)
}

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
fn gallery_display_facet_tilts_stay_under_their_ceilings() {
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
            let (folds, worst_fold) = in_face_folds(&mesh);
            let (over, worst) = facet_tilts(parts[body].brep(), &mesh);
            eprintln!(
                "{stem}#{body}: facets tilted over 12° {over} worst {worst:.1}° (folds over 12.5° {folds}, worst {worst_fold:.1}°)"
            );
            let (open, nonmanifold) = welded_defects(&mesh);
            if open + nonmanifold > 0 {
                failures.push(format!(
                    "{stem}#{body}: welded display mesh has {open} open and {nonmanifold} non-manifold edges"
                ));
            }
            if over > max_over || worst > max_worst {
                failures.push(format!(
                    "{stem}#{body}: {over} tilted facets (ceiling {max_over}), worst {worst:.1}° (ceiling {max_worst})"
                ));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
