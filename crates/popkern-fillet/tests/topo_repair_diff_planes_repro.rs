//! Task #136 first-session repro and diagnostics.
//!
//! Pyramid apex-only fillets leave cylinder faces with chord pairs
//! through an off-cylinder vertex (the cap centre or unfilleted
//! base corner).  `topology_repair`'s same-plane case lands those
//! chord pairs sit on a single planar adjacent face.  The
//! different-planes case (chord pair spans two adjacent planes)
//! is gated behind `POPKERN_TOPO_REPAIR_DIFF=1` and known to
//! regress N≥8 pyramids (max tess-nm divergence 0.12 → 1.06).
//!
//! These tests are diagnostic — they print, they don't assert
//! tight numbers — to characterise the regression's mechanism.
//! All `#[ignore]` so they don't run by default; invoke with
//! `cargo test topo_repair_diff_planes -- --ignored --nocapture`.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[path = "kp_gate_helpers.rs"]
mod helpers;
#[path = "topo_repair_diff_planes_diagnostics.rs"]
mod topo_repair_diff_planes_diagnostics;
use helpers::{incident_edges, pick_vertex_with_valence};

struct EnvGuard {
    key: &'static str,
    prior: Option<String>,
}

impl EnvGuard {
    fn set(key: &'static str, value: &str) -> Self {
        let prior = std::env::var(key).ok();
        unsafe {
            std::env::set_var(key, value);
        }
        Self { key, prior }
    }
}

impl Drop for EnvGuard {
    fn drop(&mut self) {
        unsafe {
            match &self.prior {
                Some(v) => std::env::set_var(self.key, v),
                None => std::env::remove_var(self.key),
            }
        }
    }
}

fn fillet_apex_only(n: usize, height: f64, radius: f64) -> BRep {
    let pyr = BRep::n_gonal_pyramid(n, 1.0, height, Tolerance::default()).unwrap();
    let apex = pick_vertex_with_valence(&pyr, n).expect("apex");
    let edges = incident_edges(&pyr, apex);
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(radius),
        })
        .collect();
    let opts = FilletOptions::default();
    fillet_edges(&pyr, &specs, &opts, Tolerance::default(), None)
        .unwrap()
        .brep
}

fn vertex_on_cylinder(p: &Point3, cyl: &SurfaceGeometry) -> bool {
    let SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    } = cyl
    else {
        return false;
    };
    let rel = popkern_core::geometry::Vec3::new(p.x - origin.x, p.y - origin.y, p.z - origin.z);
    let len = axis.length().max(1e-12);
    let axial = rel.dot(axis) / len;
    let axis_unit = popkern_core::geometry::Vec3::new(axis.x / len, axis.y / len, axis.z / len);
    let radial = popkern_core::geometry::Vec3::new(
        rel.x - axial * axis_unit.x,
        rel.y - axial * axis_unit.y,
        rel.z - axial * axis_unit.z,
    );
    let r = radial.length();
    (r - radius).abs() <= radius * 0.05
}

#[derive(Default)]
struct OffCylSurvey {
    cylinder_faces: usize,
    off_cyl_vertices: usize,
    same_plane_chord_pairs: usize,
    diff_plane_chord_pairs: usize,
}

fn survey_off_cyl(brep: &BRep) -> OffCylSurvey {
    let mut s = OffCylSurvey::default();
    for (fh, face) in brep.faces.iter() {
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        if !matches!(sof.geometry, SurfaceGeometry::Cylinder { .. }) {
            continue;
        }
        s.cylinder_faces += 1;
        let outer_loop = face.outer_loop;
        let Some(loop_data) = brep.loops.get(outer_loop) else {
            continue;
        };
        let start_he = loop_data.half_edge;

        // Walk loop, identify off-cyl vertices and classify each
        // chord pair as same-plane or different-planes.
        let mut hes = Vec::new();
        let mut cursor = start_he;
        let mut safety = 64;
        loop {
            if safety == 0 {
                break;
            }
            safety -= 1;
            hes.push(cursor);
            let next = brep.half_edges[cursor].next;
            cursor = next;
            if cursor == start_he {
                break;
            }
        }
        let n = hes.len();
        if n < 3 {
            continue;
        }
        for i in 0..n {
            let he_out = hes[i];
            let he_in = hes[(i + n - 1) % n];
            let v = brep.half_edges[he_out].origin;
            let p = brep.vertices[v].position;
            if vertex_on_cylinder(&p, &sof.geometry) {
                continue;
            }
            s.off_cyl_vertices += 1;
            let twin_in = brep.half_edges[he_in].twin;
            let twin_out = brep.half_edges[he_out].twin;
            let plane_face_in = brep.loops[brep.half_edges[twin_in].loop_].face;
            let plane_face_out = brep.loops[brep.half_edges[twin_out].loop_].face;
            if plane_face_in == plane_face_out {
                s.same_plane_chord_pairs += 1;
            } else {
                s.diff_plane_chord_pairs += 1;
            }
            eprintln!(
                "  cyl face {fh:?}: off-cyl v={v:?} pos=({:.4},{:.4},{:.4}) twins={plane_face_in:?}/{plane_face_out:?} {}",
                p.x,
                p.y,
                p.z,
                if plane_face_in == plane_face_out {
                    "same-plane"
                } else {
                    "diff-planes"
                }
            );
        }
    }
    s
}

fn worst_tess_nm_divergence(brep: &BRep) -> f64 {
    brep.find_tessellation_non_manifold_edges(32, 1e-9)
        .iter()
        .map(|r| r.max_divergence)
        .fold(0.0, f64::max)
}

/// What incident faces' stored surfaces does vertex `v` actually
/// lie ON, and which does it drift from?  Returns
/// `(face, surface_kind, distance_to_surface)` triples for every
/// face incident to `v`.
fn surfaces_at_vertex(brep: &BRep, v: VertexHandle) -> Vec<(FaceHandle, &'static str, f64)> {
    use popkern_core::geometry::surface::SurfaceEval;
    let p = brep.vertices[v].position;
    let mut faces = std::collections::BTreeSet::new();
    for he in brep.half_edges.values() {
        if he.origin == v {
            faces.insert(brep.loops[he.loop_].face);
        }
    }
    let mut out = Vec::new();
    for fh in faces {
        let Some(face) = brep.faces.get(fh) else {
            continue;
        };
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        let kind = match &sof.geometry {
            SurfaceGeometry::Plane { .. } => "Plane",
            SurfaceGeometry::Cylinder { .. } => "Cylinder",
            SurfaceGeometry::Sphere { .. } => "Sphere",
            SurfaceGeometry::Torus { .. } => "Torus",
            SurfaceGeometry::Cone { .. } => "Cone",
            SurfaceGeometry::Nurbs(_) => "Nurbs",
            SurfaceGeometry::Gregory(_) => "Gregory",
        };
        // Per-surface distance: just call geometry's closest-point
        // and measure.  For analytic surfaces this is exact.
        let dist = match &sof.geometry {
            SurfaceGeometry::Plane { origin, normal, .. } => {
                let rel = popkern_core::geometry::Vec3::new(
                    p.x - origin.x,
                    p.y - origin.y,
                    p.z - origin.z,
                );
                rel.dot(normal).abs() / normal.length().max(1e-12)
            }
            SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } => {
                let rel = popkern_core::geometry::Vec3::new(
                    p.x - origin.x,
                    p.y - origin.y,
                    p.z - origin.z,
                );
                let len = axis.length().max(1e-12);
                let axial = rel.dot(axis) / len;
                let au =
                    popkern_core::geometry::Vec3::new(axis.x / len, axis.y / len, axis.z / len);
                let radial = popkern_core::geometry::Vec3::new(
                    rel.x - axial * au.x,
                    rel.y - axial * au.y,
                    rel.z - axial * au.z,
                );
                (radial.length() - radius).abs()
            }
            // For non-analytic, sample the surface bbox and report
            // a coarse upper bound; we don't need exact for diag.
            other => {
                let _ = other;
                let q = sof.geometry.point_at(0.5, 0.5);
                p.distance_to(&q)
            }
        };
        out.push((fh, kind, dist));
    }
    out
}

fn run_with_flags(n: usize, height: f64, radius: f64, repair: bool, diff: bool) -> BRep {
    let _g_repair = if repair {
        Some(EnvGuard::set("POPKERN_TOPO_REPAIR", "1"))
    } else {
        None
    };
    let _g_diff = if diff {
        Some(EnvGuard::set("POPKERN_TOPO_REPAIR_DIFF", "1"))
    } else {
        None
    };
    fillet_apex_only(n, height, radius)
}

fn run_with_attempts(n: usize, height: f64, radius: f64, attempts: usize) -> BRep {
    let _g_repair = EnvGuard::set("POPKERN_TOPO_REPAIR", "1");
    let _g_diff = EnvGuard::set("POPKERN_TOPO_REPAIR_DIFF", "1");
    let v = format!("{}", attempts);
    let _g_max = EnvGuard::set("POPKERN_TOPO_REPAIR_MAX_ATTEMPTS", &v);
    fillet_apex_only(n, height, radius)
}

#[ignore]
#[test]
fn t_136_a_chord_pair_inventory() {
    eprintln!("\n=== #136 chord-pair inventory across N=4..8 pyramid apex-only ===");
    eprintln!(
        "{:>3}  {:>11}  {:>13}  {:>13}  {:>13}",
        "N", "cyl_faces", "off_cyl_v", "same_plane", "diff_planes"
    );
    eprintln!("{:-<60}", "");
    for n in 4..=8 {
        let brep = fillet_apex_only(n, 0.8, 0.05);
        let s = survey_off_cyl(&brep);
        eprintln!(
            "{:>3}  {:>11}  {:>13}  {:>13}  {:>13}",
            n,
            s.cylinder_faces,
            s.off_cyl_vertices,
            s.same_plane_chord_pairs,
            s.diff_plane_chord_pairs
        );
    }
}

#[ignore]
#[test]
fn t_136_b_diff_planes_regression() {
    eprintln!("\n=== #136 diff-planes regression: tess-nm divergence ===");
    eprintln!(
        "{:>3}  {:>16}  {:>16}  {:>16}",
        "N", "no_repair", "same_plane_only", "diff_planes_on"
    );
    eprintln!("{:-<60}", "");
    for n in 4..=8 {
        let none = run_with_flags(n, 0.8, 0.05, false, false);
        let same = run_with_flags(n, 0.8, 0.05, true, false);
        let diff = run_with_flags(n, 0.8, 0.05, true, true);
        let d_none = worst_tess_nm_divergence(&none);
        let d_same = worst_tess_nm_divergence(&same);
        let d_diff = worst_tess_nm_divergence(&diff);
        eprintln!(
            "{:>3}  {:>16.4e}  {:>16.4e}  {:>16.4e}",
            n, d_none, d_same, d_diff
        );
    }
}

/// Regression gate for the #136 chord-chain anchor-pollution
/// fix.  Asserts that enabling the diff-planes branch never
/// makes tess-nm divergence worse than the same-plane-only
/// baseline on N=4..8 pyramid apex-only fillets.
///
/// Pre-fix this would catch N=7 jumping 0.12 → 0.93 and N=8
/// 0.032 → 0.93.  Post-fix (anchor on-cyl precondition in
/// `try_repair_chord_pair_different_planes`) all values match
/// the baseline because the guard rejects every chord-chain
/// candidate.
#[test]
fn t_136_d_no_regression_from_diff_planes() {
    for n in 4..=8 {
        let baseline = run_with_flags(n, 0.8, 0.05, true, false);
        let with_diff = run_with_flags(n, 0.8, 0.05, true, true);
        let d_baseline = worst_tess_nm_divergence(&baseline);
        let d_diff = worst_tess_nm_divergence(&with_diff);
        // Allow 1% headroom for floating-point noise.
        assert!(
            d_diff <= d_baseline * 1.01 + 1e-9,
            "[N={n}] diff-planes regression: baseline={:.3e} m, diff_on={:.3e} m",
            d_baseline,
            d_diff,
        );
    }
}
