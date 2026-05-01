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

#[ignore]
#[test]
fn t_136_b3_anchor_oncyl_status() {
    // Hypothesis: pyramid apex-only fillets produce CHAINS of
    // off-cyl vertices in cylinder face loops, e.g.
    // [cyl_anchor, cap_off, base_off, base_off, cap_off, cyl_anchor].
    // Each chord-pair candidate centres on one off-cyl vertex, but
    // its `anchor_in` and `anchor_out` (the loop neighbours) may
    // themselves be OFF-cyl.  The diff-planes surgery rewrites the
    // chord edges using these anchors as endpoints — if anchors
    // aren't actually on the cylinder, the rewritten curve runs
    // from off-surface to off-surface, which is geometrically
    // wrong even before any vertex move.
    eprintln!("\n=== #136 anchor on-cyl status per candidate ===");
    for n in 4..=8 {
        let brep = fillet_apex_only(n, 0.8, 0.05);
        let mut clean = 0; // both anchors on-cyl
        let mut polluted = 0; // at least one anchor off-cyl
        for (_fh, face) in brep.faces.iter() {
            let Some(sof) = face.surface.as_ref() else {
                continue;
            };
            if !matches!(sof.geometry, SurfaceGeometry::Cylinder { .. }) {
                continue;
            }
            let outer_loop = face.outer_loop;
            let Some(loop_data) = brep.loops.get(outer_loop) else {
                continue;
            };
            let start_he = loop_data.half_edge;
            let mut hes = Vec::new();
            let mut cursor = start_he;
            let mut safety = 64;
            loop {
                if safety == 0 {
                    break;
                }
                safety -= 1;
                hes.push(cursor);
                cursor = brep.half_edges[cursor].next;
                if cursor == start_he {
                    break;
                }
            }
            let nh = hes.len();
            for i in 0..nh {
                let he_out = hes[i];
                let he_in = hes[(i + nh - 1) % nh];
                let v = brep.half_edges[he_out].origin;
                let p = brep.vertices[v].position;
                if vertex_on_cylinder(&p, &sof.geometry) {
                    continue;
                }
                let anchor_in_v = brep.half_edges[he_in].origin;
                let twin_out = brep.half_edges[he_out].twin;
                let anchor_out_v = brep.half_edges[twin_out].origin;
                let p_in = brep.vertices[anchor_in_v].position;
                let p_out = brep.vertices[anchor_out_v].position;
                let in_on = vertex_on_cylinder(&p_in, &sof.geometry);
                let out_on = vertex_on_cylinder(&p_out, &sof.geometry);
                if in_on && out_on {
                    clean += 1;
                } else {
                    polluted += 1;
                }
            }
        }
        eprintln!(
            "  N={n}: clean (both anchors on-cyl) = {clean}, polluted (>=1 anchor off-cyl) = {polluted}"
        );
    }
}

#[ignore]
#[test]
fn t_136_b2_attempt_cap_sweep() {
    eprintln!("\n=== #136 attempt-cap sweep: tess-nm vs MAX_ATTEMPTS ===");
    eprintln!(
        "{:>3}  {:>14}  {:>14}  {:>14}  {:>14}",
        "N", "attempts=1", "attempts=2", "attempts=4", "attempts=8"
    );
    eprintln!("{:-<76}", "");
    for n in 4..=8 {
        let mut row = Vec::new();
        for k in [1usize, 2, 4, 8] {
            let brep = run_with_attempts(n, 0.8, 0.05, k);
            row.push(worst_tess_nm_divergence(&brep));
        }
        eprintln!(
            "{:>3}  {:>14.4e}  {:>14.4e}  {:>14.4e}  {:>14.4e}",
            n, row[0], row[1], row[2], row[3]
        );
    }
}

#[ignore]
#[test]
fn t_136_c_diff_planes_neighbour_drift() {
    eprintln!("\n=== #136 diff-planes neighbour drift: where does v_off snap break things? ===");
    // N=7 chosen because it shows the canonical 0.12 → 0.93 regression.
    let n = 7;
    let height = 0.8;
    let radius = 0.05;

    // Run with same-plane only, find off-cyl vertices that REMAIN
    // (i.e. the diff-plane cases the same-plane branch refused).
    let pre = run_with_flags(n, height, radius, true, false);
    let pre_survey = survey_off_cyl(&pre);
    eprintln!(
        "After same-plane only (N={n}): cyl_faces={} off_cyl_v={} same_plane={} diff_planes={}",
        pre_survey.cylinder_faces,
        pre_survey.off_cyl_vertices,
        pre_survey.same_plane_chord_pairs,
        pre_survey.diff_plane_chord_pairs
    );

    // Identify the off-cyl vertices and which surfaces they drift
    // from in the pre-state.
    let mut off_cyl_vs = std::collections::BTreeSet::<VertexHandle>::new();
    for (_fh, face) in pre.faces.iter() {
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        if !matches!(sof.geometry, SurfaceGeometry::Cylinder { .. }) {
            continue;
        }
        let outer_loop = face.outer_loop;
        let Some(loop_data) = pre.loops.get(outer_loop) else {
            continue;
        };
        let start_he = loop_data.half_edge;
        let mut cursor = start_he;
        let mut safety = 64;
        loop {
            if safety == 0 {
                break;
            }
            safety -= 1;
            let v = pre.half_edges[cursor].origin;
            let p = pre.vertices[v].position;
            if !vertex_on_cylinder(&p, &sof.geometry) {
                off_cyl_vs.insert(v);
            }
            cursor = pre.half_edges[cursor].next;
            if cursor == start_he {
                break;
            }
        }
    }
    eprintln!("Off-cyl vertices to inspect: {}", off_cyl_vs.len());

    // Run with diff-planes ON, and for each of those vertices
    // (matched by approximate position since handles may differ),
    // print surface drift on every incident face — both before
    // and after.
    let post = run_with_flags(n, height, radius, true, true);
    eprintln!(
        "Post diff-planes:  worst tess-nm = {:.3e}",
        worst_tess_nm_divergence(&post)
    );

    // For each off-cyl vertex in `pre`, find the matching vertex in
    // `post` (closest position) and dump its incident-surface drift.
    for v in off_cyl_vs {
        let p_pre = pre.vertices[v].position;
        // Find closest vertex in post.
        let (v_post, p_post) = post
            .vertices
            .iter()
            .min_by(|a, b| {
                let da = a.1.position.distance_to(&p_pre);
                let db = b.1.position.distance_to(&p_pre);
                da.partial_cmp(&db).unwrap()
            })
            .map(|(h, vd)| (h, vd.position))
            .expect("non-empty post");

        let snap_dist = p_pre.distance_to(&p_post);
        eprintln!(
            "\n  v_pre={v:?} pos=({:.4},{:.4},{:.4})  →  v_post={v_post:?} pos=({:.4},{:.4},{:.4})  snap={:.4e} m",
            p_pre.x, p_pre.y, p_pre.z, p_post.x, p_post.y, p_post.z, snap_dist
        );

        let drift_pre = surfaces_at_vertex(&pre, v);
        let drift_post = surfaces_at_vertex(&post, v_post);

        let max_pre = drift_pre.iter().map(|t| t.2).fold(0.0, f64::max);
        let max_post = drift_post.iter().map(|t| t.2).fold(0.0, f64::max);
        eprintln!(
            "    incident surfaces (pre): n={}, max_drift={:.4e}",
            drift_pre.len(),
            max_pre
        );
        for (fh, kind, dist) in &drift_pre {
            eprintln!("      pre  {fh:?} {kind:>8}  drift={:.4e}", dist);
        }
        eprintln!(
            "    incident surfaces (post): n={}, max_drift={:.4e}",
            drift_post.len(),
            max_post
        );
        for (fh, kind, dist) in &drift_post {
            eprintln!("      post {fh:?} {kind:>8}  drift={:.4e}", dist);
        }
    }
}
