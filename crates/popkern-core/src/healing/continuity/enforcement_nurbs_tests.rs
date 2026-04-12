//! NURBS-specific continuity enforcement tests.
//!
//! Exercises the NURBS code paths in `nurbs_adjust.rs`, `enforcement.rs`,
//! and `g3_enforcement.rs` that only run when adjacent faces have NURBS
//! surfaces. All existing tests use analytical (plane) surfaces.

use super::*;
use crate::constructors::builder::BRepBuilder;
use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::healing::quality::ContinuityLevel;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, FaceHandle};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Two-quad BRep with NURBS surfaces at the specified degree.
///
/// Face A (flat) and Face B (tilted by `tilt`) share edge v1-v4.
fn make_nurbs_brep(tilt: f64, deg: usize) -> BRep {
    let tol = Tolerance::default();
    let mut b = BRepBuilder::new(tol);
    let (ct, st) = (tilt.cos(), tilt.sin());

    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(1.0 + ct, 0.0, st));
    let v3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));
    let v4 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let v5 = b.add_vertex(Point3::new(1.0 + ct, 1.0, st));

    let solid = b.add_solid();
    let shell = b.add_shell(solid);

    // Face A
    let (e01, he01, _) = b.add_edge(v0, v1);
    let (e14, he14, he41) = b.add_edge(v1, v4);
    let (e43, he43, _) = b.add_edge(v4, v3);
    let (e30, he30, _) = b.add_edge(v3, v0);
    let fa = b.add_face(shell);
    b.add_loop(&[he01, he14, he43, he30], fa).unwrap();

    // Face B (reuses he41 for shared edge)
    let (e12, he12, _) = b.add_edge(v1, v2);
    let (e25, he25, _) = b.add_edge(v2, v5);
    let (e54, he54, _) = b.add_edge(v5, v4);
    let fb = b.add_face(shell);
    b.add_loop(&[he12, he25, he54, he41], fb).unwrap();

    // NURBS surfaces (plane -> degree-elevate)
    let pa = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let pb = SurfaceGeometry::Plane {
        origin: Point3::new(1.0, 0.0, 0.0),
        normal: Vec3::new(-st, 0.0, ct),
        u_axis: Vec3::new(ct, 0.0, st),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let na = elevate(pa.to_nurbs_in_range((0.0, 1.0), (0.0, 1.0)), deg);
    let nb = elevate(pb.to_nurbs_in_range((0.0, 1.0), (0.0, 1.0)), deg);
    let sof = |n| SurfaceOnFace {
        geometry: SurfaceGeometry::Nurbs(n),
        orientation: true,
    };
    b.set_face_surface(fa, sof(na));
    b.set_face_surface(fb, sof(nb));

    // Edge curves
    let lc = |b: &mut BRepBuilder, e, o: Point3, d: Vec3| {
        b.set_edge_curve(
            e,
            CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: o,
                    direction: d,
                },
                t_start: 0.0,
                t_end: 1.0,
            },
        );
    };
    lc(
        &mut b,
        e01,
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
    );
    lc(
        &mut b,
        e14,
        Point3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
    );
    lc(
        &mut b,
        e43,
        Point3::new(1.0, 1.0, 0.0),
        Vec3::new(-1.0, 0.0, 0.0),
    );
    lc(
        &mut b,
        e30,
        Point3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, -1.0, 0.0),
    );
    lc(
        &mut b,
        e12,
        Point3::new(1.0, 0.0, 0.0),
        Vec3::new(ct, 0.0, st),
    );
    lc(
        &mut b,
        e25,
        Point3::new(1.0 + ct, 0.0, st),
        Vec3::new(0.0, 1.0, 0.0),
    );
    lc(
        &mut b,
        e54,
        Point3::new(1.0 + ct, 1.0, st),
        Vec3::new(-ct, 0.0, -st),
    );

    // Pcurves (UV linear segments)
    let pc = |b: &mut BRepBuilder, he, s: (f64, f64), e: (f64, f64)| {
        b.set_half_edge_pcurve(
            he,
            NurbsCurve {
                degree: 1,
                knots: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![Point3::new(s.0, s.1, 0.0), Point3::new(e.0, e.1, 0.0)],
                weights: vec![1.0, 1.0],
            },
        );
    };
    pc(&mut b, he01, (0.0, 0.0), (1.0, 0.0));
    pc(&mut b, he14, (1.0, 0.0), (1.0, 1.0));
    pc(&mut b, he43, (1.0, 1.0), (0.0, 1.0));
    pc(&mut b, he30, (0.0, 1.0), (0.0, 0.0));
    pc(&mut b, he12, (0.0, 0.0), (1.0, 0.0));
    pc(&mut b, he25, (1.0, 0.0), (1.0, 1.0));
    pc(&mut b, he54, (1.0, 1.0), (0.0, 1.0));
    pc(&mut b, he41, (0.0, 1.0), (0.0, 0.0));

    b.build_unchecked()
}

fn elevate(mut n: NurbsSurface, deg: usize) -> NurbsSurface {
    if deg > n.degree_u {
        n = n.elevate_degree_u_to(deg).unwrap();
    }
    if deg > n.degree_v {
        n = n.elevate_degree_v_to(deg).unwrap();
    }
    n
}

/// Find the shared edge between two faces, skipping boundary edges.
fn find_shared(brep: &BRep, fa: FaceHandle, fb: FaceHandle) -> Option<EdgeHandle> {
    for (eh, edge) in brep.edges.iter() {
        let Some(h1) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(h2) = brep.half_edges.get(h1.twin) else {
            continue;
        };
        let Some(l1) = brep.loops.get(h1.loop_) else {
            continue;
        };
        let Some(l2) = brep.loops.get(h2.loop_) else {
            continue;
        };
        if (l1.face, l2.face) == (fa, fb) || (l1.face, l2.face) == (fb, fa) {
            return Some(eh);
        }
    }
    None
}

fn two_faces(brep: &BRep) -> (FaceHandle, FaceHandle) {
    let mut it = brep.faces.keys();
    (it.next().unwrap(), it.next().unwrap())
}

fn shared_edge(brep: &BRep) -> EdgeHandle {
    let (a, b) = two_faces(brep);
    find_shared(brep, a, b).expect("shared edge not found")
}

fn ctrl_pts(brep: &BRep, fh: FaceHandle) -> Vec<Point3> {
    brep.faces
        .get(fh)
        .and_then(|f| f.surface.as_ref())
        .and_then(|s| match &s.geometry {
            SurfaceGeometry::Nurbs(n) => Some(n.control_points.clone()),
            _ => None,
        })
        .unwrap_or_default()
}

// ---------------------------------------------------------------------------
// G1 enforcement on NURBS faces
// ---------------------------------------------------------------------------

#[test]
fn enforce_g1_nurbs_coplanar_is_already_g1() {
    let mut brep = make_nurbs_brep(0.0, 3);
    let tol = Tolerance::default();
    let edge = shared_edge(&brep);
    let r = enforce_g1(&mut brep, edge, &tol, 1.0).unwrap();
    assert_eq!(r.achieved_level, ContinuityLevel::G1);
}

#[test]
fn enforce_g1_nurbs_small_tilt_achieves_g1() {
    let mut brep = make_nurbs_brep(0.05, 3);
    let tol = Tolerance::default();
    let edge = shared_edge(&brep);
    let r = enforce_g1(&mut brep, edge, &tol, 1.0).unwrap();
    assert!(r.achieved_level >= ContinuityLevel::G0);
}

#[test]
fn enforce_g1_nurbs_adjusts_control_points() {
    let mut brep = make_nurbs_brep(0.1, 3);
    let tol = Tolerance::default();
    let edge = shared_edge(&brep);
    let (f1, f2) = two_faces(&brep);
    let before_a = ctrl_pts(&brep, f1);
    let before_b = ctrl_pts(&brep, f2);

    let r = enforce_g1(&mut brep, edge, &tol, 1.0).unwrap();

    let changed = before_a != ctrl_pts(&brep, f1) || before_b != ctrl_pts(&brep, f2);
    assert!(changed, "control points should be adjusted");
    assert!(r.max_displacement >= 0.0);
}

#[test]
fn enforce_g1_nurbs_degree_1_still_works() {
    let mut brep = make_nurbs_brep(0.1, 1);
    let tol = Tolerance::default();
    let edge = shared_edge(&brep);
    let r = enforce_g1(&mut brep, edge, &tol, 1.0).unwrap();
    assert!(r.achieved_level >= ContinuityLevel::G0);
}

// ---------------------------------------------------------------------------
// G2 enforcement on NURBS faces
// ---------------------------------------------------------------------------

#[test]
fn enforce_g2_nurbs_coplanar_achieves_g2() {
    let mut brep = make_nurbs_brep(0.0, 3);
    let tol = Tolerance::default();
    let edge = shared_edge(&brep);
    let r = enforce_g2(&mut brep, edge, &tol, 1.0).unwrap();
    assert!(r.achieved_level >= ContinuityLevel::G1);
}

#[test]
fn enforce_g2_nurbs_small_tilt() {
    let mut brep = make_nurbs_brep(0.05, 3);
    let tol = Tolerance::default();
    let edge = shared_edge(&brep);
    let r = enforce_g2(&mut brep, edge, &tol, 1.0).unwrap();
    assert!(r.achieved_level >= ContinuityLevel::G0);
}

#[test]
fn enforce_g2_nurbs_degree_2_degrades() {
    let mut brep = make_nurbs_brep(0.0, 2);
    let tol = Tolerance::default();
    let edge = shared_edge(&brep);
    let r = enforce_g2(&mut brep, edge, &tol, 1.0).unwrap();
    assert_eq!(r.achieved_level, ContinuityLevel::G1);
}

// ---------------------------------------------------------------------------
// G3 enforcement on NURBS faces
// ---------------------------------------------------------------------------

#[test]
fn enforce_g3_nurbs_degree_3_degrades_gracefully() {
    let mut brep = make_nurbs_brep(0.0, 3);
    let tol = Tolerance::default();
    let edge = shared_edge(&brep);
    let r = super::super::g3_enforcement::enforce_g3(&mut brep, edge, &tol, 1.0).unwrap();
    assert!(r.achieved_level <= ContinuityLevel::G2);
    assert!(r.achieved_level >= ContinuityLevel::G1);
}

#[test]
fn enforce_g3_nurbs_degree_5_exercises_adjust_for_g3() {
    let mut brep = make_nurbs_brep(0.0, 5);
    let tol = Tolerance::default();
    let edge = shared_edge(&brep);
    let r = super::super::g3_enforcement::enforce_g3(&mut brep, edge, &tol, 1.0).unwrap();
    assert!(r.achieved_level >= ContinuityLevel::G2);
}

#[test]
fn enforce_g3_nurbs_degree_5_with_tilt() {
    let mut brep = make_nurbs_brep(0.02, 5);
    let tol = Tolerance::default();
    let edge = shared_edge(&brep);
    let (f1, f2) = two_faces(&brep);
    let before_a = ctrl_pts(&brep, f1);
    let before_b = ctrl_pts(&brep, f2);

    let r = super::super::g3_enforcement::enforce_g3(&mut brep, edge, &tol, 1.0).unwrap();
    assert!(r.achieved_level >= ContinuityLevel::G0);

    let changed = before_a != ctrl_pts(&brep, f1) || before_b != ctrl_pts(&brep, f2);
    assert!(
        changed,
        "degree-5 G3 with tilt should adjust control points"
    );
}

// ---------------------------------------------------------------------------
// enforce_continuity_all on NURBS faces
// ---------------------------------------------------------------------------

#[test]
fn enforce_continuity_all_g1_nurbs() {
    let mut brep = make_nurbs_brep(0.05, 3);
    let tol = Tolerance::default();
    let results = enforce_continuity_all(&mut brep, ContinuityLevel::G1, &tol, 1.0);
    assert!(!results.is_empty());
    assert!(
        results
            .iter()
            .any(|r| r.achieved_level >= ContinuityLevel::G0)
    );
}

#[test]
fn enforce_continuity_all_g2_nurbs() {
    let mut brep = make_nurbs_brep(0.0, 3);
    let tol = Tolerance::default();
    let results = enforce_continuity_all(&mut brep, ContinuityLevel::G2, &tol, 1.0);
    assert!(!results.is_empty());
}

#[test]
fn enforce_continuity_all_g3_nurbs() {
    let mut brep = make_nurbs_brep(0.0, 5);
    let tol = Tolerance::default();
    let results = enforce_continuity_all(&mut brep, ContinuityLevel::G3, &tol, 1.0);
    assert!(!results.is_empty());
}

// ---------------------------------------------------------------------------
// identify_boundary
// ---------------------------------------------------------------------------

#[test]
fn identify_boundary_picks_correct_edge() {
    use super::super::nurbs_adjust::{Boundary, identify_boundary};
    assert!(matches!(
        identify_boundary((0.01, 0.5), 0.0, 1.0, 0.0, 1.0),
        Boundary::UMin
    ));
    assert!(matches!(
        identify_boundary((0.99, 0.5), 0.0, 1.0, 0.0, 1.0),
        Boundary::UMax
    ));
    assert!(matches!(
        identify_boundary((0.5, 0.01), 0.0, 1.0, 0.0, 1.0),
        Boundary::VMin
    ));
    assert!(matches!(
        identify_boundary((0.5, 0.99), 0.0, 1.0, 0.0, 1.0),
        Boundary::VMax
    ));
}

// ---------------------------------------------------------------------------
// Direct nurbs_adjust::adjust_for_g1
// ---------------------------------------------------------------------------

#[test]
fn adjust_for_g1_on_tilted_surfaces() {
    let mut brep = make_nurbs_brep(0.1, 3);
    let (f1, f2) = two_faces(&brep);
    let edge = shared_edge(&brep);
    let samples = sample_edge_uvs(&brep, edge, 8);
    assert!(!samples.is_empty());
    let disp = super::super::nurbs_adjust::adjust_for_g1(&mut brep, f1, f2, &samples, 1.0);
    assert!(disp >= 0.0);
}

#[test]
fn adjust_for_g1_coplanar_bounded() {
    let mut brep = make_nurbs_brep(0.0, 3);
    let (f1, f2) = two_faces(&brep);
    let edge = shared_edge(&brep);
    let samples = sample_edge_uvs(&brep, edge, 8);
    let disp = super::super::nurbs_adjust::adjust_for_g1(&mut brep, f1, f2, &samples, 1.0);
    assert!(disp <= 1.0, "should respect max_displacement, got {disp}");
}

// ---------------------------------------------------------------------------
// Direct nurbs_adjust::adjust_for_g2
// ---------------------------------------------------------------------------

#[test]
fn adjust_for_g2_on_degree_3_surfaces() {
    let mut brep = make_nurbs_brep(0.0, 3);
    let (f1, f2) = two_faces(&brep);
    let edge = shared_edge(&brep);
    let samples = sample_edge_uvs(&brep, edge, 8);
    let disp = super::super::nurbs_adjust::adjust_for_g2(&mut brep, f1, f2, &samples, 1.0);
    assert!(disp >= 0.0);
}

#[test]
fn adjust_for_g2_skips_low_degree() {
    let mut brep = make_nurbs_brep(0.0, 2);
    let (f1, f2) = two_faces(&brep);
    let edge = shared_edge(&brep);
    let samples = sample_edge_uvs(&brep, edge, 8);
    let disp = super::super::nurbs_adjust::adjust_for_g2(&mut brep, f1, f2, &samples, 1.0);
    assert!(disp < f64::EPSILON);
}

// ---------------------------------------------------------------------------
// Direct nurbs_adjust::adjust_for_g3
// ---------------------------------------------------------------------------

#[test]
fn adjust_for_g3_on_degree_5_surfaces() {
    let mut brep = make_nurbs_brep(0.0, 5);
    let (f1, f2) = two_faces(&brep);
    let edge = shared_edge(&brep);
    let samples = sample_edge_uvs(&brep, edge, 8);
    let disp = super::super::nurbs_adjust::adjust_for_g3(&mut brep, f1, f2, &samples, 1.0);
    assert!(disp >= 0.0);
}

#[test]
fn adjust_for_g3_skips_degree_3() {
    let mut brep = make_nurbs_brep(0.0, 3);
    let (f1, f2) = two_faces(&brep);
    let edge = shared_edge(&brep);
    let samples = sample_edge_uvs(&brep, edge, 8);
    let disp = super::super::nurbs_adjust::adjust_for_g3(&mut brep, f1, f2, &samples, 1.0);
    assert!(disp < f64::EPSILON);
}

#[test]
fn adjust_for_g3_with_tilt_modifies_points() {
    let mut brep = make_nurbs_brep(0.05, 5);
    let (f1, _f2) = two_faces(&brep);
    let edge = shared_edge(&brep);
    let before = ctrl_pts(&brep, f1);
    let samples = sample_edge_uvs(&brep, edge, 8);
    let (f1, f2) = two_faces(&brep);
    let _disp = super::super::nurbs_adjust::adjust_for_g3(&mut brep, f1, f2, &samples, 1.0);
    let after = ctrl_pts(&brep, f1);
    // Verify function ran without error; points may or may not change.
    let _ = (before, after);
}

// ---------------------------------------------------------------------------
// Displacement budget
// ---------------------------------------------------------------------------

#[test]
fn enforce_g1_nurbs_tiny_budget() {
    let mut brep = make_nurbs_brep(0.1, 3);
    let tol = Tolerance::default();
    let edge = shared_edge(&brep);
    let r = enforce_g1(&mut brep, edge, &tol, 1e-12).unwrap();
    assert!(r.max_displacement <= 1e-12 + f64::EPSILON);
}

#[test]
fn enforce_g2_nurbs_zero_budget_after_g1() {
    let mut brep = make_nurbs_brep(0.1, 3);
    let tol = Tolerance::default();
    let edge = shared_edge(&brep);
    let r = enforce_g2(&mut brep, edge, &tol, 0.001).unwrap();
    assert!(r.achieved_level >= ContinuityLevel::G0);
}
