//! Tests for `pns_promote` — end-to-end BlendBoundary → surfaces.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::boundary::{BlendBoundary, BoundarySegment};
use super::pns_promote::promote_blend_boundary_to_pns_surfaces;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn line_curve(p0: Point3, p1: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p0, p1],
        weights: vec![1.0, 1.0],
    }
}

/// Synthesize a 6-segment BlendBoundary (n_fillets=3) with all
/// segments as FilletEnd on a z=0 plane fillet surface.
fn synth_planar_blend_for_promote(
    radius: f64,
) -> (
    BlendBoundary,
    Vec<crate::fillet_surface::FilletSurfaceResult>,
    BRep,
) {
    let n_fillets = 3;
    let mut segments = Vec::with_capacity(2 * n_fillets);
    let vertex = VertexHandle::default();

    // Plane z = 0.
    let plane = SurfaceGeometry::Plane {
        origin: pt(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };

    let mut fillets = Vec::new();
    let n_seg = 2 * n_fillets;
    for k in 0..n_seg {
        let a0 = 2.0 * std::f64::consts::PI * (k as f64) / (n_seg as f64);
        let a1 = 2.0 * std::f64::consts::PI * ((k + 1) as f64) / (n_seg as f64);
        let p0 = pt(radius * a0.cos(), radius * a0.sin(), 0.0);
        let p1 = pt(radius * a1.cos(), radius * a1.sin(), 0.0);
        let curve = line_curve(p0, p1);
        // All FilletEnd for simplicity — we need a fillet surface
        // index per segment.
        segments.push(BoundarySegment::FilletEnd {
            fillet_index: k,
            curve: curve.clone(),
        });
        let dummy = line_curve(p0, p1);
        fillets.push(crate::fillet_surface::FilletSurfaceResult {
            surface: plane.clone(),
            trim_curve_1: dummy.clone(),
            trim_curve_2: dummy,
            continuity_downgrade: None,
        });
    }

    let boundary = BlendBoundary {
        vertex,
        vertex_position: pt(0.0, 0.0, 0.4),
        segments,
        bridge_fallbacks: Vec::new(),
    };
    let brep = BRep::new(Tolerance::default());
    (boundary, fillets, brep)
}

#[test]
fn promote_produces_n_surfaces_matching_segments() {
    let (bb, fillets, brep) = synth_planar_blend_for_promote(1.5);
    let surfaces = promote_blend_boundary_to_pns_surfaces(&bb, &fillets, &brep)
        .expect("promote should succeed");
    assert_eq!(surfaces.len(), bb.segments.len());
    // Each surface should be a NurbsSurface (bi-4 plain Bezier).
    for (k, surf) in surfaces.iter().enumerate() {
        match surf {
            SurfaceGeometry::Nurbs(ns) => {
                assert_eq!(ns.degree_u, 4, "sector {k} degree_u");
                assert_eq!(ns.degree_v, 4, "sector {k} degree_v");
                assert_eq!(ns.num_u, 5);
                assert_eq!(ns.num_v, 5);
            }
            _ => panic!("sector {k}: expected NurbsSurface, got other variant"),
        }
    }
}

#[test]
fn promoted_surfaces_interpolate_pole_at_canonical_origin() {
    // Sanity: each surface's NURBS point at (u=0, v=0) = canonical
    // patch pole = the auto-derived central point.
    let (bb, fillets, brep) = synth_planar_blend_for_promote(1.5);
    let surfaces = promote_blend_boundary_to_pns_surfaces(&bb, &fillets, &brep).unwrap();

    // Derive the auto-pole the same way promote does, for the
    // expected value.
    use super::pns_adapter::blend_boundary_to_pns_ring_with_surfaces;
    use popkern_core::geometry::gregory::pole_from_boundary;
    let ring = blend_boundary_to_pns_ring_with_surfaces(&bb, &fillets, &brep).unwrap();
    let expected_pole = pole_from_boundary(&ring);

    for (k, surf) in surfaces.iter().enumerate() {
        let SurfaceGeometry::Nurbs(ns) = surf else {
            panic!("expected NurbsSurface");
        };
        let p00 = ns.point_at(0.0, 0.0);
        let diff = (p00.x - expected_pole.x).abs()
            + (p00.y - expected_pole.y).abs()
            + (p00.z - expected_pole.z).abs();
        assert!(
            diff < 1e-10,
            "sector {k} pole at (0, 0) = {p00:?}, expected {expected_pole:?}"
        );
    }
}

#[test]
fn promoted_surfaces_boundary_matches_segment_curves() {
    // Each surface's canonical v_c=4 edge (outer boundary) should
    // trace the corresponding BlendBoundary segment's NurbsCurve.
    // In the NURBS form this is the u=? line where the bi-4
    // outer-boundary row v_c=4 lives.
    //
    // Recall our canonical convention: patch[u_c][v_c] with
    // boundary at v_c=4. In the NurbsSurface conversion, the grid
    // is row-major (u_c varies within rows), so the v_c=4 edge
    // corresponds to varying u with v=1 in NURBS parameters.
    let (bb, fillets, brep) = synth_planar_blend_for_promote(1.5);
    let surfaces = promote_blend_boundary_to_pns_surfaces(&bb, &fillets, &brep).unwrap();
    for (k, surf) in surfaces.iter().enumerate() {
        let SurfaceGeometry::Nurbs(ns) = surf else {
            panic!("expected NurbsSurface");
        };
        let segment_curve = match &bb.segments[k] {
            BoundarySegment::FilletEnd { curve, .. } => curve,
            BoundarySegment::FaceEdge { curve, .. } => curve,
        };
        let (t_min, t_max) = segment_curve.parameter_range();
        // Endpoints should match.
        let seg_start = segment_curve.point_at(t_min);
        let seg_end = segment_curve.point_at(t_max);
        // At v=1 (canonical v_c=4 outer boundary edge), u varies
        // 0 to 1. u=0 and u=1 should be the sector's outer-edge
        // endpoints, matching segment start/end.
        let nurbs_u0v1 = ns.point_at(0.0, 1.0);
        let nurbs_u1v1 = ns.point_at(1.0, 1.0);
        assert!(
            (nurbs_u0v1.x - seg_start.x).abs() < 1e-10
                && (nurbs_u0v1.y - seg_start.y).abs() < 1e-10
                && (nurbs_u0v1.z - seg_start.z).abs() < 1e-10,
            "sector {k} boundary start: nurbs (0, 1) = {nurbs_u0v1:?}, segment start = {seg_start:?}"
        );
        assert!(
            (nurbs_u1v1.x - seg_end.x).abs() < 1e-10
                && (nurbs_u1v1.y - seg_end.y).abs() < 1e-10
                && (nurbs_u1v1.z - seg_end.z).abs() < 1e-10,
            "sector {k} boundary end: nurbs (1, 1) = {nurbs_u1v1:?}, segment end = {seg_end:?}"
        );
    }
}

#[test]
fn _unused_faces_stub() {
    // Touch FaceHandle to silence unused-import warnings on some
    // configurations.
    let _ = FaceHandle::default();
}
