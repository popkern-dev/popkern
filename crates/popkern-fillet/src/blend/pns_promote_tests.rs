//! Tests for `pns_promote` — end-to-end BlendBoundary → surfaces.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::boundary::{BlendBoundary, BoundarySegment};
use super::pns_promote::promote_blend_boundary_to_pns_surfaces;

#[path = "pns_promote_tests_curved_support.rs"]
mod pns_promote_tests_curved_support;
#[path = "pns_promote_tests_topology.rs"]
mod pns_promote_tests_topology;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// The promote path ships `PolyhedralSpline` (PK-200); tests evaluate
/// through its lossless NURBS image and keep every numeric assertion.
fn pns_nurbs(surface: &SurfaceGeometry) -> popkern_core::geometry::NurbsSurface {
    match surface {
        SurfaceGeometry::PolyhedralSpline(patch) => patch.to_nurbs(),
        SurfaceGeometry::Nurbs(ns) => ns.clone(),
        other => panic!("expected a PnS cap surface, got {other:?}"),
    }
}

fn line_curve(p0: Point3, p1: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p0, p1],
        weights: vec![1.0, 1.0],
    }
}

/// Synthesize a `2·n_fillets`-segment BlendBoundary with all
/// segments as FilletEnd on a z=0 plane fillet surface.
fn synth_planar_blend_for_promote(
    radius: f64,
) -> (
    BlendBoundary,
    Vec<crate::fillet_surface::FilletSurfaceResult>,
    BRep,
) {
    synth_planar_blend_for_promote_n(radius, 3)
}

fn synth_planar_blend_for_promote_n(
    radius: f64,
    n_fillets: usize,
) -> (
    BlendBoundary,
    Vec<crate::fillet_surface::FilletSurfaceResult>,
    BRep,
) {
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
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
        hex_boundary_cps_deg9: Vec::new(),
        fillet_end_arcs: Vec::new(),
    };
    let brep = BRep::new(Tolerance::default());
    (boundary, fillets, brep)
}

#[test]
fn promote_produces_n_surfaces_matching_pairs() {
    let (bb, fillets, brep) = synth_planar_blend_for_promote(1.5);
    let surfaces = promote_blend_boundary_to_pns_surfaces(&bb, &fillets, &brep)
        .expect("promote should succeed");
    // Each PnS sector covers a (arc, bridge) pair, so surface
    // count = segments.len() / 2.
    assert_eq!(surfaces.len(), bb.segments.len() / 2);
    for (k, surf) in surfaces.iter().enumerate() {
        match surf {
            _pns @ (SurfaceGeometry::PolyhedralSpline(_) | SurfaceGeometry::Nurbs(_)) => {
                let ns = pns_nurbs(_pns);
                assert_eq!(ns.degree_u, 4, "sector {k} degree_u");
                assert_eq!(ns.degree_v, 4, "sector {k} degree_v");
                assert_eq!(ns.num_u, 5);
                assert_eq!(ns.num_v, 5);
            }
            _ => panic!("sector {k}: expected a PnS cap surface, got other variant"),
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
        let ns = pns_nurbs(surf);
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
fn promoted_surfaces_boundary_matches_paired_segment_endpoints() {
    // Each surface covers a (arc, bridge) pair. The surface's
    // canonical v_c=4 boundary (v=1 in NURBS params) should go
    // from arc_start (u=0) to bridge_end (u=1).
    let (bb, fillets, brep) = synth_planar_blend_for_promote(1.5);
    let surfaces = promote_blend_boundary_to_pns_surfaces(&bb, &fillets, &brep).unwrap();
    for (k, surf) in surfaces.iter().enumerate() {
        let ns = pns_nurbs(surf);
        let seg_arc = &bb.segments[2 * k];
        let seg_bridge = &bb.segments[2 * k + 1];
        let arc_curve = match seg_arc {
            BoundarySegment::FilletEnd { curve, .. } => curve,
            BoundarySegment::FaceEdge { curve, .. } => curve,
        };
        let bridge_curve = match seg_bridge {
            BoundarySegment::FilletEnd { curve, .. } => curve,
            BoundarySegment::FaceEdge { curve, .. } => curve,
        };
        let (ta0, _) = arc_curve.parameter_range();
        let (_, tb1) = bridge_curve.parameter_range();
        let expected_start = arc_curve.point_at(ta0);
        let expected_end = bridge_curve.point_at(tb1);
        let nurbs_u0v1 = ns.point_at(0.0, 1.0);
        let nurbs_u1v1 = ns.point_at(1.0, 1.0);
        assert!(
            (nurbs_u0v1.x - expected_start.x).abs() < 1e-10
                && (nurbs_u0v1.y - expected_start.y).abs() < 1e-10
                && (nurbs_u0v1.z - expected_start.z).abs() < 1e-10,
            "sector {k} boundary start: nurbs (0, 1) = {nurbs_u0v1:?}, expected {expected_start:?}"
        );
        assert!(
            (nurbs_u1v1.x - expected_end.x).abs() < 1e-10
                && (nurbs_u1v1.y - expected_end.y).abs() < 1e-10
                && (nurbs_u1v1.z - expected_end.z).abs() < 1e-10,
            "sector {k} boundary end: nurbs (1, 1) = {nurbs_u1v1:?}, expected {expected_end:?}"
        );
    }
}

#[test]
fn cap_boundary_g1_with_planar_support() {
    // When every sector's supporting surface is the same z=0 plane
    // (the synth_planar setup), the bi-4 cap is expected to be
    // coplanar with z=0 — its outer boundary (v=1) must lie in
    // the plane and its surface normal there must be ±z. This
    // verifies the paper's eq (2) G¹ construction produces a cap
    // whose tangent plane matches the support at the boundary.
    //
    // Tested for N ∈ {3, 5, 6, 7}. N=4 is rejected by the cap
    // assembler (σ singularity); N≥5 exercises the general eq-(2)
    // path.
    for n_fillets in [3usize, 5, 6, 7] {
        let (bb, fillets, brep) = synth_planar_blend_for_promote_n(1.5, n_fillets);
        let surfaces = promote_blend_boundary_to_pns_surfaces(&bb, &fillets, &brep).unwrap();
        assert_eq!(
            surfaces.len(),
            n_fillets,
            "N={n_fillets}: expected {n_fillets} sub-face surfaces"
        );
        let sample_us = [0.0_f64, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0];
        for (k, surf) in surfaces.iter().enumerate() {
            let ns = pns_nurbs(surf);
            for &u in &sample_us {
                let p = ns.point_at(u, 1.0);
                assert!(
                    p.z.abs() < 1e-9,
                    "N={n_fillets} sector {k} u={u}: boundary point z={} should be 0",
                    p.z
                );
                let n = ns.normal_at(u, 1.0);
                let nl = n.length();
                assert!(
                    nl > 1e-9,
                    "N={n_fillets} sector {k} u={u}: degenerate normal"
                );
                let nz = (n.z / nl).abs();
                assert!(
                    (nz - 1.0).abs() < 1e-6,
                    "N={n_fillets} sector {k} u={u}: normal z-component {nz} should be ±1; \
                     n=({}, {}, {})",
                    n.x,
                    n.y,
                    n.z
                );
            }
        }
    }
}

#[test]
fn cap_boundary_g1_with_tilted_planar_support() {
    // Same shape as the planar test, but rotate the supporting
    // plane + boundary curves so the support normal is NOT aligned
    // with any axis. This catches regressions where the adapter
    // silently emits a z-axis normal regardless of the actual
    // supporting surface.
    let (bb_flat, fillets_flat, brep) = synth_planar_blend_for_promote(1.5);

    // Rotate everything about the x-axis by 30° (unit vectors):
    //   y' = cos θ · y − sin θ · z
    //   z' = sin θ · y + cos θ · z
    let theta = 30_f64.to_radians();
    let (c, s) = (theta.cos(), theta.sin());
    let rot = |p: Point3| Point3::new(p.x, c * p.y - s * p.z, s * p.y + c * p.z);
    let rot_v = |v: Vec3| Vec3::new(v.x, c * v.y - s * v.z, s * v.y + c * v.z);
    let rotate_curve = |curve: &NurbsCurve| NurbsCurve {
        degree: curve.degree,
        knots: curve.knots.clone(),
        control_points: curve.control_points.iter().copied().map(rot).collect(),
        weights: curve.weights.clone(),
    };

    let expected_normal = rot_v(Vec3::new(0.0, 0.0, 1.0));

    let tilted_plane = SurfaceGeometry::Plane {
        origin: rot(pt(0.0, 0.0, 0.0)),
        normal: expected_normal,
        u_axis: rot_v(Vec3::new(1.0, 0.0, 0.0)),
        v_axis: rot_v(Vec3::new(0.0, 1.0, 0.0)),
    };

    let segments: Vec<_> = bb_flat
        .segments
        .iter()
        .map(|s| match s {
            BoundarySegment::FilletEnd {
                fillet_index,
                curve,
            } => BoundarySegment::FilletEnd {
                fillet_index: *fillet_index,
                curve: rotate_curve(curve),
            },
            BoundarySegment::FaceEdge { face, curve } => BoundarySegment::FaceEdge {
                face: *face,
                curve: rotate_curve(curve),
            },
        })
        .collect();
    let bb = BlendBoundary {
        vertex: bb_flat.vertex,
        vertex_position: rot(bb_flat.vertex_position),
        segments,
        bridge_fallbacks: bb_flat.bridge_fallbacks.clone(),
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
        hex_boundary_cps_deg9: Vec::new(),
        fillet_end_arcs: Vec::new(),
    };
    let fillets: Vec<_> = fillets_flat
        .iter()
        .map(|f| crate::fillet_surface::FilletSurfaceResult {
            surface: tilted_plane.clone(),
            trim_curve_1: rotate_curve(&f.trim_curve_1),
            trim_curve_2: rotate_curve(&f.trim_curve_2),
            continuity_downgrade: f.continuity_downgrade.clone(),
        })
        .collect();

    let surfaces = promote_blend_boundary_to_pns_surfaces(&bb, &fillets, &brep).unwrap();

    let sample_us = [0.1_f64, 0.3, 0.5, 0.7, 0.9];
    for (k, surf) in surfaces.iter().enumerate() {
        let ns = pns_nurbs(surf);
        for &u in &sample_us {
            let n = ns.normal_at(u, 1.0);
            let nl = n.length();
            assert!(nl > 1e-9, "sector {k} u={u}: degenerate normal |n|={nl}");
            let unit = Vec3::new(n.x / nl, n.y / nl, n.z / nl);
            let dot = unit.x * expected_normal.x
                + unit.y * expected_normal.y
                + unit.z * expected_normal.z;
            assert!(
                dot.abs() > 1.0 - 1e-6,
                "sector {k} u={u}: cap normal {unit:?} not aligned with plane normal {expected_normal:?} (|dot|={})",
                dot.abs()
            );
        }
    }
}

#[test]
fn pns_env_enabled_reads_env_var() {
    use super::pns_promote::pns_env_enabled;
    // Sanity — depends on test-run env. Just check the function
    // doesn't panic and returns a bool.
    let _ = pns_env_enabled();
}

#[test]
fn _unused_faces_stub() {
    let _ = FaceHandle::default();
}
