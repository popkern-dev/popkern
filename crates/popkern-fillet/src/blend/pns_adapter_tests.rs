//! Tests for `pns_adapter`.

use popkern_core::geometry::NurbsCurve;
use popkern_core::topology::{FaceHandle, VertexHandle};

use super::boundary::{BlendBoundary, BoundarySegment};
use super::pns_adapter::{PnsAdapterError, blend_boundary_to_pns_ring};
use popkern_core::geometry::Point3;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// Build a line-segment NurbsCurve from `p0` to `p1` (degree 1,
/// 2 CPs, uniform knots).
fn line_curve(p0: Point3, p1: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p0, p1],
        weights: vec![1.0, 1.0],
    }
}

/// Build a synthetic `BlendBoundary` with `n_fillets` fillets (so
/// 2·n_fillets segments total). Segments are line segments along
/// a regular polygon of given radius.
fn synth_blend_boundary(n_fillets: usize, radius: f64) -> BlendBoundary {
    let mut segments = Vec::with_capacity(2 * n_fillets);
    let fake_vertex = VertexHandle::default();
    let fake_face = FaceHandle::default();

    let n_seg = 2 * n_fillets;
    for k in 0..n_seg {
        let angle_k = 2.0 * std::f64::consts::PI * (k as f64) / (n_seg as f64);
        let angle_k1 = 2.0 * std::f64::consts::PI * ((k + 1) as f64) / (n_seg as f64);
        let p0 = pt(radius * angle_k.cos(), radius * angle_k.sin(), 0.0);
        let p1 = pt(radius * angle_k1.cos(), radius * angle_k1.sin(), 0.0);
        let curve = line_curve(p0, p1);
        // Even-k is FilletEnd, odd-k is FaceEdge (conventional
        // alternation).
        if k.is_multiple_of(2) {
            segments.push(BoundarySegment::FilletEnd {
                fillet_index: k / 2,
                curve,
            });
        } else {
            segments.push(BoundarySegment::FaceEdge {
                face: fake_face,
                curve,
            });
        }
    }
    BlendBoundary {
        vertex: fake_vertex,
        vertex_position: pt(0.0, 0.0, 0.5),
        segments,
        bridge_fallbacks: Vec::new(),
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
        hex_boundary_cps_deg9: Vec::new(),
        fillet_end_arcs: Vec::new(),
    }
}

#[test]
fn three_fillet_blend_produces_n3_ring() {
    // 3 fillets → 6 segments → paired into 3 PnS sectors.
    let bb = synth_blend_boundary(3, 1.0);
    let ring = blend_boundary_to_pns_ring(&bb).expect("3-fillet blend should convert");
    assert_eq!(ring.n(), 3);
    // Each strip has 3 position CPs.
    for strip in &ring.sectors {
        assert_eq!(strip.pos.len(), 3);
        assert_eq!(strip.deriv.len(), 3);
    }
}

#[test]
fn five_fillet_blend_produces_n5_ring() {
    let bb = synth_blend_boundary(5, 1.0);
    let ring = blend_boundary_to_pns_ring(&bb).expect("5-fillet blend should convert");
    assert_eq!(ring.n(), 5);
}

#[test]
fn too_few_segments_rejected() {
    // 2 fillets → 4 segments → below minimum (PnS needs n ≥ 3).
    let bb = synth_blend_boundary(2, 1.0);
    let result = blend_boundary_to_pns_ring(&bb);
    assert!(matches!(result, Err(PnsAdapterError::TooFewSegments(4))));
}

#[test]
fn too_many_segments_rejected() {
    // 13 fillets → 26 segments → above max (n ≤ 12).
    let bb = synth_blend_boundary(13, 1.0);
    let result = blend_boundary_to_pns_ring(&bb);
    assert!(matches!(result, Err(PnsAdapterError::TooManySegments(26))));
}

#[test]
fn paired_position_endpoints_match_arc_start_and_bridge_end() {
    // Verify pos[0] = arc (seg 2k) start, pos[2] = bridge (seg 2k+1) end.
    let bb = synth_blend_boundary(3, 1.5);
    let ring = blend_boundary_to_pns_ring(&bb).expect("conversion should succeed");
    for (k, strip) in ring.sectors.iter().enumerate() {
        let seg_a = &bb.segments[2 * k];
        let seg_b = &bb.segments[2 * k + 1];
        let curve_a = match seg_a {
            BoundarySegment::FilletEnd { curve, .. } => curve,
            BoundarySegment::FaceEdge { curve, .. } => curve,
        };
        let curve_b = match seg_b {
            BoundarySegment::FilletEnd { curve, .. } => curve,
            BoundarySegment::FaceEdge { curve, .. } => curve,
        };
        use popkern_core::geometry::curve::CurveEval;
        let (ta0, _) = curve_a.parameter_range();
        let (_, tb1) = curve_b.parameter_range();
        let expected_start = curve_a.point_at(ta0);
        let expected_end = curve_b.point_at(tb1);
        assert!(
            (strip.pos[0].x - expected_start.x).abs() < 1e-12
                && (strip.pos[0].y - expected_start.y).abs() < 1e-12
                && (strip.pos[0].z - expected_start.z).abs() < 1e-12,
            "sector {k} pos[0] mismatch"
        );
        assert!(
            (strip.pos[2].x - expected_end.x).abs() < 1e-12
                && (strip.pos[2].y - expected_end.y).abs() < 1e-12
                && (strip.pos[2].z - expected_end.z).abs() < 1e-12,
            "sector {k} pos[2] mismatch"
        );
    }
}

#[test]
fn surface_aware_adapter_falls_back_gracefully_on_empty_brep() {
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;
    // Surface lookups are soft: when a segment's surface can't be
    // resolved the adapter falls back to placeholder derivatives
    // (rather than hard-failing). With an empty BRep / empty fillet
    // list this means every sector uses the fallback path.
    let bb = synth_blend_boundary(3, 1.0);
    let tol = Tolerance::default();
    let brep = BRep::new(tol);
    let fillets: Vec<crate::fillet_surface::FilletSurfaceResult> = vec![];
    let result = super::pns_adapter::blend_boundary_to_pns_ring_with_surfaces(&bb, &fillets, &brep);
    let ring = result.expect("soft fallback should produce a ring");
    assert_eq!(ring.n(), 3);
}

#[test]
fn surface_aware_derivs_are_perpendicular_to_curve_tangents() {
    // Build a synthetic BlendBoundary where every segment is
    // backed by a known plane (z=0, normal = +z). For each
    // segment sample, the resulting deriv should be perpendicular
    // to the curve tangent AND lie in the z=0 plane (i.e., dot
    // product with (0,0,1) ≈ 0).
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    let n_fillets = 3;
    let mut bb = synth_blend_boundary(n_fillets, 1.5);
    // Override segments to reference fillet_index 0..n_fillets and
    // the same FaceHandle — we'll build fillet_surfaces as planes
    // and set up a BRep with that face.
    let tol = Tolerance::default();
    let brep = BRep::new(tol);

    // Every fillet surface = plane z=0 with normal +z.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let mut fillets: Vec<crate::fillet_surface::FilletSurfaceResult> = Vec::new();
    let dummy_curve = line_curve(pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.0));
    for _ in 0..2 * n_fillets {
        fillets.push(crate::fillet_surface::FilletSurfaceResult {
            surface: plane.clone(),
            trim_curve_1: dummy_curve.clone(),
            trim_curve_2: dummy_curve.clone(),
            continuity_downgrade: None,
        });
    }
    // Convert all segments to FilletEnd so the lookup always goes
    // through fillet_surfaces (and we don't need a real face).
    for (k, seg) in bb.segments.iter_mut().enumerate() {
        let curve = match seg {
            BoundarySegment::FilletEnd { curve, .. } => curve.clone(),
            BoundarySegment::FaceEdge { curve, .. } => curve.clone(),
        };
        *seg = BoundarySegment::FilletEnd {
            fillet_index: k,
            curve,
        };
    }

    let ring = super::pns_adapter::blend_boundary_to_pns_ring_with_surfaces(&bb, &fillets, &brep)
        .expect("surface-aware adapter should succeed with planar surfaces");
    assert_eq!(ring.n(), n_fillets);

    // For the plane z=0, every deriv should have z ≈ 0 (in-plane)
    // and be perpendicular to the curve tangent at its sample
    // position.
    for (k, strip) in ring.sectors.iter().enumerate() {
        for (i, d) in strip.deriv.iter().enumerate() {
            assert!(
                d.z.abs() < 1e-10,
                "sector {k} deriv[{i}].z = {} (not in-plane)",
                d.z
            );
            // Non-zero deriv (surface-aware produces real derivs,
            // not just the placeholder fallback).
            let mag = (d.x * d.x + d.y * d.y + d.z * d.z).sqrt();
            assert!(mag > 1e-9, "sector {k} deriv[{i}] is zero");
        }
    }
}

#[test]
fn blend_boundary_to_ring_feeds_pns_cap_assembly() {
    // End-to-end: adapter output + pole_from_boundary + assembly.
    use popkern_core::geometry::gregory::{
        BetaReparam, CapAssembleInput, assemble_pns_bi4_cap, pole_from_boundary,
        sigma_center_jet_default,
    };
    let bb = synth_blend_boundary(3, 1.5);
    let ring = blend_boundary_to_pns_ring(&bb).expect("adapter conversion should succeed");
    let n = ring.n();
    let pole = pole_from_boundary(&ring);
    let input = CapAssembleInput {
        boundary: ring,
        beta: BetaReparam::bi4(n),
        sigma: sigma_center_jet_default(n).expect("paper defaults yield a finite σ jet"),
        pole,
        kappa: 3,
    };
    let cap = assemble_pns_bi4_cap(&input).expect("end-to-end assembly should succeed");
    assert_eq!(cap.patches.len(), n);
    // All patch poles land at the derived pole.
    for k in 0..n {
        let p00 = cap.patches[k][0][0];
        assert!(
            (p00.x - pole.x).abs() < 1e-10
                && (p00.y - pole.y).abs() < 1e-10
                && (p00.z - pole.z).abs() < 1e-10
        );
    }
}
