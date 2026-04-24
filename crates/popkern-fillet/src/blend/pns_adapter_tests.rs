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
    }
}

#[test]
fn three_fillet_blend_produces_n6_ring() {
    // 3 fillets → 6 segments → valence-6 PnS cap.
    let bb = synth_blend_boundary(3, 1.0);
    let ring = blend_boundary_to_pns_ring(&bb).expect("3-fillet blend should convert");
    assert_eq!(ring.n(), 6);
    // Each strip has 3 position CPs.
    for strip in &ring.sectors {
        assert_eq!(strip.pos.len(), 3);
        assert_eq!(strip.deriv.len(), 3);
    }
}

#[test]
fn four_fillet_blend_produces_n8_ring() {
    let bb = synth_blend_boundary(4, 1.0);
    let ring = blend_boundary_to_pns_ring(&bb).expect("4-fillet blend should convert");
    assert_eq!(ring.n(), 8);
}

#[test]
fn too_few_segments_rejected() {
    // 2 fillets → 4 segments → below minimum.
    let bb = synth_blend_boundary(2, 1.0);
    let result = blend_boundary_to_pns_ring(&bb);
    assert!(matches!(result, Err(PnsAdapterError::TooFewSegments(4))));
}

#[test]
fn too_many_segments_rejected() {
    // 7 fillets → 14 segments → above max.
    let bb = synth_blend_boundary(7, 1.0);
    let result = blend_boundary_to_pns_ring(&bb);
    assert!(matches!(result, Err(PnsAdapterError::TooManySegments(14))));
}

#[test]
fn position_endpoints_interpolate_segment_curve() {
    // Verify pos[0] and pos[2] match the segment curve's endpoints.
    let bb = synth_blend_boundary(3, 1.5);
    let ring = blend_boundary_to_pns_ring(&bb).expect("conversion should succeed");
    for (k, strip) in ring.sectors.iter().enumerate() {
        let segment = &bb.segments[k];
        let curve = match segment {
            BoundarySegment::FilletEnd { curve, .. } => curve,
            BoundarySegment::FaceEdge { curve, .. } => curve,
        };
        use popkern_core::geometry::curve::CurveEval;
        let (t0, t1) = curve.parameter_range();
        let expected_start = curve.point_at(t0);
        let expected_end = curve.point_at(t1);
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
        sigma: sigma_center_jet_default(n),
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
