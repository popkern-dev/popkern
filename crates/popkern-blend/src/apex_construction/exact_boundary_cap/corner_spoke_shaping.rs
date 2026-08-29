//! Corner-compatible spoke shaping for exact-boundary apex caps.
//!
//! Nested inside `corner` (shares its private types); extracted to
//! hold the parent under the 500-line module budget.

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;

use super::{
    GENERATED_ANCHOR_RECONCILIATION_TOLERANCE, GeneratedInnerSpoke,
    PreparedExactBoundarySectorBoundary,
};

/// Bend each generated spoke so its ANCHOR-side tangent lies in the
/// cap tangent plane at the anchor — the plane spanned by the real
/// outer-boundary tangents meeting there (PK-125).
///
/// A straight anchor→pole spoke leaves the anchor OUT of that plane
/// (the pole sits above the cap), so the solved sector surface tears
/// away from the support tangent plane exactly at the seam ends
/// (measured pyramid4_uniform: cap-vs-cylinder unit-normal deficit
/// 4.8e-1 at s = 0.05/0.95 on every fillet seam, chart-independent).
/// At a setback corner the meeting boundary tangents both lie in the
/// common support tangent plane (the fillet arcs end on the support
/// tangency lines), so an in-plane spoke start makes the cap tangent
/// plane at the corner EQUAL the support's — G1 at the seam ends by
/// construction. The spoke stays a point-anchored cubic ending
/// straight into the pole; anchors and pole are unchanged, so spoke
/// G0 reconciliation and the pole meeting point are unaffected.
pub(super) fn shape_spokes_for_corner_compatibility(
    sectors: &mut [PreparedExactBoundarySectorBoundary],
) {
    // Collect (endpoint, inward tangent) of every non-degenerate
    // prepared boundary curve end.
    let mut boundary_end_tangents: Vec<(Point3, popkern_core::geometry::Vec3)> = Vec::new();
    for sector in sectors.iter() {
        for boundary in [&sector.fillet_boundary, &sector.face_boundary] {
            let curve = &boundary.curve;
            let start = curve.point_at(0.0);
            let end = curve.point_at(1.0);
            if start.distance_to(&end) <= 1.0e-6 {
                continue;
            }
            let eps = 1.0e-3;
            if let Some(t_start) = (curve.point_at(eps) - start).normalized() {
                boundary_end_tangents.push((start, t_start));
            }
            if let Some(t_end) = (curve.point_at(1.0 - eps) - end).normalized() {
                boundary_end_tangents.push((end, t_end));
            }
        }
    }
    for sector in sectors.iter_mut() {
        for spoke in [&mut sector.fillet_inner_spoke, &mut sector.face_inner_spoke] {
            shape_one_spoke(spoke, &boundary_end_tangents);
        }
    }
}

fn shape_one_spoke(
    spoke: &mut GeneratedInnerSpoke,
    boundary_end_tangents: &[(Point3, popkern_core::geometry::Vec3)],
) {
    let anchor = spoke.anchor;
    let chord = spoke.pole - anchor;
    let Some(chord_dir) = chord.normalized() else {
        return;
    };
    let length = chord.length();
    // The two nearest distinct boundary tangents anchored at this
    // spoke's foot span the cap tangent plane there.
    let mut local: Vec<popkern_core::geometry::Vec3> = boundary_end_tangents
        .iter()
        .filter(|(point, _)| {
            point.distance_to(&anchor) <= GENERATED_ANCHOR_RECONCILIATION_TOLERANCE
        })
        .map(|(_, tangent)| *tangent)
        .collect();
    local.dedup_by(|a, b| a.cross(b).length() < 1.0e-6);
    let [t1, t2, ..] = local.as_slice() else {
        return;
    };
    let Some(normal) = t1.cross(t2).normalized() else {
        return;
    };
    let in_plane = chord - normal * chord.dot(&normal);
    let Some(start_dir) = in_plane.normalized() else {
        return;
    };
    // Keep the straight spoke when the chord is already in-plane.
    if start_dir.cross(&chord_dir).length() < 1.0e-9 {
        return;
    }
    let handle = length / 3.0;
    spoke.curve = NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![
            anchor,
            anchor + start_dir * handle,
            spoke.pole - chord_dir * handle,
            spoke.pole,
        ],
        weights: vec![1.0; 4],
    };
}
