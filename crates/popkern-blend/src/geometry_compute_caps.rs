//! §8.3 blend cap endpoint extraction.
//!
//! Given a computed [`CenterCurveResult`] and the pair of endpoint
//! radii, produces the two [`CapGeometry`] values a fillet operation
//! needs for its start and end endpoints. Pure helpers extracted
//! from `geometry_compute.rs` to keep that file under the 500-line
//! structural cap; no BRep or topology interaction.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::{Point3, Vec3};

use crate::cap::{self, CapEndpoint, CapGeometry};
use crate::center_curve::CenterCurveResult;
use crate::types::CapKind;

/// Build §8.3 cap geometries at both endpoints of a blend from the
/// computed center curve + footprints.
///
/// Returns `(None, None)` for [`CapKind::Open`] and
/// [`CapKind::UserDefined`] (the latter is a [BEYOND] placeholder
/// that currently falls back to open). All other kinds produce
/// concrete cap geometry at both endpoints — the cap is fully
/// determined by the endpoint contact data regardless of whether
/// that endpoint is a "free boundary" in the BRep sense. The final
/// topology insertion pass decides which caps to actually insert.
#[allow(clippy::too_many_arguments)]
pub(crate) fn build_edge_caps(
    cap_kind: CapKind,
    center_result: &CenterCurveResult,
    radius_at_start: f64,
    radius_at_end: f64,
    vertex_start_pos: Point3,
    vertex_end_pos: Point3,
    trim_3d_1: &NurbsCurve,
    trim_3d_2: &NurbsCurve,
) -> (Option<CapGeometry>, Option<CapGeometry>) {
    if matches!(cap_kind, CapKind::Open | CapKind::UserDefined) {
        return (None, None);
    }

    // Contacts and centre come from each curve's point at the
    // vertex's station, not from the curve endpoints: trims and the
    // centre curve may be extended past the edge (on an open shell
    // they overshoot by a full unit at both ends, which put a free
    // end's flat cap 1.0 off its own boundary, PK-530). The station
    // is the curve's closest point to the vertex.
    let at_vertex = |curve: &NurbsCurve, vertex: Point3| -> Point3 {
        let t = curve.closest_point(&vertex).unwrap_or_else(|| {
            let (lo, hi) = curve.parameter_range();
            if curve.point_at(lo).distance_to(&vertex) <= curve.point_at(hi).distance_to(&vertex) {
                lo
            } else {
                hi
            }
        });
        curve.point_at(t)
    };
    let start_cap = build_endpoint_cap(
        cap_kind,
        center_result,
        at_vertex(trim_3d_1, vertex_start_pos),
        at_vertex(trim_3d_2, vertex_start_pos),
        radius_at_start,
        (vertex_start_pos, vertex_end_pos),
    );
    let end_cap = build_endpoint_cap(
        cap_kind,
        center_result,
        at_vertex(trim_3d_1, vertex_end_pos),
        at_vertex(trim_3d_2, vertex_end_pos),
        radius_at_end,
        (vertex_end_pos, vertex_start_pos),
    );
    (start_cap, end_cap)
}

/// Build one endpoint's cap from the contact points and the center
/// curve. `ends` is `(this endpoint's vertex, the other end's
/// vertex)`: the centre is the centre curve's point at this vertex's
/// station, and the outward spine tangent points away from the other
/// end, whatever the curve's parameter direction.
fn build_endpoint_cap(
    cap_kind: CapKind,
    center_result: &CenterCurveResult,
    contact_1: Point3,
    contact_2: Point3,
    radius: f64,
    ends: (Point3, Point3),
) -> Option<CapGeometry> {
    let center_curve = &center_result.center_curve;
    let (t_lo, t_hi) = center_curve.parameter_range();
    let (vertex, other) = ends;
    let t = center_curve.closest_point(&vertex).unwrap_or_else(|| {
        if center_curve.point_at(t_lo).distance_to(&vertex)
            <= center_curve.point_at(t_hi).distance_to(&vertex)
        {
            t_lo
        } else {
            t_hi
        }
    });
    let center = center_curve.point_at(t);
    let tangent = center_curve.derivative_at(t, 1);
    let normalized = tangent.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    // Outward: away from the blend interior, i.e. from the other end
    // toward this one.
    let spine_tangent_outward = if normalized.dot(&(vertex - other)) >= 0.0 {
        normalized
    } else {
        normalized * -1.0
    };

    // Support normals: the rolling ball touches each support with
    // the outward surface normal pointing from the ball centre
    // toward the contact point. So `(contact - center) / r` is the
    // unit outward normal of that support at the contact.
    let inv_r = 1.0 / radius.max(1e-30);
    let normal_1 = Vec3::new(
        (contact_1.x - center.x) * inv_r,
        (contact_1.y - center.y) * inv_r,
        (contact_1.z - center.z) * inv_r,
    )
    .normalized()
    .unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let normal_2 = Vec3::new(
        (contact_2.x - center.x) * inv_r,
        (contact_2.y - center.y) * inv_r,
        (contact_2.z - center.z) * inv_r,
    )
    .normalized()
    .unwrap_or(Vec3::new(1.0, 0.0, 0.0));

    let endpoint = CapEndpoint {
        contact_1,
        contact_2,
        normal_1,
        normal_2,
        center,
        spine_tangent_outward,
        radius,
    };

    cap::build_cap(cap_kind, &endpoint)
}
