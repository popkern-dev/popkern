//! §8.3 blend cap endpoint extraction.
//!
//! Given a computed [`CenterCurveResult`] and the pair of endpoint
//! radii, produces the two [`CapGeometry`] values a fillet operation
//! needs for its start and end endpoints. Pure helpers extracted
//! from `geometry_compute.rs` to keep that file under the 500-line
//! structural cap; no BRep or topology interaction.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::CurveEval;

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
pub(crate) fn build_edge_caps(
    cap_kind: CapKind,
    center_result: &CenterCurveResult,
    radius_at_start: f64,
    radius_at_end: f64,
) -> (Option<CapGeometry>, Option<CapGeometry>) {
    if matches!(cap_kind, CapKind::Open | CapKind::UserDefined) {
        return (None, None);
    }

    let fp1 = &center_result.footprint_1;
    let fp2 = &center_result.footprint_2;
    if fp1.is_empty() || fp2.is_empty() {
        return (None, None);
    }

    let start_cap = build_endpoint_cap(
        cap_kind,
        center_result,
        fp1[0],
        fp2[0],
        radius_at_start,
        true,
    );
    let last = fp1.len() - 1;
    let last2 = fp2.len() - 1;
    let end_cap = build_endpoint_cap(
        cap_kind,
        center_result,
        fp1[last],
        fp2[last2],
        radius_at_end,
        false,
    );
    (start_cap, end_cap)
}

/// Build one endpoint's cap from the contact points and the center
/// curve. `is_start = true` means the cap is at parameter `t = 0`
/// (the spine's start); `false` means `t = 1` (the spine's end).
/// The spine tangent is reversed at the start so both caps have an
/// "outward" orientation pointing away from the blend interior.
fn build_endpoint_cap(
    cap_kind: CapKind,
    center_result: &CenterCurveResult,
    contact_1: popkern_core::geometry::Point3,
    contact_2: popkern_core::geometry::Point3,
    radius: f64,
    is_start: bool,
) -> Option<CapGeometry> {
    let center_curve = &center_result.center_curve;
    let (t_lo, t_hi) = center_curve.parameter_range();
    let t = if is_start { t_lo } else { t_hi };
    let center = center_curve.point_at(t);
    let tangent = center_curve.derivative_at(t, 1);
    let normalized = tangent.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    // Outward tangent: points away from the blend interior. At the
    // spine start (`t = 0`) the natural tangent points into the
    // blend, so we reverse; at `t = 1` it already points outward.
    let spine_tangent_outward = if is_start {
        normalized * -1.0
    } else {
        normalized
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
