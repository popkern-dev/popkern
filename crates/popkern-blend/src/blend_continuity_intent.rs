//! Seam intent classification for the continuity evidence (PK-592):
//! which live seams are gated for G1/G2 and which are designed creases.

use popkern_core::topology::{BRep, FaceHandle};

use super::{BlendSeamKind, ContinuityGates};

/// Smallest variation of the normal angle along a seam that marks a
/// strip end (degrees). A strip's normal sweeps its cross-section profile
/// along an end, so the crease against the end face varies (measured:
/// flat pyramid ends 5.4°, steep pyramid ends 9.2°); a rail degraded
/// uniformly keeps a constant crease (the PK-498 row-tilt control: ~0°).
const END_PROFILE_SWEEP_DEG: f64 = 1.0;

/// Serde default for a record without a best-sample reading.
pub(super) fn no_tangency() -> f64 {
    f64::INFINITY
}

/// Whether a seam is gated for G1/G2 (intended smooth).
///
/// Blend-to-support seams are classified structurally. A seam is
/// intended smooth, and gated at its worst sample, unless it is a
/// designed crease. It is one when it never runs tangent anywhere along
/// its length and, in addition, one of these holds:
///
/// - it creases past the angle bound (a cube strip end is 90°);
/// - its blend side is planar (a plane-plane fragment is never a
///   rolling-ball rail);
/// - its crease varies along it by [`END_PROFILE_SWEEP_DEG`] or more,
///   the strip normal sweeping its profile along an end.
///
/// A rail is built tangent to its contact face (0.00° at every sample).
/// A rail degraded along its length still touches tangency, and a rail
/// degraded uniformly keeps a constant sub-bound crease on a curved
/// blend, so both stay gated and fail. This replaces the bare 30° bound
/// on these seams, under which every flat pyramid read G0 on designed
/// strip ends at 14–30° (PK-592).
///
/// Blend-to-blend seams keep the angle bound.
pub(super) fn smooth_intent(
    kind: BlendSeamKind,
    worst: &crate::cap_seam_measure::SeamWorst,
    blend_side_planar: bool,
    gates: &ContinuityGates,
) -> bool {
    let within_bound = worst.angle_deg <= gates.smooth_intent_max_angle_deg;
    match kind {
        BlendSeamKind::BlendToSupport if worst.samples > 0 => {
            if worst.best_g1 <= gates.g1_normal {
                return true;
            }
            let best_angle_deg = worst.best_g1.min(1.0).asin().to_degrees();
            let sweep = worst.angle_deg - best_angle_deg;
            within_bound && !blend_side_planar && sweep < END_PROFILE_SWEEP_DEG
        }
        _ => within_bound,
    }
}

/// Whether `face` carries a plane.
pub(super) fn planar(brep: &BRep, face: FaceHandle) -> bool {
    matches!(
        brep.faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry),
        Some(popkern_core::geometry::surface::SurfaceGeometry::Plane { .. })
    )
}
