//! Seam intent classification for the continuity evidence (PK-592):
//! which live seams are gated for G1/G2 and which are designed creases.

use std::collections::BTreeSet;

use popkern_core::topology::{BRep, FaceHandle};

use crate::cap_seam_measure::EdgeSeamFaces;

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

/// Seam intent over one blend: the classification rules plus the
/// blend's corner-cap interior (PK-609).
pub(super) struct SeamIntent<'a> {
    brep: &'a BRep,
    gates: &'a ContinuityGates,
    /// Blend faces every one of whose edge neighbours is a blend face:
    /// corner-cap pieces (KP20 fan pieces, fitted, direct-apex and
    /// vertex-blend sectors). Fillet strips and chamfer faces border
    /// their supports, so they are never interior.
    interior: BTreeSet<FaceHandle>,
}

impl<'a> SeamIntent<'a> {
    pub(super) fn new(
        brep: &'a BRep,
        blend: &BTreeSet<FaceHandle>,
        gates: &'a ContinuityGates,
    ) -> Self {
        let interior = blend
            .iter()
            .copied()
            .filter(|&f| {
                edge_neighbours(brep, f).is_some_and(|ns| ns.iter().all(|n| blend.contains(n)))
            })
            .collect();
        Self {
            brep,
            gates,
            interior,
        }
    }

    /// Whether the seam is gated for G1/G2 (intended smooth).
    ///
    /// A seam between two sectors of one corner cap is smooth by
    /// construction: no cap is designed with an internal crease. It is
    /// gated at any angle — the angle bound let a direct-apex cube cap
    /// whose sector seam certifies a normal sine of 0.9998 read as a
    /// designed crease and the blend as G1 (PK-609). Two planar pieces
    /// (a chamfer's corner patch) stay designed creases. Every other
    /// seam follows [`smooth_intent`].
    pub(super) fn gated(
        &self,
        kind: BlendSeamKind,
        worst: &crate::cap_seam_measure::SeamWorst,
        adjacency: &EdgeSeamFaces,
        blend_face: FaceHandle,
    ) -> bool {
        let (a, b) = (adjacency.face_a, adjacency.face_b);
        if kind == BlendSeamKind::BlendToBlend
            && self.interior.contains(&a)
            && self.interior.contains(&b)
            && !(planar(self.brep, a) && planar(self.brep, b))
        {
            return true;
        }
        smooth_intent(kind, worst, planar(self.brep, blend_face), self.gates)
    }
}

/// Every face sharing an edge with `face` (all loops), or `None` when
/// a loop is broken or an edge has no second side.
fn edge_neighbours(brep: &BRep, face: FaceHandle) -> Option<Vec<FaceHandle>> {
    let data = brep.faces.get(face)?;
    let mut out = Vec::new();
    for &lh in std::iter::once(&data.outer_loop).chain(&data.inner_loops) {
        let start = brep.loops.get(lh)?.half_edge;
        let mut cursor = start;
        for _ in 0..=brep.half_edges.len() {
            let he = brep.half_edges.get(cursor)?;
            let twin = brep.half_edges.get(he.twin)?;
            out.push(brep.loops.get(twin.loop_)?.face);
            cursor = he.next;
            if cursor == start {
                break;
            }
        }
    }
    Some(out)
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
