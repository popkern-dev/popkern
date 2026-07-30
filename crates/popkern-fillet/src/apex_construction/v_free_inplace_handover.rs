//! In-place V-free support-rewrite handover completion (PK-214).
//!
//! On the rational-Gregory direct-apex route the in-place rewrite
//! serves a body whose NEW shell (caps, strips, rewritten support
//! quads) is mirror-wound against the untouched far faces, with the
//! original surface fragments surviving as a GLUED ghost complex:
//! surfaceless far rectangles twinning the far planes, retired
//! fillet-strip bands, corner D-segments still holding the cylinder
//! far-end twins, zero-area corner trigons, and bigon-held diagonal
//! chords inside the production plane loops. The winding
//! reconciliation (v_free_shell_orientation) is the designed repair
//! for exactly this parallel-dead-pair class — but its 2-coloring
//! only sees the anatomy once both halves of each coincident pair
//! are DEAD, and in the served body one half of every pair is still
//! twinned into a ghost. This stage retires the ghost complex so the
//! reconciliation, ring rebuilds, and coincident welds have their
//! designed substrate.
//!
//! Fallback boundary: the caller (the serving-shape gate) applies
//! this on a CLONE and accepts only a defect-free census; anything
//! partial falls through to the recorded aggregate redo unchanged.

use std::collections::BTreeSet;

use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, TopologyRole};

use super::v_free_shell_consolidation::remove_face;

/// Retire the glued mirror-ghost complex: seed on construction
/// residuals (surfaceless membranes, degenerate surfaced slivers),
/// then admit faces — regardless of recorded role — whose every side
/// either twins the complex, is already dead, or twins a live side
/// that has a coincident DEAD counterpart elsewhere (the far-plane
/// terminal-arc class: the cylinder far end twins a corner
/// D-segment while the plane's true arc side lies dead). Returns the
/// number of faces removed.
pub(crate) fn retire_mirror_ghost_complex(brep: &mut BRep) -> usize {
    let mut complex: BTreeSet<FaceHandle> = brep
        .faces
        .iter()
        .filter(|(_, face)| {
            face.role == TopologyRole::ConstructionResidual
                && (face.surface.is_none() || face_is_degenerate_sliver(brep, face.outer_loop))
        })
        .map(|(handle, _)| handle)
        .collect();
    if complex.is_empty() {
        return 0;
    }
    // Grow by the strict admission rule until stable. The rule is
    // deliberately conservative: a face with ANY live twin lacking a
    // dead coincident counterpart stays (production planes twin the
    // ghosts too, but always carry live production seams as well).
    loop {
        let mut grew = false;
        for (candidate, _) in brep.faces.iter().collect::<Vec<_>>() {
            if complex.contains(&candidate) {
                continue;
            }
            let sides = loop_sides(brep, candidate);
            if sides.is_empty() {
                continue;
            }
            let mut touches_complex = false;
            let mut admissible = true;
            for side in &sides {
                match twin_face_of(brep, *side) {
                    None => {}
                    Some(face) if complex.contains(&face) => touches_complex = true,
                    Some(_) => {
                        let Some(twin) = brep.half_edges.get(*side).map(|d| d.twin) else {
                            admissible = false;
                            break;
                        };
                        if !has_dead_coincident_counterpart(brep, twin, &sides) {
                            admissible = false;
                            break;
                        }
                    }
                }
            }
            if touches_complex && admissible {
                complex.insert(candidate);
                grew = true;
            }
        }
        if !grew {
            break;
        }
    }
    let removed = complex.len();
    for face in &complex {
        remove_face(brep, *face);
    }
    // Displaced records: scaffolding carries half-edges whose loop_
    // points at a retired loop while their next-chains run elsewhere
    // (the bigon far twins). They are unreachable by remove_face's
    // loop walk; sweep them by dead-loop membership.
    let displaced: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, data)| !brep.loops.contains_key(data.loop_))
        .map(|(handle, _)| handle)
        .collect();
    for side in displaced {
        let Some(data) = brep.half_edges.remove(side) else {
            continue;
        };
        if brep.half_edges.contains_key(data.twin) {
            brep.boundary_edges.insert(data.edge);
            if let Some(edge) = brep.edges.get_mut(data.edge)
                && edge.half_edge == side
            {
                edge.half_edge = data.twin;
            }
        } else if brep.edges.remove(data.edge).is_some() {
            brep.boundary_edges.remove(&data.edge);
        }
    }
    removed
}

/// Remap the exact-boundary sector ownership records after the
/// handover completion: the retired trigons/bigons were recorded as
/// the cap's outer support sides, and the welds re-twinned the cap
/// boundary directly onto the production support faces. Every seam
/// is re-derived from its live edge and the (stable) cap face;
/// returns false when any seam fails to resolve — the caller then
/// declines the completion so ownership records never go stale.
pub(crate) fn remap_exact_boundary_sector_records(
    brep: &BRep,
    sectors: &mut [super::exact_boundary_cap::sector::ExactBoundaryCapSectorInput],
    report_sectors: &mut [crate::types::ExactBoundaryReportSector],
) -> bool {
    for sector in sectors.iter_mut() {
        for side_index in 0..2 {
            let edge = if side_index == 0 {
                sector.fillet_end_edge
            } else if let Some(face_edge) = sector.face_edge {
                face_edge
            } else {
                // Welded shared corner (PK-219): no bridge side to
                // remap for this sector.
                continue;
            };
            let Some(edge_data) = brep.edges.get(edge) else {
                return false;
            };
            let a = edge_data.half_edge;
            let Some(b) = brep.half_edges.get(a).map(|d| d.twin) else {
                return false;
            };
            if !brep.half_edges.contains_key(b) {
                return false;
            }
            let face_of = |handle: HalfEdgeHandle| -> Option<FaceHandle> {
                brep.loops
                    .get(brep.half_edges.get(handle)?.loop_)
                    .map(|loop_data| loop_data.face)
            };
            let (Some(face_a), Some(face_b)) = (face_of(a), face_of(b)) else {
                return false;
            };
            let (cap_side, support_side, support_face) = if face_a == sector.cap_face {
                (a, b, face_b)
            } else if face_b == sector.cap_face {
                (b, a, face_a)
            } else {
                return false;
            };
            let seam = if side_index == 0 {
                sector.fillet_support_face = support_face;
                &mut sector.fillet_seam
            } else {
                sector.bridge_support_face = Some(support_face);
                let Some(bridge_seam) = sector.bridge_seam.as_mut() else {
                    return false;
                };
                bridge_seam
            };
            seam.cap_half_edge = cap_side;
            seam.support_half_edge = support_side;
            seam.support_face = support_face;
            seam.cap_pcurve_present = brep
                .half_edges
                .get(cap_side)
                .is_some_and(|d| d.pcurve.is_some());
            seam.support_pcurve_present = brep
                .half_edges
                .get(support_side)
                .is_some_and(|d| d.pcurve.is_some());
        }
        // The user-facing report carries its own sector copies —
        // keep them identical to the remapped ownership facts.
        if let Some(report_sector) = report_sectors
            .iter_mut()
            .find(|r| r.sector_index == sector.sector_index)
        {
            report_sector.fillet_support_face = sector.fillet_support_face;
            report_sector.bridge_support_face = sector.bridge_support_face;
            report_sector.cap_half_edges = [
                Some(sector.fillet_seam.cap_half_edge),
                sector.bridge_seam.as_ref().map(|seam| seam.cap_half_edge),
            ];
            report_sector.support_half_edges = [
                Some(sector.fillet_seam.support_half_edge),
                sector
                    .bridge_seam
                    .as_ref()
                    .map(|seam| seam.support_half_edge),
            ];
            report_sector.cap_pcurves_present = [
                sector.fillet_seam.cap_pcurve_present,
                sector
                    .bridge_seam
                    .as_ref()
                    .is_some_and(|seam| seam.cap_pcurve_present),
            ];
            report_sector.support_pcurves_present = [
                sector.fillet_seam.support_pcurve_present,
                sector
                    .bridge_seam
                    .as_ref()
                    .is_some_and(|seam| seam.support_pcurve_present),
            ];
        }
    }
    true
}

/// A surfaced residual whose loop is geometrically a zero-area
/// sliver: every vertex lies on the segment between the loop's two
/// extreme points (the collinear corner trigons pinned at a retired
/// corner vertex).
fn face_is_degenerate_sliver(brep: &BRep, outer_loop: popkern_core::topology::LoopHandle) -> bool {
    let Some(loop_data) = brep.loops.get(outer_loop) else {
        return false;
    };
    let mut positions = Vec::new();
    let start = loop_data.half_edge;
    let mut cursor = start;
    for _ in 0..brep.half_edges.len() + 1 {
        let Some(data) = brep.half_edges.get(cursor) else {
            return false;
        };
        let Some(vertex) = brep.vertices.get(data.origin) else {
            return false;
        };
        positions.push(vertex.position);
        cursor = data.next;
        if cursor == start {
            break;
        }
    }
    if positions.len() < 3 {
        return true;
    }
    let mut best = (0usize, 1usize, 0.0f64);
    for i in 0..positions.len() {
        for j in (i + 1)..positions.len() {
            let d = positions[i].distance_to(&positions[j]);
            if d > best.2 {
                best = (i, j, d);
            }
        }
    }
    if best.2 <= 0.0 {
        return true;
    }
    let (a, b) = (positions[best.0], positions[best.1]);
    let axis = b - a;
    let len2 = axis.dot(&axis);
    let tol = brep.tolerance.linear().max(1.0e-9) * 10.0;
    positions.iter().all(|p| {
        let t = ((*p - a).dot(&axis) / len2).clamp(0.0, 1.0);
        let foot = a + axis * t;
        p.distance_to(&foot) <= tol
    })
}

fn loop_sides(brep: &BRep, face: FaceHandle) -> Vec<HalfEdgeHandle> {
    let Some(face_data) = brep.faces.get(face) else {
        return Vec::new();
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return Vec::new();
    };
    let start = loop_data.half_edge;
    let mut sides = Vec::new();
    let mut cursor = start;
    for _ in 0..brep.half_edges.len() + 1 {
        if !brep.half_edges.contains_key(cursor) {
            return Vec::new();
        }
        sides.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            break;
        }
    }
    sides
}

fn twin_face_of(brep: &BRep, side: HalfEdgeHandle) -> Option<FaceHandle> {
    let twin = brep.half_edges.get(side)?.twin;
    let twin_data = brep.half_edges.get(twin)?;
    brep.loops.get(twin_data.loop_).map(|l| l.face)
}

/// True when some DEAD half-edge elsewhere in the body spans the
/// same two positions as `side` (either direction) within the
/// model's linear tolerance — the evidence that deleting `side`'s
/// face leaves its live twin weldable instead of orphaned.
fn has_dead_coincident_counterpart(
    brep: &BRep,
    side: HalfEdgeHandle,
    own_sides: &[HalfEdgeHandle],
) -> bool {
    let Some((a, b)) = span_of(brep, side) else {
        return false;
    };
    let tol = brep.tolerance.linear().max(1.0e-9) * 10.0;
    brep.half_edges.iter().any(|(handle, data)| {
        if handle == side || own_sides.contains(&handle) {
            return false;
        }
        if brep.half_edges.contains_key(data.twin) {
            return false;
        }
        let Some((c, d)) = span_of(brep, handle) else {
            return false;
        };
        (a.distance_to(&c) <= tol && b.distance_to(&d) <= tol)
            || (a.distance_to(&d) <= tol && b.distance_to(&c) <= tol)
    })
}

fn span_of(
    brep: &BRep,
    side: HalfEdgeHandle,
) -> Option<(
    popkern_core::geometry::Point3,
    popkern_core::geometry::Point3,
)> {
    let data = brep.half_edges.get(side)?;
    let a = brep.vertices.get(data.origin)?.position;
    let next = brep.half_edges.get(data.next)?;
    let b = brep.vertices.get(next.origin)?.position;
    Some((a, b))
}
