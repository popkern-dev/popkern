//! PK-163 — material-adding attach for THROUGH gaps.
//!
//! A through gap's strip end arc lies on an end face's SURFACE but
//! over the void (outside its trim): the gap region is open at that
//! end and the closed-pocket end split declines. The commercial
//! shape adds the blend material through the opening and EXTENDS the
//! end cap down to the arc:
//!
//! 1. the doomed-region flood treats the end face as a BARRIER (the
//!    region is bounded there by the pre-existing notch-profile
//!    edges instead of a split edge);
//! 2. after the region merges into one face `D`, `D` is `mef`-split
//!    at each through end between the two end contact vertices with
//!    the section end ARC as the new edge — the sliver piece that
//!    emerges (bounded by arc + notch profile) is exactly the added
//!    material's end cross-section;
//! 3. the sliver is kef-merged into the end face across the profile
//!    edges (coplanar extension; leftover same-face profile edges
//!    dissolve as spur slits), leaving the arc as the strip × end
//!    face boundary and `D`'s loop as the four strip borders.
//!
//! Declines are structured (`AttachDeclined`) and the caller's
//! transactional clone contract is unchanged.

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use crate::error::BlendError;
use crate::face_blend_attach_geometry::decline;

/// One through end awaiting the cap extension.
pub(crate) struct ThroughEnd {
    /// The end face whose trim the arc escaped (the flood barrier).
    pub end_face: FaceHandle,
    /// The section end arc, contact-A → contact-B.
    pub arc: Vec<Point3>,
}

/// Extend one through end: mef `strip_region` along the arc, merge
/// the sliver into the end face, dissolve the profile leftovers.
/// Returns the surviving strip-region face (the mef may leave the
/// sliver on the original handle, which then merges into the cap).
pub(crate) fn extend_through_end(
    work: &mut BRep,
    strip_region: FaceHandle,
    end: &ThroughEnd,
    tolerance: &Tolerance,
) -> Result<FaceHandle, BlendError> {
    let tol = tolerance.linear().max(1.0e-9) * 100.0;
    let first = *end
        .arc
        .first()
        .ok_or_else(|| decline("cap extension", "end arc is empty".into()))?;
    let last = *end
        .arc
        .last()
        .ok_or_else(|| decline("cap extension", "end arc is empty".into()))?;
    // The two contact vertices in the region's loop.
    let find_he = |work: &BRep, p: Point3| -> Option<HalfEdgeHandle> {
        let l = work.loops.get(work.faces.get(strip_region)?.outer_loop)?;
        let start = l.half_edge;
        let mut cur = start;
        for _ in 0..work.half_edges.len() + 1 {
            let he = work.half_edges.get(cur)?;
            if work.vertices.get(he.origin)?.position.distance_to(&p) < tol {
                return Some(cur);
            }
            cur = he.next;
            if cur == start {
                return None;
            }
        }
        None
    };
    let he_a = find_he(work, first).ok_or_else(|| {
        decline(
            "cap extension",
            format!("no region-loop vertex at the arc start {first:?}"),
        )
    })?;
    let he_b = find_he(work, last).ok_or_else(|| {
        decline(
            "cap extension",
            format!("no region-loop vertex at the arc end {last:?}"),
        )
    })?;
    let mef = work
        .mef(he_a, he_b)
        .map_err(|e| decline("cap extension", format!("mef along the end arc: {e}")))?;
    // Author the arc curve on the new edge (primary direction:
    // origin(he_a) = arc start).
    let poly = crate::face_blend_attach_geometry::uniform_polyline_nurbs(&end.arc);
    let (t0, t1) = {
        use popkern_core::geometry::CurveEval;
        poly.parameter_range()
    };
    if let Some(e) = work.edges.get_mut(mef.edge) {
        e.curve = Some(popkern_core::geometry::CurveOnEdge {
            geometry: popkern_core::geometry::CurveGeometry::Nurbs(poly),
            t_start: t0,
            t_end: t1,
        });
    }
    // The sliver = whichever of the two pieces borders the end face.
    let borders_end_face = |work: &BRep, f: FaceHandle| -> bool {
        let Some(fd) = work.faces.get(f) else {
            return false;
        };
        let Some(l) = work.loops.get(fd.outer_loop) else {
            return false;
        };
        let start = l.half_edge;
        let mut cur = start;
        for _ in 0..work.half_edges.len() + 1 {
            let Some(he) = work.half_edges.get(cur) else {
                return false;
            };
            let Some(tw) = work.half_edges.get(he.twin) else {
                return false;
            };
            if work
                .loops
                .get(tw.loop_)
                .map(|l| l.face == end.end_face)
                .unwrap_or(false)
            {
                return true;
            }
            cur = he.next;
            if cur == start {
                return false;
            }
        }
        false
    };
    let (sliver, main) = if borders_end_face(work, mef.face) {
        (mef.face, strip_region)
    } else if borders_end_face(work, strip_region) {
        (strip_region, mef.face)
    } else {
        return Err(decline(
            "cap extension",
            "neither mef piece borders the end face".into(),
        ));
    };
    // Coplanar merge: kef the first profile edge (survivor = the end
    // face side), then dissolve leftover same-face profile spurs.
    let profile_he = {
        let fd = work
            .faces
            .get(end.end_face)
            .ok_or_else(|| decline("cap extension", "end face vanished".into()))?;
        let mut found = None;
        let mut loops = vec![fd.outer_loop];
        loops.extend(fd.inner_loops.iter().copied());
        'outer: for lh in loops {
            let Some(l) = work.loops.get(lh) else {
                continue;
            };
            let start = l.half_edge;
            let mut cur = start;
            for _ in 0..work.half_edges.len() + 1 {
                let Some(he) = work.half_edges.get(cur) else {
                    break;
                };
                let Some(tw) = work.half_edges.get(he.twin) else {
                    break;
                };
                if work
                    .loops
                    .get(tw.loop_)
                    .map(|l| l.face == sliver)
                    .unwrap_or(false)
                {
                    found = Some(cur);
                    break 'outer;
                }
                cur = he.next;
                if cur == start {
                    break;
                }
            }
        }
        found.ok_or_else(|| decline("cap extension", "no profile edge to merge across".into()))?
    };
    work.merge_faces_across_edge(profile_he)
        .map_err(|e| decline("cap extension", format!("profile merge: {e}")))?;
    // Tip-dissolve the leftover profile spurs on the end face.
    while let Some(fd) = work.faces.get(end.end_face) {
        let mut tip: Option<HalfEdgeHandle> = None;
        let mut loops = vec![fd.outer_loop];
        loops.extend(fd.inner_loops.iter().copied());
        'outer: for lh in loops {
            let Some(l) = work.loops.get(lh) else {
                continue;
            };
            let start = l.half_edge;
            let mut cur = start;
            for _ in 0..work.half_edges.len() + 1 {
                let Some(he) = work.half_edges.get(cur) else {
                    break;
                };
                if he.next == he.twin {
                    tip = Some(cur);
                    break 'outer;
                }
                cur = he.next;
                if cur == start {
                    break;
                }
            }
        }
        let Some(he) = tip else { break };
        let tip_vertex = {
            let Some(hed) = work.half_edges.get(he) else {
                break;
            };
            let Some(tw) = work.half_edges.get(hed.twin) else {
                break;
            };
            tw.origin
        };
        work.dissolve_same_loop_bridge(he)
            .map_err(|e| decline("cap extension", format!("profile spur dissolve: {e}")))?;
        if work.half_edges.iter().all(|(_, d)| d.origin != tip_vertex) {
            work.vertices.remove(tip_vertex);
        }
    }
    // Refresh the end face's pcurves (its boundary gained the arc).
    let _ = work.recompute_face_pcurves(end.end_face);
    Ok(main)
}
