//! Perpendicular anchor position selection for Phase A splitting.

use super::*;
use popkern_core::topology::LoopHandle;

/// Compute the 3D positions of the two perpendicular anchors for
/// a single face. Read-only — this function never mutates `brep`.
/// Called on the Phase A snapshot (un-fragmented BRep), so anchor
/// positions are always computed on the *original* full-length
/// perpendicular edges rather than whatever sub-edge an earlier
/// batch member already split off.
///
/// Walks every loop on the face (outer + inner). The fillet edge
/// can sit on an inner loop when the face has a hole — for example,
/// a cube top face with a cylindrical drilling — so walking only
/// the outer loop would miss the incident perpendicular edges.
/// Place anchor positions on the face boundary at a fillet edge's
/// start and end vertices.
///
/// **Commercial-grade path** (per #136): when this face's stored
/// surface is a plane and we have the fillet's analytic trim
/// curve (`trim_3d_for_face` — already the cyl × plane analytic
/// intersection produced by `compute_trim_curves`), find where
/// that curve crosses each adjacent BRep edge.  The crossing is
/// by construction ON the cylinder, so the resulting trim edges
/// have on-surface endpoints and the cylinder face's outer loop
/// never contains chord-line edges through the cylinder's interior.
///
/// **Fallback path**: walk distance `radius` along the adjacent
/// BRep edge from the fillet edge's vertex.  This is correct for
/// perpendicular faces (cube), but for non-perpendicular faces
/// (pyramid slants) it lands the anchor off the cylinder by
/// `(1 − cos θ)·r`.  Used when no trim curve is provided or the
/// geometry isn't a cyl × plane combination.
#[allow(clippy::too_many_arguments)]
pub(super) fn compute_perp_anchor_positions_for_face(
    brep: &BRep,
    face: FaceHandle,
    fillet_edge: EdgeHandle,
    v_start: VertexHandle,
    v_end: VertexHandle,
    radius_at_start: f64,
    radius_at_end: f64,
    trim_3d_for_face: Option<&NurbsCurve>,
    anchor_on_trim: bool,
    bevel_plane: Option<(Point3, popkern_core::geometry::Vec3)>,
    _tolerance: &Tolerance,
) -> Result<(Point3, Point3), String> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("face {face:?} not found"))?;

    // Use the precomputed analytic trim curve when the face is
    // planar (cyl × plane case).  trim_3d_for_face is already on
    // the cylinder; intersecting it with the adjacent BRep edge
    // gives an on-cyl anchor position.
    let analytic_curve: Option<&NurbsCurve> = match (
        trim_3d_for_face,
        face_data.surface.as_ref().map(|s| &s.geometry),
    ) {
        (Some(c), Some(SurfaceGeometry::Plane { .. })) => Some(c),
        _ => None,
    };

    let mut loop_handles: Vec<LoopHandle> = Vec::with_capacity(1 + face_data.inner_loops.len());
    loop_handles.push(face_data.outer_loop);
    loop_handles.extend(face_data.inner_loops.iter().copied());

    // Slot quality: 2 = exact trim crossing, 1 = perp-walk/bevel.
    // A later edge's CROSSING must outrank an earlier edge's walk
    // placement — slot order follows the loop walk, and on an
    // interrupted rim the crossing edge is not always first
    // (measured: 3 of 4 chamfer anchor slots kept their ~2e-2
    // perp-walk noise because the walk edge came first).
    let mut perp_start: Option<(Point3, u8)> = None;
    let mut perp_end: Option<(Point3, u8)> = None;

    let max_iters = brep.half_edges.len() + 1;
    for loop_handle in loop_handles {
        let start_he = brep
            .loops
            .get(loop_handle)
            .ok_or_else(|| format!("loop {loop_handle:?} not found"))?
            .half_edge;

        // Collect the loop once so the anchor placement can walk
        // COLLINEAR CONTINUATIONS in both directions (a boolean may
        // pre-split a support rim mid-span; the anchor chain must
        // cross those junction vertices).
        let mut loop_hes: Vec<popkern_core::topology::HalfEdgeHandle> = Vec::new();
        {
            let mut cursor = start_he;
            let mut iters = 0;
            loop {
                loop_hes.push(cursor);
                cursor = brep
                    .half_edges
                    .get(cursor)
                    .ok_or_else(|| format!("half-edge {cursor:?} not found"))?
                    .next;
                iters += 1;
                if cursor == start_he || iters >= max_iters {
                    break;
                }
            }
        }
        let seg_of = |heh: popkern_core::topology::HalfEdgeHandle| -> Option<(Point3, Point3)> {
            let he = brep.half_edges.get(heh)?;
            let next_he = brep.half_edges.get(he.next)?;
            let p0 = brep.vertices.get(he.origin)?.position;
            let p1 = brep.vertices.get(next_he.origin)?.position;
            Some((p0, p1))
        };

        let mut cursor = start_he;
        let mut iters = 0;
        let mut loop_index = 0usize;
        loop {
            let he = brep
                .half_edges
                .get(cursor)
                .ok_or_else(|| format!("half-edge {cursor:?} not found"))?;
            let next_he = brep
                .half_edges
                .get(he.next)
                .ok_or_else(|| "next half-edge not found".to_string())?;

            if he.edge != fillet_edge {
                let origin = he.origin;
                let dest = next_he.origin;
                let p0 = brep.vertices[origin].position;
                let p1 = brep.vertices[dest].position;
                let edge_len = p0.distance_to(&p1);
                let along = |t: f64| -> Point3 {
                    Point3::new(
                        p0.x + t * (p1.x - p0.x),
                        p0.y + t * (p1.y - p0.y),
                        p0.z + t * (p1.z - p0.z),
                    )
                };

                // Place an anchor for `vertex` (either v_start or
                // v_end) on this BRep edge, given that `vertex`
                // sits at edge endpoint `origin` (forward) or
                // `dest` (reversed).  Returns the anchor position
                // or `None` if the edge isn't suitable.
                let trace = std::env::var("POPKERN_ANCHOR_TRACE").is_ok();
                let place_anchor = |radius: f64, forward: bool, slot: &mut Option<(Point3, u8)>| {
                    let slot_quality = slot.map(|(_, q)| q).unwrap_or(0);
                    if slot_quality >= 2 || edge_len <= 1e-15 {
                        return;
                    }
                    // Compute perp-walk t — the canonical legacy
                    // position.  See `analytic_curve` doc above:
                    // the on-cyl path computes a different (true
                    // on-surface) anchor for non-perpendicular
                    // faces, but downstream Phase B/C/apex stages
                    // currently assume perp-walk positions.
                    // Promoting to true on-cyl anchors requires
                    // updating those stages — out of scope for
                    // this slice (#136 multi-session work).
                    let t_perp = if forward {
                        (radius / edge_len).min(0.999)
                    } else {
                        1.0 - (radius / edge_len).min(0.999)
                    };
                    let perp_accept = if forward {
                        t_perp > 0.001
                    } else {
                        (0.001..0.999).contains(&t_perp)
                    };

                    // TRUE trim-crossing anchor (the deferred #136
                    // promotion): when the trim curve crosses this
                    // boundary edge, the anchor is that crossing —
                    // the perp-walk position (arc-length `radius`
                    // from the vertex) is where an ARC cross-section
                    // would land and disagrees with the G2 sweep's
                    // contact line (PK-095: uniform 0.06 seam between
                    // every strip and its faces on the pentagon).
                    // Perp-walk math stays as the fallback when the
                    // trim doesn't cross this edge (or wasn't
                    // computed).
                    if anchor_on_trim
                        && let Some(curve) = analytic_curve
                        && let Some((s, _)) =
                            crate::construct::cyl_plane_intersect::intersect_curve_with_segment(
                                curve, p0, p1,
                            )
                        && (0.001..0.999).contains(&s)
                    {
                        if trace {
                            let delta = (s - t_perp).abs();
                            eprintln!(
                                "[anchor] face={face:?} trim-crossing s={s:.6} vs perp t={t_perp:.6} (delta={delta:.4e}) — using trim crossing"
                            );
                        }
                        *slot = Some((along(s), 2));
                        return;
                    }
                    if slot_quality >= 1 {
                        return;
                    }
                    // PLANAR-BEVEL anchor (chamfers): the strip's
                    // true end boundary is where the BEVEL PLANE
                    // crosses this perpendicular edge. On orthogonal
                    // corners this equals the perp-walk position
                    // exactly; on skewed corners the perp-walk
                    // anchor sits OFF the bevel plane, the end
                    // transition gets merged into the planar strip
                    // face anyway, and the face no longer contains
                    // its own boundary (measured: a tilted
                    // boolean-cut bar's mitre corner 0.46 off the
                    // recognized plane — the mesher and OCC's sewing
                    // both cracked along both rails). Capped at 8×
                    // the chamfer distance so a near-parallel edge
                    // falls back to perp-walk instead of a grotesque
                    // mitre.
                    if let Some((po, pn)) = bevel_plane {
                        let d = p1 - p0;
                        let denom = pn.dot(&d);
                        if denom.abs() > 1e-12 {
                            let s = pn.dot(&(po - p0)) / denom;
                            let dist_from_vertex = if forward {
                                s * edge_len
                            } else {
                                (1.0 - s) * edge_len
                            };
                            if (0.001..0.999).contains(&s)
                                && dist_from_vertex > 1e-9
                                && dist_from_vertex <= radius * 8.0
                            {
                                if trace {
                                    eprintln!(
                                        "[anchor] face={face:?} bevel-plane s={s:.6} vs perp t={t_perp:.6} — using plane crossing"
                                    );
                                }
                                *slot = Some((along(s), 1));
                                return;
                            }
                        }
                    }
                    // COLLINEAR-CHAIN CONTINUATION: when the anchor
                    // distance overshoots this edge, the true anchor
                    // lies on a collinear continuation — a prior
                    // boolean may have split the support rim mid-span
                    // (measured: an L-bracket union splits the base
                    // top rim 2.0 from the corner; radius 4 clamped
                    // the anchor onto the junction vertex, the strip
                    // side trim slanted to it, the terminal stub
                    // settle declined on the off-cylinder anchor, and
                    // the stale corner stub shipped with deviation
                    // exactly r(√2−1) — check-9 catch, PK-171). Walk
                    // the loop while edges stay collinear, trying the
                    // trim crossing on each piece and spending the
                    // remaining radius.
                    if radius >= edge_len * 0.999 {
                        let dir0 = (p1 - p0) * (1.0 / edge_len);
                        let mut remaining = radius - edge_len;
                        let n = loop_hes.len();
                        for step in 1..n {
                            let j = if forward {
                                (loop_index + step) % n
                            } else {
                                (loop_index + n - step) % n
                            };
                            let Some((q0, q1)) = seg_of(loop_hes[j]) else {
                                break;
                            };
                            // Orient the piece along the walk.
                            let (a, b) = if forward { (q0, q1) } else { (q1, q0) };
                            let len = a.distance_to(&b);
                            if len <= 1e-15 {
                                continue;
                            }
                            let dir = (b - a) * (1.0 / len);
                            let walk_dir = if forward { dir0 } else { dir0 * -1.0 };
                            if dir.dot(&walk_dir) < 1.0 - 1e-9 {
                                break; // chain ends at a genuine corner
                            }
                            if anchor_on_trim
                                && let Some(curve) = analytic_curve
                                && let Some((s, _)) =
                                    crate::construct::cyl_plane_intersect::intersect_curve_with_segment(
                                        curve, q0, q1,
                                    )
                                && (0.001..0.999).contains(&s)
                            {
                                *slot = Some((
                                    Point3::new(
                                        q0.x + s * (q1.x - q0.x),
                                        q0.y + s * (q1.y - q0.y),
                                        q0.z + s * (q1.z - q0.z),
                                    ),
                                    2,
                                ));
                                return;
                            }
                            if remaining <= len * 0.999 {
                                let f = (remaining / len).max(0.001);
                                *slot = Some((
                                    Point3::new(
                                        a.x + f * (b.x - a.x),
                                        a.y + f * (b.y - a.y),
                                        a.z + f * (b.z - a.z),
                                    ),
                                    1,
                                ));
                                return;
                            }
                            remaining -= len;
                        }
                    }
                    if perp_accept {
                        *slot = Some((along(t_perp), 1));
                    }
                };

                if origin == v_start {
                    place_anchor(radius_at_start, true, &mut perp_start);
                } else if dest == v_start {
                    place_anchor(radius_at_start, false, &mut perp_start);
                }
                if origin == v_end {
                    place_anchor(radius_at_end, true, &mut perp_end);
                } else if dest == v_end {
                    place_anchor(radius_at_end, false, &mut perp_end);
                }
            }

            cursor = he.next;
            iters += 1;
            loop_index += 1;
            if cursor == start_he || iters >= max_iters {
                break;
            }
        }
    }

    let start_pos = perp_start
        .map(|(p, _)| p)
        .ok_or_else(|| format!("no perpendicular edge found at v_start on face {face:?}"))?;
    let end_pos = perp_end
        .map(|(p, _)| p)
        .ok_or_else(|| format!("no perpendicular edge found at v_end on face {face:?}"))?;
    Ok((start_pos, end_pos))
}
