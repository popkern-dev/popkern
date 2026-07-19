//! Multi-chord face splits for Phase B (Round 5r+).
//!
//! When a single original face carries multiple fillet trim chords
//! whose 3D lines intersect in the face interior (e.g. all four
//! chords on a cube face filleted on every edge), a naive
//! sequential `split_face` approach fails on the second chord
//! because its anchor pair has been separated onto two different
//! sub-faces by the first chord's split.
//!
//! This module handles that case with a lazy "insert crossing
//! vertex, then split in sub-chord segments" algorithm:
//!
//! 1. Precompute all pairwise segment intersections among the
//!    trim chords on this face.
//! 2. Apply trims in sequence. For each trim, before its own
//!    `split_face`, walk its precomputed intersection list and
//!    — for every *earlier* trim whose chord crossed this one —
//!    locate the matching sub-edge of the earlier trim in the
//!    current BRep and split it at the crossing point. The new
//!    vertex becomes an intermediate anchor on the current trim.
//! 3. Apply `split_face` on each sub-chord segment in turn. With
//!    the intermediate vertices now on the earlier trims' sub-face
//!    boundaries, each sub-chord has both endpoints on one
//!    sub-face and the split succeeds.
//! 4. Pick the LONGEST sub-chord's consumed face as the "true"
//!    strip for this trim. On a cube's inset-square decomposition
//!    the longest sub-chord is the middle chord between the two
//!    inset corners, and its consumed half is the rectangular
//!    fillet strip. Shorter outer sub-chords produce corner
//!    regions that become sphere vertex blends later.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use super::helpers;
use super::phase_split::{FaceTrim, TrimSplitResult, find_face_containing_vertices};
use crate::error::BlendFailure;

/// Round 5r+: given two candidate descendant faces `f_a`/`f_b`,
/// pick the one that descends from `desired_parent`.
///
/// "Descends" is checked by surface identity: sub-faces created
/// by `split_face` inherit the surface of the original face, so
/// if `desired_parent`'s surface in the original BRep matches
/// `f_a`'s current surface, `f_a` is the descendant we want.
/// Falls back to `f_a` if neither matches.
pub(super) fn choose_side_face(
    result: &BRep,
    original: &BRep,
    desired_parent: FaceHandle,
    f_a: FaceHandle,
    f_b: FaceHandle,
) -> FaceHandle {
    if f_a == desired_parent {
        return f_a;
    }
    if f_b == desired_parent {
        return f_b;
    }
    let Some(parent_face) = original.faces.get(desired_parent) else {
        return f_a;
    };
    let Some(parent_surf) = parent_face.surface.as_ref() else {
        return f_a;
    };
    let surf_matches = |fh: FaceHandle| -> bool {
        let Some(face) = result.faces.get(fh) else {
            return false;
        };
        let Some(sof) = face.surface.as_ref() else {
            return false;
        };
        match (&parent_surf.geometry, &sof.geometry) {
            (
                SurfaceGeometry::Plane {
                    origin: o1,
                    normal: n1,
                    ..
                },
                SurfaceGeometry::Plane {
                    origin: o2,
                    normal: n2,
                    ..
                },
            ) => {
                let dot_n = n1.x * n2.x + n1.y * n2.y + n1.z * n2.z;
                if dot_n.abs() < 0.999 {
                    return false;
                }
                let dx = o2.x - o1.x;
                let dy = o2.y - o1.y;
                let dz = o2.z - o1.z;
                let off = dx * n1.x + dy * n1.y + dz * n1.z;
                off.abs() < 1e-6
            }
            _ => false,
        }
    };
    if surf_matches(f_a) {
        f_a
    } else if surf_matches(f_b) {
        f_b
    } else {
        f_a
    }
}

/// Compute the intersection of two 3D line segments.
///
/// Returns `(s, t, point)` where `s ∈ (0,1)` is the parameter
/// along segment `a → b` and `t ∈ (0,1)` the parameter along
/// `c → d`. Returns `None` for parallel, skew, or endpoint-only
/// intersections.
pub(super) fn intersect_segments(
    a: Point3,
    b: Point3,
    c: Point3,
    d: Point3,
) -> Option<(f64, f64, Point3)> {
    let d1 = b - a;
    let d2 = d - c;
    let cross = d1.cross(&d2);
    let cross_sq = cross.dot(&cross);
    if cross_sq < 1e-18 {
        return None;
    }
    let diff = c - a;
    let s = diff.cross(&d2).dot(&cross) / cross_sq;
    let t = diff.cross(&d1).dot(&cross) / cross_sq;
    let interior = |x: f64| x > 1e-4 && x < 1.0 - 1e-4;
    if !interior(s) || !interior(t) {
        return None;
    }
    let pa = Point3::new(a.x + s * d1.x, a.y + s * d1.y, a.z + s * d1.z);
    let pb = Point3::new(c.x + t * d2.x, c.y + t * d2.y, c.z + t * d2.z);
    let gap_sq = (pa.x - pb.x).powi(2) + (pa.y - pb.y).powi(2) + (pa.z - pb.z).powi(2);
    if gap_sq > 1e-10 {
        return None;
    }
    Some((s, t, pa))
}

/// Find the sub-edge in `sub_edges` whose line segment contains
/// `point` (within tolerance), returning the edge and its
/// parameter along the edge in `(0.001, 0.999)`.
fn find_sub_edge_containing_point(
    brep: &BRep,
    sub_edges: &[EdgeHandle],
    point: Point3,
) -> Option<(EdgeHandle, f64)> {
    for &eh in sub_edges {
        let edata = brep.edges.get(eh)?;
        let he = brep.half_edges.get(edata.half_edge)?;
        let twin = brep.half_edges.get(he.twin)?;
        let v0 = brep.vertices[he.origin].position;
        let v1 = brep.vertices[twin.origin].position;
        let dx = v1.x - v0.x;
        let dy = v1.y - v0.y;
        let dz = v1.z - v0.z;
        let len_sq = dx * dx + dy * dy + dz * dz;
        if len_sq < 1e-18 {
            continue;
        }
        let t = ((point.x - v0.x) * dx + (point.y - v0.y) * dy + (point.z - v0.z) * dz) / len_sq;
        if t <= 1e-4 || t >= 1.0 - 1e-4 {
            continue;
        }
        let px = v0.x + t * dx;
        let py = v0.y + t * dy;
        let pz = v0.z + t * dz;
        let dist_sq = (point.x - px).powi(2) + (point.y - py).powi(2) + (point.z - pz).powi(2);
        if dist_sq < 1e-10 {
            return Some((eh, t));
        }
    }
    None
}

/// Replace `old` with `[new_a, new_b]` in every `Vec<EdgeHandle>`
/// in `per_trim_sub_edges`. Used after `split_edge` to keep the
/// per-trim sub-edge bookkeeping in sync with the BRep.
fn remap_split(
    per_trim_sub_edges: &mut [Vec<EdgeHandle>],
    old: EdgeHandle,
    new_a: EdgeHandle,
    new_b: EdgeHandle,
) {
    for list in per_trim_sub_edges.iter_mut() {
        if let Some(pos) = list.iter().position(|&e| e == old) {
            list.splice(pos..=pos, [new_a, new_b]);
        }
    }
}

/// Outcome of applying all trims on one original face via
/// `apply_multi_chord_for_face`.
pub(super) struct MultiChordOutcome {
    pub trim_results: Vec<TrimSplitResult>,
    pub failures: Vec<(EdgeHandle, BlendFailure)>,
    /// All new descendant edges produced by `split_edge` calls
    /// inside the multi-chord algorithm (splitting an earlier
    /// trim's sub-chord at a crossing point). The caller folds
    /// these into its fillet-edge descendants map, so Phase C
    /// sees them as valid descendants when identifying consumed
    /// faces.
    pub extra_edge_splits: Vec<(EdgeHandle, EdgeHandle, EdgeHandle)>,
}

/// Apply all of a face's trims, honouring chord-chord
/// intersections. This is the Round 5r+ replacement for the
/// per-face inner loop in Phase B.
///
/// `fillet_edges_by_trim[k]` is the canonical fillet edge for
/// `trims[k]`. `descendants` is the per-fillet-edge descendants
/// map — used here in read-only fashion to identify the strip via
/// `helpers::identify_consumed_face_any`.
pub(super) fn apply_multi_chord_for_face(
    brep: &mut BRep,
    original_face: FaceHandle,
    trims: &[FaceTrim],
    fillet_edges_by_trim: &[EdgeHandle],
    descendants: &std::collections::HashMap<EdgeHandle, Vec<EdgeHandle>>,
    tolerance: &Tolerance,
) -> MultiChordOutcome {
    let mut trim_results: Vec<TrimSplitResult> = Vec::new();
    let mut failures: Vec<(EdgeHandle, BlendFailure)> = Vec::new();
    let mut extra_edge_splits: Vec<(EdgeHandle, EdgeHandle, EdgeHandle)> = Vec::new();

    let n = trims.len();
    if n == 0 {
        return MultiChordOutcome {
            trim_results,
            failures,
            extra_edge_splits,
        };
    }

    // Precompute pairwise chord intersections.
    let segments: Vec<(Point3, Point3)> = trims
        .iter()
        .map(|t| {
            (
                brep.vertices[t.anchors.0].position,
                brep.vertices[t.anchors.1].position,
            )
        })
        .collect();

    let mut intersections: Vec<Vec<(f64, Point3, usize)>> = vec![Vec::new(); n];
    for i in 0..n {
        for j in (i + 1)..n {
            let (a, b) = segments[i];
            let (c, d) = segments[j];
            if let Some((s, t, pt)) = intersect_segments(a, b, c, d) {
                intersections[i].push((s, pt, j));
                intersections[j].push((t, pt, i));
            }
        }
    }
    for list in intersections.iter_mut() {
        list.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    }

    let mut per_trim_sub_edges: Vec<Vec<EdgeHandle>> = vec![Vec::new(); n];
    let mut sub_faces: Vec<FaceHandle> = vec![original_face];

    for i in 0..n {
        let trim = &trims[i];
        let fillet_edge = fillet_edges_by_trim[i];
        let fillet_desc: &[EdgeHandle] = descendants
            .get(&fillet_edge)
            .map(|v| v.as_slice())
            .unwrap_or(std::slice::from_ref(&fillet_edge));

        // Step 1: for every precomputed intersection with an
        // earlier trim `j < i`, split the earlier trim's current
        // sub-edge at the crossing point to create an
        // intermediate vertex.
        let mut intermediate_vertices: Vec<VertexHandle> = Vec::new();
        for &(_, pt, j) in &intersections[i] {
            if j >= i {
                continue;
            }
            let Some((eh, t_on)) = find_sub_edge_containing_point(brep, &per_trim_sub_edges[j], pt)
            else {
                continue;
            };
            let t_clamped = t_on.clamp(0.001, 0.999);
            match brep.split_edge(eh, t_clamped) {
                Ok(split) => {
                    extra_edge_splits.push((eh, split.edge_a, split.edge_b));
                    remap_split(&mut per_trim_sub_edges, eh, split.edge_a, split.edge_b);
                    intermediate_vertices.push(split.vertex);
                }
                Err(_) => {
                    // Fall through — the sub-chord below will
                    // still be attempted without this crossing
                    // vertex and likely fail cleanly.
                }
            }
        }

        // Step 2: build the ordered vertex sequence for this
        // trim and apply one `split_face` per sub-chord.
        let mut endpoints: Vec<VertexHandle> = Vec::with_capacity(intermediate_vertices.len() + 2);
        endpoints.push(trim.anchors.0);
        endpoints.extend(intermediate_vertices.iter().copied());
        endpoints.push(trim.anchors.1);

        let mut produced_edges: Vec<EdgeHandle> = Vec::new();
        // Candidates: (length, consumed, pcurve, trim_3d).
        // We pick the LONGEST sub-chord's consumed as the true
        // strip for this fillet edge; shorter outer sub-chords
        // produce corner regions.
        let mut strip_candidates: Vec<(f64, FaceHandle, NurbsCurve, NurbsCurve)> = Vec::new();
        let mut sub_chord_failed = false;

        for k in 0..endpoints.len() - 1 {
            let start_v = endpoints[k];
            let end_v = endpoints[k + 1];
            let target = find_face_containing_vertices(brep, &sub_faces, start_v, end_v);
            let sub_midpoint = {
                let p_a = brep.vertices[start_v].position;
                let p_b = brep.vertices[end_v].position;
                Some(Point3::new(
                    (p_a.x + p_b.x) * 0.5,
                    (p_a.y + p_b.y) * 0.5,
                    (p_a.z + p_b.z) * 0.5,
                ))
            };
            let Some(target) = target else {
                failures.push((
                    fillet_edge,
                    BlendFailure::InvalidTopology {
                        reason: format!(
                            "multi-chord sub-chord anchors not found on any sub-face of {original_face:?}"
                        ),
                        location: sub_midpoint,
                    },
                ));
                sub_chord_failed = true;
                break;
            };

            let p_a = brep.vertices[start_v].position;
            let p_b = brep.vertices[end_v].position;
            // #136 slice 6 / PK-139: the sub-chord's 3-D trim is the
            // SUB-SPAN of the fillet's actual trim curve, not the
            // straight chord between the anchors — a by-points
            // variable law's rail bulges between its breakpoints and
            // the chord line cut straight through it (measured: rail
            // Line at offset 0.05 under a 0.12 bulge apex, 40 welded
            // opens). Straight trims still canonicalize back to
            // analytic lines inside `split_face`. Falls back to the
            // chord when the anchors do not invert onto the trim
            // (corner sub-chords off the trim span).
            let t3d = trim_subcurve_between(&trim.trim_3d, &p_a, &p_b, tolerance)
                .unwrap_or_else(|| helpers::line_nurbs(&p_a, &p_b));
            let pc = match helpers::pcurve_for_face_along(brep, target, &t3d) {
                Ok(pc) => pc,
                Err(e) => {
                    failures.push((
                        fillet_edge,
                        BlendFailure::InvalidTopology {
                            reason: e,
                            location: sub_midpoint,
                        },
                    ));
                    sub_chord_failed = true;
                    break;
                }
            };

            let split = match brep.split_face(target, &pc, &t3d, tolerance) {
                Ok(s) => s,
                Err(e) => {
                    failures.push((
                        fillet_edge,
                        BlendFailure::InvalidTopology {
                            reason: format!("face split failed: {e}"),
                            location: sub_midpoint,
                        },
                    ));
                    sub_chord_failed = true;
                    break;
                }
            };
            produced_edges.extend(split.new_edges.iter().copied());

            // identify_consumed_face_any matches any descendant
            // sub-edge. For the full-strip case this returns the
            // strip. For outer sub-chords on a multi-fillet corner
            // it returns a corner region (which also contains a
            // fillet-edge descendant via the corner's trim).
            let (kept, consumed) = match helpers::identify_consumed_face_any(
                brep,
                split.face_a,
                split.face_b,
                fillet_desc,
            ) {
                Ok(kc) => kc,
                Err(_) => {
                    // Neither half contains any descendant —
                    // this sub-chord is interior. Both halves
                    // are valid sub-faces going forward, but
                    // neither is a strip candidate.
                    if let Some(pos) = sub_faces.iter().position(|&f| f == target) {
                        sub_faces[pos] = split.face_a;
                        sub_faces.push(split.face_b);
                    }
                    continue;
                }
            };

            if let Some(pos) = sub_faces.iter().position(|&f| f == target) {
                sub_faces[pos] = kept;
                sub_faces.push(consumed);
            }
            let length = p_a.distance_to(&p_b);
            strip_candidates.push((length, consumed, pc, t3d));
        }

        if sub_chord_failed {
            per_trim_sub_edges[i] = produced_edges;
            continue;
        }
        per_trim_sub_edges[i] = produced_edges.clone();

        // Pick the LONGEST sub-chord's consumed face as the true
        // strip for this fillet edge.
        strip_candidates.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(std::cmp::Ordering::Equal));
        let Some((_, consumed, pcurve, trim_3d)) = strip_candidates.into_iter().next() else {
            failures.push((
                fillet_edge,
                BlendFailure::InvalidTopology {
                    reason: "multi-chord trim produced no sub-face containing the fillet edge"
                        .to_string(),
                    location: crate::geometry_compute::edge_midpoint_by_handle(brep, fillet_edge),
                },
            ));
            continue;
        };

        trim_results.push(TrimSplitResult {
            edge_idx: trim.edge_idx,
            side: trim.side,
            consumed,
            trim_edges: produced_edges,
            pcurve,
            trim_3d,
        });
    }

    MultiChordOutcome {
        trim_results,
        failures,
        extra_edge_splits,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intersect_crossing_segments() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 1.0, 0.0);
        let c = Point3::new(0.0, 1.0, 0.0);
        let d = Point3::new(1.0, 0.0, 0.0);
        let result = intersect_segments(a, b, c, d).unwrap();
        assert!((result.0 - 0.5).abs() < 1e-9);
        assert!((result.1 - 0.5).abs() < 1e-9);
        assert!((result.2.x - 0.5).abs() < 1e-9);
        assert!((result.2.y - 0.5).abs() < 1e-9);
    }

    #[test]
    fn parallel_segments_no_intersection() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let c = Point3::new(0.0, 1.0, 0.0);
        let d = Point3::new(1.0, 1.0, 0.0);
        assert!(intersect_segments(a, b, c, d).is_none());
    }

    #[test]
    fn endpoint_touching_segments_no_intersection() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let c = Point3::new(1.0, 0.0, 0.0);
        let d = Point3::new(1.0, 1.0, 0.0);
        assert!(intersect_segments(a, b, c, d).is_none());
    }

    #[test]
    fn skew_segments_no_intersection() {
        let a = Point3::new(0.0, 0.0, 0.0);
        let b = Point3::new(1.0, 0.0, 0.0);
        let c = Point3::new(0.5, -1.0, 1.0);
        let d = Point3::new(0.5, 1.0, 1.0);
        assert!(intersect_segments(a, b, c, d).is_none());
    }
}

/// Sub-span of the fillet's actual 3-D trim curve between two anchor
/// positions, as a uniform degree-1 polyline. `None` when either
/// anchor does not lie on the trim (corner sub-chords) or the span is
/// effectively straight (the chord is already exact).
fn trim_subcurve_between(
    trim: &NurbsCurve,
    p_a: &Point3,
    p_b: &Point3,
    tolerance: &Tolerance,
) -> Option<NurbsCurve> {
    use popkern_core::geometry::CurveEval;
    use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
    let (t0, t1) = trim.parameter_range();
    let on_trim = CurveOnEdge {
        geometry: CurveGeometry::Nurbs(trim.clone()),
        t_start: t0,
        t_end: t1,
    };
    let gate = tolerance.linear().max(1.0e-9) * 10.0;
    let ta =
        popkern_core::geometry::curve_inversion::invert_curve_at_point(&on_trim, *p_a, tolerance)
            .ok()?;
    let tb =
        popkern_core::geometry::curve_inversion::invert_curve_at_point(&on_trim, *p_b, tolerance)
            .ok()?;
    if trim.point_at(ta).distance_to(p_a) > gate || trim.point_at(tb).distance_to(p_b) > gate {
        return None;
    }
    let (lo, hi) = if ta <= tb { (ta, tb) } else { (tb, ta) };
    if (hi - lo).abs() < 1.0e-12 {
        return None;
    }
    let n = 17;
    let mut pts = Vec::with_capacity(n);
    for k in 0..n {
        let t = lo + (hi - lo) * k as f64 / (n - 1) as f64;
        pts.push(trim.point_at(t));
    }
    if ta > tb {
        pts.reverse();
    }
    // Endpoints must be EXACTLY the anchors (split_face matches them
    // to vertices).
    pts[0] = *p_a;
    let last = pts.len() - 1;
    pts[last] = *p_b;
    // Straight span: let the chord path (and its analytic
    // canonicalization) serve it.
    let chord = *p_b - *p_a;
    let chord_len = chord.length();
    if chord_len > 1.0e-12 {
        let dir = chord * (1.0 / chord_len);
        let max_dev = pts
            .iter()
            .map(|p| {
                let d = *p - *p_a;
                (d - dir * d.dot(&dir)).length()
            })
            .fold(0.0_f64, f64::max);
        if max_dev <= gate {
            return None;
        }
    }
    let n_pts = pts.len();
    let mut knots = Vec::with_capacity(n_pts + 2);
    knots.push(0.0);
    for k in 0..n_pts {
        knots.push(k as f64 / (n_pts - 1) as f64);
    }
    knots.push(1.0);
    Some(NurbsCurve {
        degree: 1,
        knots,
        control_points: pts,
        weights: vec![1.0; n_pts],
    })
}
