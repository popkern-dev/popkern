//! Low-valence PnS rim materialization: loop/ring alignment, rim
//! curve + vertex authoring, analytic cap pcurve authoring, and the
//! measured sub-face/patch assignment. Split from `pns_promote.rs`
//! (the promote/dispatch entry points), which re-exports what its
//! callers use. See each function's doc for the measured provenance —
//! this family is the PK-192/PK-203 campaign's construction record.

use popkern_core::error::CoreError;
use popkern_core::topology::{BRep, FaceHandle};

use super::boundary::BlendBoundary;

#[allow(unused_imports)]
use popkern_core::geometry::CurveGeometry;

/// If the face has exactly `n` half-edges and `n ∈ {3, 4}` (too
/// few for the poke's normal `N ∈ [5, 16]` split path), split each
/// HE at the arc/bridge junction (from `boundary`) to reach `2n`
/// half-edges. Splitting at the junction — rather than a naive
/// `t = 0.5` — keeps the sub-face loop aligned with the cap's
/// bi-2 strip midpoint (strip[1] = arc/bridge junction point),
/// so the cap's `uv = (0.5, 1)` boundary midpoint matches the
/// sub-face's mid-split vertex position.
pub(crate) fn normalize_low_valence_face_for_poke(
    brep: &mut BRep,
    face: FaceHandle,
    expected_n: usize,
    boundary: &BlendBoundary,
) -> Result<(), CoreError> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::topology::EdgeHandle;
    let outer_loop = brep
        .faces
        .get(face)
        .ok_or_else(|| CoreError::InvalidHandle("face not in BRep".into()))?
        .outer_loop;
    let loop_start = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| CoreError::InvalidHandle("loop not in BRep".into()))?
        .half_edge;
    let mut count = 0usize;
    let mut he_edges: Vec<EdgeHandle> = Vec::new();
    let mut cursor = loop_start;
    loop {
        count += 1;
        he_edges.push(brep.half_edges[cursor].edge);
        cursor = brep.half_edges[cursor].next;
        if cursor == loop_start {
            break;
        }
        if count > 32 {
            return Err(CoreError::PreconditionFailed(
                "face loop exceeds 32 half-edges — likely unterminated".into(),
            ));
        }
    }
    if count != expected_n || !(3..=4).contains(&expected_n) {
        return Ok(());
    }
    let junctions: Vec<popkern_core::geometry::Point3> = (0..expected_n)
        .map(|k| {
            let arc = match &boundary.segments[2 * k] {
                super::boundary::BoundarySegment::FilletEnd { curve, .. } => curve,
                super::boundary::BoundarySegment::FaceEdge { curve, .. } => curve,
            };
            let (_, t_max) = arc.parameter_range();
            arc.point_at(t_max)
        })
        .collect();
    for edge in he_edges {
        let curve_on_edge = brep
            .edges
            .get(edge)
            .and_then(|e| e.curve.as_ref())
            .ok_or_else(|| CoreError::InvalidHandle("edge has no curve".into()))?
            .clone();
        let he_curve = &curve_on_edge.geometry;
        // The EDGE's stored window, not the geometry's intrinsic
        // parameter_range(): a Line's intrinsic range is (−inf, inf),
        // which disarmed the endpoint-exclusion guards below and let
        // a junction that projects exactly onto an edge endpoint
        // split the edge AT its end — the split's empty half then
        // panicked in `parameter_range` (PK-294, flagship boss
        // plate).
        let (t_min, t_max) = (
            curve_on_edge.t_start.min(curve_on_edge.t_end),
            curve_on_edge.t_start.max(curve_on_edge.t_end),
        );
        let mut best: Option<(f64, f64)> = None;
        for jp in &junctions {
            if let Some(t) = he_curve.closest_point(jp) {
                let band = ((t_max - t_min) * 1e-6).max(1e-9);
                if !t.is_finite() || t <= t_min + band || t >= t_max - band {
                    continue;
                }
                let d = he_curve.point_at(t).distance_to(jp);
                if best.map(|(_, bd)| d < bd).unwrap_or(true) {
                    best = Some((t, d));
                }
            }
        }
        let t_split = best.map(|(t, _)| t).unwrap_or(0.5 * (t_min + t_max));
        brep.split_edge(edge, t_split)?;
    }
    Ok(())
}

/// Walk `face`'s outer loop from its stored half-edge (the same
/// enumeration the poke uses) and pair each half-edge with the
/// boundary segment that LIES ON the surface across it.
///
/// Pairing is semantic, not proximity: the measured separation
/// between a segment and its own support surface is ~1e-5, against
/// ~1e-2 for every other surface at the corner — an earlier
/// endpoint-proximity matcher mis-paired arcs across strips and made
/// the rim worse. Returns `None` (nothing consumed) unless every
/// half-edge pairs with a distinct segment and consecutive loop
/// edges carry ring-adjacent segments.
pub(crate) fn pair_loop_half_edges_to_segments(
    brep: &BRep,
    face: FaceHandle,
    boundary: &BlendBoundary,
) -> Option<(Vec<popkern_core::topology::HalfEdgeHandle>, Vec<usize>)> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;
    /// True pairs measure ~1e-5 (station drift included); false
    /// pairs ~1e-2.
    const LIES_ON: f64 = 1.0e-3;
    let outer_loop = brep.faces.get(face)?.outer_loop;
    let start = brep.loops.get(outer_loop)?.half_edge;
    let mut loop_hes = Vec::new();
    let mut cursor = start;
    for _ in 0..32 {
        loop_hes.push(cursor);
        cursor = brep.half_edges.get(cursor)?.next;
        if cursor == start {
            break;
        }
    }
    let n_seg = boundary.segments.len();
    if loop_hes.len() != n_seg {
        return None;
    }
    let mut assignment = Vec::with_capacity(n_seg);
    let mut used = vec![false; n_seg];
    for &he in &loop_hes {
        let across_surface = brep
            .half_edges
            .get(he)
            .and_then(|hd| brep.half_edges.get(hd.twin))
            .and_then(|t| brep.loops.get(t.loop_))
            .and_then(|l| brep.faces.get(l.face))
            .and_then(|f| f.surface.as_ref())?;
        let mut best: Option<(usize, f64)> = None;
        for (k, seg_used) in used.iter().enumerate() {
            if *seg_used {
                continue;
            }
            let curve = pns_segment_curve(boundary, k);
            let (t0, t1) = curve.parameter_range();
            let mut worst = 0.0_f64;
            for s in 0..=8 {
                let p = curve.point_at(t0 + (t1 - t0) * f64::from(s) / 8.0);
                let Some((u, v)) = across_surface.geometry.closest_point(&p) else {
                    worst = f64::MAX;
                    break;
                };
                worst = worst.max(across_surface.geometry.point_at(u, v).distance_to(&p));
            }
            if best.map(|(_, bd)| worst < bd).unwrap_or(true) {
                best = Some((k, worst));
            }
        }
        let (k, d) = best?;
        if d > LIES_ON {
            return None;
        }
        assignment.push(k);
        used[k] = true;
    }
    for i in 0..assignment.len() {
        let a = assignment[i] as isize;
        let b = assignment[(i + 1) % assignment.len()] as isize;
        let diff = (a - b).rem_euclid(n_seg as isize);
        if diff != 1 && diff != n_seg as isize - 1 {
            return None;
        }
    }
    Some((loop_hes, assignment))
}

/// The segment's curve regardless of variant.
pub(crate) fn pns_segment_curve(
    boundary: &BlendBoundary,
    k: usize,
) -> &popkern_core::geometry::nurbs_curve::NurbsCurve {
    match &boundary.segments[k] {
        super::boundary::BoundarySegment::FilletEnd { curve, .. } => curve,
        super::boundary::BoundarySegment::FaceEdge { curve, .. } => curve,
    }
}

/// Reverse a NURBS curve's direction (same image, opposite
/// parameterization). Mirror of `popkern_core`'s private
/// `reverse_nurbs_curve`.
pub(crate) fn reverse_pns_segment_curve(
    nurbs: &popkern_core::geometry::nurbs_curve::NurbsCurve,
) -> popkern_core::geometry::nurbs_curve::NurbsCurve {
    use popkern_core::geometry::curve::CurveEval;
    let (t_min, t_max) = nurbs.parameter_range();
    let mut knots = nurbs
        .knots
        .iter()
        .rev()
        .map(|&k| t_min + t_max - k)
        .collect::<Vec<_>>();
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut control_points = nurbs.control_points.clone();
    control_points.reverse();
    let mut weights = nurbs.weights.clone();
    weights.reverse();
    popkern_core::geometry::nurbs_curve::NurbsCurve {
        degree: nurbs.degree,
        knots,
        control_points,
        weights,
    }
}

/// Rebuild the boundary in the committed loop's own traversal order:
/// segment `i` of the result is the segment across loop half-edge
/// `i`, its curve oriented along the loop direction, starting from
/// the loop's stored half-edge (which the poke also starts from).
///
/// Sector parity is part of the contract: the poke pairs loop
/// half-edges `(2k, 2k+1)` into sector `k` and the cap adapter pairs
/// segments the same way, so segment 0 must be a `FilletEnd`. When
/// the loop's stored half-edge sits on a bridge instead, the loop
/// entry is rotated forward by one half-edge — an entry-pointer
/// update, not a topology change.
pub(crate) fn align_boundary_to_loop(
    brep: &mut BRep,
    face: FaceHandle,
    boundary: &BlendBoundary,
) -> Option<BlendBoundary> {
    let (mut loop_hes, mut assignment) = pair_loop_half_edges_to_segments(brep, face, boundary)?;
    // Rotate the loop entry so it starts on a FilletEnd (even index).
    if assignment.first().is_some_and(|&k| k % 2 == 1) {
        let outer_loop = brep.faces.get(face)?.outer_loop;
        let next = brep.half_edges.get(*loop_hes.first()?)?.next;
        brep.loops[outer_loop].half_edge = next;
        loop_hes.rotate_left(1);
        assignment.rotate_left(1);
    }
    if assignment.first().is_some_and(|&k| k % 2 == 1) {
        return None;
    }
    let n_seg = boundary.segments.len();
    let reversed = {
        let a = assignment[0] as isize;
        let b = assignment[1] as isize;
        (a - b).rem_euclid(n_seg as isize) == 1
    };
    let segments = assignment
        .iter()
        .map(|&k| {
            let curve = if reversed {
                reverse_pns_segment_curve(pns_segment_curve(boundary, k))
            } else {
                pns_segment_curve(boundary, k).clone()
            };
            match &boundary.segments[k] {
                super::boundary::BoundarySegment::FilletEnd { fillet_index, .. } => {
                    super::boundary::BoundarySegment::FilletEnd {
                        fillet_index: *fillet_index,
                        curve,
                    }
                }
                super::boundary::BoundarySegment::FaceEdge { face, .. } => {
                    super::boundary::BoundarySegment::FaceEdge { face: *face, curve }
                }
            }
        })
        .collect();
    let mut aligned = boundary.clone();
    aligned.segments = segments;
    Some(aligned)
}
/// Materialize the PnS rim: rewrite the blend face's outer-loop edge
/// curves and vertex positions from the loop-ALIGNED boundary ring.
///
/// At the low valences the committed corner loop is a mid-pipeline
/// scaffold: every edge carries a stale `Line` and consecutive edge
/// curves do not even meet (measured 2026-07-26 on the n=3 mixed
/// pyramid: 4.632e-2 curve-off-support, vertex positions that chain
/// nowhere). The incumbent direct-apex route rebuilds this rim inside
/// its own materialization step, which is why its finished bodies
/// measure clean; PnS previously authored only cap surfaces on top of
/// the scaffold, and the force-accept probe showed the resulting body
/// keeps the full defect. The boundary is exact — each segment lies
/// on its support surface at ~1e-5 or better — so it is the authority
/// for the rim.
///
/// Takes the output of [`align_boundary_to_loop`]: segment `i`
/// corresponds to loop half-edge `i` and is oriented along the loop,
/// so junction `i` (start of segment `i`) is the position of half-
/// edge `i`'s origin vertex. All-or-nothing: if the aligned segments
/// do not chain end-to-start the rewrite is abandoned before touching
/// anything, and the acceptance gate measures the scaffold as before.
pub(crate) fn materialize_pns_rim_from_boundary(
    brep: &mut BRep,
    face: FaceHandle,
    boundary: &BlendBoundary,
) {
    use popkern_core::geometry::CurveGeometry;
    use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};

    /// Aligned segments must chain within this slack; the ring closes
    /// to ~1e-12, so anything looser means the alignment is wrong.
    const CHAIN_SLACK: f64 = 1.0e-6;

    let Some(outer_loop) = brep.faces.get(face).map(|f| f.outer_loop) else {
        return;
    };
    let Some(start) = brep.loops.get(outer_loop).map(|l| l.half_edge) else {
        return;
    };
    let mut loop_hes = Vec::new();
    let mut cursor = start;
    for _ in 0..32 {
        loop_hes.push(cursor);
        let Some(hd) = brep.half_edges.get(cursor) else {
            return;
        };
        cursor = hd.next;
        if cursor == start {
            break;
        }
    }
    let n_seg = boundary.segments.len();
    if loop_hes.len() != n_seg {
        return;
    }
    // Verify the aligned ring chains before mutating anything.
    for i in 0..n_seg {
        let cur = pns_segment_curve(boundary, i);
        let next = pns_segment_curve(boundary, (i + 1) % n_seg);
        let (_, ct1) = cur.parameter_range();
        let (nt0, _) = next.parameter_range();
        if cur.point_at(ct1).distance_to(&next.point_at(nt0)) > CHAIN_SLACK {
            return;
        }
    }
    for (i, &he) in loop_hes.iter().enumerate() {
        let curve = pns_segment_curve(boundary, i).clone();
        let (t0, t1) = curve.parameter_range();
        let junction = curve.point_at(t0);
        let Some(hd) = brep.half_edges.get(he) else {
            continue;
        };
        let (edge, origin) = (hd.edge, hd.origin);
        // CurveOnEdge's t_start belongs to the edge's CANONICAL
        // half-edge origin, which may be the twin of `he`.
        let canonical_is_he = brep.edges.get(edge).map(|e| e.half_edge) == Some(he);
        let (t_start, t_end) = if canonical_is_he { (t0, t1) } else { (t1, t0) };
        if let Some(edge_data) = brep.edges.get_mut(edge) {
            edge_data.curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Nurbs(curve),
                t_start,
                t_end,
            });
        }
        brep.vertices[origin].position = junction;
    }
    // The moved vertices are shared with neighbouring faces whose
    // edge curves still end at the old positions — pull those curve
    // endpoints to the new vertices so every incident loop stays
    // closed. The rim curves just written end exactly at the moved
    // vertices, so this pass leaves them untouched. The threshold is
    // a gap-detection floor (same floor the pipeline's own coherence
    // pass uses), not a geometric tolerance.
    let _ = brep.repair_edge_curve_endpoint_coherence(1.0e-8);
}
