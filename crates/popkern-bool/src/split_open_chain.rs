//! Open-chain pre-pass for boolean face splitting: chains OPEN
//! intersection fragments that dead-end at shared interior points into
//! single boundary-to-boundary polylines and splits the face with each
//! combined trim. Extracted from `split.rs` (PK-304).

use std::collections::HashSet;

use popkern_core::{
    tolerance::Tolerance,
    topology::{BRep, FaceHandle},
};

use crate::error::BoolDiagnostics;
use crate::intersect::IntersectionRecord;

use super::FragmentMap;
use super::helpers::{find_fragment_for_curve, split_one_face};

use popkern_core::geometry::curve::CurveEval as _;

/// Reverse a NURBS curve (control points, weights, mirrored knots).
pub(super) fn reversed_nurbs(
    c: &popkern_core::geometry::NurbsCurve,
) -> popkern_core::geometry::NurbsCurve {
    let (lo, hi) = c.parameter_range();
    let mut knots: Vec<f64> = c.knots.iter().map(|&k| lo + hi - k).collect();
    knots.reverse();
    popkern_core::geometry::NurbsCurve {
        degree: c.degree,
        knots,
        control_points: c.control_points.iter().rev().copied().collect(),
        weights: c.weights.iter().rev().copied().collect(),
    }
}

/// Chain OPEN intersection fragments that dead-end at shared interior
/// points into single boundary-to-boundary polylines, and split the
/// face with each combined curve. Scoped to degree-1 (line) segments
/// — the bounded-cutter crease case; curved dead-ends stay on the
/// per-curve path. Returns the record indices consumed here.
#[allow(clippy::too_many_arguments)]
pub(super) fn open_chain_pre_pass(
    brep: &mut BRep,
    original_face: FaceHandle,
    entries: &[(usize, bool)],
    records: &[IntersectionRecord],
    consumed: &HashSet<usize>,
    tolerance: &Tolerance,
    fragments: &mut FragmentMap,
) -> Result<HashSet<usize>, BoolDiagnostics> {
    let mut chained: HashSet<usize> = HashSet::new();
    // Curved records join the chain pool since PK-304 — but ONLY on
    // cell-resolvable faces (doubly-periodic, non-winding ring): a
    // marched torus arc can be delivered as two fragments dead-ending
    // at the seed junction mid-face (measured: face 16v1's two pieces
    // each failed the vertex match, the face stayed unsplit, and
    // classification flood-fill leaked through the unpartitioned
    // disc). On winding fillet-band tori the same relaxation chained
    // genuinely distinct curved trims into a corner-rounding polyline
    // and regressed the R8 fillet boolean parity fixture (measured
    // 2026-08-21), so those keep the historic degree-1-only pool.
    // Curved members are polyline-sampled below, not endpoint-joined.
    let curved_ok =
        crate::clip::chart_cell::face_chart_cell(brep, original_face, tolerance).is_some();
    let candidates: Vec<(usize, bool)> = entries
        .iter()
        .copied()
        .filter(|(ri, _)| {
            !consumed.contains(ri) && (curved_ok || records[*ri].curve.curve_3d.degree == 1)
        })
        .collect();
    if candidates.len() < 2 {
        return Ok(chained);
    }
    let ends: Vec<(
        popkern_core::geometry::Point3,
        popkern_core::geometry::Point3,
    )> = candidates
        .iter()
        .map(|&(ri, _)| {
            let c = &records[ri].curve.curve_3d;
            let (t0, t1) = c.parameter_range();
            (c.point_at(t0), c.point_at(t1))
        })
        .collect();
    let gate = tolerance.linear() * 100.0;
    let near = |a: &popkern_core::geometry::Point3, b: &popkern_core::geometry::Point3| {
        a.distance_to(b) < gate
    };
    // A junction is "interior" when it sits on no boundary VERTEX
    // and no boundary EDGE of the face. Mid-edge chain ENDS are
    // boundary-reaching — the split stage's T-junction pre-split
    // mints their vertices (measured: a shaped cutter's face took a
    // three-segment ∏ of trims whose outer ends landed mid-edge on
    // its top edge; the vertex-only test called them interior and
    // the chain was rejected, failing the whole boolean). Snapshot
    // the boundary up front — the face is mutated below.
    let boundary_segments: Vec<(
        popkern_core::geometry::Point3,
        popkern_core::geometry::Point3,
    )> = {
        let Some(face) = brep.faces.get(original_face) else {
            return Ok(chained);
        };
        let mut loops = vec![face.outer_loop];
        loops.extend(face.inner_loops.iter().copied());
        let mut segs = Vec::new();
        for lh in loops {
            let Some(loop_data) = brep.loops.get(lh) else {
                continue;
            };
            let start = loop_data.half_edge;
            let mut he = start;
            loop {
                let a = brep.vertices[brep.half_edges[he].origin].position;
                let next = brep.half_edges[he].next;
                let b = brep.vertices[brep.half_edges[next].origin].position;
                // CURVED edges contribute sampled sub-segments, not
                // the vertex chord: a torus quadrant rim's chord sits
                // a sagitta (~0.1) off the arc, so real
                // boundary-reaching chain ends read interior and the
                // chain was rejected (PK-304 — face 16v1's two-piece
                // arc never partitioned and classification flood-fill
                // leaked through the unsplit disc).
                let edge = &brep.edges[brep.half_edges[he].edge];
                if popkern_core::trace::trace_enabled("face-split") {
                    eprintln!(
                        "    [chain-bnd] {original_face:?} edge {:?} curve={} a=({:.3},{:.3},{:.3}) b=({:.3},{:.3},{:.3})",
                        brep.half_edges[he].edge,
                        edge.curve.is_some(),
                        a.x,
                        a.y,
                        a.z,
                        b.x,
                        b.y,
                        b.z
                    );
                }
                match edge.curve.as_ref() {
                    Some(c)
                        if !matches!(
                            c.geometry,
                            popkern_core::geometry::curve::CurveGeometry::Line { .. }
                        ) =>
                    {
                        use popkern_core::geometry::curve::CurveEval as _;
                        // 64 spans: a quarter-arc rim's 16-span
                        // sagitta (~1.2e-3) still exceeded the 1e-4
                        // on-boundary gate and boundary heads read
                        // interior; 64 spans sit at ~7e-5.
                        let n = 64;
                        let mut prev = a;
                        for i in 1..=n {
                            let t = c.t_start + (c.t_end - c.t_start) * (i as f64) / (n as f64);
                            let q = c.geometry.point_at(t);
                            // The half-edge may run against the curve;
                            // endpoints still bound the span, so the
                            // sampled polyline covers the same arc
                            // either way.
                            segs.push((prev, q));
                            prev = q;
                        }
                        segs.push((prev, b));
                    }
                    _ => segs.push((a, b)),
                }
                he = next;
                if he == start {
                    break;
                }
            }
        }
        segs
    };
    let interior_point = |p: &popkern_core::geometry::Point3| -> bool {
        !boundary_segments.iter().any(|(a, b)| {
            if near(a, p) || near(b, p) {
                return true;
            }
            let ab = *b - *a;
            let len_sq = ab.length_squared();
            if len_sq < 1e-24 {
                return false;
            }
            let t = ((*p - *a).dot(&ab) / len_sq).clamp(0.0, 1.0);
            let closest = *a + ab * t;
            closest.distance_to(p) < gate
        })
    };

    let mut used = vec![false; candidates.len()];
    'seeds: for seed in 0..candidates.len() {
        if used[seed] {
            continue;
        }
        // Try both orientations of the seed: the junction may sit at
        // either of its ends.
        for seed_rev in [false, true] {
            let mut chain: Vec<(usize, bool)> = vec![(seed, seed_rev)];
            while let Some(&(tail_idx, tail_rev)) = chain.last() {
                let tip = if tail_rev {
                    ends[tail_idx].0
                } else {
                    ends[tail_idx].1
                };
                if !interior_point(&tip) {
                    break;
                }
                let mut extended = false;
                for k in 0..candidates.len() {
                    if used[k] || k == seed || chain.iter().any(|&(i, _)| i == k) {
                        continue;
                    }
                    if near(&ends[k].0, &tip) {
                        chain.push((k, false));
                        extended = true;
                        break;
                    }
                    if near(&ends[k].1, &tip) {
                        chain.push((k, true));
                        extended = true;
                        break;
                    }
                }
                if !extended {
                    break;
                }
            }
            if chain.len() < 2 {
                continue;
            }
            // Both free ends must be non-interior (reach the boundary).
            let (head_idx, head_rev) = chain[0];
            let head = if head_rev {
                ends[head_idx].1
            } else {
                ends[head_idx].0
            };
            let Some(&(tail_idx, tail_rev)) = chain.last() else {
                continue; // unreachable: the chain is seeded non-empty
            };
            let tail = if tail_rev {
                ends[tail_idx].0
            } else {
                ends[tail_idx].1
            };
            if popkern_core::trace::trace_enabled("face-split") {
                eprintln!(
                    "  [open-chain] {original_face:?} len {} head ({:.3},{:.3},{:.3}) int={} tail ({:.3},{:.3},{:.3}) int={}",
                    chain.len(),
                    head.x,
                    head.y,
                    head.z,
                    interior_point(&head),
                    tail.x,
                    tail.y,
                    tail.z,
                    interior_point(&tail)
                );
            }
            if interior_point(&head) || interior_point(&tail) {
                continue;
            }
            // Concatenate the chain into degree-1 multi-span curves in 3D
            // and both UV charts.
            let mut pts3: Vec<popkern_core::geometry::Point3> = Vec::new();
            // THIS face's UV chart: each record contributes its own
            // side (mixed sides are legal — self-splits pass per-entry
            // flags; two-body splits pass a constant flag).
            let mut pts_uv_face: Vec<popkern_core::geometry::Point3> = Vec::new();
            for (pos, &(k, rev)) in chain.iter().enumerate() {
                let (ri, use_uv1) = candidates[k];
                let c = &records[ri].curve;
                let (t0, t1) = c.curve_3d.parameter_range();
                let (ta, tb) = if rev { (t1, t0) } else { (t0, t1) };
                // Straight members contribute their endpoints (exact);
                // curved members are sampled densely so the joined
                // polyline stays within the boolean-era accuracy class
                // (24 spans over a sub-unit arc ≈ 1e-4 sagitta).
                let spans = if c.curve_3d.degree == 1 { 1 } else { 24 };
                let start = if pos == 0 { 0 } else { 1 };
                let uv_curve = if use_uv1 { &c.curve_uv1 } else { &c.curve_uv2 };
                for i in start..=spans {
                    let t = ta + (tb - ta) * (i as f64) / (spans as f64);
                    pts3.push(c.curve_3d.point_at(t));
                    pts_uv_face.push(uv_curve.point_at(t));
                }
            }
            let polyline = |pts: &[popkern_core::geometry::Point3]| {
                let n = pts.len();
                let mut knots = vec![0.0, 0.0];
                for i in 1..n - 1 {
                    knots.push(i as f64 / (n - 1) as f64);
                }
                knots.push(1.0);
                knots.push(1.0);
                popkern_core::geometry::NurbsCurve {
                    degree: 1,
                    knots,
                    control_points: pts.to_vec(),
                    weights: vec![1.0; n],
                }
            };
            // Both slots carry this face's chart: split_one_face only
            // reads the side we hand it, and the first-control-point
            // disambiguation in the presplit helper then resolves to
            // the same chart either way.
            let uv_face = polyline(&pts_uv_face);
            let combined = popkern_core::geometry::curve_on_surface::CurveOnSurface {
                curve_3d: polyline(&pts3),
                curve_uv1: uv_face.clone(),
                curve_uv2: uv_face,
                is_closed: false,
            };
            let curve_uv = &combined.curve_uv1;
            let target_face = match find_fragment_for_curve(
                brep,
                original_face,
                curve_uv,
                fragments,
                tolerance,
            ) {
                Some(frag) => frag,
                None => continue,
            };
            if !brep.faces.contains_key(target_face) {
                continue;
            }
            split_one_face(
                brep,
                target_face,
                &combined.curve_3d,
                curve_uv,
                &combined,
                tolerance,
                records,
                fragments,
            )?;
            for &(k, _) in &chain {
                used[k] = true;
                chained.insert(candidates[k].0);
            }
            if popkern_core::trace::trace_enabled("boolean") {
                eprintln!(
                    "  [open-chain] face {original_face:?}: chained {} segments into one trim",
                    chain.len()
                );
            }
            continue 'seeds;
        }
    }
    Ok(chained)
}
