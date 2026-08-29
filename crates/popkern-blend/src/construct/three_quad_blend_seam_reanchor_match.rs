//! Structural matching of live seam edges to blend-boundary
//! segments by cyclic ring alignment.
//!
//! The mitred ring can sit farther from the boundary than a whole
//! segment chord (hexagon6_variable: 1.5e-2 offsets on 1.4e-2
//! chords, exact per-face distance ties, and no reliable twin
//! surface parity), so neither identity nor orientation may rely on
//! local proximity. Both rings are cycles, and cycles align
//! globally:
//!
//! * the **live cycle** is pure topology — consecutive seam edges
//!   share a ring vertex;
//! * the **segment cycle** is machine-exact — adjacent boundary
//!   segments share junction endpoints (measured ~1e-16 on
//!   hexagon6_variable);
//! * the correspondence is one of `2 × 2n` rotation/flip
//!   hypotheses, chosen by the global sum of live-vertex-to-junction
//!   distances. Local ties cannot flip a global alignment.

use super::*;

pub(super) struct SeamEdge {
    pub half_edge: HalfEdgeHandle,
    pub edge: EdgeHandle,
    pub v_origin: VertexHandle,
    pub v_dest: VertexHandle,
}

/// Junction-coherence outcome alongside the jobs: the worst
/// distance between the two segment endpoints merged into one
/// junction (adjacent segments share junctions by construction, so
/// this measures the boundary's own internal coherence).
pub(super) struct MatchedRing {
    pub jobs: Vec<ReanchorJob>,
    pub worst_junction_spread: f64,
}

/// One junction of the segment cycle: cycle position `i` joins
/// segment `seg_before` (consuming its end `end_before`) to
/// `seg_after` (consuming `end_after`), at `point`.
struct SegmentJunction {
    seg_before: usize,
    end_before: usize,
    seg_after: usize,
    end_after: usize,
    point: Point3,
}

pub(super) fn match_ring_to_segments(
    brep: &BRep,
    boundary: &BlendBoundary,
    sub_faces: &[FaceHandle],
) -> Result<MatchedRing, &'static str> {
    let endpoints: Vec<[Point3; 2]> = boundary
        .segments
        .iter()
        .map(|segment| {
            let curve = segment_curve(segment);
            let (t0, t1) = curve.parameter_range();
            [curve.point_at(t0), curve.point_at(t1)]
        })
        .collect();
    let m = endpoints.len();
    let (seg_order, junctions, worst_junction_spread) = segment_cycle(&endpoints)?;
    let (live_edges, walk_forward) = live_seam_cycle(brep, sub_faces)?;
    if live_edges.len() != m {
        return Err("ring_edge_count_mismatch");
    }

    // Live junction positions: live_junction[i] is the vertex shared
    // by live edge i and live edge i+1 (cyclically).
    let mut live_junctions = Vec::with_capacity(m);
    for i in 0..m {
        let edge = &live_edges[i];
        let v = if walk_forward[i] {
            edge.v_dest
        } else {
            edge.v_origin
        };
        let Some(vd) = brep.vertices.get(v) else {
            return Err("ring_vertex_missing");
        };
        live_junctions.push(vd.position);
    }

    // Global alignment: map live edge i to cycle position
    // (r + i) mod m (flip = false) or (r - i) mod m (flip = true);
    // the junction AFTER live edge i then sits at cycle position
    // p(i) (flip = false) or p(i) - 1 (flip = true). Cost is the sum
    // of live-junction-to-segment-junction distances.
    let mut best: Option<(usize, bool, f64)> = None;
    for flip in [false, true] {
        for r in 0..m {
            let mut cost = 0.0_f64;
            for (i, live_pos) in live_junctions.iter().enumerate() {
                let j = junction_position_after(i, r, flip, m);
                cost += live_pos.distance_to(&junctions[j].point);
            }
            if best.is_none_or(|(_, _, bc)| cost < bc) {
                best = Some((r, flip, cost));
            }
        }
    }
    let Some((r, flip, cost)) = best else {
        return Err("ring_alignment_failed");
    };
    // Sanity bound — average junction displacement past a quarter of
    // the boundary perimeter means the live ring is not this
    // boundary's ring; anchoring to it would be a guess.
    let perimeter: f64 = endpoints
        .iter()
        .map(|eps| eps[0].distance_to(&eps[1]))
        .sum();
    if cost / m as f64 > 0.25 * perimeter {
        return Err("match_distance_implausible");
    }
    if crate::direct_apex_trace::cap_g0_enabled() {
        eprintln!(
            "[direct-apex-trace cap-g0] reanchor ring alignment: r={r} flip={flip} \
             mean_gap={:.3e} junction_spread={worst_junction_spread:.3e}",
            cost / m as f64
        );
    }

    let mut jobs = Vec::with_capacity(m);
    for (i, edge) in live_edges.iter().enumerate() {
        let p = cycle_position(i, r, flip, m);
        let segment_index = seg_order[p];
        let j_after = &junctions[junction_position_after(i, r, flip, m)];
        let j_before = &junctions[junction_position_after((i + m - 1) % m, r, flip, m)];
        // Walk direction runs before-junction -> after-junction; map
        // back to the half-edge's own origin/dest.
        let (origin_junction, dest_junction) = if walk_forward[i] {
            (j_before, j_after)
        } else {
            (j_after, j_before)
        };
        let end_at = |junction: &SegmentJunction| -> Result<usize, &'static str> {
            if junction.seg_before == segment_index {
                Ok(junction.end_before)
            } else if junction.seg_after == segment_index {
                Ok(junction.end_after)
            } else {
                Err("ring_order_mismatch")
            }
        };
        let end_origin = end_at(origin_junction)?;
        let end_dest = end_at(dest_junction)?;
        if end_origin == end_dest {
            return Err("segment_end_reused");
        }
        jobs.push(ReanchorJob {
            half_edge: edge.half_edge,
            edge: edge.edge,
            segment_index,
            reversed: end_origin == 1,
            v_origin: edge.v_origin,
            v_dest: edge.v_dest,
            target_origin: origin_junction.point,
            target_dest: dest_junction.point,
        });
    }
    Ok(MatchedRing {
        jobs,
        worst_junction_spread,
    })
}

fn cycle_position(i: usize, r: usize, flip: bool, m: usize) -> usize {
    if flip {
        (r + m - (i % m)) % m
    } else {
        (r + i) % m
    }
}

/// Cycle position of the junction that follows live edge `i` in walk
/// order. Junction `p` sits between cycle positions `p` and `p + 1`,
/// so under a flip the "after" junction is the one BEFORE the mapped
/// position.
fn junction_position_after(i: usize, r: usize, flip: bool, m: usize) -> usize {
    let p = cycle_position(i, r, flip, m);
    if flip { (p + m - 1) % m } else { p }
}

/// Order the boundary segments into their closed cycle by pairing
/// segment ends globally (nearest cross-segment pairs first, each
/// end consumed once — exact on a coherent boundary). Returns the
/// cycle order, the junction between consecutive cycle positions
/// (junction `i` joins `order[i]` to `order[i + 1]`), and the worst
/// end-pair distance.
#[allow(clippy::type_complexity)]
fn segment_cycle(
    endpoints: &[[Point3; 2]],
) -> Result<(Vec<usize>, Vec<SegmentJunction>, f64), &'static str> {
    let m = endpoints.len();
    let mut pair_candidates: Vec<(usize, usize, usize, usize, f64)> = Vec::new();
    for sa in 0..m {
        for sb in (sa + 1)..m {
            for ia in 0..2 {
                for ib in 0..2 {
                    let d = endpoints[sa][ia].distance_to(&endpoints[sb][ib]);
                    pair_candidates.push((sa, ia, sb, ib, d));
                }
            }
        }
    }
    pair_candidates.sort_by(|a, b| a.4.partial_cmp(&b.4).unwrap_or(std::cmp::Ordering::Equal));
    type EndPartner = Option<(usize, usize, Point3, f64)>;
    let mut end_partner: Vec<[EndPartner; 2]> = vec![[None, None]; m];
    let mut paired = 0usize;
    for (sa, ia, sb, ib, d) in pair_candidates {
        if paired == m {
            break;
        }
        if end_partner[sa][ia].is_some() || end_partner[sb][ib].is_some() {
            continue;
        }
        let junction = average_point(&[endpoints[sa][ia], endpoints[sb][ib]]);
        end_partner[sa][ia] = Some((sb, ib, junction, d));
        end_partner[sb][ib] = Some((sa, ia, junction, d));
        paired += 1;
    }
    if paired != m {
        return Err("junction_pairing_incomplete");
    }
    let worst_spread = end_partner
        .iter()
        .flatten()
        .flatten()
        .map(|(_, _, _, d)| *d)
        .fold(0.0_f64, f64::max);

    // Walk the cycle: enter each segment through one end, exit the
    // other; the exit end's partner is the next segment.
    let mut seg_order = Vec::with_capacity(m);
    let mut junctions = Vec::with_capacity(m);
    let mut seg = 0usize;
    let mut entry_end = 0usize;
    let mut visited = vec![false; m];
    for _ in 0..m {
        if visited[seg] {
            return Err("segment_cycle_open");
        }
        visited[seg] = true;
        seg_order.push(seg);
        let exit_end = 1 - entry_end;
        let Some((next_seg, next_end, junction, _)) = end_partner[seg][exit_end] else {
            return Err("junction_pairing_incomplete");
        };
        junctions.push(SegmentJunction {
            seg_before: seg,
            end_before: exit_end,
            seg_after: next_seg,
            end_after: next_end,
            point: junction,
        });
        seg = next_seg;
        entry_end = next_end;
    }
    if seg != 0 || entry_end != 0 {
        return Err("segment_cycle_open");
    }
    if crate::direct_apex_trace::cap_g0_enabled() {
        eprintln!(
            "[direct-apex-trace cap-g0] reanchor segment cycle: order={seg_order:?} \
             spread={worst_spread:.3e}"
        );
    }
    Ok((seg_order, junctions, worst_spread))
}

/// Walk the live seam ring: collect every cap outer-loop half-edge
/// whose twin is on a support face, then chain them through their
/// shared vertices into one closed cycle. `walk_forward[i]` records
/// whether edge `i`'s origin->dest direction agrees with the walk.
fn live_seam_cycle(
    brep: &BRep,
    sub_faces: &[FaceHandle],
) -> Result<(Vec<SeamEdge>, Vec<bool>), &'static str> {
    let cap_faces: std::collections::BTreeSet<FaceHandle> = sub_faces.iter().copied().collect();
    let mut pool: Vec<SeamEdge> = Vec::new();
    for &sub_face in sub_faces {
        let Some(loop_data) = brep
            .faces
            .get(sub_face)
            .and_then(|f| brep.loops.get(f.outer_loop))
        else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut safety = 64usize;
        loop {
            if safety == 0 {
                break;
            }
            safety -= 1;
            let Some(he) = brep.half_edges.get(cursor) else {
                break;
            };
            let twin_face = brep
                .half_edges
                .get(he.twin)
                .and_then(|twin| brep.loops.get(twin.loop_))
                .map(|l| l.face);
            if twin_face.is_some_and(|f| !cap_faces.contains(&f))
                && let Some(next) = brep.half_edges.get(he.next)
            {
                pool.push(SeamEdge {
                    half_edge: cursor,
                    edge: he.edge,
                    v_origin: he.origin,
                    v_dest: next.origin,
                });
            }
            cursor = he.next;
            if cursor == start {
                break;
            }
        }
    }
    if pool.is_empty() {
        return Err("ring_edge_count_mismatch");
    }
    // Chain: at each ring vertex exactly two pool edges meet.
    let mut used = vec![false; pool.len()];
    let mut order = Vec::with_capacity(pool.len());
    let mut walk_forward = Vec::with_capacity(pool.len());
    let mut current = 0usize;
    let mut forward = true;
    for step in 0..pool.len() {
        used[current] = true;
        order.push(current);
        walk_forward.push(forward);
        let far_vertex = if forward {
            pool[current].v_dest
        } else {
            pool[current].v_origin
        };
        let mut next: Option<(usize, bool)> = None;
        for (idx, candidate) in pool.iter().enumerate() {
            if used[idx] {
                continue;
            }
            let dir = if candidate.v_origin == far_vertex {
                true
            } else if candidate.v_dest == far_vertex {
                false
            } else {
                continue;
            };
            if next.is_some() {
                return Err("ring_vertex_valence");
            }
            next = Some((idx, dir));
        }
        match next {
            Some((idx, dir)) => {
                current = idx;
                forward = dir;
            }
            None => {
                // Only the LAST edge may find no unused successor,
                // and it must close back to the start vertex.
                let start_vertex = if walk_forward[0] {
                    pool[order[0]].v_origin
                } else {
                    pool[order[0]].v_dest
                };
                if step + 1 != pool.len() || far_vertex != start_vertex {
                    return Err("ring_not_closed");
                }
            }
        }
    }
    let ordered: Vec<SeamEdge> = order
        .into_iter()
        .map(|idx| SeamEdge {
            half_edge: pool[idx].half_edge,
            edge: pool[idx].edge,
            v_origin: pool[idx].v_origin,
            v_dest: pool[idx].v_dest,
        })
        .collect();
    Ok((ordered, walk_forward))
}
