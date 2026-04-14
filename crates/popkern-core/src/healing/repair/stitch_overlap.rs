//! Partial-overlap detection between free boundary edges.
//!
//! Split out of `stitch.rs` to keep the parent module under the
//! 500-line structural cap. Operates on the same sampled-point
//! representation used by [`super::find_matching_edges`].

use crate::topology::{BRep, EdgeHandle};

use super::sample_edge_points;

/// Detect partial overlap between two free boundary edges.
///
/// Returns the parameter ranges on each edge where they overlap,
/// or `None` if there is no partial overlap (full match or no match).
pub fn detect_partial_overlap(
    brep: &BRep,
    edge_a: EdgeHandle,
    edge_b: EdgeHandle,
    max_gap: f64,
) -> Option<((f64, f64), (f64, f64))> {
    let pts_a = sample_edge_points(brep, edge_a, 16)?;
    let pts_b = sample_edge_points(brep, edge_b, 16)?;
    if pts_a.len() < 2 || pts_b.len() < 2 {
        return None;
    }
    // For each point on edge_a, check closest distance to edge_b.
    let in_range: Vec<bool> = pts_a
        .iter()
        .map(|pa| {
            pts_b
                .iter()
                .map(|pb| pa.distance_to(pb))
                .fold(f64::MAX, f64::min)
                <= max_gap
        })
        .collect();
    // Find the longest contiguous run of matches.
    let (mut best_start, mut best_len) = (0, 0_usize);
    let (mut cur_start, mut cur_len) = (0, 0_usize);
    for (i, &m) in in_range.iter().enumerate() {
        if m {
            if cur_len == 0 {
                cur_start = i;
            }
            cur_len += 1;
        } else {
            if cur_len > best_len {
                best_start = cur_start;
                best_len = cur_len;
            }
            cur_len = 0;
        }
    }
    if cur_len > best_len {
        best_start = cur_start;
        best_len = cur_len;
    }
    if best_len == 0 || best_len == pts_a.len() {
        return None;
    }
    // Map sample indices to parameter ranges.
    let na = (pts_a.len() - 1) as f64;
    let ca = brep.edges.get(edge_a)?.curve.as_ref()?;
    let span_a = ca.t_end - ca.t_start;
    let t_start_a = ca.t_start + (best_start as f64 / na) * span_a;
    let t_end_a = ca.t_start + ((best_start + best_len - 1) as f64 / na) * span_a;
    // Find which edge_b samples are closest to the overlapping edge_a points.
    let overlap_a = &pts_a[best_start..(best_start + best_len)];
    let (mut min_b, mut max_b) = (pts_b.len(), 0_usize);
    for pa in overlap_a {
        let idx = pts_b
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                pa.distance_to(a)
                    .partial_cmp(&pa.distance_to(b))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(i, _)| i)?;
        if idx < min_b {
            min_b = idx;
        }
        if idx > max_b {
            max_b = idx;
        }
    }
    let nb = (pts_b.len() - 1) as f64;
    let cb = brep.edges.get(edge_b)?.curve.as_ref()?;
    let span_b = cb.t_end - cb.t_start;
    let t_start_b = cb.t_start + (min_b as f64 / nb) * span_b;
    let t_end_b = cb.t_start + (max_b as f64 / nb) * span_b;
    Some(((t_start_a, t_end_a), (t_start_b, t_end_b)))
}
