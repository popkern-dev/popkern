//! Measured boundary-row targeting for KP20 support-cell reconciliation.

use std::collections::HashMap;

use popkern_core::error::CoreError;
use popkern_core::topology::{BRep, FaceHandle};

use super::{BoundarySide, FullBoundaryTrim, rectangular_binonic_on_face, trim_gap};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct BoundaryRowGroupKey {
    face: FaceHandle,
    side: BoundarySide,
}

#[derive(Debug, Clone)]
struct BoundaryRowGroupScore {
    worst_gap: f64,
    trims: Vec<FullBoundaryTrim>,
}

pub(super) fn worst_boundary_row_group_trims(
    brep: &BRep,
    trims: &[FullBoundaryTrim],
) -> Result<Vec<FullBoundaryTrim>, CoreError> {
    let mut groups: HashMap<BoundaryRowGroupKey, BoundaryRowGroupScore> = HashMap::new();
    for &trim in trims {
        let patch_a = rectangular_binonic_on_face(brep, trim.face_a)?;
        let patch_b = rectangular_binonic_on_face(brep, trim.face_b)?;
        let gap = trim_gap(&patch_a, &patch_b, trim);
        push_boundary_row_group(&mut groups, trim.face_a, trim.side_a, trim, gap);
        push_boundary_row_group(&mut groups, trim.face_b, trim.side_b, trim, gap);
    }
    // Full deterministic ordering: HashMap value order must not
    // break (worst_gap, len) ties — the group key is the final
    // tie-break so two runs always pick the same group.
    Ok(groups
        .into_iter()
        .max_by(|(key_a, a), (key_b, b)| {
            a.worst_gap
                .total_cmp(&b.worst_gap)
                .then_with(|| a.trims.len().cmp(&b.trims.len()))
                .then_with(|| format!("{key_a:?}").cmp(&format!("{key_b:?}")))
        })
        .map(|(_, group)| group.trims)
        .unwrap_or_default())
}

fn push_boundary_row_group(
    groups: &mut HashMap<BoundaryRowGroupKey, BoundaryRowGroupScore>,
    face: FaceHandle,
    side: BoundarySide,
    trim: FullBoundaryTrim,
    gap: f64,
) {
    let score = groups
        .entry(BoundaryRowGroupKey { face, side })
        .or_insert_with(|| BoundaryRowGroupScore {
            worst_gap: 0.0,
            trims: Vec::new(),
        });
    score.worst_gap = score.worst_gap.max(gap);
    if !score
        .trims
        .iter()
        .any(|existing| same_trim(*existing, trim))
    {
        score.trims.push(trim);
    }
}

fn same_trim(a: FullBoundaryTrim, b: FullBoundaryTrim) -> bool {
    a.edge == b.edge
        && a.half_edge_a == b.half_edge_a
        && a.half_edge_b == b.half_edge_b
        && a.face_a == b.face_a
        && a.side_a == b.side_a
        && a.face_b == b.face_b
        && a.side_b == b.side_b
        && (a.start_a - b.start_a).abs() <= 1.0e-12
        && (a.end_a - b.end_a).abs() <= 1.0e-12
        && (a.start_b - b.start_b).abs() <= 1.0e-12
        && (a.end_b - b.end_b).abs() <= 1.0e-12
}
