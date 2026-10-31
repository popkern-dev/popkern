//! Branches of one intersection meeting (plan §4.3a, the branch rule).
//!
//! Two records of one face pair are two branches of that pair's
//! intersection, or two pieces of one branch cut at a hit, a seam or a
//! band lattice. Where the faces cross transversally the intersection
//! is locally one smooth curve, so two branches meet only at a singular
//! point, and the SSI's own tangency analysis is the evidence: it ends
//! a branch there. Two records ending at one point are pieces of one
//! branch or branches the SSI stopped at a branch point; a record
//! whose interior passes within tolerance of another record of the same
//! pair was marched through, the SSI seeing no branch point there, and
//! is a second branch passing close.

use super::candidates::{CandidateVertex, Incidence};
use super::resolve::Resolution;
use crate::intersect::IntersectionRecord;

/// Whether records `ra` and `rb` are distinct records of one face pair.
pub(crate) fn same_pair(records: &[IntersectionRecord], ra: usize, rb: usize) -> bool {
    match (records.get(ra), records.get(rb)) {
        (Some(a), Some(b)) => ra != rb && a.face_a == b.face_a && a.face_b == b.face_b,
        _ => false,
    }
}

/// Groups where a record's interior meets another record of the same
/// face pair: what the branch rule would keep apart.
pub(crate) fn interior_branch_meets(
    cands: &[CandidateVertex],
    res: &Resolution,
    records: &[IntersectionRecord],
) -> usize {
    res.groups
        .iter()
        .filter(|g| {
            let on: Vec<(usize, bool)> = g
                .members
                .iter()
                .flat_map(|&m| cands[m].incidence.iter())
                .filter_map(|i| match *i {
                    Incidence::Record { record, end, .. } => Some((record, end.is_none())),
                    _ => None,
                })
                .collect();
            on.iter().any(|&(ra, interior)| {
                interior && on.iter().any(|&(rb, _)| same_pair(records, ra, rb))
            })
        })
        .count()
}
