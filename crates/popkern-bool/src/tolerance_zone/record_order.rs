//! The record-order rule (plan §4.3a): a vertex never reverses the
//! order of points along a record. Points s₁ < s₂ < s₃ of one open
//! record, s₁ and s₃ one vertex and s₂ another, would make the record
//! run V → W → V: an edge out to W and back that the record never
//! traced. Separation along the record bounds how far apart two members
//! of one vertex may be, but not what lies between them, so the rule
//! looks at the whole partition, not at a pair.
//!
//! A candidate joining a group is kept out when a point of another
//! group lies strictly between it and a member along a record they
//! share. A candidate that must start a vertex of its own inside
//! another vertex's span along a record (it is incompatible with that
//! vertex by another rule) is a conflict no join choice in canonical
//! order avoids; the resolution records it for the exact partition
//! search to settle.
//!
//! Closed records have no order to keep: their parameter wraps.
//!
//! Measured (bool, examples, mesh, blend; table on): 2788/2788 pass,
//! and across 1449 shadowed booleans the rule keeps no pair apart and
//! records no conflict. It guards an order the canonical partition
//! already kept there; its fixture shows what it catches.

use super::candidates::{CandidateVertex, Incidence};
use super::resolve::VertexGroup;
use crate::intersect::IntersectionRecord;

/// The open records `c` lies on, each with its parameter.
fn on_open_records<'a>(
    c: &'a CandidateVertex,
    records: &'a [IntersectionRecord],
) -> impl Iterator<Item = (usize, f64)> + 'a {
    c.incidence.iter().filter_map(move |i| match *i {
        Incidence::Record { record, s, .. }
            if records.get(record).is_some_and(|r| !r.curve.is_closed) =>
        {
            Some((record, s))
        }
        _ => None,
    })
}

/// The span `group`'s members cover along `record`, `None` when none
/// lies on it.
fn span(
    cands: &[CandidateVertex],
    records: &[IntersectionRecord],
    group: &VertexGroup,
    record: usize,
) -> Option<(f64, f64)> {
    group
        .members
        .iter()
        .flat_map(|&m| on_open_records(&cands[m], records))
        .filter(|&(r, _)| r == record)
        .fold(None, |acc, (_, s)| {
            Some(acc.map_or((s, s), |(lo, hi): (f64, f64)| (lo.min(s), hi.max(s))))
        })
}

/// A point of a group other than the one `c` joins (`joining`, `None`
/// when `c` starts a vertex of its own) that puts `c` out of order
/// along a record they share: it lies strictly between `c` and a member
/// of the group `c` joins, or `c` lies strictly between two members of
/// its group. Returns the point and the record.
pub(crate) fn reversal(
    cands: &[CandidateVertex],
    records: &[IntersectionRecord],
    groups: &[VertexGroup],
    joining: Option<usize>,
    c: usize,
) -> Option<(usize, usize)> {
    for (r, s) in on_open_records(&cands[c], records) {
        let (lo, hi) = joining
            .and_then(|g| span(cands, records, &groups[g], r))
            .map_or((s, s), |(lo, hi)| (lo.min(s), hi.max(s)));
        for (h, other) in groups.iter().enumerate() {
            if Some(h) == joining {
                continue;
            }
            if let Some((a, b)) = span(cands, records, other, r)
                && s > a
                && s < b
            {
                return Some((other.members[0], r));
            }
            for &w in &other.members {
                if on_open_records(&cands[w], records).any(|(rw, sw)| rw == r && sw > lo && sw < hi)
                {
                    return Some((w, r));
                }
            }
        }
    }
    None
}
