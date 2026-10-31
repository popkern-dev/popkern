//! Consequences of a resolution for the records (plan §4.3a): a record
//! whose two ends resolved into one vertex is a micro edge when its
//! whole extent stays within that vertex's tolerance (it collapses into
//! the vertex), and a loop through the vertex otherwise. A closed
//! record's ends meet by construction and are neither.
//!
//! **Stage.** Measured in shadow beside the graph's point-contact guard
//! (`insert_intersection_curves_at`), which skips a record whose ends
//! share a node and whose extent stays within 500 tolerances (capped by
//! the model size): the table decides by the vertex's own tolerance
//! instead, and the counts show which records the two disagree on
//! before either decides the graph.

use popkern_core::geometry::curve::CurveEval;

use super::candidates::{CandidateVertex, Incidence, SourceKind};
use super::resolve::Resolution;
use crate::intersect::IntersectionRecord;

/// Samples along a record for its extent from the vertex.
const EXTENT_SAMPLES: usize = 8;

/// Each record's (start, end) candidate, `None` where it has none.
pub(crate) fn record_end_candidates(
    cands: &[CandidateVertex],
    records: usize,
) -> Vec<[Option<usize>; 2]> {
    let mut ends: Vec<[Option<usize>; 2]> = vec![[None; 2]; records];
    for (i, c) in cands.iter().enumerate() {
        if c.source != SourceKind::RecordEndpoint {
            continue;
        }
        if let Some(&Incidence::Record {
            record,
            end: Some(end),
            ..
        }) = c.incidence.first()
            && let Some(slot) = ends.get_mut(record)
        {
            slot[usize::from(end)] = Some(i);
        }
    }
    ends
}

/// A record whose ends resolved into one vertex.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct CollapsedEnds {
    /// The record.
    pub record: usize,
    /// The farthest the record strays from the vertex.
    pub extent: f64,
    /// The vertex's tolerance.
    pub tolerance: f64,
}

impl CollapsedEnds {
    /// Whether the whole record lies within the vertex: a micro edge,
    /// which collapses into it.
    pub(crate) fn is_micro_edge(&self) -> bool {
        self.extent <= self.tolerance
    }
}

/// Every open record of `records` whose two ends resolved into one
/// vertex.
pub(crate) fn collapsed_ends(
    cands: &[CandidateVertex],
    res: &Resolution,
    records: &[IntersectionRecord],
) -> Vec<CollapsedEnds> {
    let mut out = Vec::new();
    for (ri, [s, e]) in record_end_candidates(cands, records.len())
        .into_iter()
        .enumerate()
    {
        let (Some(s), Some(e)) = (s, e) else {
            continue;
        };
        if records[ri].curve.is_closed {
            continue;
        }
        let group = res.group_of[s];
        if group != res.group_of[e] {
            continue;
        }
        let vertex = &res.groups[group];
        let anchor = cands[vertex.anchor].position;
        let curve = &records[ri].curve.curve_3d;
        let (t0, t1) = curve.parameter_range();
        let extent = (0..=EXTENT_SAMPLES)
            .map(|k| {
                let t = t0 + (t1 - t0) * k as f64 / EXTENT_SAMPLES as f64;
                curve.point_at(t).distance_to(&anchor)
            })
            .fold(0.0_f64, f64::max);
        out.push(CollapsedEnds {
            record: ri,
            extent,
            tolerance: vertex.tolerance,
        });
    }
    out
}
