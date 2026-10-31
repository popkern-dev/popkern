//! A boolean's candidate-vertex table: its candidates, their interior
//! incidences and crossing zones, and their resolution, built once for
//! the shadow comparison and for graph node identity.

use popkern_core::geometry::Point3;
use popkern_core::topology::BRep;

use super::candidates::{CandidateVertex, Incidence, SourceKind, collect_candidates};
use super::resolve::{
    Resolution, ResolutionContext, ResolutionDecline, ResolutionOptions, resolve,
};
use crate::intersect::{IntersectionRecord, IntersectionResult};

/// Multiple of the linear tolerance a vertex may grow to without a
/// crossing angle explaining it: the graph's own record-end adoption
/// radius, so the table never identifies what the graph would not.
const GROWTH_CAP_FACTOR: f64 = 20.0;

/// Where a record's ends lie, as the table decided them.
#[derive(Debug, Clone, Copy)]
pub(crate) struct RecordEnds {
    /// The start's vertex: the anchor of the group it resolved into.
    pub start: Point3,
    /// The end's vertex.
    pub end: Point3,
    /// The larger of the two vertices' tolerances: what decides a
    /// record whose ends share a node (`consequences::collapse_of`).
    pub tolerance: f64,
}

/// The candidates and their resolution.
pub(crate) struct VertexTable {
    /// The candidates, sorted by key.
    pub cands: Vec<CandidateVertex>,
    /// Their resolution.
    pub resolution: Result<Resolution, ResolutionDecline>,
}

/// Whether `POPKERN_VERTEX_TABLE=1` asks the table to decide the graph's
/// record-end vertices (plan Stage 2: the table replaces the graph's
/// snap-then-adopt identity once its shadow agrees).
pub(crate) fn table_enabled() -> bool {
    std::env::var_os("POPKERN_VERTEX_TABLE").is_some_and(|v| v == "1")
}

/// Build the table of a boolean between `a` and `b` whose graph inserts
/// `records`.
pub(crate) fn build(
    (a, b): (&BRep, &BRep),
    ix: &IntersectionResult,
    records: &[IntersectionRecord],
    linear: f64,
) -> VertexTable {
    let cap = GROWTH_CAP_FACTOR * linear;
    let mut cands = collect_candidates(a, b, ix, records, linear, cap);
    super::incidence::add_interior_incidences(&mut cands, [a, b], records, linear);
    super::incidence::add_crossing_zones(&mut cands, [a, b], cap);
    let ctx = ResolutionContext {
        operands: [a, b],
        records,
    };
    let resolution = resolve(&cands, &ctx, ResolutionOptions { growth_cap: cap });
    VertexTable { cands, resolution }
}

impl VertexTable {
    /// Each record's end vertices: the anchors of the groups its ends
    /// resolved into. `None` when the resolution declined or a record
    /// end has no candidate.
    pub(crate) fn record_ends(&self, records: usize) -> Option<Vec<RecordEnds>> {
        let res = self.resolution.as_ref().ok()?;
        let mut ends: Vec<[Option<(Point3, f64)>; 2]> = vec![[None; 2]; records];
        for (i, c) in self.cands.iter().enumerate() {
            if c.source != SourceKind::RecordEndpoint {
                continue;
            }
            let Some(&Incidence::Record {
                record,
                end: Some(end),
                ..
            }) = c.incidence.first()
            else {
                continue;
            };
            let group = &res.groups[res.group_of[i]];
            let anchor = self.cands[group.anchor].position;
            *ends.get_mut(record)?.get_mut(usize::from(end))? = Some((anchor, group.tolerance));
        }
        ends.into_iter()
            .map(|[s, e]| {
                let ((start, ts), (end, te)) = (s?, e?);
                Some(RecordEnds {
                    start,
                    end,
                    tolerance: ts.max(te),
                })
            })
            .collect()
    }
}
