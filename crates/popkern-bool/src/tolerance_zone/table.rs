//! A boolean's candidate-vertex table: its candidates, their interior
//! incidences and crossing zones, and their resolution, built once for
//! the shadow comparison and for graph node identity.

use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, EdgeHandle};

use super::candidates::{CandidateVertex, Incidence, SourceKind, collect_candidates};
use super::resolve::{
    Resolution, ResolutionContext, ResolutionDecline, ResolutionOptions, resolve,
};
use crate::intersect::{IntersectionRecord, IntersectionResult};

/// Multiple of the linear tolerance a vertex may grow to without a
/// crossing angle explaining it: the graph's own record-end adoption
/// radius, so the table never identifies what the graph would not.
const GROWTH_CAP_FACTOR: f64 = 20.0;

/// The largest component searched exactly for its fewest vertices
/// (plan §4.3a: components are small in practice; the search over 12
/// candidates visits at most the Bell number B₁₂ ≈ 4.2e6 partitions,
/// pruned by the best found and the lower bound).
const EXACT_PARTITION_LIMIT: usize = 12;

/// Where a record's vertices lie, as the table decided them.
#[derive(Debug, Clone)]
pub(crate) struct RecordEnds {
    /// The start's vertex: the anchor of the group it resolved into.
    pub start: Point3,
    /// The end's vertex.
    pub end: Point3,
    /// The larger of the two vertices' tolerances: what decides a
    /// record whose ends share a node (`consequences::collapse_of`).
    pub tolerance: f64,
    /// The record's T-junctions: each group holding an interior
    /// incidence on the record, as (record parameter, the group's
    /// anchor), by parameter. The record's edge is split there at
    /// insertion (a pave on the record, plan §4.3).
    pub paves: Vec<(f64, Point3)>,
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
    let resolution = resolve(
        &cands,
        &ctx,
        ResolutionOptions {
            growth_cap: cap,
            exact_limit: EXACT_PARTITION_LIMIT,
        },
    );
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
        let mut paves = self.record_paves(res, records);
        ends.into_iter()
            .zip(paves.iter_mut())
            .map(|([s, e], paves)| {
                let ((start, ts), (end, te)) = (s?, e?);
                Some(RecordEnds {
                    start,
                    end,
                    tolerance: ts.max(te),
                    paves: std::mem::take(paves),
                })
            })
            .collect()
    }

    /// Each record's interior paves, one per group however many of its
    /// candidates fall on the record, at the parameter of the first. A
    /// group holding one of the record's own ends is that end, not a
    /// pave.
    fn record_paves(&self, res: &Resolution, records: usize) -> Vec<Vec<(f64, Point3)>> {
        let mut seen = std::collections::BTreeSet::new();
        for (record, ends) in super::consequences::record_end_candidates(&self.cands, records)
            .into_iter()
            .enumerate()
        {
            for c in ends.into_iter().flatten() {
                seen.insert((record, res.group_of[c]));
            }
        }
        let mut paves: Vec<Vec<(f64, Point3)>> = vec![Vec::new(); records];
        for (i, c) in self.cands.iter().enumerate() {
            let group = res.group_of[i];
            for inc in &c.incidence {
                if let &Incidence::Record {
                    record,
                    s,
                    end: None,
                } = inc
                    && record < records
                    && seen.insert((record, group))
                {
                    let anchor = self.cands[res.groups[group].anchor].position;
                    paves[record].push((s, anchor));
                }
            }
        }
        for p in &mut paves {
            p.sort_by(|a, b| a.0.total_cmp(&b.0));
        }
        paves
    }

    /// Every operand edge a record node lies on mid-span: (operand,
    /// edge, the group's anchor) for each group holding a record
    /// incidence (an end, or a pave the record is split at) and an
    /// incidence on the edge, one per group and edge. The edge is split
    /// there (a pave on the operand edge, plan §4.3).
    pub(crate) fn operand_edge_paves(&self) -> Vec<(usize, EdgeHandle, Point3)> {
        let Ok(res) = self.resolution.as_ref() else {
            return Vec::new();
        };
        let record_groups: std::collections::BTreeSet<usize> = (0..self.cands.len())
            .filter(|&i| {
                self.cands[i]
                    .incidence
                    .iter()
                    .any(|inc| matches!(inc, Incidence::Record { .. }))
            })
            .map(|i| res.group_of[i])
            .collect();
        let mut seen = std::collections::BTreeSet::new();
        let mut out = Vec::new();
        for (i, c) in self.cands.iter().enumerate() {
            let group = res.group_of[i];
            if !record_groups.contains(&group) {
                continue;
            }
            for inc in &c.incidence {
                if let &Incidence::Edge { operand, edge, .. } = inc
                    && seen.insert((group, operand, edge))
                {
                    out.push((operand, edge, self.cands[res.groups[group].anchor].position));
                }
            }
        }
        out
    }

    /// Whether some record is a common block with an operand edge: a
    /// same-domain overlap (plan §4.4), where the graph's own splitting
    /// still decides the junctions until the join does.
    pub(crate) fn has_common_blocks(
        &self,
        operands: [&BRep; 2],
        records: &[IntersectionRecord],
        linear: f64,
    ) -> bool {
        self.resolution.as_ref().is_ok_and(|res| {
            !super::common_blocks::common_block_records(&self.cands, res, operands, records, linear)
                .is_empty()
        })
    }
}
