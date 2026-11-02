//! A boolean's candidate-vertex table: its candidates, their interior
//! incidences and crossing zones, and their resolution, built once for
//! the shadow comparison and for graph node identity.

use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, EdgeHandle};

use super::candidates::{CandidateVertex, Incidence, SourceKind, collect_candidates};
use super::resolve::{
    Resolution, ResolutionContext, ResolutionDecline, ResolutionOptions, resolve,
};
use crate::intersect::IntersectionRecord;
use crate::intersect_edge_face::EdgeFaceHit;

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

/// A record, or a piece of one, lying along an operand edge or a piece
/// of it between two resolved vertices: a common block (plan §4.3a,
/// from E-E compatibility). One edge of the result: the graph keeps
/// the record's piece and drops the edge's from the face both lie on.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct RecordEdgeBlock {
    /// The record.
    pub record: usize,
    /// The edge's operand (0 = A).
    pub operand: usize,
    /// The edge.
    pub edge: EdgeHandle,
    /// The two vertices the block runs between, as their anchors.
    pub ends: [Point3; 2],
}

/// The candidates and their resolution.
pub(crate) struct VertexTable {
    /// The candidates, sorted by key.
    pub cands: Vec<CandidateVertex>,
    /// Their resolution.
    pub resolution: Result<Resolution, ResolutionDecline>,
}

/// Whether the table decides the graph's record-end vertices: always,
/// unless `POPKERN_VERTEX_TABLE=0` (the kill switch) restores the
/// graph's own snap-then-adopt identity. A table that declines falls
/// back to that identity for its boolean, recorded as a
/// `VertexTableFallback` shadow line and a decline artifact.
pub(crate) fn table_enabled() -> bool {
    std::env::var_os("POPKERN_VERTEX_TABLE").is_none_or(|v| v != "0")
}

/// Build the table of a boolean between `a` and `b` whose graph inserts
/// `records`, with `hits` its interfering edge-face hits
/// (`candidates::interfering_hits`).
pub(crate) fn build(
    (a, b): (&BRep, &BRep),
    hits: &[EdgeFaceHit],
    records: &[IntersectionRecord],
    linear: f64,
) -> VertexTable {
    let cap = GROWTH_CAP_FACTOR * linear;
    let mut cands = collect_candidates(a, b, hits, records, linear, cap);
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
    if resolution.is_err() {
        emit_decline_artifact(&cands, &resolution, records);
    }
    VertexTable { cands, resolution }
}

/// With `POPKERN_COMPLEX_ARTIFACT=<dir>`, write a declined resolution's
/// candidates and decline to `<dir>` as JSON (`resolution_record`), so
/// the decline can be inspected without rerunning the boolean.
fn emit_decline_artifact(
    cands: &[CandidateVertex],
    resolution: &Result<Resolution, ResolutionDecline>,
    records: &[IntersectionRecord],
) {
    use std::sync::atomic::{AtomicUsize, Ordering};
    static SEQUENCE: AtomicUsize = AtomicUsize::new(0);
    let Some(dir) = std::env::var_os("POPKERN_COMPLEX_ARTIFACT") else {
        return;
    };
    let n = SEQUENCE.fetch_add(1, Ordering::Relaxed);
    let path = std::path::Path::new(&dir).join(format!(
        "vertex-table-decline-{}-{n}.json",
        std::process::id()
    ));
    let body = super::resolution_record::resolution_json(cands, resolution, records, &[]);
    if std::fs::write(&path, body.to_string()).is_ok() {
        eprintln!("[vertex table] decline artifact: {}", path.display());
    }
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

    /// Every common block of a record with an operand edge: whole
    /// records lying along an edge (`common_blocks`) and record pieces
    /// along pieces of an edge (`pave_blocks`). Empty when the
    /// resolution declined.
    pub(crate) fn record_edge_blocks(
        &self,
        operands: [&BRep; 2],
        records: &[IntersectionRecord],
        linear: f64,
    ) -> Vec<RecordEdgeBlock> {
        use super::pave_blocks::BlockOwner;
        let Ok(res) = self.resolution.as_ref() else {
            return Vec::new();
        };
        let anchor = |g: usize| self.cands[res.groups[g].anchor].position;
        let mut out = Vec::new();
        let ends = super::consequences::record_end_candidates(&self.cands, records.len());
        for (record, operand, edge) in
            super::common_blocks::common_block_records(&self.cands, res, operands, records, linear)
        {
            if let [Some(s), Some(e)] = ends[record] {
                out.push(RecordEdgeBlock {
                    record,
                    operand,
                    edge,
                    ends: [anchor(res.group_of[s]), anchor(res.group_of[e])],
                });
            }
        }
        let blocks = super::pave_blocks::pave_blocks(&self.cands, res, operands, records, linear);
        for (i, j) in super::pave_blocks::common_block_pairs(&blocks, res, operands, records) {
            let (p, q) = (&blocks[i], &blocks[j]);
            let (record, operand, bits) = match (p.owner, q.owner) {
                (BlockOwner::Record(r), BlockOwner::Edge(o, e))
                | (BlockOwner::Edge(o, e), BlockOwner::Record(r)) => (r, o, e),
                _ => continue,
            };
            // Paired blocks join the same two vertices.
            let [Some(ga), Some(gb)] = p.ends else {
                continue;
            };
            let Some(edge) = operands[operand]
                .edges
                .keys()
                .find(|eh| slotmap::Key::data(eh).as_ffi() == bits)
            else {
                continue;
            };
            out.push(RecordEdgeBlock {
                record,
                operand,
                edge,
                ends: [anchor(ga), anchor(gb)],
            });
        }
        out
    }
}
