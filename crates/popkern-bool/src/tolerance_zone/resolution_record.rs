//! Resolution records (plan §4.3a output): for every vertex the
//! resolution committed, its members, anchor and growth, the rule
//! behind every split that kept a candidate out of it, and the
//! consequences it carries, so a decline or a surprising vertex can be
//! explained afterwards without rerunning the boolean.

use serde_json::{Value, json};

use super::candidates::CandidateVertex;
use super::consequences::{Collapse, OpenIntersectionLoop};
use super::pave_blocks::{BlockOwner, PaveBlock};
use super::resolve::{Incompatibility, Resolution, ResolutionDecline};
use crate::intersect::IntersectionRecord;

/// A consequence a vertex carries.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum VertexConsequence {
    /// A member lies in another edge's or record's interior: a pave.
    TJunction,
    /// A record's two ends resolved into the vertex.
    CollapsedEnds {
        /// The record.
        record: usize,
        /// What the record became.
        kind: Collapse,
    },
    /// A record end resolved into the vertex onto no boundary.
    OpenEnd(OpenIntersectionLoop),
    /// A pave block from this vertex lies within it and the next.
    MicroBlock {
        /// The block's owner.
        owner: BlockOwner,
    },
    /// A member started the vertex inside another's span along a
    /// record.
    OrderConflict {
        /// The record.
        record: usize,
    },
}

/// The record of one committed vertex.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct VertexRecord {
    /// The anchor candidate.
    pub anchor: usize,
    /// Every member, the anchor first.
    pub members: Vec<usize>,
    /// The vertex tolerance.
    pub tolerance: f64,
    /// How far the tolerance grew past the anchor's radius.
    pub growth: f64,
    /// Every split involving a member: (member, other, rule).
    pub splits: Vec<(usize, usize, Incompatibility)>,
    /// The consequences the vertex carries.
    pub consequences: Vec<VertexConsequence>,
}

/// The record of every committed vertex, by group.
pub(crate) fn vertex_records(
    cands: &[CandidateVertex],
    res: &Resolution,
    records: &[IntersectionRecord],
    blocks: &[PaveBlock],
) -> Vec<VertexRecord> {
    let mut out: Vec<VertexRecord> = res
        .groups
        .iter()
        .map(|g| VertexRecord {
            anchor: g.anchor,
            members: g.members.clone(),
            tolerance: g.tolerance,
            growth: g.tolerance - cands[g.anchor].radius,
            splits: Vec::new(),
            consequences: Vec::new(),
        })
        .collect();
    for &(m, c, why) in &res.separations {
        out[res.group_of[m]].splits.push((m, c, why));
        out[res.group_of[c]].splits.push((c, m, why));
    }
    for (g, group) in res.groups.iter().enumerate() {
        if group
            .members
            .iter()
            .any(|&m| super::incidence::is_t_junction(&cands[m]))
        {
            out[g].consequences.push(VertexConsequence::TJunction);
        }
    }
    let ends = super::consequences::record_end_candidates(cands, records.len());
    for c in super::consequences::collapsed_ends(cands, res, records) {
        if let Some(s) = ends[c.record][0] {
            out[res.group_of[s]]
                .consequences
                .push(VertexConsequence::CollapsedEnds {
                    record: c.record,
                    kind: c.kind(),
                });
        }
    }
    for open in super::consequences::open_loops(cands, res) {
        if let Some(i) = ends[open.record][usize::from(open.end)] {
            out[res.group_of[i]]
                .consequences
                .push(VertexConsequence::OpenEnd(open));
        }
    }
    for b in blocks.iter().filter(|b| b.shrunk.is_none()) {
        if let Some(g) = b.ends[0] {
            out[g]
                .consequences
                .push(VertexConsequence::MicroBlock { owner: b.owner });
        }
    }
    for &(_, c, record) in &res.order_conflicts {
        out[res.group_of[c]]
            .consequences
            .push(VertexConsequence::OrderConflict { record });
    }
    out
}

fn point(p: popkern_core::geometry::Point3) -> Value {
    json!([p.x, p.y, p.z])
}

/// The candidates as JSON: position, radius, zone, source, incidence.
fn candidates_json(cands: &[CandidateVertex]) -> Value {
    Value::Array(
        cands
            .iter()
            .map(|c| {
                json!({
                    "position": point(c.position),
                    "radius": c.radius,
                    "zone": c.zone,
                    "source": format!("{:?}", c.source),
                    "incidence": c.incidence.iter().map(|i| format!("{i:?}")).collect::<Vec<_>>(),
                })
            })
            .collect(),
    )
}

/// A resolution as JSON: the candidates and, when it committed, every
/// vertex record; when it declined, the decline.
pub(crate) fn resolution_json(
    cands: &[CandidateVertex],
    resolution: &Result<Resolution, ResolutionDecline>,
    records: &[IntersectionRecord],
    blocks: &[PaveBlock],
) -> Value {
    let outcome = match resolution {
        Ok(res) => json!({
            "vertices": vertex_records(cands, res, records, blocks)
                .iter()
                .map(|v| json!({
                    "anchor": v.anchor,
                    "position": point(cands[v.anchor].position),
                    "members": v.members,
                    "tolerance": v.tolerance,
                    "growth": v.growth,
                    "splits": v.splits.iter().map(|s| format!("{s:?}")).collect::<Vec<_>>(),
                    "consequences": v.consequences.iter().map(|c| format!("{c:?}")).collect::<Vec<_>>(),
                }))
                .collect::<Vec<_>>(),
        }),
        Err(d) => json!({ "decline": format!("{d:?}") }),
    };
    json!({ "candidates": candidates_json(cands), "resolution": outcome })
}
