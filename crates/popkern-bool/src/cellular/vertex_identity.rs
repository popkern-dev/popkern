//! What the vertex table decides for the cellular graph (plan §4.3a):
//! where operand vertices and record ends sit, the paves splitting
//! records and operand edges, and the record/edge common blocks.
//!
//! The table is built before any graph entity, from the operands, the
//! processed records and the interfering edge-face hits, so every node
//! the graph inserts for a resolved vertex sits at its anchor. A
//! declined table, or `POPKERN_VERTEX_TABLE=0`, leaves the graph its
//! own snap-then-adopt identity (`graph_identity`).

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

use crate::intersect::{IntersectionRecord, IntersectionResult};
use crate::intersect_edge_face::EdgeFaceHit;
use crate::tolerance_zone::table::{RecordEdgeBlock, RecordEnds};

/// The table's decisions, or none of them under graph identity.
#[derive(Default)]
pub(super) struct VertexIdentity {
    /// The edge-face hits the graph splits operand edges at.
    pub hits: Vec<EdgeFaceHit>,
    /// Each record's ends and paves; `None` under graph identity.
    pub ends: Option<Vec<RecordEnds>>,
    /// Where each operand's table-held vertices sit.
    pub vertex_anchors: [HashMap<VertexHandle, Point3>; 2],
    /// The operand edges' paves: (operand, edge, anchor).
    pub operand_paves: Vec<(usize, EdgeHandle, Point3)>,
    /// Records, or pieces, along operand edges (partial-overlap scenes).
    pub common_blocks: Vec<RecordEdgeBlock>,
    /// Each decided vertex's anchor and grown tolerance.
    pub reaches: Vec<(Point3, f64)>,
}

impl VertexIdentity {
    /// Whether the graph decides vertices itself: the table is off or
    /// declined, and the graph's split passes and snaps run.
    pub(super) fn graph_identity(&self) -> bool {
        self.ends.is_none()
    }
}

/// Build the table over `records` and decide the graph's vertices,
/// returning the records as the graph inserts them: a closed record
/// whose closure no other entity needs closes at one of its paves
/// (`tolerance_zone::closure`), and the table is rebuilt over the moved
/// records, so its ends, paves and blocks describe them. Under graph
/// identity the records close at an edge-face hit on them instead
/// (`record_rotation`), with the kill switch from the start.
pub(super) fn decide(
    (a, b): (&BRep, &BRep),
    ix: &IntersectionResult,
    mut records: Vec<IntersectionRecord>,
    tolerance: &popkern_core::tolerance::Tolerance,
    partial_overlaps: bool,
) -> (Vec<IntersectionRecord>, VertexIdentity) {
    use crate::tolerance_zone::{candidates, table};
    let linear = tolerance.linear();
    if !table::table_enabled() {
        let identity = VertexIdentity {
            hits: ix.edge_face_hits.clone(),
            ..VertexIdentity::default()
        };
        return (records, identity);
    }
    // The operand edges split only at the hits that are interferences
    // (on their face), the table's own candidates.
    let mut hits = candidates::interfering_hits((a, b), ix, &records, linear);
    let mut table = table::build((a, b), &hits, &records, linear);
    let moves = table.closure_moves(&records);
    if moves.iter().any(Option::is_some) {
        for (record, mv) in records.iter_mut().zip(&moves) {
            let moved =
                mv.and_then(|m| super::record_rotation::rotate_closed_record_at(record, m.s, m.at));
            if let Some(moved) = moved {
                *record = moved;
            }
        }
        crate::complex_shadow::subsumed_pass_fired(
            "closure_moved",
            moves.iter().flatten().count(),
            0.0,
        );
        hits = candidates::interfering_hits((a, b), ix, &records, linear);
        table = table::build((a, b), &hits, &records, linear);
    }
    let Some(ends) = table.record_ends(records.len()) else {
        crate::complex_shadow::append_log_line(&format!(
            "VertexTableFallback {:?}",
            table.resolution.as_ref().err()
        ));
        let records = records
            .iter()
            .map(|r| {
                super::record_rotation::rotate_closed_record_to_hits(
                    r,
                    &ix.edge_face_hits,
                    tolerance,
                )
            })
            .collect();
        let identity = VertexIdentity {
            hits,
            ..VertexIdentity::default()
        };
        return (records, identity);
    };
    let identity = VertexIdentity {
        vertex_anchors: table.operand_vertex_anchors(),
        operand_paves: table.operand_edge_paves(),
        common_blocks: if partial_overlaps {
            table.record_edge_blocks([a, b], &records, linear)
        } else {
            Vec::new()
        },
        ends: Some(ends),
        reaches: table.vertex_reaches(),
        hits,
    };
    (records, identity)
}

#[cfg(test)]
#[path = "vertex_identity_tests.rs"]
mod tests;
