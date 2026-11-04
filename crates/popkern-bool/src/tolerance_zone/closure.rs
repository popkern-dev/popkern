//! Where a closed record closes (plan §4.3a: the closure is a vertex
//! like any other, decided by the table).
//!
//! A closed record's two ends meet at a point its tracer happened to
//! start from. When nothing else lies there, the point is no vertex of
//! the result: a node there splits the record, and every operand edge
//! it lies on, for nothing (a ring flush on a shaft, rotated, closed
//! its rim records off the seam and split both rims there). The record
//! closes instead at one of its paves, a vertex the result has anyway.
//!
//! A closure group is free when every member is the end of a closed
//! record, its faces' and edge incidences allowed: no hit, operand vertex or open
//! record end needs a vertex there. Coincident closed records (a cone's
//! base circle traced from two face pairs) share their closure group,
//! so they move together, to a pave group every one of them holds; a
//! group with no common pave stays where it is.
//!
//! A record traced round to an end the table cannot make one vertex
//! with its start (a torn closure) closes at a pave too, which welds
//! its two ends into one.

use std::collections::BTreeSet;

use popkern_core::geometry::Point3;

use super::candidates::{Incidence, SourceKind};
use super::resolve::Resolution;
use super::table::VertexTable;
use crate::intersect::IntersectionRecord;

/// Where a record closes instead of its own closure: the parameter of
/// the pave on it, and the pave's anchor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ClosureMove {
    /// The record parameter the pave sits at.
    pub s: f64,
    /// The pave's vertex.
    pub at: Point3,
}

impl VertexTable {
    /// Each record's closure move: `Some` for a closed record whose
    /// closure group is free and whose closure group's records share a
    /// pave group, or whose two ends resolved apart and which holds a
    /// pave. Empty when the resolution declined.
    pub(crate) fn closure_moves(&self, records: &[IntersectionRecord]) -> Vec<Option<ClosureMove>> {
        let mut out = vec![None; records.len()];
        let Ok(res) = self.resolution.as_ref() else {
            return out;
        };
        let ends = super::consequences::record_end_candidates(&self.cands, records.len());
        let closed = |r: usize| records.get(r).is_some_and(|rec| rec.curve.is_closed);
        let mut done = BTreeSet::new();
        for (r, [s, e]) in ends.iter().enumerate() {
            let (Some(s), Some(e)) = (*s, *e) else {
                continue;
            };
            let (g, ge) = (res.group_of[s], res.group_of[e]);
            if !closed(r) {
                continue;
            }
            if g != ge {
                // A torn closure: the ends lie apart beyond what one
                // vertex holds (measured: a plane chord of a torus tube
                // traced round from 5.8e-8 beside its seam crossing,
                // ending 1.2e-7 from where it began), so the record
                // would leave the graph open there. It closes at a pave
                // instead, both ends welded onto it.
                let own = BTreeSet::from([r]);
                if let Some(p) = self.common_pave(res, &own, &[g, ge]) {
                    out[r] = self.move_to(res, p, r);
                }
                continue;
            }
            if !done.insert(g) {
                continue;
            }
            // The records with an end in the group, each closed.
            let group = &res.groups[g];
            let owners: BTreeSet<usize> = group
                .members
                .iter()
                .flat_map(|&m| &self.cands[m].incidence)
                .filter_map(|inc| match *inc {
                    Incidence::Record {
                        record,
                        end: Some(_),
                        ..
                    } => Some(record),
                    _ => None,
                })
                .collect();
            let free = owners.iter().all(|&o| closed(o))
                && group.members.iter().all(|&m| {
                    let c = &self.cands[m];
                    c.source == SourceKind::RecordEndpoint
                        && c.incidence.iter().all(|inc| {
                            matches!(
                                inc,
                                Incidence::Record { end: Some(_), .. }
                                    | Incidence::Edge { .. }
                                    | Incidence::Face { .. }
                            )
                        })
                });
            if !free {
                continue;
            }
            let Some(p) = self.common_pave(res, &owners, &[g]) else {
                continue;
            };
            for &o in &owners {
                if let Some(slot) = out.get_mut(o) {
                    *slot = self.move_to(res, p, o);
                }
            }
        }
        out
    }

    /// The lowest pave group every one of `owners` holds, other than
    /// the groups in `exclude`.
    fn common_pave(
        &self,
        res: &Resolution,
        owners: &BTreeSet<usize>,
        exclude: &[usize],
    ) -> Option<usize> {
        let mut held: Option<BTreeSet<usize>> = None;
        for &o in owners {
            let mine: BTreeSet<usize> = self
                .cands
                .iter()
                .enumerate()
                .filter(|(_, c)| {
                    c.incidence.iter().any(|inc| {
                        matches!(*inc, Incidence::Record { record, end: None, .. } if record == o)
                    })
                })
                .map(|(i, _)| res.group_of[i])
                .filter(|p| !exclude.contains(p))
                .collect();
            held = Some(match held {
                None => mine,
                Some(h) => h.intersection(&mine).copied().collect(),
            });
        }
        held?.into_iter().next()
    }

    /// Record `o` closing at pave group `p`: the parameter of its first
    /// candidate there, and the group's anchor.
    fn move_to(&self, res: &Resolution, p: usize, o: usize) -> Option<ClosureMove> {
        let s = self.cands.iter().enumerate().find_map(|(i, c)| {
            (res.group_of[i] == p)
                .then(|| {
                    c.incidence.iter().find_map(|inc| match *inc {
                        Incidence::Record {
                            record,
                            s,
                            end: None,
                        } if record == o => Some(s),
                        _ => None,
                    })
                })
                .flatten()
        })?;
        Some(ClosureMove {
            s,
            at: self.cands[res.groups[p].anchor].position,
        })
    }
}
