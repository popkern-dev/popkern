//! Structural invariants of the cell complex, reported as typed
//! findings (never a panic, never a bare string).
//!
//! Checked at the end of every construction stage. The stage that
//! built a complex with findings declines; the findings and the
//! complex artifact (`artifact.rs`) go into the diagnostics.

use super::ids::{EdgeId, FaceId, FinId, LoopId, ShellId};
use super::model::{CellComplex, SideIndex};

/// A violated invariant of the complex.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) enum ComplexFinding {
    /// A fin, loop, face, edge or vertex id points outside its arena.
    DanglingId {
        /// What held the bad id.
        what: String,
    },
    /// A loop's fins do not close (the walk from `first` never returns
    /// to it) or a fin's `prev` disagrees with its predecessor.
    OpenLoop {
        /// The loop.
        loop_id: LoopId,
    },
    /// Consecutive fins of a loop do not meet at a vertex.
    LoopGap {
        /// The fin whose end is not its successor's start.
        fin: FinId,
    },
    /// A fin of a loop names another loop or face.
    LoopMembership {
        /// The fin.
        fin: FinId,
    },
    /// A fin is missing from its edge's radial list, listed twice, or
    /// listed on another edge.
    RadialMembership {
        /// The edge.
        edge: EdgeId,
        /// The fin.
        fin: FinId,
    },
    /// A face has no loop, or its first loop is not outer.
    FaceLoops {
        /// The face.
        face: FaceId,
    },
    /// An edge has exactly one fin: a free edge, which no closed
    /// operand produces.
    FreeEdge {
        /// The edge.
        edge: EdgeId,
    },
    /// Shells exist but a face side belongs to none (or to one whose
    /// side list does not contain it).
    UnassignedSide {
        /// The face.
        face: FaceId,
        /// The side.
        side: SideIndex,
    },
    /// Membership is recorded on some shells but not on this one.
    MissingMembership {
        /// The shell.
        shell: ShellId,
    },
}

/// Check every structural invariant; an empty result means the complex
/// is well formed.
pub(crate) fn validate_complex(c: &CellComplex) -> Vec<ComplexFinding> {
    let mut out = Vec::new();
    check_ids(c, &mut out);
    if !out.is_empty() {
        return out;
    }
    check_loops(c, &mut out);
    check_radial(c, &mut out);
    check_faces(c, &mut out);
    check_shells(c, &mut out);
    out
}

fn dangling(out: &mut Vec<ComplexFinding>, what: String) {
    out.push(ComplexFinding::DanglingId { what });
}

fn check_ids(c: &CellComplex, out: &mut Vec<ComplexFinding>) {
    let (nv, ne, nf, nl, nfa) = (
        c.vertices.len(),
        c.edges.len(),
        c.fins.len(),
        c.loops.len(),
        c.faces.len(),
    );
    for (i, e) in c.edges.iter().enumerate() {
        if e.start.index() >= nv || e.end.index() >= nv {
            dangling(out, format!("edge {i} names a missing vertex"));
        }
        if e.radial.iter().any(|f| f.index() >= nf) {
            dangling(out, format!("edge {i} lists a missing fin"));
        }
    }
    for (i, f) in c.fins.iter().enumerate() {
        if f.edge.index() >= ne
            || f.face.index() >= nfa
            || f.loop_.index() >= nl
            || f.next.index() >= nf
            || f.prev.index() >= nf
        {
            dangling(out, format!("fin {i} names a missing entity"));
        }
    }
    for (i, l) in c.loops.iter().enumerate() {
        if l.face.index() >= nfa || l.first.index() >= nf {
            dangling(out, format!("loop {i} names a missing entity"));
        }
    }
    for (i, face) in c.faces.iter().enumerate() {
        if face.loops.iter().any(|l| l.index() >= nl) {
            dangling(out, format!("face {i} lists a missing loop"));
        }
        if face
            .sides
            .iter()
            .any(|s| s.shell.is_some_and(|sh| sh.index() >= c.shells.len()))
        {
            dangling(out, format!("face {i} names a missing shell"));
        }
    }
}

fn check_loops(c: &CellComplex, out: &mut Vec<ComplexFinding>) {
    for (li, l) in c.loops.iter().enumerate() {
        let loop_id = LoopId::from_index(li);
        let fins = c.loop_fins(loop_id);
        let closes = fins
            .last()
            .and_then(|last| c.fins.get(last.index()))
            .is_some_and(|f| f.next == l.first);
        if !closes {
            out.push(ComplexFinding::OpenLoop { loop_id });
            continue;
        }
        for (k, &fin) in fins.iter().enumerate() {
            let f = &c.fins[fin.index()];
            if f.loop_ != loop_id || f.face != l.face {
                out.push(ComplexFinding::LoopMembership { fin });
            }
            let next = fins[(k + 1) % fins.len()];
            if c.fins[next.index()].prev != fin {
                out.push(ComplexFinding::OpenLoop { loop_id });
            }
            if c.fin_end(fin) != c.fin_start(next) {
                out.push(ComplexFinding::LoopGap { fin });
            }
        }
    }
}

fn check_radial(c: &CellComplex, out: &mut Vec<ComplexFinding>) {
    let mut seen = vec![0usize; c.fins.len()];
    for (ei, e) in c.edges.iter().enumerate() {
        let edge = EdgeId::from_index(ei);
        for &fin in &e.radial {
            seen[fin.index()] += 1;
            if c.fins[fin.index()].edge != edge {
                out.push(ComplexFinding::RadialMembership { edge, fin });
            }
        }
        // A sheet's boundary is a free edge by design (its one fin a
        // lamina's); any other single-fin edge is a defect.
        if e.radial.len() == 1 && !e.frontier && !c.is_free_lamina_edge(edge) {
            out.push(ComplexFinding::FreeEdge { edge });
        }
    }
    for (fi, count) in seen.iter().enumerate() {
        if *count != 1 {
            out.push(ComplexFinding::RadialMembership {
                edge: c.fins[fi].edge,
                fin: FinId::from_index(fi),
            });
        }
    }
}

fn check_faces(c: &CellComplex, out: &mut Vec<ComplexFinding>) {
    for (fi, face) in c.faces.iter().enumerate() {
        let outer_first = face.loops.first().is_some_and(|l| c.loops[l.index()].outer);
        if !outer_first {
            out.push(ComplexFinding::FaceLoops {
                face: FaceId::from_index(fi),
            });
        }
    }
}

fn check_shells(c: &CellComplex, out: &mut Vec<ComplexFinding>) {
    if c.shells.is_empty() {
        return;
    }
    for (fi, face) in c.faces.iter().enumerate() {
        let face_id = FaceId::from_index(fi);
        for side in [SideIndex::Ahead, SideIndex::Behind] {
            let listed = face.sides[side as usize].shell.is_some_and(|sh| {
                c.shells[sh.index()]
                    .sides
                    .iter()
                    .any(|&(f, s)| f == face_id && s == side)
            });
            if !listed {
                out.push(ComplexFinding::UnassignedSide {
                    face: face_id,
                    side,
                });
            }
        }
    }
    let any_membership = c.shells.iter().any(|s| s.membership.is_some());
    if any_membership {
        for (si, s) in c.shells.iter().enumerate() {
            if s.membership.is_none() {
                out.push(ComplexFinding::MissingMembership {
                    shell: ShellId::from_index(si),
                });
            }
        }
    }
}
