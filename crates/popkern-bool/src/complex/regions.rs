//! Shells and region membership (PK-974).
//!
//! **Shells.** Around each edge, the wedge between two radially
//! consecutive fins is one region's slice of space: the side of the
//! first fin's face that faces the wedge and the side of the second
//! fin's face that faces it bound the same region. Joining face sides
//! through every wedge of every edge gives the shells — connected
//! boundary components of the regions (Weiler 1986, radial edge
//! traversal). Shells meet only through edges, so two regions touching
//! at a vertex stay separate.
//!
//! **Membership.** A face owned by operand k fixes bit k on both its
//! sides: the side behind k's outward normal is in k, the side ahead is
//! not. Crossing a face k does NOT own leaves bit k unchanged, so bit k
//! propagates between the shells on either side of such a face. A set
//! of shells joined that way and touching no face of k (a region inside
//! k's interior with no contact, say) takes bit k from ONE point
//! classification — OCCT's non-BRep interference — instead of a ray per
//! fragment. Two determinations that disagree are a
//! `MembershipConflict`: subdivision, orientation or a phantom fragment
//! is wrong, and the stage declines.

use super::ids::{FaceId, ShellId};
use super::model::{CellComplex, ComplexShell, SideIndex};
use super::operand_set::OperandSet;
use super::radial::RadialSides;

/// Why shell or membership assembly declined.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum RegionDecline {
    /// Both sides of a face joined into one shell: the face is not a
    /// separating wall (a lamina or an orientation fault).
    SidesJoined {
        /// The face.
        face: FaceId,
    },
    /// Two faces adjacent to one region disagree on an operand.
    MembershipConflict {
        /// The shell.
        shell: ShellId,
        /// The operand.
        operand: usize,
        /// The face side that first fixed the bit, and the face side
        /// that contradicts it.
        faces: [(FaceId, SideIndex); 2],
    },
    /// A shell's membership in an operand could not be decided (no
    /// owned face and the point classification failed).
    MembershipUndetermined {
        /// The shell.
        shell: ShellId,
        /// The operand.
        operand: usize,
    },
}

/// Disjoint-set forest.
pub(crate) struct Forest(Vec<usize>);

impl Forest {
    /// `n` singleton classes.
    pub(crate) fn new(n: usize) -> Self {
        Self((0..n).collect())
    }
    /// The representative of `x`'s class.
    pub(crate) fn find(&mut self, mut x: usize) -> usize {
        while self.0[x] != x {
            self.0[x] = self.0[self.0[x]];
            x = self.0[x];
        }
        x
    }
    /// Join the classes of `a` and `b`.
    pub(crate) fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            // Lower root wins: deterministic class representatives.
            let (lo, hi) = if ra < rb { (ra, rb) } else { (rb, ra) };
            self.0[hi] = lo;
        }
    }
}

fn side_slot(face: FaceId, side: SideIndex) -> usize {
    face.index() * 2 + side as usize
}

/// Assemble shells from the radial order. Fills `c.shells` and every
/// face side's `shell`.
pub(crate) fn assemble_shells(
    c: &mut CellComplex,
    sides: &RadialSides,
) -> Result<(), RegionDecline> {
    let mut forest = Forest::new(c.faces.len() * 2);
    let trace = popkern_core::trace::trace_enabled("complex-radial");
    for (ei, e) in c.edges.iter().enumerate() {
        let m = e.radial.len();
        // Round a sheet's free edge one side of the sheet meets the
        // other (stage 9).
        if m == 1 && c.is_free_lamina_edge(super::ids::EdgeId::from_index(ei)) {
            let face = c.fins[e.radial[0].index()].face;
            forest.union(
                side_slot(face, SideIndex::Ahead),
                side_slot(face, SideIndex::Behind),
            );
            continue;
        }
        if m < 2 {
            continue;
        }
        if trace {
            let order: Vec<(usize, bool)> = e
                .radial
                .iter()
                .map(|f| {
                    (
                        c.fins[f.index()].face.index(),
                        sides.ahead_faces_next[f.index()],
                    )
                })
                .collect();
            eprintln!("[complex-radial] edge {ei}: (face, ahead_faces_next) {order:?}");
        }
        for i in 0..m {
            let a = e.radial[i];
            let b = e.radial[(i + 1) % m];
            // a's side facing the wedge after a; b's side facing the
            // wedge before b.
            let a_side = if sides.ahead_faces_next[a.index()] {
                SideIndex::Ahead
            } else {
                SideIndex::Behind
            };
            let b_side = if sides.ahead_faces_next[b.index()] {
                SideIndex::Behind
            } else {
                SideIndex::Ahead
            };
            let fa = c.fins[a.index()].face;
            let fb = c.fins[b.index()].face;
            forest.union(side_slot(fa, a_side), side_slot(fb, b_side));
        }
    }
    let mut shell_of_root = std::collections::BTreeMap::new();
    c.shells.clear();
    for fi in 0..c.faces.len() {
        let face = FaceId::from_index(fi);
        // A lamina's sides join where the sheet ends inside a region.
        if !c.is_lamina(face)
            && forest.find(side_slot(face, SideIndex::Ahead))
                == forest.find(side_slot(face, SideIndex::Behind))
        {
            return Err(RegionDecline::SidesJoined { face });
        }
        for side in [SideIndex::Ahead, SideIndex::Behind] {
            let root = forest.find(side_slot(face, side));
            let shell = *shell_of_root.entry(root).or_insert_with(|| {
                c.shells.push(ComplexShell::default());
                ShellId::from_index(c.shells.len() - 1)
            });
            c.shells[shell.index()].sides.push((face, side));
            c.faces[fi].sides[side as usize].shell = Some(shell);
        }
    }
    Ok(())
}

/// Bit `k` a face side fixes directly, if the face is owned by `k`.
fn owned_bit(c: &CellComplex, face: FaceId, side: SideIndex, k: usize) -> Option<bool> {
    let owner = c.faces[face.index()]
        .owners
        .iter()
        .find(|o| o.operand == k)?;
    // Behind an owner's outward normal is inside it.
    let behind_owner = (side == SideIndex::Behind) == owner.same_sense;
    Some(behind_owner)
}

/// Assign every shell's membership. `classify(shell, operand)` answers
/// for a shell set no face of the operand touches: `Some(true)` inside,
/// `Some(false)` outside, `None` undecidable. Returns how many such
/// classifications ran.
pub(crate) fn assign_membership(
    c: &mut CellComplex,
    mut classify: impl FnMut(&CellComplex, ShellId, usize) -> Option<bool>,
) -> Result<usize, RegionDecline> {
    let ns = c.shells.len();
    let mut bits: Vec<OperandSet> = vec![OperandSet::new(); ns];
    let mut classified = 0;
    for k in 0..c.operand_count {
        // A sheet bounds no volume: no region lies inside it.
        if c.lamina_operands.contains(k) {
            continue;
        }
        // Join shells across faces k does not own.
        let mut forest = Forest::new(ns);
        for face in &c.faces {
            if face.owners.iter().any(|o| o.operand == k) {
                continue;
            }
            if let (Some(a), Some(b)) = (face.sides[0].shell, face.sides[1].shell) {
                forest.union(a.index(), b.index());
            }
        }
        // Direct determinations per class.
        let mut class_bit: std::collections::BTreeMap<usize, (bool, Option<(FaceId, SideIndex)>)> =
            std::collections::BTreeMap::new();
        for (si, shell) in c.shells.iter().enumerate() {
            for &(face, side) in &shell.sides {
                let Some(bit) = owned_bit(c, face, side, k) else {
                    continue;
                };
                let root = forest.find(si);
                match class_bit.get(&root) {
                    Some(&(prev, first)) if prev != bit => {
                        return Err(RegionDecline::MembershipConflict {
                            shell: ShellId::from_index(si),
                            operand: k,
                            faces: [first.unwrap_or((face, side)), (face, side)],
                        });
                    }
                    Some(_) => {}
                    None => {
                        class_bit.insert(root, (bit, Some((face, side))));
                    }
                }
            }
        }
        for (si, set) in bits.iter_mut().enumerate() {
            let root = forest.find(si);
            let bit = match class_bit.get(&root) {
                Some(&(bit, _)) => bit,
                None => {
                    classified += 1;
                    let bit = classify(c, ShellId::from_index(root), k).ok_or(
                        RegionDecline::MembershipUndetermined {
                            shell: ShellId::from_index(root),
                            operand: k,
                        },
                    )?;
                    class_bit.insert(root, (bit, None));
                    bit
                }
            };
            if bit {
                set.insert(k);
            }
        }
    }
    for (shell, b) in c.shells.iter_mut().zip(bits) {
        shell.membership = Some(b);
    }
    Ok(classified)
}
