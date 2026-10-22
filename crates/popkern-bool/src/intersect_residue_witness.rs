//! Closure witness for trim residues (PK-473): a residue survives the
//! Stage 1 filters only when it can close a gap the ordinary records
//! leave open on its owner face.
//!
//! A trim residue lies on the other operand's surface EXTENSION,
//! beyond that operand's trim, so it is never part of the other
//! operand's boundary: the regions on its two sides classify alike,
//! and in exact arithmetic it separates nothing the ordinary records
//! do not already separate. It carries information only where the
//! ordinary records fail to close: an ordinary record that ends in the
//! face interior with no other ordinary record and no face boundary at
//! that end (a graze the intersector missed) lets the fragment walk
//! leak one region into the next, and a residue continuing from that
//! end stands in for the missing boundary. So a residue is kept when
//! its CHAIN (residues sharing ends) touches such a dangling ordinary
//! end, and retired otherwise.
//!
//! This replaces the residue-swarm cap (PK-472's containment, which
//! dropped every residue on a face holding more than 16). Measured
//! 2026-10-21: the centrifugal impeller's flow-cavity floor held 835
//! residues against 294 ordinary records and the polyhedral dice's pip
//! faces 53..137 each, every one of them inside an outline the
//! ordinary records already closed; with the cap removed the impeller
//! declined after 285 s and a dice pip cut declined, and a rule keyed
//! on residue ends meeting ordinary-curve interiors broke the dice's
//! small pip faces. With every residue dropped the whole workspace
//! suite stays green, so no fixture today reaches the dangling class;
//! the witness keeps the mechanism for it without the swarm.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval as _;
use popkern_core::topology::FaceHandle;

use super::intersect::IntersectionRecord;

/// The ends of a record's 3-D curve.
fn ends(record: &IntersectionRecord) -> [Point3; 2] {
    let (t0, t1) = record.curve.curve_3d.parameter_range();
    [
        record.curve.curve_3d.point_at(t0),
        record.curve.curve_3d.point_at(t1),
    ]
}

/// One owner face's witness outcome, for diagnostics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ResidueWitness {
    /// `true` when the owner face belongs to operand A.
    pub(crate) owner_is_a: bool,
    /// The owner face.
    pub(crate) face: FaceHandle,
    /// Ordinary record ends on the face with nothing at them.
    pub(crate) dangling_ends: usize,
    /// Residues kept (their chain reaches a dangling end).
    pub(crate) kept: usize,
    /// Residues retired.
    pub(crate) retired: usize,
}

/// Retire every residue whose chain reaches no dangling ordinary end
/// on its owner face. Ends meet within `gate` (the graph's
/// junction-adoption distance, so a meeting here is a shared node
/// there); `on_owner_boundary(owner_is_a, face, point)` says whether a
/// point lies on that face's trim boundary. Returns one row per owner
/// face that held residues, in first-seen record order.
pub(crate) fn retire_unwitnessed_residues(
    records: &mut Vec<IntersectionRecord>,
    gate: f64,
    on_owner_boundary: impl Fn(bool, FaceHandle, &Point3) -> bool,
) -> Vec<ResidueWitness> {
    let owner = |r: &IntersectionRecord| -> Option<(bool, FaceHandle)> {
        r.residue_owner_a
            .map(|a| (a, if a { r.face_a } else { r.face_b }))
    };
    let mut faces: Vec<(bool, FaceHandle)> = Vec::new();
    for r in records.iter() {
        if let Some(key) = owner(r)
            && !faces.contains(&key)
        {
            faces.push(key);
        }
    }
    let record_ends: Vec<[Point3; 2]> = records.iter().map(ends).collect();
    let mut keep = vec![true; records.len()];
    let mut rows = Vec::with_capacity(faces.len());
    for (owner_is_a, face) in faces {
        let on_face = |r: &IntersectionRecord| {
            if owner_is_a {
                r.face_a == face
            } else {
                r.face_b == face
            }
        };
        let ordinary: Vec<usize> = (0..records.len())
            .filter(|&i| records[i].residue_owner_a.is_none() && on_face(&records[i]))
            .collect();
        let residues: Vec<usize> = (0..records.len())
            .filter(|&i| owner(&records[i]) == Some((owner_is_a, face)))
            .collect();
        let dangling = dangling_ends(&ordinary, &record_ends, gate, |p| {
            on_owner_boundary(owner_is_a, face, p)
        });
        let live = witnessed_chains(&residues, &record_ends, &dangling, gate);
        let kept = live.iter().filter(|&&l| l).count();
        for (&i, &l) in residues.iter().zip(&live) {
            keep[i] = l;
        }
        rows.push(ResidueWitness {
            owner_is_a,
            face,
            dangling_ends: dangling.len(),
            kept,
            retired: residues.len() - kept,
        });
    }
    let mut it = keep.into_iter();
    records.retain(|_| it.next().unwrap_or(true));
    rows
}

/// Ends of the `ordinary` records met by no other ordinary end and not
/// on the face boundary. A closed record has no ends.
fn dangling_ends(
    ordinary: &[usize],
    record_ends: &[[Point3; 2]],
    gate: f64,
    on_boundary: impl Fn(&Point3) -> bool,
) -> Vec<Point3> {
    let mut out = Vec::new();
    for &i in ordinary {
        let [p0, p1] = record_ends[i];
        if p0.distance_to(&p1) < gate {
            continue;
        }
        for (k, p) in [p0, p1].iter().enumerate() {
            let met = ordinary.iter().any(|&j| {
                record_ends[j]
                    .iter()
                    .enumerate()
                    .any(|(m, q)| (j, m) != (i, k) && q.distance_to(p) < gate)
            });
            if !met && !on_boundary(p) {
                out.push(*p);
            }
        }
    }
    out
}

/// Per residue (in `residues` order): whether its chain — residues
/// joined end to end within `gate` — touches a dangling end.
fn witnessed_chains(
    residues: &[usize],
    record_ends: &[[Point3; 2]],
    dangling: &[Point3],
    gate: f64,
) -> Vec<bool> {
    let n = residues.len();
    if dangling.is_empty() {
        return vec![false; n];
    }
    let mut parent: Vec<usize> = (0..n).collect();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    for x in 0..n {
        for y in x + 1..n {
            let (ex, ey) = (&record_ends[residues[x]], &record_ends[residues[y]]);
            if ex
                .iter()
                .any(|p| ey.iter().any(|q| p.distance_to(q) < gate))
            {
                let (rx, ry) = (root(&mut parent, x), root(&mut parent, y));
                parent[rx] = ry;
            }
        }
    }
    let mut live_root = vec![false; n];
    for x in 0..n {
        let touches = record_ends[residues[x]]
            .iter()
            .any(|p| dangling.iter().any(|q| q.distance_to(p) < gate));
        if touches {
            let rx = root(&mut parent, x);
            live_root[rx] = true;
        }
    }
    (0..n).map(|x| live_root[root(&mut parent, x)]).collect()
}

#[cfg(test)]
#[path = "intersect_residue_witness_tests.rs"]
mod tests;
