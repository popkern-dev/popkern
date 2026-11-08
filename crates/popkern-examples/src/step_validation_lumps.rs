//! The bodies a STEP export re-reads as under the round-trip's import
//! options (PK-1115: lumps that touch without sharing are separate).

use std::collections::{HashMap, HashSet};

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, RegionKind, ShellHandle, VertexHandle};
use popkern_io::StepNonManifold;

use super::StepRoundTripValidationOptions;

/// The bodies `breps` re-read as: one each, except that under the
/// default import a body touching itself (a vertex or edge its pieces
/// share) re-reads as its manifold pieces, one body per piece (PK-1115).
pub(super) fn expected_import_bodies(
    breps: &[BRep],
    options: &StepRoundTripValidationOptions,
) -> usize {
    if options.import.non_manifold == StepNonManifold::General {
        return breps.len();
    }
    breps
        .iter()
        .map(|body| {
            if body.contacts().is_empty() {
                separate_lumps(body)
            } else {
                body.split_to_manifold().map_or(1, |pieces| pieces.len())
            }
        })
        .sum()
}

/// The lumps of `body` the import separates: face sets joined through a
/// shared edge or vertex, those in different shells of one region (a
/// solid and its cavity) one lump. Face sets of one shell sharing no
/// entity (parts touching only where their vertices coincide) are
/// separate lumps; the outside joins nothing.
fn separate_lumps(body: &BRep) -> usize {
    let faces: Vec<FaceHandle> = body.faces.keys().collect();
    let index: HashMap<FaceHandle, usize> =
        faces.iter().enumerate().map(|(i, &f)| (f, i)).collect();
    let mut parent: Vec<usize> = (0..faces.len()).collect();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    let mut join = |a: usize, b: usize| {
        let (a, b) = (root(&mut parent, a), root(&mut parent, b));
        parent[a] = b;
    };
    let mut by_vertex: HashMap<VertexHandle, usize> = HashMap::new();
    let mut by_edge: HashMap<EdgeHandle, usize> = HashMap::new();
    let mut uses: Vec<(usize, usize)> = Vec::new();
    for fin in body.half_edges.values() {
        let Some(&i) = body.loops.get(fin.loop_).and_then(|l| index.get(&l.face)) else {
            continue;
        };
        uses.push((i, *by_vertex.entry(fin.origin).or_insert(i)));
        uses.push((i, *by_edge.entry(fin.edge).or_insert(i)));
    }
    for l in body.loops.values() {
        if let (Some(v), Some(&i)) = (l.vertex, index.get(&l.face)) {
            uses.push((i, *by_vertex.entry(v).or_insert(i)));
        }
    }
    for (i, j) in uses {
        join(i, j);
    }
    let mut region_shells: HashMap<_, Vec<(ShellHandle, usize)>> = HashMap::new();
    for (&f, &i) in &index {
        for side in [body.faces[f].back, body.faces[f].front] {
            let r = body.shell_region(side);
            if let Some(r) = r.filter(|&r| body.region_kind(r) == Some(RegionKind::Solid)) {
                region_shells.entry(r).or_default().push((side, i));
            }
        }
    }
    for uses in region_shells.values() {
        if let Some(&(shell, i)) = uses.first() {
            for &(other, j) in uses {
                if other != shell {
                    join(i, j);
                }
            }
        }
    }
    let roots: HashSet<usize> = (0..faces.len()).map(|i| root(&mut parent, i)).collect();
    roots.len().max(1)
}
