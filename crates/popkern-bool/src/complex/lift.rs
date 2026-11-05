//! Lifting the chart walk's fragments into the cell complex (PK-972).
//!
//! The intersection graph already shares vertex and edge identity
//! between the two operands: a record (intersection curve) is ONE
//! graph half-edge pair used by up to four fragments, and operand
//! edges are pairs used by the fragments of their two faces. The lift
//! keeps that identity — one complex edge per graph edge key
//! (`min(half-edge, twin)`), one fin per fragment loop entry — and adds
//! what the graph lacks:
//!
//! - **Same-domain faces.** A coincident patch is walked once on each
//!   operand's face; the two fragments bound the same edge set. They
//!   become ONE complex face owned by both operands, so the coincident
//!   pair can never form a zero-thickness wedge (Gursoz 1990 §4.6:
//!   overlapping faces detected by connectivity and merged).
//! - **Material sense.** Fragment loops wind counter-clockwise about
//!   their source face's folded normal (the assembler keeps the
//!   source orientation flag for an unflipped fragment); whether that
//!   normal points out of the material comes from the operand's
//!   material sense, decided once per closed component.

use std::collections::{BTreeMap, HashMap};

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_core::validation::MaterialSense;

use super::ids::{EdgeId, FaceId, VertexId};
use super::model::{CellComplex, FaceOwner};
use super::spans::unify_spans;
use crate::merge::graph::{FaceRef, IntersectionGraph, Operand};
use crate::merge::walk::Fragment;

/// Why the lift declined.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum LiftDecline {
    /// A fragment's source face has no surface.
    MissingSurface {
        /// The fragment.
        fragment: usize,
    },
    /// An operand component's material sense could not be decided.
    UnknownMaterialSense {
        /// Operand index.
        operand: usize,
        /// A face of the component.
        face: FaceHandle,
    },
    /// A fragment's outer loop is empty after removing degenerate
    /// (zero-length, curveless) edges.
    EmptyOuterLoop {
        /// The fragment.
        fragment: usize,
    },
}

/// A lifted complex and its correspondence with the fragments.
#[derive(Debug, Clone)]
pub(crate) struct Lift {
    /// The complex (faces, loops, fins, edges, vertices; no shells).
    pub complex: CellComplex,
    /// Operand A's folded-normal convention (folded normals point out
    /// of its material): the result keeps it.
    pub output_folded_outward: bool,
}

/// Operand index of a graph operand tag.
pub(crate) fn operand_index(op: Operand) -> usize {
    op.index()
}

/// The source BRep of a face reference.
pub(crate) fn operand_brep<'a>(operands: &[&'a BRep], face: FaceRef) -> &'a BRep {
    operands[operand_index(face.operand)]
}

/// Whether a graph edge is degenerate residue: coincident endpoint
/// nodes and no curve with extent (the assembler's rule — a genuine
/// closed rim carries its circle).
pub(crate) fn is_degenerate_key(graph: &IntersectionGraph, key: usize) -> bool {
    use popkern_core::geometry::curve::CurveEval;
    let twin = graph.half_edges[key].twin;
    let na = graph.half_edges[key].origin;
    let nb = graph.half_edges[twin].origin;
    let (pa, pb) = (graph.nodes[na].position, graph.nodes[nb].position);
    if pa.distance_to(&pb) >= 1.0e-9 {
        return false;
    }
    let curve = graph.half_edges[key]
        .curve
        .as_ref()
        .or(graph.half_edges[twin].curve.as_ref());
    match curve {
        None => true,
        Some(c) => {
            let pm = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
            pm.distance_to(&pa) < 1.0e-9
        }
    }
}

/// Closed components of `brep` (faces reachable over shared edges),
/// each sorted, in order of their lowest face.
fn components(brep: &BRep) -> Vec<Vec<FaceHandle>> {
    let mut seen = std::collections::HashSet::new();
    let mut faces: Vec<FaceHandle> = brep.faces.keys().collect();
    faces.sort();
    let mut out = Vec::new();
    for seed in faces {
        if seen.contains(&seed) {
            continue;
        }
        let mut component = Vec::new();
        let mut stack = vec![seed];
        while let Some(f) = stack.pop() {
            if !seen.insert(f) {
                continue;
            }
            component.push(f);
            let Some(face) = brep.faces.get(f) else {
                continue;
            };
            for &l in std::iter::once(&face.outer_loop).chain(face.inner_loops.iter()) {
                let Some(start) = brep.loops.get(l).map(|l| l.half_edge) else {
                    continue;
                };
                let mut he = start;
                for _ in 0..brep.half_edges.len() {
                    let Some(h) = brep.half_edges.get(he) else {
                        break;
                    };
                    if let Some(tw) = h.mate().and_then(|m| brep.half_edges.get(m))
                        && let Some(tl) = brep.loops.get(tw.loop_)
                        && !seen.contains(&tl.face)
                    {
                        stack.push(tl.face);
                    }
                    he = h.next;
                    if he == start {
                        break;
                    }
                }
            }
        }
        component.sort();
        out.push(component);
    }
    out
}

/// Material sense of every face of `brep`: whether its folded normal
/// points out of the material, decided once per closed component.
///
/// A component is decided by point classifications just off its faces
/// (`probe_votes`). The sign of the folded-normal volume integral
/// (`component_material_sense`) is the fallback: it samples a trim
/// classification on a grid over every face (11.8 of the 12 ms a
/// 48-gon prism union spent in the complex), and it is the WRONG test
/// for a cavity shell — correctly oriented, its folded normals point
/// into the cavity and its volume integral is negative. A
/// multi-component operand takes its first decisive probe and falls
/// back to the volume sign only where no face votes. A
/// single-component operand takes the probes only when two faces' votes
/// agree; otherwise the volume sign decides.
fn material_senses(
    brep: &BRep,
    operand: usize,
    boxes: &[(
        FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )],
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<HashMap<FaceHandle, bool>, LiftDecline> {
    let comps = components(brep);
    let mut out = HashMap::new();
    for component in &comps {
        let Some(&seed) = component.first() else {
            continue;
        };
        let by_volume = || {
            brep.component_material_sense(seed)
                .map(|s| matches!(s, MaterialSense::FoldedOutward))
        };
        let votes = probe_votes(brep, component, boxes, tolerance);
        let outward = match (votes.as_slice(), comps.len()) {
            ([a, b], _) if a == b => Some(*a),
            ([a, ..], n) if n > 1 => Some(*a),
            _ => None,
        }
        .or_else(by_volume)
        .ok_or(LiftDecline::UnknownMaterialSense {
            operand,
            face: seed,
        })?;
        for &f in component {
            out.insert(f, outward);
        }
    }
    Ok(out)
}

/// Whether the folded normals of `component` point out of `brep`'s
/// material, read off its first faces: at most two votes. A face votes
/// only when it separates material at its sample point — the point just
/// off it along its folded normal and the point just off it against
/// classify one Inside, one Outside. Across a wall thinner than the
/// offset both read the far side alike, and that face does not vote
/// (a one-sided probe read a thin plate's caps as folding outward when
/// they fold in).
fn probe_votes(
    brep: &BRep,
    component: &[FaceHandle],
    boxes: &[(
        FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )],
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Vec<bool> {
    let mut votes = Vec::with_capacity(2);
    let scale = boxes
        .iter()
        .map(|(_, b)| (b.max - b.min).length())
        .fold(0.0_f64, f64::max);
    let eps = (scale * 1.0e-5).max(tolerance.linear() * 50.0);
    for &fh in component.iter().take(6) {
        let Some(surf) = brep.faces.get(fh).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let sample = crate::classify::pick_sample_point(brep, fh);
        let Some((u, v)) = surf.geometry.closest_point(&sample) else {
            continue;
        };
        let Some(n) = surf.folded_normal(u, v) else {
            continue;
        };
        let at = surf.geometry.point_at(u, v);
        let side = |p| crate::classify::point_in_solid_final(&p, brep, boxes, tolerance).ok();
        use crate::error::PointClassification::{Inside, Outside};
        match (side(at + n * eps), side(at - n * eps)) {
            (Some(Outside), Some(Inside)) => votes.push(true),
            (Some(Inside), Some(Outside)) => votes.push(false),
            _ => continue,
        }
        if votes.len() == 2 {
            break;
        }
    }
    votes
}

/// Outward normals of two coincident faces agree at the fragment's
/// first boundary node.
fn same_outward(
    operands: &[&BRep],
    senses: &[HashMap<FaceHandle, bool>],
    a: FaceRef,
    b: FaceRef,
    at: popkern_core::geometry::Point3,
) -> Option<bool> {
    let outward = |f: FaceRef| {
        let face = operand_brep(operands, f).faces.get(f.face)?;
        let surf = face.surface.as_ref()?;
        let (u, v) = surf.geometry.closest_point(&at)?;
        let n = surf.folded_normal(u, v)?;
        let out = *senses[operand_index(f.operand)].get(&f.face)?;
        Some(if out { n } else { n * -1.0 })
    };
    Some(outward(a)?.dot(&outward(b)?) > 0.0)
}

/// Lift `fragments` of `graph` into a complex. `overlaps` are the
/// coincident face pairs (operand A face, operand B face); only
/// fragments of such pairs are candidates for same-domain merging.
pub(crate) fn lift_fragments(
    graph: &IntersectionGraph,
    fragments: &[Fragment],
    operands: &[&BRep],
    boxes: &[&[(
        FaceHandle,
        popkern_core::geometry::intersection::surface_subdivision::Aabb,
    )]],
    overlaps: &[crate::merge::graph::FacePair],
    tol: &popkern_core::tolerance::Tolerance,
) -> Result<Lift, LiftDecline> {
    let tolerance = tol.linear();
    let reps = unify_spans(graph, fragments, tolerance);
    // A sheet has no material: its faces' sense is a convention the
    // complex never reads as material (lamina faces join no volume).
    let is_sheet = |brep: &BRep| brep.body_kind() == popkern_core::topology::BodyKind::Sheet;
    let senses = operands
        .iter()
        .zip(boxes)
        .enumerate()
        .map(|(k, (brep, boxes))| {
            if is_sheet(brep) {
                Ok(brep.faces.keys().map(|f| (f, true)).collect())
            } else {
                material_senses(brep, k, boxes, tol)
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    let output_folded_outward = operands[0]
        .faces
        .keys()
        .min()
        .and_then(|f| senses[0].get(&f).copied())
        .unwrap_or(true);

    // Fragments filling their own face's holes bound no material: they
    // neither become complex faces nor own a partner's (the selection
    // never keeps them).
    let hole_fillers =
        super::hole_phantoms::hole_filling_fragments(graph, fragments, operands, &reps);
    let mut groups = super::same_domain::same_domain_groups(graph, fragments, &reps, overlaps);
    for g in &mut groups {
        g.retain(|&k| !hole_fillers[k]);
    }
    groups.retain(|g| !g.is_empty());

    let mut complex = CellComplex::new(operands.len());
    for (k, brep) in operands.iter().enumerate() {
        if is_sheet(brep) {
            complex.lamina_operands.insert(k);
        }
    }
    let mut vertex_of_node: BTreeMap<usize, VertexId> = BTreeMap::new();
    let mut edge_of_key: BTreeMap<usize, EdgeId> = BTreeMap::new();
    let mut face_of_fragment = vec![FaceId::from_index(usize::MAX); fragments.len()];

    for group in &groups {
        let primary = group[0];
        let frag = &fragments[primary];
        let pref = frag.face;
        let pi = operand_index(pref.operand);
        let source = operand_brep(operands, pref);
        if source
            .faces
            .get(pref.face)
            .and_then(|f| f.surface.as_ref())
            .is_none()
        {
            return Err(LiftDecline::MissingSurface { fragment: primary });
        }
        let folded_is_outward =
            *senses[pi]
                .get(&pref.face)
                .ok_or(LiftDecline::UnknownMaterialSense {
                    operand: pi,
                    face: pref.face,
                })?;
        let anchor = frag
            .outer_loop
            .first()
            .map(|&h| graph.nodes[graph.half_edges[h].origin].position);
        let mut owners = vec![FaceOwner {
            operand: pi,
            face: pref.face,
            same_sense: true,
        }];
        for &other in &group[1..] {
            let oref = fragments[other].face;
            let same = anchor
                .and_then(|p| same_outward(operands, &senses, pref, oref, p))
                .unwrap_or(true);
            owners.push(FaceOwner {
                operand: operand_index(oref.operand),
                face: oref.face,
                same_sense: same,
            });
        }
        let face = complex.add_face(owners, folded_is_outward, group.clone());
        for &k in group {
            face_of_fragment[k] = face;
        }
        for (li, lp) in std::iter::once(&frag.outer_loop)
            .chain(frag.inner_loops.iter())
            .enumerate()
        {
            let mut uses = Vec::with_capacity(lp.len());
            for &h in lp {
                let raw = h.min(graph.half_edges[h].twin);
                let Some(&rep) = reps.get(&raw) else {
                    continue;
                };
                let key = rep.key;
                let twin = graph.half_edges[key].twin;
                let edge = match edge_of_key.get(&key) {
                    Some(&e) => e,
                    None => {
                        let mut vertex = |node: usize, c: &mut CellComplex| {
                            *vertex_of_node.entry(node).or_insert_with(|| {
                                c.add_vertex(graph.nodes[node].position, tolerance, Some(node))
                            })
                        };
                        let start = vertex(graph.half_edges[key].origin, &mut complex);
                        let end = vertex(
                            graph.half_edges[graph.half_edges[key].twin].origin,
                            &mut complex,
                        );
                        let curve = graph.half_edges[key]
                            .curve
                            .clone()
                            .or_else(|| graph.half_edges[twin].curve.clone());
                        let e = complex.add_edge(start, end, curve, tolerance, Some(key));
                        edge_of_key.insert(key, e);
                        e
                    }
                };
                uses.push((edge, (h == raw) == rep.same_dir, Some(h)));
            }
            if uses.is_empty() {
                if li == 0 {
                    return Err(LiftDecline::EmptyOuterLoop { fragment: primary });
                }
                continue;
            }
            complex.add_loop(face, &uses, li == 0);
        }
    }

    Ok(Lift {
        complex,
        output_folded_outward,
    })
}

#[cfg(test)]
#[path = "lift_tests.rs"]
mod tests;
