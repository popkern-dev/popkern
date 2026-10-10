//! Open tangent chains (PK-764, PK-735): runs of fillet edges joined
//! end to end on one tangent line, built whole by the chain surgery
//! (`tangent_chain_surgery`). Each of the chain's two free ends must
//! be of a kind the surgery builds (`ChainEnd`):
//! - VANISHING (`docs/architecture/fillet-vanishing-ends.md`): where
//!   the supports turn tangent — an equal-radius tee's singular points
//!   — the strip runs out to zero width, both rails ending on the end
//!   vertex;
//! - TANGENT CAP (`docs/architecture/fillet-open-tangent-chains.md`):
//!   the rim continues tangent and unfilleted past the end vertex, and
//!   the strip stops on its section under a planar cap.
//!
//! Chains with any other end stay with the per-edge pipeline.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{BRep, VertexHandle};

use super::EdgeGeometry;
use super::tangent_chain_cap::ChainEnd;

/// Whether `rail` ends at `p` (either parameter end within `gate`).
fn rail_ends_at(rail: &NurbsCurve, p: Point3, gate: f64) -> bool {
    let (a, b) = rail.parameter_range();
    [rail.point_at(a), rail.point_at(b)]
        .into_iter()
        .any(|q| q.distance_to(&p) <= gate)
}

/// The gate a vanishing end is measured with: ten times the model's
/// linear tolerance (contacts converge onto the tangent point exactly;
/// the band covers the fitted rail's end).
fn vanishing_gate(brep: &BRep) -> f64 {
    brep.tolerance.linear().max(1e-9) * 10.0
}

/// Whether `geometry`'s strip vanishes at `vertex`: both rails end on it.
fn vanishes_at(brep: &BRep, geometry: &EdgeGeometry, vertex: VertexHandle) -> bool {
    let Some(p) = brep.vertices.get(vertex).map(|v| v.position) else {
        return false;
    };
    let gate = vanishing_gate(brep);
    rail_ends_at(&geometry.trim_3d_1, p, gate) && rail_ends_at(&geometry.trim_3d_2, p, gate)
}

/// Whether the rim continues TANGENT past `vertex` on an unfilleted
/// edge: another edge there, not among `fillet_edges`, shares the
/// filleted edge's tangent line (valence 2 or 3; the surgery checks the
/// rest of the cap's preconditions).
fn continues_tangent(
    brep: &BRep,
    edge: popkern_core::topology::EdgeHandle,
    vertex: VertexHandle,
    fillet_edges: &[popkern_core::topology::EdgeHandle],
) -> bool {
    let out: Vec<_> = brep
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    (2..=3).contains(&out.len())
        && out.iter().any(|&other| {
            other != edge
                && !fillet_edges.contains(&other)
                && super::g1_junction::tangent_at_vertex(brep, edge, other, vertex)
        })
}

/// Whether `vertex` is a sharp corner the strip runs out on: valence
/// 3, the other two edges unfilleted and on one PLANAR end face.
fn runs_out_on_end_face(
    brep: &BRep,
    edge: popkern_core::topology::EdgeHandle,
    vertex: VertexHandle,
    fillet_edges: &[popkern_core::topology::EdgeHandle],
) -> bool {
    let out: Vec<_> = brep
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .collect();
    let others: Vec<_> = out.iter().filter(|he| he.edge != edge).collect();
    if out.len() != 3
        || others.len() != 2
        || others.iter().any(|he| fillet_edges.contains(&he.edge))
    {
        return false;
    }
    // The face between the two others: the one both border.
    let faces = |he: &popkern_core::topology::HalfEdge| {
        [
            brep.loops[he.loop_].face,
            brep.loops[brep.half_edges[he.twin].loop_].face,
        ]
    };
    let (f0, f1) = (faces(others[0]), faces(others[1]));
    f0.iter().filter(|f| f1.contains(f)).any(|&f| {
        matches!(
            brep.faces[f].surface.as_ref().map(|s| &s.geometry),
            Some(popkern_core::geometry::surface::SurfaceGeometry::Plane { .. })
        ) && !super::g1_junction::tangent_at_vertex(brep, edge, others[0].edge, vertex)
            && !super::g1_junction::tangent_at_vertex(brep, edge, others[1].edge, vertex)
    })
}

/// How `geometry`'s strip ends at its free end `vertex`, if the chain
/// surgery builds that end.
fn end_kind(
    brep: &BRep,
    geometry: &EdgeGeometry,
    vertex: VertexHandle,
    fillet_edges: &[popkern_core::topology::EdgeHandle],
) -> Option<ChainEnd> {
    if vanishes_at(brep, geometry, vertex) {
        Some(ChainEnd::Vanishing)
    } else if continues_tangent(brep, geometry.analysis.edge, vertex, fillet_edges) {
        Some(ChainEnd::TangentCap)
    } else if runs_out_on_end_face(brep, geometry.analysis.edge, vertex, fillet_edges) {
        Some(ChainEnd::EndFace)
    } else {
        None
    }
}

/// Mitre corners among the fillet edges: valence-3 vertices where
/// exactly two of them meet NOT tangent, the third edge unfilleted,
/// and neither strip vanishes there (an equal-radius tee's singular
/// point joins two saddle chains at an angle on a valence-3 vertex: a
/// vanishing end, not a mitre).
fn mitre_corners(
    brep: &BRep,
    edge_geometries: &[EdgeGeometry],
) -> std::collections::BTreeSet<VertexHandle> {
    let mut by_vertex: std::collections::BTreeMap<VertexHandle, Vec<&EdgeGeometry>> =
        Default::default();
    for g in edge_geometries {
        if g.analysis.vertex_start == g.analysis.vertex_end {
            continue;
        }
        for v in [g.analysis.vertex_start, g.analysis.vertex_end] {
            by_vertex.entry(v).or_default().push(g);
        }
    }
    by_vertex
        .into_iter()
        .filter(|(v, gs)| {
            gs.len() == 2
                && brep
                    .half_edges
                    .values()
                    .filter(|he| he.origin == *v)
                    .count()
                    == 3
                && !gs.iter().any(|g| vanishes_at(brep, g, *v))
                && !super::g1_junction::tangent_at_vertex(
                    brep,
                    gs[0].analysis.edge,
                    gs[1].analysis.edge,
                    *v,
                )
        })
        .map(|(v, _)| v)
        .collect()
}

/// One open tangent chain: its members (indices into the edge
/// geometries, ascending), the half-edge leaving its start vertex along
/// its first member, and its start and end kinds.
type OpenChain = (
    Vec<usize>,
    popkern_core::topology::HalfEdgeHandle,
    (ChainEnd, ChainEnd),
);

/// The OPEN tangent chains among `edge_geometries`: groups of fillet
/// edges joined end to end on one tangent line at every interior
/// junction (PK-733's junction class), whose two free ends are both of
/// a kind the surgery builds — an equal-radius tee's saddle
/// half-ellipse S− → seam crossing → S+ (vanishing), part of a rounded
/// rim (tangent caps), or one edge alone.
pub(super) fn open_tangent_chains(brep: &BRep, edge_geometries: &[EdgeGeometry]) -> Vec<OpenChain> {
    let tangent_junctions = super::g1_junction::junction_vertices(
        brep,
        edge_geometries.iter().map(|g| {
            (
                g.analysis.edge,
                g.analysis.vertex_start,
                g.analysis.vertex_end,
            )
        }),
    );
    // Links join members: tangent junctions, and mitre corners (two
    // fillet edges meeting at an angle, the third edge sharp, PK-768).
    let mut junctions = tangent_junctions.clone();
    junctions.extend(mitre_corners(brep, edge_geometries));
    let fillet_edges: Vec<_> = edge_geometries.iter().map(|g| g.analysis.edge).collect();
    let ends = |g: &EdgeGeometry| [g.analysis.vertex_start, g.analysis.vertex_end];
    let kind = |g: &EdgeGeometry, v: VertexHandle| end_kind(brep, g, v, &fillet_edges);
    // Members: open edges each of whose ends is a junction or a free
    // end the surgery builds.
    let member: Vec<bool> = edge_geometries
        .iter()
        .map(|g| {
            g.analysis.vertex_start != g.analysis.vertex_end
                && ends(g)
                    .iter()
                    .all(|&v| junctions.contains(&v) || kind(g, v).is_some())
                && ends(g).iter().any(|&v| !junctions.contains(&v))
        })
        .collect();
    let mut seen = vec![false; edge_geometries.len()];
    let mut chains = Vec::new();
    for start in 0..edge_geometries.len() {
        if seen[start] || !member[start] {
            continue;
        }
        // Grow through junctions only (a vanishing vertex where two
        // chains meet at an angle is not a junction).
        let mut group = vec![start];
        seen[start] = true;
        let mut grew = true;
        while grew {
            grew = false;
            for j in 0..edge_geometries.len() {
                if seen[j]
                    || !(member[j]
                        || ends(&edge_geometries[j])
                            .iter()
                            .all(|v| junctions.contains(v)))
                {
                    continue;
                }
                let joined = group.iter().any(|&g| {
                    ends(&edge_geometries[g])
                        .iter()
                        .any(|v| junctions.contains(v) && ends(&edge_geometries[j]).contains(v))
                });
                if joined {
                    seen[j] = true;
                    group.push(j);
                    grew = true;
                }
            }
        }
        // An open path: exactly two free (non-junction) ends, each of a
        // kind the surgery builds, on different members or one
        // member's two ends.
        let free: Vec<(usize, VertexHandle)> = group
            .iter()
            .flat_map(|&g| ends(&edge_geometries[g]).into_iter().map(move |v| (g, v)))
            .filter(|(_, v)| !junctions.contains(v))
            .collect();
        let kinds: Vec<Option<ChainEnd>> = free
            .iter()
            .map(|&(g, v)| kind(&edge_geometries[g], v))
            .collect();
        let (2, [Some(k0), Some(k1)]) = (free.len(), kinds.as_slice()) else {
            continue;
        };
        // The surgery takes a chain only for what the per-edge pipeline
        // cannot build — a tangent junction, a tangent cap or a
        // vanishing end. Lone edges and mitre pairs between end faces
        // are ordinary fillets: the per-edge corner and terminal
        // machinery owns them.
        let tangent_inside = group.iter().any(|&g| {
            ends(&edge_geometries[g])
                .iter()
                .any(|v| tangent_junctions.contains(v))
        });
        if !tangent_inside && *k0 == ChainEnd::EndFace && *k1 == ChainEnd::EndFace {
            continue;
        }
        let (g0, v0) = free[0];
        let e0 = edge_geometries[g0].analysis.edge;
        let Some(h) = brep.edges.get(e0).map(|e| e.half_edge) else {
            continue;
        };
        let h0 = if brep.half_edges[h].origin == v0 {
            h
        } else {
            brep.half_edges[h].twin
        };
        group.sort_unstable();
        chains.push((group, h0, (*k0, *k1)));
    }
    chains
}

/// Build every open tangent chain (`apply_open_tangent_chain`),
/// recording its members as handled and its strips and supports as
/// protected; a declined chain's failure joins `failures`.
pub(super) fn apply_open_chains(
    result: &mut BRep,
    edge_geometries: &[EdgeGeometry],
    tolerance: &popkern_core::tolerance::Tolerance,
    chain_members: &mut std::collections::HashSet<usize>,
    protected: &mut Vec<popkern_core::topology::FaceHandle>,
    failures: &mut Vec<(
        popkern_core::topology::EdgeHandle,
        crate::error::BlendFailure,
    )>,
) {
    for (chain, h0, kinds) in open_tangent_chains(result, edge_geometries) {
        if chain.iter().any(|i| chain_members.contains(i)) {
            continue;
        }
        match super::tangent_chain_surgery::apply_open_tangent_chain(
            result,
            edge_geometries,
            &chain,
            h0,
            kinds,
            tolerance,
        ) {
            Ok(strips) => {
                protected.extend(strips);
                for &i in &chain {
                    protected.push(edge_geometries[i].analysis.face_1);
                    protected.push(edge_geometries[i].analysis.face_2);
                }
            }
            Err(failure) => {
                if popkern_core::trace::trace_enabled("fillet") {
                    eprintln!(
                        "  [tangent-chain] open chain {kinds:?} declined {chain:?}: {failure:?}"
                    );
                }
                failures.push((edge_geometries[chain[0]].analysis.edge, failure));
            }
        }
        chain_members.extend(chain);
    }
}
