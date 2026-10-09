//! Vanishing strip ends (PK-764, `docs/architecture/fillet-vanishing-ends.md`):
//! where the supports turn tangent — an equal-radius tee's singular
//! points — a fillet strip runs out to zero width, both rails ending on
//! the edge's end vertex. Such strips are built as whole open tangent
//! chains, their tips meeting at the point with no corner blend.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{BRep, VertexHandle};

use super::EdgeGeometry;

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

/// The OPEN tangent chains among `edge_geometries` whose two ends
/// vanish: groups of fillet edges joined end to end on one tangent
/// line at every interior junction (PK-733's junction class), whose
/// two free ends are vanishing ends — an equal-radius tee's saddle
/// half-ellipse S− → seam crossing → S+. Each chain with the half-edge
/// leaving its start vertex along the first member.
pub(super) fn vanishing_tangent_chains(
    brep: &BRep,
    edge_geometries: &[EdgeGeometry],
) -> Vec<(Vec<usize>, popkern_core::topology::HalfEdgeHandle)> {
    let junctions = super::g1_junction::junction_vertices(
        brep,
        edge_geometries.iter().map(|g| {
            (
                g.analysis.edge,
                g.analysis.vertex_start,
                g.analysis.vertex_end,
            )
        }),
    );
    let ends = |g: &EdgeGeometry| [g.analysis.vertex_start, g.analysis.vertex_end];
    // Members: open edges each of whose ends is a junction or vanishes.
    let member: Vec<bool> = edge_geometries
        .iter()
        .map(|g| {
            g.analysis.vertex_start != g.analysis.vertex_end
                && ends(g)
                    .iter()
                    .all(|&v| junctions.contains(&v) || vanishes_at(brep, g, v))
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
        // An open path: exactly two free (non-junction) ends, both
        // vanishing, on different members or one member's two ends.
        let free: Vec<(usize, VertexHandle)> = group
            .iter()
            .flat_map(|&g| ends(&edge_geometries[g]).into_iter().map(move |v| (g, v)))
            .filter(|(_, v)| !junctions.contains(v))
            .collect();
        if free.len() != 2
            || !free
                .iter()
                .all(|&(g, v)| vanishes_at(brep, &edge_geometries[g], v))
        {
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
        chains.push((group, h0));
    }
    chains
}

/// Build every open tangent chain whose ends vanish (`apply_open_tangent_chain`),
/// recording its members as handled and its strips and supports as
/// protected; a declined chain's failure joins `failures`.
pub(super) fn apply_vanishing_chains(
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
    for (chain, h0) in vanishing_tangent_chains(result, edge_geometries) {
        if chain.iter().any(|i| chain_members.contains(i)) {
            continue;
        }
        match super::tangent_chain_surgery::apply_open_tangent_chain(
            result,
            edge_geometries,
            &chain,
            h0,
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
                    eprintln!("  [tangent-chain] open chain declined {chain:?}: {failure:?}");
                }
                failures.push((edge_geometries[chain[0]].analysis.edge, failure));
            }
        }
        chain_members.extend(chain);
    }
}
