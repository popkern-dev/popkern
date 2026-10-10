//! Tangent-chain junctions (PK-733): two filleted edges meeting at a
//! vertex with one tangent line there. Their strips meet on the
//! junction's cross-section — the ball at the vertex is the same ball
//! for both — so the junction is a continuation, never a corner.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

/// Two unit tangents closer than this (sine) are one tangent line.
/// Fitted intersection edges carry ~1e-6 of tangent noise at their
/// split points (measured: an r40 × r20 tee saddle split at a seam
/// crossing, sine 1.2e-6); a genuine corner is orders above 1e-4, and
/// a near-tangent one is still caught by the chain construction's
/// section-agreement check.
const TANGENT_SINE: f64 = 1e-4;

/// The unit tangent of `edge` at its endpoint `vertex` (along the
/// stored curve, or the chord for an edge without one).
pub(crate) fn edge_tangent(brep: &BRep, edge: EdgeHandle, vertex: VertexHandle) -> Option<Vec3> {
    let data = brep.edges.get(edge)?;
    let he = brep.half_edges.get(data.half_edge)?;
    let twin = brep.half_edges.get(he.twin)?;
    let (from, to) = (he.origin, twin.origin);
    if from == to {
        return None; // a closed edge has no junction end
    }
    let tangent = match data.curve.as_ref() {
        Some(curve) => {
            let t = if vertex == from {
                curve.t_start
            } else {
                curve.t_end
            };
            curve.geometry.tangent_at(t)
        }
        None => brep.vertices.get(to)?.position - brep.vertices.get(from)?.position,
    };
    tangent.normalized()
}

/// Whether filleted edges `a` and `b` share one tangent line at
/// `vertex` (either sense: one runs into the vertex, one out).
pub(crate) fn tangent_at_vertex(
    brep: &BRep,
    a: EdgeHandle,
    b: EdgeHandle,
    vertex: VertexHandle,
) -> bool {
    match (edge_tangent(brep, a, vertex), edge_tangent(brep, b, vertex)) {
        (Some(ta), Some(tb)) => ta.cross(&tb).length() < TANGENT_SINE,
        _ => false,
    }
}

/// The tangent-chain junctions among the fillet edges `ends`
/// (`(edge, start, end)` per edge): vertices where exactly two distinct
/// fillet edges meet on one tangent line. They are continuations, so
/// the corner machinery must not see them as corners.
pub(crate) fn junction_vertices(
    brep: &BRep,
    ends: impl Iterator<Item = (EdgeHandle, VertexHandle, VertexHandle)>,
) -> std::collections::BTreeSet<VertexHandle> {
    let mut by_vertex: std::collections::BTreeMap<VertexHandle, Vec<EdgeHandle>> =
        std::collections::BTreeMap::new();
    for (edge, start, end) in ends {
        by_vertex.entry(start).or_default().push(edge);
        if end != start {
            by_vertex.entry(end).or_default().push(edge);
        }
    }
    by_vertex
        .into_iter()
        .filter(|(v, edges)| {
            edges.len() == 2
                && edges[0] != edges[1]
                && tangent_at_vertex(brep, edges[0], edges[1], *v)
        })
        .map(|(v, _)| v)
        .collect()
}

/// The CLOSED tangent chains among `edge_geometries`: groups of fillet
/// edges whose every endpoint is a tangent-chain junction, joined into
/// cycles (each junction shared by exactly two members). Indices in
/// `edge_geometries`, one group per chain, ascending.
pub(crate) fn closed_tangent_chains(
    brep: &BRep,
    edge_geometries: &[super::EdgeGeometry],
) -> Vec<Vec<usize>> {
    let junctions = junction_vertices(
        brep,
        edge_geometries.iter().map(|eg| {
            (
                eg.analysis.edge,
                eg.analysis.vertex_start,
                eg.analysis.vertex_end,
            )
        }),
    );
    let open: Vec<bool> = edge_geometries
        .iter()
        .map(|eg| {
            eg.analysis.vertex_start != eg.analysis.vertex_end
                && junctions.contains(&eg.analysis.vertex_start)
                && junctions.contains(&eg.analysis.vertex_end)
        })
        .collect();
    let mut seen = vec![false; edge_geometries.len()];
    let mut chains = Vec::new();
    for start in 0..edge_geometries.len() {
        if seen[start] || !open[start] {
            continue;
        }
        let mut group = vec![start];
        seen[start] = true;
        let mut grew = true;
        while grew {
            grew = false;
            for j in 0..edge_geometries.len() {
                if seen[j] || !open[j] {
                    continue;
                }
                let (a, b) = (
                    edge_geometries[j].analysis.vertex_start,
                    edge_geometries[j].analysis.vertex_end,
                );
                let touches = group.iter().any(|&g| {
                    let (c, d) = (
                        edge_geometries[g].analysis.vertex_start,
                        edge_geometries[g].analysis.vertex_end,
                    );
                    a == c || a == d || b == c || b == d
                });
                if touches {
                    seen[j] = true;
                    group.push(j);
                    grew = true;
                }
            }
        }
        group.sort_unstable();
        // Closed: every member endpoint is shared by exactly two members.
        let closed = group.iter().all(|&g| {
            [
                edge_geometries[g].analysis.vertex_start,
                edge_geometries[g].analysis.vertex_end,
            ]
            .iter()
            .all(|v| {
                group
                    .iter()
                    .filter(|&&o| {
                        edge_geometries[o].analysis.vertex_start == *v
                            || edge_geometries[o].analysis.vertex_end == *v
                    })
                    .count()
                    == 2
            })
        });
        if closed && group.len() >= 2 {
            chains.push(group);
        }
    }
    chains
}

/// At a tangent-chain junction `vertex` between fillet edge `fillet`
/// and its filleted neighbour `other`, the anchor is the section
/// contact on `trim`'s face — not a corner walk `r` along the tangent
/// neighbour, which lands off the section entirely (measured: (27, −15)
/// for the section (27, −12) on a rounded-rectangle rim). It is the
/// EXACT section's contact (the ball of `radius` in the plane normal to
/// the chain at the junction, touching both supports), so both strips
/// meeting there land on one point: each strip's own trim stops short
/// or overshoots differently (measured, PK-770: (25.580, 15.493) and
/// (25.493, 15.580) for the exact (25.536, 15.536) at a split corner
/// arc). The trim's nearest point seeds the solve and is the fallback.
/// `None` when the neighbour does not continue `fillet` tangentially.
pub(crate) fn junction_section_contact(
    brep: &BRep,
    fillet: EdgeHandle,
    other: EdgeHandle,
    vertex: VertexHandle,
    trim: &popkern_core::geometry::NurbsCurve,
    radius: f64,
) -> Option<popkern_core::geometry::Point3> {
    if !tangent_at_vertex(brep, fillet, other, vertex) {
        return None;
    }
    let p = brep.vertices.get(vertex)?.position;
    let s = trim.closest_point(&p)?;
    let seed = trim.point_at(s);
    Some(exact_junction_contact(brep, fillet, vertex, seed, radius).unwrap_or(seed))
}

/// The exact section contact nearest `seed` at `vertex` (see
/// [`junction_section_contact`]).
fn exact_junction_contact(
    brep: &BRep,
    fillet: EdgeHandle,
    vertex: VertexHandle,
    seed: popkern_core::geometry::Point3,
    radius: f64,
) -> Option<popkern_core::geometry::Point3> {
    use popkern_core::geometry::surface::SurfaceEval;
    let p = brep.vertices.get(vertex)?.position;
    let he = brep.half_edges.get(brep.edges.get(fillet)?.half_edge)?;
    let surface = |h: popkern_core::topology::HalfEdgeHandle| {
        let l = brep.half_edges.get(h)?.loop_;
        Some(
            brep.faces
                .get(brep.loops.get(l)?.face)?
                .surface
                .as_ref()?
                .geometry
                .clone(),
        )
    };
    let (sp, sq) = (surface(brep.edges[fillet].half_edge)?, surface(he.twin)?);
    let n = edge_tangent(brep, fillet, vertex)?;
    // Ball-centre guess: from the seed, `radius` along the seed's
    // support normal, on the side whose distance to the other
    // support is nearest `radius` (convex and concave alike).
    let on = |s: &popkern_core::geometry::surface::SurfaceGeometry| {
        s.closest_point(&seed)
            .map(|(u, v)| s.point_at(u, v).distance_to(&seed))
            .unwrap_or(f64::INFINITY)
    };
    let (own, far) = if on(&sp) <= on(&sq) {
        (&sp, &sq)
    } else {
        (&sq, &sp)
    };
    let (u, v) = own.closest_point(&seed)?;
    let normal = own.normal_at(u, v).normalized()?;
    let miss = |c: popkern_core::geometry::Point3| {
        far.closest_point(&c)
            .map(|(u, v)| (far.point_at(u, v).distance_to(&c) - radius).abs())
            .unwrap_or(f64::INFINITY)
    };
    let guess = [seed + normal * radius, seed + normal * (-radius)]
        .into_iter()
        .min_by(|a, b| miss(*a).total_cmp(&miss(*b)))?;
    let (_, cp, cq) = super::tangent_chain_geometry::exact_section(&sp, &sq, p, n, radius, guess)?;
    let contact = if cp.distance_to(&seed) <= cq.distance_to(&seed) {
        cp
    } else {
        cq
    };
    // A gross jump is a different ball, not fit slack: keep the seed.
    (contact.distance_to(&seed) <= 0.25 * radius).then_some(contact)
}

/// Build every closed tangent chain (`apply_tangent_chain`), recording
/// its members as handled and its strips and supports as protected; a
/// declined chain's failure joins `failures`.
pub(crate) fn apply_closed_chains(
    result: &mut BRep,
    edge_geometries: &[super::EdgeGeometry],
    tolerance: &popkern_core::tolerance::Tolerance,
    chain_members: &mut std::collections::HashSet<usize>,
    protected: &mut Vec<popkern_core::topology::FaceHandle>,
    failures: &mut Vec<(
        popkern_core::topology::EdgeHandle,
        crate::error::BlendFailure,
    )>,
) {
    for chain in super::g1_junction::closed_tangent_chains(result, edge_geometries) {
        match super::tangent_chain_surgery::apply_tangent_chain(
            result,
            edge_geometries,
            &chain,
            tolerance,
        ) {
            Ok(strips) => {
                if popkern_core::trace::trace_enabled("fillet") {
                    eprintln!(
                        "  [tangent-chain] built {} strip(s) for {chain:?}",
                        strips.len()
                    );
                }
                protected.extend(strips);
                for &i in &chain {
                    protected.push(edge_geometries[i].analysis.face_1);
                    protected.push(edge_geometries[i].analysis.face_2);
                }
            }
            Err(failure) => {
                if popkern_core::trace::trace_enabled("fillet") {
                    eprintln!("  [tangent-chain] declined {chain:?}: {failure:?}");
                }
                failures.push((edge_geometries[chain[0]].analysis.edge, failure));
            }
        }
        chain_members.extend(chain);
    }
}
