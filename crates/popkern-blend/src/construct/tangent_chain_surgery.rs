//! Closed tangent-chain fillets (PK-733): a loop of filleted edges
//! joined end to end on one tangent line at every junction — a
//! rounded-rectangle rim, a cylinder tee's saddle split at seam
//! crossings — built by direct topology construction.
//!
//! **Why not the per-edge pipeline.** Its anchors, face trims and
//! corner phase are built for edge ENDS: a junction there read as a
//! corner, the anchor walked `r` along the tangent neighbour (off the
//! true section) and the face splits wired chords between distant
//! junctions. A closed tangent chain has no ends: its topology is
//! fully determined, so it is constructed whole.
//!
//! **Construction** (edges `e₁…e_k` in chain order, `h_i` the
//! half-edge of `e_i` running along the chain, `P_i` its face and `Q_i`
//! the twin's; junction `v_i` between `e_i` and `e_{i+1}`):
//!
//! - every junction has valence 3: its third edge `t_i` lies on ONE
//!   side (a crease between `P_i`/`P_{i+1}`, or a seam of that side's
//!   face); the other side is the same face across the junction;
//! - the ball at `v_i` touches each side at its section contact (the
//!   strip's trim point nearest `v_i`, agreed by both strips); `t_i`
//!   must pass through its side's contact — `v_i` moves there and
//!   `t_i` shortens; the other side gets a new vertex;
//! - each support loop has `h_i` replaced in place by the trim edge
//!   between consecutive contacts (exact arcs when the trim is a
//!   circle); each strip is a quad: trim on `P`, section arc, trim on
//!   `Q`, previous section arc; the section arcs are the ball's
//!   circle in the plane normal to the chain at the junction, shared
//!   by the two strips that meet there.
//!
//! **Invariants.** Support faces, third edges and junction vertices
//! keep their handles; each strip carries its edge's fillet surface;
//! topology-reference and manifold checks run here, the pipeline's
//! curve-surface gate after. **Fallback.** Built on a clone: any
//! decline leaves the input untouched with a typed [`BlendFailure`].
//! **Tolerance.** Contact seeds agree to max(100 × linear, r / 4);
//! third edges meet their contact within 100 × linear; every position
//! is an exact point (solved section, snapped onto the third edge).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceOnFace;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, Edge, Face, FaceHandle, HalfEdgeHandle, Loop, LoopHandle, Vertex,
};

use super::EdgeGeometry;
use super::tangent_chain_cap::{Cap, ChainEnd, prepare_cap, wire_cap};
use super::tangent_chain_junction::strip_rail;
use crate::analysis::Convexity;
use crate::error::BlendFailure;

use super::tangent_chain_geometry::{
    ball_center, blank, chain_tangent, chain_tangent_at_origin, decline, insert_edge, orient_strip,
    section_arc, shorten_third, substitute, trim_piece,
};

/// Order the chain `members` (indices into `edge_geometries`) and
/// build it. On success the strips are returned and `brep` holds the
/// filleted chain; on failure `brep` is untouched.
///
/// # Errors
///
/// [`BlendFailure::InvalidTopology`] naming why the chain is outside
/// the class or could not be built exactly.
pub(crate) fn apply_tangent_chain(
    brep: &mut BRep,
    edge_geometries: &[EdgeGeometry],
    members: &[usize],
    tolerance: &Tolerance,
) -> Result<Vec<FaceHandle>, BlendFailure> {
    let mut work = brep.clone();
    let first = edge_geometries[members[0]].analysis.edge;
    let h0 = work.edges[first].half_edge;
    let strips = build(&mut work, edge_geometries, members, tolerance, h0, None)?;
    *brep = work;
    Ok(strips)
}

/// Build every closed tangent chain (`apply_tangent_chain`), recording
/// its members as handled and its strips and supports as protected; a
/// declined chain's failure joins `failures`.
pub(crate) fn apply_closed_chains(
    result: &mut BRep,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
    chain_members: &mut std::collections::HashSet<usize>,
    protected: &mut Vec<FaceHandle>,
    failures: &mut Vec<(popkern_core::topology::EdgeHandle, BlendFailure)>,
) {
    for chain in super::g1_junction::closed_tangent_chains(result, edge_geometries) {
        match apply_tangent_chain(result, edge_geometries, &chain, tolerance) {
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

/// Build an OPEN tangent chain: `h0` is the first member's half-edge
/// leaving the chain's start vertex, `ends` its start and end kinds.
/// The junctions are built as in a closed chain. At a VANISHING end
/// (PK-764, `docs/architecture/fillet-vanishing-ends.md`) the strip
/// runs out to the end vertex itself (both rails end there, no section
/// arc); at a TANGENT CAP (PK-735,
/// `docs/architecture/fillet-open-tangent-chains.md`) it stops on its
/// section, closed by a planar cap face.
///
/// # Errors
///
/// As [`apply_tangent_chain`].
pub(crate) fn apply_open_tangent_chain(
    brep: &mut BRep,
    edge_geometries: &[EdgeGeometry],
    members: &[usize],
    h0: HalfEdgeHandle,
    ends: (ChainEnd, ChainEnd),
    tolerance: &Tolerance,
) -> Result<Vec<FaceHandle>, BlendFailure> {
    let mut work = brep.clone();
    let strips = build(
        &mut work,
        edge_geometries,
        members,
        tolerance,
        h0,
        Some(ends),
    )?;
    *brep = work;
    Ok(strips)
}

fn build(
    brep: &mut BRep,
    edge_geometries: &[EdgeGeometry],
    members: &[usize],
    tolerance: &Tolerance,
    h0: HalfEdgeHandle,
    ends: Option<(ChainEnd, ChainEnd)>,
) -> Result<Vec<FaceHandle>, BlendFailure> {
    let closed = ends.is_none();
    let at = brep
        .vertices
        .get(brep.half_edges[h0].origin)
        .map(|v| v.position)
        .unwrap_or(Point3::new(0.0, 0.0, 0.0));
    let links = super::tangent_chain_junction::chain_links(
        brep,
        edge_geometries,
        members,
        h0,
        closed,
        at,
        tolerance,
    )?;
    let k = links.len();
    let band = tolerance.linear().max(1e-9) * 100.0;
    // Junction i joins link i to link i + 1 (cyclically when closed).
    let nj = if closed { k } else { k - 1 };
    // An open chain's end vertices: the strips' tips.
    let start_vertex = brep.half_edges[links[0].h].origin;
    let end_vertex = brep.half_edges[brep.half_edges[links[k - 1].h].twin].origin;

    // Junction analysis: v_i = dest(h_i), third edge, section contacts.
    let mut thirds = Vec::with_capacity(nj);
    let mut contacts: Vec<(Point3, Point3)> = Vec::with_capacity(nj);
    for i in 0..nj {
        let (third, contact) = super::tangent_chain_junction::junction_section(
            brep,
            &links[i],
            &links[(i + 1) % k],
            band,
        )?;
        thirds.push(third);
        contacts.push(contact);
    }

    if popkern_core::trace::trace_enabled("tangent-chain") {
        let surfaces: Vec<_> = links.iter().map(|l| &l.eg.surface).collect();
        super::tangent_chain_trace::trace_sections(&surfaces, &contacts);
    }
    // Vertices: the third edge's side reuses v_i (moved), the other
    // side is new.
    let mut vp_of = Vec::with_capacity(nj);
    let mut vq_of = Vec::with_capacity(nj);
    for i in 0..nj {
        let v = brep.half_edges[brep.half_edges[links[i].h].twin].origin;
        let (cp, cq) = contacts[i];
        let fresh = |brep: &mut BRep, position: Point3| {
            brep.vertices.insert(Vertex {
                position,
                outgoing: HalfEdgeHandle::default(),
                name: None,
                tolerance: None,
            })
        };
        if thirds[i].on_p {
            brep.vertices[v].position = cp;
            vp_of.push(v);
            vq_of.push(fresh(brep, cq));
        } else {
            brep.vertices[v].position = cq;
            vq_of.push(v);
            vp_of.push(fresh(brep, cp));
        }
        if let Some(edge) = thirds[i].edge {
            shorten_third(brep, edge, v);
        }
    }

    // Edges and half-edges.
    let shell = brep.faces[brep.loops[brep.half_edges[links[0].h].loop_].face].shell;
    // Tangent caps at the open ends.
    let cap_at = |brep: &mut BRep, kind: Option<ChainEnd>, at_start: bool| {
        let link = if at_start { &links[0] } else { &links[k - 1] };
        match kind {
            Some(ChainEnd::TangentCap) => {
                prepare_cap(brep, link, at_start, shell, band, tolerance).map(Some)
            }
            Some(ChainEnd::EndFace) => super::tangent_chain_face_end::prepare_face_end(
                brep, link, at_start, band, tolerance,
            )
            .map(Some),
            _ => Ok(None),
        }
    };
    let start_cap: Option<Cap> = cap_at(brep, ends.map(|e| e.0), true)?;
    let end_cap: Option<Cap> = cap_at(brep, ends.map(|e| e.1), false)?;
    let mut x_edges = Vec::with_capacity(nj);
    let mut centers = Vec::with_capacity(nj);
    for i in 0..nj {
        let (cp, cq) = contacts[i];
        let tangent = chain_tangent(brep, links[i].h)
            .ok_or_else(|| decline("junction tangent undefined".into(), cp))?;
        let radius = links[i].eg.radius_at_end;
        let support_p = brep.faces[brep.loops[brep.half_edges[links[i].h].loop_].face]
            .surface
            .as_ref()
            .map(|s| s.geometry.clone())
            .ok_or_else(|| decline("support without a surface".into(), cp))?;
        let center = ball_center(&support_p, cp, cq, tangent, radius).ok_or_else(|| {
            decline(
                format!("no ball of radius {radius} touches both contacts"),
                cp,
            )
        })?;
        centers.push(center);
        let curve = section_arc(cp, cq, center, tangent, radius).ok_or_else(|| {
            decline(
                format!("no section arc of radius {radius} spans the contacts"),
                cp,
            )
        })?;
        x_edges.push(brep.insert_edge_tracked(Edge {
            half_edge: HalfEdgeHandle::default(),
            curve: Some(curve),
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        }));
    }
    let mut strips = Vec::with_capacity(k);
    let mut supports = Vec::with_capacity(k);
    // Half-edges of X_i: (in strip i: Q→P, in strip i+1: P→Q).
    let mut x_hes: Vec<(HalfEdgeHandle, HalfEdgeHandle)> = Vec::with_capacity(nj);
    for i in 0..nj {
        let a = brep.half_edges.insert(blank(vq_of[i], x_edges[i]));
        let b = brep.half_edges.insert(blank(vp_of[i], x_edges[i]));
        brep.half_edges[a].twin = b;
        brep.half_edges[b].twin = a;
        brep.edges[x_edges[i]].half_edge = b; // curve runs P → Q
        x_hes.push((a, b));
    }
    for (i, link) in links.iter().enumerate() {
        // The section before link i and after it: absent at an open
        // chain's ends, where the strip runs out to the end vertex.
        let before = if closed || i > 0 {
            Some((i + k - 1) % k)
        } else {
            None
        };
        let after = if closed || i + 1 < k { Some(i) } else { None };
        let h = link.h;
        let tw = brep.half_edges[h].twin;
        // A cap end stops on the cap's section; a vanishing end runs
        // out to the end vertex.
        let start_end = start_cap.filter(|_| i == 0);
        let end_end = end_cap.filter(|_| i + 1 == k);
        let (pa, qa, x_before) = match (before, start_end) {
            (Some(j), _) => (vp_of[j], vq_of[j], Some(x_hes[j].1)),
            (None, Some(cap)) => (cap.vp, cap.vq, Some(cap.strip_x())),
            (None, None) => (start_vertex, start_vertex, None),
        };
        let (pb, qb, x_after) = match (after, end_end) {
            (Some(j), _) => (vp_of[j], vq_of[j], Some(x_hes[j].0)),
            (None, Some(cap)) => (cap.vp, cap.vq, Some(cap.strip_x())),
            (None, None) => (end_vertex, end_vertex, None),
        };
        // The arcs are oriented at their START point, so take the chain
        // direction there (at the end of a quarter arc it is
        // perpendicular and picked the 270° complement — measured: three
        // corner cylinders at 3× their area).
        let dir = chain_tangent_at_origin(brep, h).unwrap_or(Vec3::new(1.0, 0.0, 0.0));
        let (a_p, b_p) = (brep.vertices[pa].position, brep.vertices[pb].position);
        let (a_q, b_q) = (brep.vertices[qa].position, brep.vertices[qb].position);
        let tp = strip_rail(&link.eg.surface, a_p, b_p, band)
            .unwrap_or_else(|| trim_piece(link.p_trim, a_p, b_p, dir));
        let tq = strip_rail(&link.eg.surface, a_q, b_q, band)
            .unwrap_or_else(|| trim_piece(link.q_trim, a_q, b_q, dir));
        let e_p = insert_edge(brep, tp);
        let e_q = insert_edge(brep, tq);
        // Support half-edges replace h (P) and its twin (Q) in place.
        let sp = brep.half_edges.insert(blank(pa, e_p));
        substitute(brep, h, sp);
        let sq = brep.half_edges.insert(blank(qb, e_q));
        substitute(brep, tw, sq);
        // The strip.
        let face = brep.faces.insert(Face {
            outer_loop: LoopHandle::default(),
            inner_loops: vec![],
            shell,
            surface: Some(SurfaceOnFace {
                geometry: link.eg.surface.clone(),
                orientation: true,
            }),
            name: None,
            tolerance: None,
            role: Default::default(),
        });
        let lp = brep.loops.insert(Loop {
            half_edge: HalfEdgeHandle::default(),
            face,
        });
        brep.faces[face].outer_loop = lp;
        let fp = brep.half_edges.insert(blank(pb, e_p));
        let fq = brep.half_edges.insert(blank(qa, e_q));
        brep.half_edges[sp].twin = fp;
        brep.half_edges[fp].twin = sp;
        brep.half_edges[sq].twin = fq;
        brep.half_edges[fq].twin = sq;
        brep.edges[e_p].half_edge = sp; // P side runs pa → pb with the curve
        brep.edges[e_q].half_edge = fq; // strip side runs qa → qb with the curve
        // Strip loop: P trim back, section X_prev (P→Q), Q trim, X_i
        // (Q→P); a vanishing end has no section — its trims meet there.
        let mut ring = vec![fp];
        ring.extend(x_before);
        ring.push(fq);
        ring.extend(x_after);
        for (j, &he) in ring.iter().enumerate() {
            let nx = ring[(j + 1) % ring.len()];
            brep.half_edges[he].next = nx;
            brep.half_edges[nx].prev = he;
            brep.half_edges[he].loop_ = lp;
        }
        brep.loops[lp].half_edge = fp;
        strips.push(face);
        supports.push((sp, sq));
        // The old chain edge goes.
        brep.remove_half_edge_reanchored(h);
        brep.remove_half_edge_reanchored(tw);
        brep.edges.remove(link.eg.analysis.edge);
    }
    for (cap, i) in [(start_cap, 0), (end_cap, k - 1)] {
        if let Some(cap) = cap {
            wire_cap(brep, &cap, supports[i].0, supports[i].1)?;
        }
    }
    let caps: Vec<Cap> = start_cap.into_iter().chain(end_cap).collect();
    // Vertices whose outgoing half-edge may have gone: junctions,
    // vanishing tips, cap ends and contacts.
    let mut tips = Vec::new();
    if let Some((first, last)) = ends {
        if first == ChainEnd::Vanishing {
            tips.push(start_vertex);
        }
        if last == ChainEnd::Vanishing {
            tips.push(end_vertex);
        }
    }
    for c in &caps {
        if !c.end_face {
            tips.push(c.v);
        }
        tips.extend([c.vp, c.vq]);
    }
    for v in vp_of.iter().chain(vq_of.iter()).chain(tips.iter()) {
        if (!brep.half_edges.contains_key(brep.vertices[*v].outgoing)
            || brep.half_edges[brep.vertices[*v].outgoing].origin != *v)
            && let Some((he, _)) = brep.half_edges.iter().find(|(_, he)| he.origin == *v)
        {
            brep.vertices[*v].outgoing = he;
        }
    }
    // Orientation first (pcurve recompute reads the face sense), then
    // pcurves for every face the chain touched, recomputed per face so
    // each loop stays chart-coherent on periodic supports (individually
    // authored arc pcurves landed a period off their loop-mates and the
    // corner cylinders read unorientable — measured in OCC).
    for (i, &face) in strips.iter().enumerate() {
        let convex = matches!(links[i].eg.analysis.convexity, Convexity::Convex);
        // Oriented from a section on the strip's own boundary: its
        // following junction, else the one before it, else its cap.
        let cap_section = |cap: Option<Cap>| cap.map(|c| (c.center, brep.half_edges[c.x_pq].edge));
        let section = if i < nj {
            Some((centers[i], x_edges[i]))
        } else if i > 0 {
            Some((centers[i - 1], x_edges[i - 1]))
        } else {
            None
        }
        .or_else(|| cap_section(start_cap.filter(|_| i == 0)))
        .or_else(|| cap_section(end_cap.filter(|_| i + 1 == k)));
        let Some((center, x)) = section else {
            return Err(decline("a one-edge open chain has no section".into(), at));
        };
        let Some(mid) = brep.edges[x]
            .curve
            .as_ref()
            .map(|c| c.geometry.point_at(0.5 * (c.t_start + c.t_end)))
        else {
            return Err(decline("a strip section has no curve".into(), at));
        };
        orient_strip(brep, face, convex, center, mid);
    }
    let cap_edges: Vec<_> = caps.iter().flat_map(|c| c.edges(brep)).collect();
    let mut faces: Vec<FaceHandle> = strips.clone();
    faces.extend(caps.iter().map(|c| c.face));
    for he in brep.half_edges.values() {
        let f = brep.loops[he.loop_].face;
        let new_edge = strips.contains(&f)
            || thirds.iter().any(|t| t.edge == Some(he.edge))
            || cap_edges.contains(&he.edge);
        if new_edge && !faces.contains(&f) {
            faces.push(f);
        }
    }
    for link in &links {
        for f in [link.eg.analysis.face_1, link.eg.analysis.face_2] {
            if !faces.contains(&f) {
                faces.push(f);
            }
        }
    }
    for face in faces {
        brep.recompute_face_pcurves(face)
            .map_err(|e| decline(format!("pcurves of {face:?}: {e}"), at))?;
    }
    for cap in &caps {
        super::tangent_chain_cap::check_cap_winding(brep, cap)?;
    }
    if popkern_core::trace::trace_enabled("tangent-chain") {
        super::tangent_chain_trace::trace_edge_ends(brep, &strips);
    }
    brep.validate_topology_references()
        .and_then(|()| brep.validate_manifold_edges())
        .map_err(|e| decline(format!("built chain fails validation: {e}"), at))?;
    Ok(strips)
}
