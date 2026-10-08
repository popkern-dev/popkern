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
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, Edge, Face, FaceHandle, HalfEdgeHandle, Loop, LoopHandle, Vertex,
};

use super::EdgeGeometry;
use crate::analysis::Convexity;
use crate::error::BlendFailure;

use super::tangent_chain_geometry::{
    Circle, ball_center, chain_tangent, chain_tangent_at_origin, exact_section, nearest,
    orient_strip, section_arc, segment_distance, shorten_third, strip_rail, trim_piece,
};
use super::tangent_chain_geometry::{
    Link, Third, blank, decline, insert_edge, substitute, verified_circle,
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
    let strips = build(&mut work, edge_geometries, members, tolerance)?;
    *brep = work;
    Ok(strips)
}

fn build(
    brep: &mut BRep,
    edge_geometries: &[EdgeGeometry],
    members: &[usize],
    tolerance: &Tolerance,
) -> Result<Vec<FaceHandle>, BlendFailure> {
    let first = edge_geometries[members[0]].analysis.edge;
    let at = brep
        .vertices
        .get(edge_geometries[members[0]].analysis.vertex_start)
        .map(|v| v.position)
        .unwrap_or(Point3::new(0.0, 0.0, 0.0));
    // Chain order: follow half-edges head to tail.
    let mut links: Vec<Link<'_>> = Vec::with_capacity(members.len());
    let mut h = brep.edges[first].half_edge;
    for _ in 0..members.len() {
        let edge = brep.half_edges[h].edge;
        let eg = members
            .iter()
            .map(|&i| &edge_geometries[i])
            .find(|g| g.analysis.edge == edge)
            .ok_or_else(|| decline(format!("edge {edge:?} left the chain"), at))?;
        let p = brep.loops[brep.half_edges[h].loop_].face;
        let trims = |face: FaceHandle| -> Option<(&NurbsCurve, Circle)> {
            if face == eg.analysis.face_1 {
                Some((&eg.trim_3d_1, eg.trim_circle_1))
            } else if face == eg.analysis.face_2 {
                Some((&eg.trim_3d_2, eg.trim_circle_2))
            } else {
                None
            }
        };
        let q = brep.loops[brep.half_edges[brep.half_edges[h].twin].loop_].face;
        let (Some(p_trim), Some(q_trim)) = (trims(p), trims(q)) else {
            return Err(decline(
                format!("edge {edge:?} faces are not its supports"),
                at,
            ));
        };
        // An exact trim circle is trusted only where it lies on its
        // support (measured: a pocket-floor rim's wall circle carried
        // the cylinder origin's height, z = 0, for contacts at z = 10).
        let support = |face: FaceHandle| {
            brep.faces[face]
                .surface
                .as_ref()
                .map(|s| s.geometry.clone())
        };
        let p_trim = (
            p_trim.0,
            verified_circle(p_trim.1, support(p).as_ref(), tolerance),
        );
        let q_trim = (
            q_trim.0,
            verified_circle(q_trim.1, support(q).as_ref(), tolerance),
        );
        if matches!(eg.analysis.convexity, Convexity::Mixed) {
            return Err(decline(format!("edge {edge:?} changes convexity"), at));
        }
        links.push(Link {
            eg,
            h,
            p_trim,
            q_trim,
        });
        // The next chain edge leaves this half-edge's destination.
        let v = brep.half_edges[brep.half_edges[h].twin].origin;
        let next = members
            .iter()
            .map(|&i| edge_geometries[i].analysis.edge)
            .filter(|e| *e != edge)
            .find_map(|e| {
                let he = brep.edges[e].half_edge;
                let tw = brep.half_edges[he].twin;
                if brep.half_edges[he].origin == v {
                    Some(he)
                } else if brep.half_edges[tw].origin == v {
                    Some(tw)
                } else {
                    None
                }
            })
            .ok_or_else(|| decline("the chain does not close".into(), at))?;
        h = next;
    }
    if h != links[0].h && brep.half_edges[h].twin != links[0].h {
        return Err(decline("the chain does not close on itself".into(), at));
    }
    let k = links.len();
    let band = tolerance.linear().max(1e-9) * 100.0;

    // Junction analysis: v_i = dest(h_i), third edge, section contacts.
    let mut thirds = Vec::with_capacity(k);
    let mut contacts: Vec<(Point3, Point3)> = Vec::with_capacity(k);
    for i in 0..k {
        let (a, b) = (&links[i], &links[(i + 1) % k]);
        let hi = &brep.half_edges[a.h];
        let v = brep.half_edges[hi.twin].origin;
        let vp = brep.vertices[v].position;
        let valence = brep.half_edges.values().filter(|he| he.origin == v).count();
        if valence != 3 {
            return Err(decline(
                format!("junction {v:?} has valence {valence}, not 3"),
                vp,
            ));
        }
        let next_p = hi.next;
        let next_q = brep.half_edges[brep.half_edges[b.h].twin].next;
        let third = if next_p != b.h {
            Third {
                edge: brep.half_edges[next_p].edge,
                on_p: true,
            }
        } else if next_q != hi.twin {
            Third {
                edge: brep.half_edges[next_q].edge,
                on_p: false,
            }
        } else {
            return Err(decline(format!("junction {v:?} has no third edge"), vp));
        };
        // The fitted contacts seed the exact solve; they must at least
        // agree to a fraction of the radius (a gross mismatch is a
        // different strip geometry, not fit slack).
        let agree = |t1: (&NurbsCurve, Circle),
                     t2: (&NurbsCurve, Circle),
                     r: f64|
         -> Result<Point3, BlendFailure> {
            let c1 = nearest(t1, vp);
            let c2 = nearest(t2, vp);
            if c1.distance_to(&c2) > band.max(0.25 * r) {
                return Err(decline(
                    format!(
                        "strips disagree on the section at {v:?} by {:.3e} ({c1:?} vs {c2:?})",
                        c1.distance_to(&c2)
                    ),
                    vp,
                ));
            }
            Ok(Point3::new(
                0.5 * (c1.x + c2.x),
                0.5 * (c1.y + c2.y),
                0.5 * (c1.z + c2.z),
            ))
        };
        let cp0 = agree(a.p_trim, b.p_trim, a.eg.effective_radius)?;
        let cq0 = agree(a.q_trim, b.q_trim, a.eg.effective_radius)?;
        // The exact section, seeded by the fitted contacts' ball.
        let surf = |he: HalfEdgeHandle| {
            brep.faces[brep.loops[brep.half_edges[he].loop_].face]
                .surface
                .as_ref()
                .map(|s| s.geometry.clone())
        };
        let (Some(sp), Some(sq), Some(tangent)) =
            (surf(a.h), surf(hi.twin), chain_tangent(brep, a.h))
        else {
            return Err(decline(
                format!("junction {v:?} lacks supports or a tangent"),
                vp,
            ));
        };
        let r = a.eg.radius_at_end;
        let seed = ball_center(&sp, cp0, cq0, tangent, r).unwrap_or(vp);
        let Some((_, cp, cq)) = exact_section(&sp, &sq, vp, tangent, r, seed) else {
            return Err(decline(
                format!("no ball of radius {r} sits on both supports at {v:?}"),
                vp,
            ));
        };
        // The third edge must reach its side's contact: it is kept and
        // shortened to it.
        let on_third = if third.on_p { cp } else { cq };
        let t_curve = brep.edges[third.edge].curve.clone();
        let off = match &t_curve {
            Some(c) => c
                .geometry
                .closest_point(&on_third)
                .map(|s| c.geometry.point_at(s).distance_to(&on_third))
                .unwrap_or(f64::INFINITY),
            None => segment_distance(brep, third.edge, on_third),
        };
        if off > band {
            return Err(decline(
                format!(
                    "third edge {:?} misses the section contact by {off:.3e} (contacts {cp:?} / {cq:?})",
                    third.edge
                ),
                vp,
            ));
        }
        // The kept third edge carries the vertex: land the contact
        // exactly on it (the solved foot sits within the band, ~1e-6
        // off a seam line — an OpenLoop at validation otherwise).
        let snapped = t_curve
            .as_ref()
            .and_then(|c| {
                c.geometry
                    .closest_point(&on_third)
                    .map(|s| c.geometry.point_at(s))
            })
            .unwrap_or(on_third);
        let (cp, cq) = if third.on_p {
            (snapped, cq)
        } else {
            (cp, snapped)
        };
        thirds.push(third);
        contacts.push((cp, cq));
    }

    if popkern_core::trace::trace_enabled("tangent-chain") {
        let surfaces: Vec<_> = links.iter().map(|l| &l.eg.surface).collect();
        super::tangent_chain_trace::trace_sections(&surfaces, &contacts);
    }
    // Vertices: the third edge's side reuses v_i (moved), the other
    // side is new.
    let mut vp_of = Vec::with_capacity(k);
    let mut vq_of = Vec::with_capacity(k);
    for i in 0..k {
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
        shorten_third(brep, thirds[i].edge, v);
    }

    // Edges and half-edges.
    let shell = brep.faces[brep.loops[brep.half_edges[links[0].h].loop_].face].shell;
    let mut x_edges = Vec::with_capacity(k);
    let mut centers = Vec::with_capacity(k);
    for i in 0..k {
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
    // Half-edges of X_i: (in strip i: Q→P, in strip i+1: P→Q).
    let mut x_hes: Vec<(HalfEdgeHandle, HalfEdgeHandle)> = Vec::with_capacity(k);
    for i in 0..k {
        let a = brep.half_edges.insert(blank(vq_of[i], x_edges[i]));
        let b = brep.half_edges.insert(blank(vp_of[i], x_edges[i]));
        brep.half_edges[a].twin = b;
        brep.half_edges[b].twin = a;
        brep.edges[x_edges[i]].half_edge = b; // curve runs P → Q
        x_hes.push((a, b));
    }
    for i in 0..k {
        let prev = (i + k - 1) % k;
        let link = &links[i];
        let h = link.h;
        let tw = brep.half_edges[h].twin;
        let (pa, pb) = (vp_of[prev], vp_of[i]);
        let (qa, qb) = (vq_of[prev], vq_of[i]);
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
        // Strip loop: P trim back, section X_prev (P→Q), Q trim, X_i (Q→P).
        let ring = [fp, x_hes[prev].1, fq, x_hes[i].0];
        for (j, &he) in ring.iter().enumerate() {
            let nx = ring[(j + 1) % 4];
            brep.half_edges[he].next = nx;
            brep.half_edges[nx].prev = he;
            brep.half_edges[he].loop_ = lp;
        }
        brep.loops[lp].half_edge = fp;
        strips.push(face);
        // The old chain edge goes.
        brep.remove_half_edge_reanchored(h);
        brep.remove_half_edge_reanchored(tw);
        brep.edges.remove(link.eg.analysis.edge);
    }
    for v in vp_of.iter().chain(vq_of.iter()) {
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
        let mid = brep.edges[x_edges[i]]
            .curve
            .as_ref()
            .map(|c| c.geometry.point_at(0.5 * (c.t_start + c.t_end)))
            .unwrap_or(contacts[i].0);
        orient_strip(brep, face, convex, centers[i], mid);
    }
    let mut faces: Vec<FaceHandle> = strips.clone();
    for he in brep.half_edges.values() {
        let f = brep.loops[he.loop_].face;
        let new_edge = strips.contains(&f) || thirds.iter().any(|t| t.edge == he.edge);
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
    if popkern_core::trace::trace_enabled("tangent-chain") {
        super::tangent_chain_trace::trace_edge_ends(brep, &strips);
    }
    brep.validate_topology_references()
        .and_then(|()| brep.validate_manifold_edges())
        .map_err(|e| decline(format!("built chain fails validation: {e}"), at))?;
    Ok(strips)
}
