//! Splice a planned fillet overflow (`strip_overflow`) into the applied
//! body, with Euler operators only (PK-1189; design
//! `docs/architecture/fillet-overflow-onto-neighbour-strip.md`).
//!
//! The apply phase leaves the support face F with the rolling strip's
//! foot (an inner loop) crossing the convex strip's rail (its outer
//! loop). The splice:
//! - splits the rail and the foot where they cross (P1, P2) and the
//!   rolling strip's wall edge at the overflow's ends (Q1, Q2);
//! - cuts the rolling strip T along its sections P1Q1 and P2Q2, and the
//!   convex strip S along the overflow's contact curve P1P2;
//! - glues foot to rail in F, which leaves the lens between the rail's
//!   middle and the foot's short arc as its own face;
//! - merges the lens, S's cap and T's overflow pieces into one face,
//!   which takes the overflow surface.
//!
//! Its wall edges are re-authored onto the overflow's wall contact. It
//! runs on a copy, committed only when the copy validates and every
//! edge of the new face lies on both its faces.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use super::strip_overflow::OverflowPlan;
use super::strip_overflow_topology::*;

/// Splice `plan` into `brep`; on any failure `brep` is unchanged and the
/// reason is returned.
pub(crate) fn splice_overflow(
    brep: &mut BRep,
    plan: &OverflowPlan,
    tolerance: &Tolerance,
) -> Result<(), String> {
    let mut b = brep.clone();
    let reach = (tolerance.linear() * 10.0).max(1e-9);
    let f = support_face(&b, plan, reach)?;
    let neighbours = face_neighbours(&b, f);
    // The strip faces beside F through P1 (a strip may already be in
    // pieces: the one holding the crossing is the one cut).
    let face_with = |surface: &SurfaceGeometry| -> Result<FaceHandle, String> {
        let found: Vec<FaceHandle> = neighbours
            .iter()
            .copied()
            .filter(|&n| surface_of(&b, n).is_some_and(|s| same_surface(s, surface, reach)))
            .filter(|&n| edge_through(&b, n, f, &plan.p1, reach))
            .collect();
        match found.as_slice() {
            [one] => Ok(*one),
            _ => {
                if popkern_core::trace::trace_enabled("interference") {
                    for face in faces_with_surface(&b, &plan.support_surface) {
                        let kinds: Vec<_> = face_neighbours(&b, face)
                            .into_iter()
                            .map(|n| surface_of(&b, n).map(|s| s.kind_name()))
                            .collect();
                        eprintln!("[overflow] support piece {face:?} (F {f:?}) beside {kinds:?}");
                    }
                    for face in faces_with_surface(&b, &plan.rolling_surface) {
                        let kinds: Vec<_> = face_neighbours(&b, face)
                            .into_iter()
                            .map(|n| (n, surface_of(&b, n).map(|s| s.kind_name())))
                            .collect();
                        eprintln!("[overflow] rolling piece {face:?} beside {kinds:?}");
                    }
                }
                Err(format!(
                    "{} faces beside F through P1 carry the strip's surface",
                    found.len()
                ))
            }
        }
    };
    let t = face_with(&plan.rolling_surface)?;
    let s = face_with(&plan.onto_surface)?;
    let w = plan.face_w;
    // Each surface's pieces (a strip may already be in pieces, and an
    // earlier splice splits F).
    let f_set = faces_with_surface(&b, &plan.support_surface);
    let t_set = faces_with_surface(&b, &plan.rolling_surface);
    let s_set = faces_with_surface(&b, &plan.onto_surface);
    let _ = (t, s);
    let p1_rail = split_at(&mut b, &f_set, &s_set, plan.p1, reach, tolerance)
        .map_err(|e| format!("P1 on the rail: {e}"))?;
    let p2_rail = split_at(&mut b, &f_set, &s_set, plan.p2, reach, tolerance)
        .map_err(|e| format!("P2 on the rail: {e}"))?;
    let p1_foot = split_at(&mut b, &f_set, &t_set, plan.p1, reach, tolerance)
        .map_err(|e| format!("P1 on the foot: {e}"))?;
    let p2_foot = split_at(&mut b, &f_set, &t_set, plan.p2, reach, tolerance)
        .map_err(|e| format!("P2 on the foot: {e}"))?;
    let q1 = split_at(&mut b, &t_set, &[w], plan.q1, reach, tolerance)
        .map_err(|e| format!("Q1 on the wall: {e}"))?;
    let q2 = split_at(&mut b, &t_set, &[w], plan.q2, reach, tolerance)
        .map_err(|e| format!("Q2 on the wall: {e}"))?;
    let rail_middle = edge_between(&b, p1_rail, p2_rail).ok_or("no rail edge between P1 and P2")?;
    // The rolling strip's sections and the contact curve on S.
    let section_1 = make_edge_in(
        &mut b,
        &plan.rolling_surface,
        p1_foot,
        q1,
        &plan.section_start(),
    )?;
    let _section_2 = make_edge_in(
        &mut b,
        &plan.rolling_surface,
        p2_foot,
        q2,
        &plan.section_end(),
    )?;
    let s_faces: Vec<FaceHandle> = faces_with_surface(&b, &plan.onto_surface);
    make_edge_among(&mut b, &s_faces, p1_rail, p2_rail, &plan.curve_on_onto())?;
    // Glue foot to rail in F: a bridge joins the foot's loop to the
    // rail's, an edge at P2 splits off the lens, both contracted.
    let he_rail = half_edge_from(&b, f, p1_rail).ok_or("P1 not on F's rail")?;
    let he_foot = half_edge_from(&b, f, p1_foot).ok_or("P1 not on F's foot")?;
    let (outer, inner) = if b.faces[f].outer_loop == b.half_edges[he_rail].loop_ {
        (he_rail, he_foot)
    } else {
        (he_foot, he_rail)
    };
    // Foot and rail in two loops of F (the foot an inner loop) are joined
    // by a bridge; in one loop (an earlier overflow already joined them)
    // the bridge splits F, the foot cutting a piece of it off.
    let bridge_1 = if b.half_edges[outer].loop_ == b.half_edges[inner].loop_ {
        b.mef(outer, inner)
            .map_err(|e| format!("mef at P1: {e}"))?
            .half_edge
    } else {
        b.mekr(outer, inner)
            .map_err(|e| format!("mekr: {e}"))?
            .half_edge
    };
    let f_faces = faces_with_surface(&b, &plan.support_surface);
    let (he_p2_rail, he_p2_foot) = f_faces
        .iter()
        .find_map(|&face| {
            let r = half_edge_from(&b, face, p2_rail)?;
            let ft = half_edge_from(&b, face, p2_foot)?;
            (b.half_edges[r].loop_ == b.half_edges[ft].loop_).then_some((r, ft))
        })
        .ok_or("P2's rail and foot are not on one loop of F")?;
    let bridge_2 = b
        .mef(he_p2_rail, he_p2_foot)
        .map_err(|e| format!("mef at P2: {e}"))?;
    contract_onto_rail(&mut b, bridge_1, p1_rail)?;
    contract_onto_rail(&mut b, bridge_2.half_edge, p2_rail)?;
    // The faces the overflow replaces: the lens (beside the rail's
    // middle in F's split), S's cap (across it) and T's pieces not
    // bordering F's main part.
    let lens = side_face(&b, rail_middle, |face| {
        face != s && !s_faces.contains(&face)
    })
    .ok_or("no lens beside the rail's middle")?;
    let cap = side_face(&b, rail_middle, |face| s_faces.contains(&face) || face == s)
        .ok_or("no cap of S beside the rail's middle")?;
    let f_pieces: Vec<FaceHandle> = faces_with_surface(&b, &plan.support_surface)
        .into_iter()
        .filter(|&face| face != lens)
        .collect();
    let t_main = side_face(&b, section_1, |face| {
        face_neighbours(&b, face)
            .iter()
            .any(|n| f_pieces.contains(n))
    })
    .ok_or("no main piece of the rolling strip")?;
    let main_f = *f_pieces
        .iter()
        .find(|&&face| face_neighbours(&b, face).contains(&t_main))
        .ok_or("no piece of F beside the rolling strip")?;
    let mut doomed: Vec<FaceHandle> = faces_with_surface(&b, &plan.rolling_surface)
        .into_iter()
        .filter(|&face| face != t_main && face_neighbours(&b, face).contains(&lens))
        .collect();
    doomed.push(lens);
    doomed.push(cap);
    let omega = merge_faces(&mut b, &doomed)?;
    // The rolling strip's remaining pieces meeting across its seam (the
    // seam outside the overflow, so the sections cut it in three) are one
    // face again.
    loop {
        let pieces: Vec<FaceHandle> = faces_with_surface(&b, &plan.rolling_surface)
            .into_iter()
            .filter(|&face| face != omega)
            .collect();
        let seam = pieces.iter().find_map(|&face| {
            face_half_edges(&b, face)
                .into_iter()
                .find(|&h| mate_face(&b, h).is_some_and(|m| m != face && pieces.contains(&m)))
        });
        let Some(h) = seam else { break };
        b.kef(h).map_err(|e| format!("kef across the seam: {e}"))?;
    }
    let t_main = faces_with_surface(&b, &plan.rolling_surface)
        .into_iter()
        .find(|&face| face != omega && face_neighbours(&b, face).contains(&omega))
        .ok_or("no rolling strip beside the overflow")?;
    // The merged face takes the overflow surface, wound like the rolling
    // strip relative to its chart.
    let orientation = matching_orientation(&b, t_main, plan, section_1)?;
    if let Some(face) = b.faces.get_mut(omega) {
        face.surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Nurbs(plan.surface.clone()),
            orientation,
        });
    }
    reauthor_wall(&mut b, omega, w, &plan.curve_on_wall(), tolerance)?;
    let touched = [omega, t_main, w, main_f]
        .into_iter()
        .chain(faces_with_surface(&b, &plan.onto_surface))
        .chain(faces_with_surface(&b, &plan.support_surface))
        .chain(faces_with_surface(&b, &plan.rolling_surface))
        .collect::<Vec<_>>();
    // Every pcurve on the touched faces is authored afresh: a split edge
    // keeps its piece of the old pcurve, which the recompute left in
    // place (a foot arc split by a second overflow kept a pcurve 5e-5
    // off its edge on the F piece between the two).
    for &face in &touched {
        for h in face_half_edges(&b, face) {
            if let Some(he) = b.half_edges.get_mut(h) {
                he.pcurve = None;
            }
            b.face_use_pcurves
                .remove(&popkern_core::topology::FaceUseTrimKey { face, half_edge: h });
        }
    }
    b.recompute_pcurves_for(&touched)
        .map_err(|e| format!("pcurves: {e}"))?;
    for edge in face_edges(&b, omega) {
        let dev = b.measure_edge_deviation(edge);
        if dev > reach {
            return Err(format!(
                "overflow edge {edge:?} lies {dev:.3e} off its faces"
            ));
        }
    }
    b.validate().map_err(|e| format!("validate: {e}"))?;
    *brep = b;
    Ok(())
}

/// Merge `faces` into one by killing the edges between them, then the
/// spurs left inside; returns the survivor.
fn merge_faces(b: &mut BRep, faces: &[FaceHandle]) -> Result<FaceHandle, String> {
    let mut alive: Vec<FaceHandle> = faces.to_vec();
    alive.sort();
    alive.dedup();
    loop {
        let shared = alive.iter().find_map(|&face| {
            face_half_edges(b, face)
                .into_iter()
                .find(|&h| mate_face(b, h).is_some_and(|m| m != face && alive.contains(&m)))
        });
        let Some(h) = shared else { break };
        let killed = b.kef(h).map_err(|e| format!("kef: {e}"))?;
        alive.retain(|&face| face != killed.removed_face);
    }
    let [survivor] = alive.as_slice() else {
        return Err(format!("{} overflow faces do not connect", alive.len()));
    };
    let survivor = *survivor;
    // Spurs: edges with the survivor on both sides ending at a vertex of
    // valence one.
    for _ in 0..64 {
        let spur = face_half_edges(b, survivor).into_iter().find(|&h| {
            mate_face(b, h) == Some(survivor)
                && valence(b, b.half_edges[b.half_edges[h].next].origin) == 1
        });
        let Some(h) = spur else { break };
        b.kev(h).map_err(|e| format!("kev: {e}"))?;
    }
    if face_half_edges(b, survivor)
        .into_iter()
        .any(|h| mate_face(b, h) == Some(survivor))
    {
        return Err("an edge inside the overflow face survives".into());
    }
    Ok(survivor)
}

/// The orientation flag for the overflow face: the rolling strip's, when
/// their chart normals agree along the section they share.
fn matching_orientation(
    b: &BRep,
    t_main: FaceHandle,
    plan: &OverflowPlan,
    section: EdgeHandle,
) -> Result<bool, String> {
    let t_flag = b
        .faces
        .get(t_main)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.orientation)
        .ok_or("rolling strip has no surface")?;
    let c = b
        .edges
        .get(section)
        .and_then(|e| e.curve.as_ref())
        .ok_or("section has no curve")?;
    let p = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
    let n_t = {
        let uv = (&plan.rolling_surface as &dyn SurfaceEval)
            .closest_point(&p)
            .ok_or("no chart point on the rolling strip")?;
        plan.rolling_surface.normal_at(uv.0, uv.1)
    };
    let overflow = SurfaceGeometry::Nurbs(plan.surface.clone());
    let n_o = {
        let uv = (&overflow as &dyn SurfaceEval)
            .closest_point(&p)
            .ok_or("no chart point on the overflow")?;
        overflow.normal_at(uv.0, uv.1)
    };
    Ok(if n_t.dot(&n_o) >= 0.0 {
        t_flag
    } else {
        !t_flag
    })
}

/// Re-author the edges between the overflow face and W onto the
/// overflow's wall contact, moving their interior vertices onto it (and
/// the straight edges meeting those vertices with them).
fn reauthor_wall(
    b: &mut BRep,
    omega: FaceHandle,
    w: FaceHandle,
    contact: &NurbsCurve,
    tolerance: &Tolerance,
) -> Result<(), String> {
    let window = contact.parameter_range();
    let param_of = |p: &Point3| {
        popkern_core::geometry::curve::closest_point_in_window(contact, p, window.0, window.1)
    };
    let wall_edges: Vec<HalfEdgeHandle> = face_half_edges(b, omega)
        .into_iter()
        .filter(|&h| mate_face(b, h) == Some(w))
        .collect();
    if wall_edges.is_empty() {
        return Err("the overflow face does not border W".into());
    }
    // Interior vertices of the wall chain: shared by two wall edges.
    let mut ends: Vec<VertexHandle> = Vec::new();
    for &h in &wall_edges {
        ends.push(b.half_edges[h].origin);
        ends.push(b.half_edges[b.half_edges[h].next].origin);
    }
    let interior: Vec<VertexHandle> = ends
        .iter()
        .copied()
        .filter(|v| ends.iter().filter(|x| *x == v).count() == 2)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    for v in interior {
        let old = b.vertices[v].position;
        let new = contact.point_at(param_of(&old));
        b.vertices[v].position = new;
        // Straight edges meeting the moved vertex follow it.
        let outgoing: Vec<HalfEdgeHandle> = b
            .half_edges
            .iter()
            .filter(|(_, he)| he.origin == v)
            .map(|(h, _)| h)
            .collect();
        for h in outgoing {
            let e = b.half_edges[h].edge;
            if wall_edges.iter().any(|&x| b.half_edges[x].edge == e) {
                continue;
            }
            let other = b.vertices[b.half_edges[b.half_edges[h].next].origin].position;
            let primary = b.edges[e].half_edge;
            let (start, end) = if b.half_edges[primary].origin == v {
                (new, other)
            } else {
                (other, new)
            };
            let Some(edge) = b.edges.get_mut(e) else {
                continue;
            };
            let straight = match edge.curve.as_ref().map(|c| &c.geometry) {
                Some(CurveGeometry::Line { .. }) => true,
                Some(CurveGeometry::Nurbs(n)) => is_straight(n, tolerance),
                _ => false,
            };
            match straight {
                true => {
                    edge.curve = Some(
                        CurveOnEdge::line_segment(start, end, tolerance.linear())
                            .ok_or("a moved straight edge collapsed")?,
                    );
                }
                false => return Err("a curved edge meets a moved wall vertex".into()),
            }
        }
    }
    for h in wall_edges {
        let e = b.half_edges[h].edge;
        let primary = b.edges[e].half_edge;
        let a = b.vertices[b.half_edges[primary].origin].position;
        let z = b.vertices[b.half_edges[b.half_edges[primary].next].origin].position;
        let (ta, tz) = (param_of(&a), param_of(&z));
        if contact.point_at(ta).distance_to(&a) > (tolerance.linear() * 10.0).max(1e-9) {
            return Err("a wall edge's end lies off the overflow's wall contact".into());
        }
        let geometry = if ta <= tz {
            CurveOnEdge {
                geometry: CurveGeometry::Nurbs(contact.clone()),
                t_start: ta,
                t_end: tz,
            }
        } else {
            let rev = reversed(contact);
            let (r0, r1) = rev.parameter_range();
            CurveOnEdge {
                geometry: CurveGeometry::Nurbs(rev),
                t_start: r0 + r1 - ta,
                t_end: r0 + r1 - tz,
            }
        };
        if let Some(edge) = b.edges.get_mut(e) {
            edge.curve = Some(geometry);
        }
    }
    Ok(())
}

/// Splice every planned overflow; a pair that does not splice declines
/// the operation as the interference it is (`FilletInterference`), the
/// input unchanged.
pub(crate) fn splice_all(
    brep: &mut BRep,
    plans: &[OverflowPlan],
    tolerance: &Tolerance,
) -> Result<(), crate::error::BlendError> {
    for plan in plans {
        if let Err(reason) = splice_overflow(brep, plan, tolerance) {
            if popkern_core::trace::trace_enabled("interference") {
                eprintln!(
                    "[overflow] pair ({},{}) did not splice: {reason}",
                    plan.rolling, plan.onto
                );
            }
            return Err(crate::error::BlendError::FilletInterference {
                fillet_a: plan.rolling.min(plan.onto),
                fillet_b: plan.rolling.max(plan.onto),
            });
        }
    }
    Ok(())
}
