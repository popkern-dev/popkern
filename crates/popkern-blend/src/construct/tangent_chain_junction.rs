//! Shared pieces of the tangent-chain constructions — closed rings
//! (PK-733, `tangent_chain_surgery`) and open chains whose ends vanish
//! (PK-764, `tangent_chain_open`): the chain's links in order, and one
//! junction's exact section. Split from `tangent_chain_surgery` for
//! the line cap.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use super::EdgeGeometry;
use super::tangent_chain_geometry::{
    Circle, Link, Third, ball_center, chain_tangent, decline, exact_section, nearest, reversed,
    segment_distance, verified_circle,
};
use crate::analysis::Convexity;
use crate::error::BlendFailure;

/// The chain's links in order, following half-edges head to tail from
/// `h0`. A closed chain must return to `h0`; an open one ends at its
/// last member (`members.len()` links either way).
pub(super) fn chain_links<'a>(
    brep: &BRep,
    edge_geometries: &'a [EdgeGeometry],
    members: &[usize],
    h0: HalfEdgeHandle,
    closed: bool,
    at: Point3,
    tolerance: &Tolerance,
) -> Result<Vec<Link<'a>>, BlendFailure> {
    let mut links: Vec<Link<'a>> = Vec::with_capacity(members.len());
    let mut h = h0;
    for step in 0..members.len() {
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
            verified_circle(p_trim.0, p_trim.1, support(p).as_ref(), tolerance),
        );
        let q_trim = (
            q_trim.0,
            verified_circle(q_trim.0, q_trim.1, support(q).as_ref(), tolerance),
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
        // The next chain edge leaves this half-edge's destination; an
        // open chain ends at its last member.
        if !closed && step + 1 == members.len() {
            break;
        }
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
    if closed && h != links[0].h && brep.half_edges[h].twin != links[0].h {
        return Err(decline("the chain does not close on itself".into(), at));
    }
    Ok(links)
}

/// The junction between links `a` and `b` (`v = dest(a.h)`): its third
/// edge, the exact section contacts on `a`'s P and Q supports (the
/// third edge's side snapped onto it: the third edge is kept and
/// shortened to its contact), and the section's ball centre and normal.
/// When the section at V misses the third edge, the junction slides
/// along the chain to the station where its contact lies on it
/// (`tangent_chain_station`, PK-736).
pub(super) fn junction_section(
    brep: &BRep,
    a: &Link<'_>,
    b: &Link<'_>,
    band: f64,
) -> Result<JunctionSection, BlendFailure> {
    let hi = &brep.half_edges[a.h];
    let v = brep.half_edges[hi.twin].origin;
    let vp = brep.vertices[v].position;
    let valence = brep.half_edges.values().filter(|he| he.origin == v).count();
    if valence != 2 && valence != 3 {
        return Err(decline(
            format!("junction {v:?} has valence {valence}, not 2 or 3"),
            vp,
        ));
    }
    let next_p = hi.next;
    let next_q = brep.half_edges[brep.half_edges[b.h].twin].next;
    let third = if next_p != b.h {
        Third {
            edge: Some(brep.half_edges[next_p].edge),
            on_p: true,
        }
    } else if next_q != hi.twin {
        Third {
            edge: Some(brep.half_edges[next_q].edge),
            on_p: false,
        }
    } else if valence == 2 {
        Third {
            edge: None,
            on_p: true,
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
    let (Some(sp), Some(sq), Some(tangent)) = (surf(a.h), surf(hi.twin), chain_tangent(brep, a.h))
    else {
        return Err(decline(
            format!("junction {v:?} lacks supports or a tangent"),
            vp,
        ));
    };
    let r = a.eg.radius_at_end;
    let seed = ball_center(&sp, cp0, cq0, tangent, r).unwrap_or(vp);
    let Some((center, cp, cq)) = exact_section(&sp, &sq, vp, tangent, r, seed) else {
        return Err(decline(
            format!("no ball of radius {r} sits on both supports at {v:?}"),
            vp,
        ));
    };
    // The third edge must reach its side's contact: it is kept and
    // shortened to it.
    let Some(third_edge) = third.edge else {
        return Ok((third, (cp, cq), (center, tangent)));
    };
    let on_third = if third.on_p { cp } else { cq };
    let t_curve = brep.edges[third_edge].curve.clone();
    let off = match &t_curve {
        Some(c) => c
            .geometry
            .closest_point(&on_third)
            .map(|s| c.geometry.point_at(s).distance_to(&on_third))
            .unwrap_or(f64::INFINITY),
        None => segment_distance(brep, third_edge, on_third),
    };
    let (cp, cq, center, normal) = if off <= band {
        (cp, cq, center, tangent)
    } else if let Some(st) = super::tangent_chain_station::station_on_third(
        brep,
        (hi.edge, brep.half_edges[b.h].edge),
        v,
        (&sp, &sq),
        third_edge,
        third.on_p,
        r,
        (center, tangent),
        band,
    ) {
        (st.cp, st.cq, st.center, st.normal)
    } else {
        return Err(decline(
            format!(
                "third edge {:?} misses the section contact by {off:.3e} (contacts {cp:?} / {cq:?})",
                third_edge
            ),
            vp,
        ));
    };
    let on_third = if third.on_p { cp } else { cq };
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
    Ok((third, (cp, cq), (center, normal)))
}

/// A junction's third edge, contacts `(P, Q)` and section `(centre,
/// normal)`.
pub(super) type JunctionSection = (Third, (Point3, Point3), (Point3, Vec3));

/// A swept strip's own boundary rail between `a` and `b`: the
/// iso-curve of the strip surface both points lie on. The walker's
/// strip is exact at its rails (the sections lie on it to 1e-9) while
/// the separately fitted trims strayed 0.03–0.14 off their supports
/// at r = 6 on a cylinder tee's saddle — so for NURBS strips the rail
/// is the trim. `None` for analytic strips (their trims are exact) or
/// when the points do not share a boundary iso-line within `band`.
pub(super) fn strip_rail(
    surface: &SurfaceGeometry,
    a: Point3,
    b: Point3,
    band: f64,
) -> Option<CurveOnEdge> {
    let SurfaceGeometry::Nurbs(ns) = surface else {
        return None;
    };
    let (ua, va) = surface.closest_point(&a)?;
    let (ub, vb) = surface.closest_point(&b)?;
    if surface.point_at(ua, va).distance_to(&a) > band
        || surface.point_at(ub, vb).distance_to(&b) > band
    {
        return None;
    }
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    let eps_u = (u1 - u0).abs() * 1e-7;
    let eps_v = (v1 - v0).abs() * 1e-7;
    let at_bound = |x: f64, lo: f64, hi: f64, eps: f64| {
        if (x - lo).abs() <= eps {
            Some(lo)
        } else if (x - hi).abs() <= eps {
            Some(hi)
        } else {
            None
        }
    };
    // A VANISHING end (PK-764): the point is the strip's collapsed end
    // row, where every v is the same point, so its projection's v is
    // arbitrary — it takes the other end's boundary v (the rail runs
    // along that boundary into the tip). Before, the rail fell back to
    // the fitted trim, 1.9e-3 off its support.
    let collapsed = |u: f64| {
        surface
            .point_at(u, v0)
            .distance_to(&surface.point_at(u, v1))
            <= band
    };
    let tip = |u: f64| at_bound(u, u0, u1, eps_u).is_some_and(collapsed);
    let (va, vb) = match (tip(ua), tip(ub)) {
        (true, false) => (at_bound(vb, v0, v1, eps_v).unwrap_or(va), vb),
        (false, true) => (va, at_bound(va, v0, v1, eps_v).unwrap_or(vb)),
        _ => (va, vb),
    };
    // A point a hair inside the chart (a junction slid off its strip's
    // stations, PK-736: v = 0.999999) is on the boundary when the
    // boundary's own point there is within `band` of it.
    let snap_v = |u: f64, v: f64, p: Point3| {
        at_bound(v, v0, v1, eps_v)
            .or_else(|| {
                [v0, v1]
                    .into_iter()
                    .find(|&b| surface.point_at(u, b).distance_to(&p) <= band)
            })
            .unwrap_or(v)
    };
    let (va, vb) = (snap_v(ua, va, a), snap_v(ub, vb, b));
    let (curve, mut sa, mut sb) = if let (Some(fa), Some(fb)) =
        (at_bound(va, v0, v1, eps_v), at_bound(vb, v0, v1, eps_v))
        && fa == fb
    {
        (ns.iso_curve_at_v(fa), ua, ub)
    } else if let (Some(fa), Some(fb)) = (at_bound(ua, u0, u1, eps_u), at_bound(ub, u0, u1, eps_u))
        && fa == fb
    {
        (ns.iso_curve_at_u(fa), va, vb)
    } else {
        return None;
    };
    // The window ends from the rail's own projection: the surface
    // projection's parameters left the rail ends ~1e-6 off the
    // junction vertices (an OpenLoop at validation).
    sa = curve.closest_point(&a).unwrap_or(sa);
    sb = curve.closest_point(&b).unwrap_or(sb);
    if sa <= sb {
        Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(curve),
            t_start: sa,
            t_end: sb,
        })
    } else {
        let (lo, hi) = curve.parameter_range();
        Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(reversed(&curve)),
            t_start: lo + hi - sa,
            t_end: lo + hi - sb,
        })
    }
}
