//! Mitre corners in open tangent chains (PK-768,
//! `docs/architecture/fillet-open-tangent-chains.md`): a chain's end
//! link meets ANOTHER fillet edge at a corner V whose third edge `g`
//! stays sharp. The per-edge pipeline builds this corner as its classic
//! mitre: the two strips meet along their intersection curve M and each
//! runs out on the face it does not share with the other. Here it is a
//! junction of the chain surgery whose section is M.
//!
//! **Geometry (exact).** Faces at V: the shared support B (between the
//! two filleted edges), A (the end link's other support) and C (the
//! other fillet's), with `g = A ∩ C`. For two cylinder strips of one
//! radius their axes meet at O, the ball centre at the corner:
//! - `m1` — both strips' rails on B pass through O's foot on B;
//! - `m2` — the end link's rail on A and the other strip's rail on C
//!   both reach `g` there (each rail is O's foot on its face run along
//!   its axis to the other face);
//! - M is the strip's section by the plane through O, m1 and m2 (an
//!   exact ellipse; a circle when the plane is normal to the axis), and
//!   lies on both cylinders.
//!
//! **Topology.** `g` is split at m2 and its corner piece V..m2 cut away;
//! m1 is a new vertex on B; V goes. Strip `i` ends on M and strip `i+1`
//! starts on it, exactly as at a tangent junction.
//!
//! **Declines** (typed, the surgery rolls back): faces that are not
//! planes, strips that are not cylinders of one radius, axes that do
//! not meet, rails that do not reach `g` at one point, or an M that
//! leaves the second cylinder.

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::tangent_chain_geometry::{Link, decline};
use crate::error::BlendFailure;

/// A solved mitre between links `a` (ending at V) and `b` (leaving it).
pub(super) struct Mitre {
    /// The sharp third edge `g`, split at `on_g`.
    pub(super) g: EdgeHandle,
    /// m2: the outer rails' common point on `g`.
    pub(super) on_g: Point3,
    /// m1: the shared face's rail crossing.
    pub(super) on_shared: Point3,
    /// Whether link `a`'s P side is the `g` side (its P contact is m2).
    pub(super) p_is_g: bool,
    /// M.
    pub(super) curve: CurveOnEdge,
    /// Whether M's curve runs from `a`'s P contact to its Q contact.
    pub(super) starts_at_p: bool,
    /// The ball centre at the corner (both axes meet there).
    pub(super) center: Point3,
}

fn plane(brep: &BRep, face: FaceHandle) -> Option<(Point3, Vec3)> {
    match brep.faces.get(face)?.surface.as_ref()?.geometry {
        SurfaceGeometry::Plane { origin, normal, .. } => Some((origin, normal.normalized()?)),
        _ => None,
    }
}

fn cylinder(s: &SurfaceGeometry) -> Option<(Point3, Vec3, f64)> {
    match *s {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => Some((origin, axis.normalized()?, radius)),
        _ => None,
    }
}

/// The foot of `p` on the plane `(o, n)`.
fn foot(p: Point3, (o, n): (Point3, Vec3)) -> Point3 {
    p - n * (p - o).dot(&n)
}

/// Where the line `p + t·d` meets the plane `(o, n)`.
fn line_plane(p: Point3, d: Vec3, (o, n): (Point3, Vec3)) -> Option<Point3> {
    let den = d.dot(&n);
    (den.abs() > 1e-9).then(|| p + d * ((o - p).dot(&n) / den))
}

/// Solve the mitre at `v = dest(a.h)`, where link `b` leaves.
///
/// # Errors
///
/// A typed decline naming the unmet precondition.
pub(super) fn mitre_section(
    brep: &BRep,
    a: &Link<'_>,
    b: &Link<'_>,
    band: f64,
) -> Result<Mitre, BlendFailure> {
    let face = |he: HalfEdgeHandle| brep.loops[brep.half_edges[he].loop_].face;
    let ha = &brep.half_edges[a.h];
    let v = brep.half_edges[ha.twin].origin;
    let vp = brep.vertices[v].position;
    let (pa, qa) = (face(a.h), face(ha.twin));
    let (pb, qb) = (face(b.h), face(brep.half_edges[b.h].twin));
    let shared = [pa, qa]
        .into_iter()
        .filter(|f| *f == pb || *f == qb)
        .collect::<Vec<_>>();
    let [shared] = shared.as_slice() else {
        return Err(decline(
            format!("mitre at {v:?}: the strips do not share exactly one support"),
            vp,
        ));
    };
    let shared = *shared;
    let fa = if pa == shared { qa } else { pa };
    let fb = if pb == shared { qb } else { pb };
    // The third edge: the one at V that is neither link's edge.
    let edges: Vec<EdgeHandle> = brep
        .half_edges
        .values()
        .filter(|he| he.origin == v)
        .map(|he| he.edge)
        .collect();
    let g = match edges
        .iter()
        .copied()
        .filter(|e| *e != ha.edge && *e != brep.half_edges[b.h].edge)
        .collect::<Vec<_>>()
        .as_slice()
    {
        [g] if edges.len() == 3 => *g,
        _ => {
            return Err(decline(
                format!("mitre at {v:?} has valence {}, not 3", edges.len()),
                vp,
            ));
        }
    };
    let (Some(plane_a), Some(plane_b), Some(plane_c)) =
        (plane(brep, fa), plane(brep, shared), plane(brep, fb))
    else {
        return Err(decline(format!("mitre at {v:?} on curved faces"), vp));
    };
    let (Some((oa, da, ra)), Some((ob, db, rb))) =
        (cylinder(&a.eg.surface), cylinder(&b.eg.surface))
    else {
        return Err(decline(
            format!("mitre at {v:?} of non-cylinder strips"),
            vp,
        ));
    };
    if (ra - rb).abs() > band {
        return Err(decline(
            format!("mitre at {v:?} of unequal radii {ra} / {rb}"),
            vp,
        ));
    }
    // O: where the axes meet (closest points of the two lines).
    let w = oa - ob;
    let (aa, ab, bb) = (da.dot(&da), da.dot(&db), db.dot(&db));
    let (dw_a, dw_b) = (da.dot(&w), db.dot(&w));
    let den = aa * bb - ab * ab;
    if den.abs() < 1e-12 {
        return Err(decline(format!("mitre at {v:?}: parallel strip axes"), vp));
    }
    let s = (ab * dw_b - bb * dw_a) / den;
    let t = (aa * dw_b - ab * dw_a) / den;
    let (on_a, on_b) = (oa + da * s, ob + db * t);
    if on_a.distance_to(&on_b) > band {
        return Err(decline(
            format!(
                "mitre at {v:?}: strip axes miss by {:.3e}",
                on_a.distance_to(&on_b)
            ),
            vp,
        ));
    }
    let center = on_a + (on_b - on_a) * 0.5;
    let m1 = foot(center, plane_b);
    let m2a = line_plane(foot(center, plane_a), da, plane_c);
    let m2b = line_plane(foot(center, plane_c), db, plane_a);
    let (Some(m2a), Some(m2b)) = (m2a, m2b) else {
        return Err(decline(format!("mitre at {v:?}: a rail misses g"), vp));
    };
    let on_g_curve = brep.edges[g].curve.as_ref().and_then(|c| {
        let s = c.geometry.closest_point(&m2a)?;
        Some(c.geometry.point_at(s).distance_to(&m2a))
    });
    if m2a.distance_to(&m2b) > band || on_g_curve.is_none_or(|d| d > band) {
        return Err(decline(
            format!(
                "mitre at {v:?}: rails reach g {:.3e} apart, {:.3e} off it",
                m2a.distance_to(&m2b),
                on_g_curve.unwrap_or(f64::INFINITY)
            ),
            vp,
        ));
    }
    let m2 = m2a;
    // M: the end link's cylinder cut by the plane through O, m1, m2.
    let Some(n) = (m1 - center).cross(&(m2 - center)).normalized() else {
        return Err(decline(
            format!("mitre at {v:?}: degenerate mitre plane"),
            vp,
        ));
    };
    let p_is_g = pa != shared;
    let (p1, p2) = if p_is_g { (m2, m1) } else { (m1, m2) };
    let cos = da.dot(&n);
    let section = if cos.abs() >= 1.0 - 1e-9 {
        super::strip_end_arcs::strip_end_section_curves::verified_circle_section(
            &oa, &da, ra, center, p1, p2, &center, &n,
        )
    } else {
        super::strip_end_arcs::strip_end_section_curves::verified_ellipse_section(
            &oa, &da, ra, center, p1, p2, &center, &n, cos,
        )
    };
    let Some((curve, _, _)) = section else {
        return Err(decline(format!("mitre at {v:?}: no verified section"), vp));
    };
    // M lies on the second cylinder too.
    let mid = curve.geometry.point_at(0.5 * (curve.t_start + curve.t_end));
    let rel = mid - ob;
    let off = ((rel - db * rel.dot(&db)).length() - rb).abs();
    if off > band {
        return Err(decline(
            format!("mitre at {v:?}: M leaves the second strip by {off:.3e}"),
            vp,
        ));
    }
    // Whether the curve runs from a's P contact to its Q contact.
    let starts_at_p = curve.geometry.point_at(curve.t_start).distance_to(&p1)
        <= curve.geometry.point_at(curve.t_end).distance_to(&p1);
    Ok(Mitre {
        g,
        on_g: m2,
        on_shared: m1,
        p_is_g,
        curve,
        starts_at_p,
        center,
    })
}

/// The mitre's two contact vertices, P side first (for link `a`), and
/// the corner piece of `g` (V..m2): `g` is split at m2 and m1 is a new
/// vertex on the shared face.
///
/// # Errors
///
/// A typed decline when `g` cannot be split at m2.
pub(super) fn mitre_vertices(
    brep: &mut BRep,
    mitre: &Mitre,
    v: popkern_core::topology::VertexHandle,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<
    (
        popkern_core::topology::VertexHandle,
        popkern_core::topology::VertexHandle,
        EdgeHandle,
    ),
    BlendFailure,
> {
    let split = brep
        .split_edge_at_point(mitre.g, mitre.on_g, tolerance)
        .map_err(|e| decline(format!("mitre: g split at m2: {e}"), mitre.on_g))?;
    brep.vertices[split.vertex].position = mitre.on_g;
    let has_v = |brep: &BRep, e: EdgeHandle| {
        let he = &brep.half_edges[brep.edges[e].half_edge];
        he.origin == v || brep.half_edges[he.twin].origin == v
    };
    let corner = if has_v(brep, split.edge_a) {
        split.edge_a
    } else {
        split.edge_b
    };
    let m1 = brep.vertices.insert(popkern_core::topology::Vertex {
        position: mitre.on_shared,
        outgoing: HalfEdgeHandle::default(),
        name: None,
        tolerance: None,
    });
    Ok(if mitre.p_is_g {
        (split.vertex, m1, corner)
    } else {
        (m1, split.vertex, corner)
    })
}

/// Cut the corner away once both strips are wired: each strip's support
/// half-edge on the `g` side touches m2 next to a half of the corner
/// piece; both halves leave their loops, and the corner piece and V go.
///
/// # Errors
///
/// A typed decline when a support loop does not run through the corner
/// piece at m2.
pub(super) fn cut_mitre_corner(
    brep: &mut BRep,
    corner: EdgeHandle,
    v: popkern_core::topology::VertexHandle,
    m2: popkern_core::topology::VertexHandle,
    supports: [HalfEdgeHandle; 2],
) -> Result<(), BlendFailure> {
    let at = brep.vertices[m2].position;
    for s in supports {
        let removed = if brep.half_edges[s].origin == m2 {
            brep.half_edges[s].prev
        } else {
            brep.half_edges[s].next
        };
        if brep.half_edges[removed].edge != corner {
            return Err(decline(
                "mitre: support loop does not run through g's corner piece".into(),
                at,
            ));
        }
        let (pr, nx) = (brep.half_edges[removed].prev, brep.half_edges[removed].next);
        brep.half_edges[pr].next = nx;
        brep.half_edges[nx].prev = pr;
        let lp = brep.half_edges[removed].loop_;
        if brep.loops[lp].half_edge == removed {
            brep.loops[lp].half_edge = nx;
        }
        brep.remove_half_edge_reanchored(removed);
    }
    brep.edges.remove(corner);
    if brep.half_edges.values().any(|he| he.origin == v) {
        return Err(decline("mitre: the corner vertex keeps edges".into(), at));
    }
    brep.vertices.remove(v);
    Ok(())
}
