//! Trim-rail proximity helpers for the anchor-split stage anchor placement
//! (PK-185/PK-186 reflex hole-corner gate). Split from
//! `perp_anchor_positions.rs` (file-size limit).

use popkern_core::geometry::{NurbsCurve, Point3};

/// Minimum distance from `point` to the trim curve, through the
/// refined foot of [`closest_point_on_curve`]: a bare sample grid
/// over-reports by up to half the sample spacing, which on a trim
/// extended far past the edge ends exceeds the gates it feeds
/// (PK-730).
pub(super) fn distance_to_curve(curve: &NurbsCurve, point: &Point3) -> f64 {
    closest_point_on_curve(curve, point).distance_to(point)
}

/// Closest point on the trim curve to `point` (coarse sample plus
/// ternary refinement — the rail is smooth and unimodal near the
/// foot for the class this serves).
pub(super) fn closest_point_on_curve(curve: &NurbsCurve, point: &Point3) -> Point3 {
    use popkern_core::geometry::CurveEval;
    let (t0, t1) = curve.parameter_range();
    let mut bt = t0;
    let mut bd = f64::INFINITY;
    for k in 0..=64 {
        let t = t0 + (t1 - t0) * f64::from(k) / 64.0;
        let d = curve.point_at(t).distance_to(point);
        if d < bd {
            bd = d;
            bt = t;
        }
    }
    let step = (t1 - t0) / 64.0;
    let (mut lo, mut hi) = ((bt - step).max(t0), (bt + step).min(t1));
    for _ in 0..40 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if curve.point_at(m1).distance_to(point) < curve.point_at(m2).distance_to(point) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    curve.point_at(0.5 * (lo + hi))
}

/// Where two support rails on one face meet near `seed` (PK-740: the
/// mitre point of two fillets sharing a reflex corner). Alternating
/// projection converges at once for the orthogonal boss-base corner
/// and linearly for any transversal pair; `None` when the rails do
/// not meet within `gate`.
pub(super) fn mitre_point(
    own: &NurbsCurve,
    neighbour: &NurbsCurve,
    seed: Point3,
    gate: f64,
) -> Option<Point3> {
    let mut p = closest_point_on_curve(own, &seed);
    for _ in 0..64 {
        let q = closest_point_on_curve(neighbour, &p);
        let next = closest_point_on_curve(own, &q);
        let moved = next.distance_to(&p);
        p = next;
        if moved < gate * 1e-3 {
            break;
        }
    }
    let q = closest_point_on_curve(neighbour, &p);
    (p.distance_to(&q) <= gate)
        .then(|| Point3::new(0.5 * (p.x + q.x), 0.5 * (p.y + q.y), 0.5 * (p.z + q.z)))
}

/// What [`reflex_mitre`] needs about one reflex pair corner.
pub(super) struct ReflexMitreQuery<'a> {
    pub(super) brep: &'a popkern_core::topology::BRep,
    pub(super) blend_edges: &'a crate::construct::blend_edge_set::BlendEdgeSet,
    /// The support both strips share (reflex at the corner).
    pub(super) support: popkern_core::topology::FaceHandle,
    pub(super) fillet_edge: popkern_core::topology::EdgeHandle,
    pub(super) neighbour_edge: popkern_core::topology::EdgeHandle,
    pub(super) own_rail: &'a NurbsCurve,
    pub(super) neighbour_rail: &'a NurbsCurve,
    /// The own rail's end nearest the corner.
    pub(super) seed: Point3,
    pub(super) radius: f64,
    pub(super) gate: f64,
}

/// The mitre point of two fillets sharing a reflex support corner
/// (PK-740): where their support rails meet. Rails that reach each
/// other give it directly; otherwise (rails trimmed at their edges'
/// ends, as on a curved housing) it is the support foot of the ball
/// tangent to all three faces — support, own wall, neighbour's wall —
/// solved by Newton on the three offset distances, each face's offset
/// side read from the strips' own rail ends.
pub(super) fn reflex_mitre(q: &ReflexMitreQuery) -> Option<Point3> {
    if let Some(m) = mitre_point(q.own_rail, q.neighbour_rail, q.seed, q.gate) {
        return Some(m);
    }
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    let surface = |f| {
        q.brep
            .faces
            .get(f)
            .and_then(|d| d.surface.as_ref())
            .map(|s| s.geometry.clone())
    };
    let tr = popkern_core::trace::trace_enabled("anchor");
    if tr {
        eprintln!("[anchor] reflex mitre: rails do not meet; solving the corner ball");
    }
    let own_wall = q.blend_edges.other_face(q.fillet_edge, q.support)?;
    let nb_wall = q.blend_edges.other_face(q.neighbour_edge, q.support)?;
    let own_wall_rail = q.blend_edges.rail(q.fillet_edge, own_wall)?;
    let nb_wall_rail = q.blend_edges.rail(q.neighbour_edge, nb_wall)?;
    let s_surf = surface(q.support)?;
    // Foot and unit normal on a surface.
    let foot = |g: &popkern_core::geometry::surface::SurfaceGeometry, p: Point3| {
        let (u, v) = project_point_to_uv(g, &p)?;
        Some((g.point_at(u, v), g.normal_at(u, v).normalized()?))
    };
    // One strip's end-ball centre, with the offset signs it fixes on
    // the support and on its wall.
    let end_ball = |support_end: Point3, wall_rail: &NurbsCurve, wall: &_| {
        let wall_end = closest_point_on_curve(wall_rail, &support_end);
        let (fs, ns) = foot(&s_surf, support_end)?;
        let best = [1.0, -1.0]
            .into_iter()
            .map(|s| (s, fs + ns * (s * q.radius)))
            .min_by(|a, b| {
                let da = (a.1.distance_to(&wall_end) - q.radius).abs();
                let db = (b.1.distance_to(&wall_end) - q.radius).abs();
                da.total_cmp(&db)
            })?;
        let (fw, nw) = foot(wall, wall_end)?;
        let s_wall = if (best.1 - fw).dot(&nw) >= 0.0 {
            1.0
        } else {
            -1.0
        };
        Some((best.1, best.0, s_wall))
    };
    let own_surf = surface(own_wall)?;
    let nb_surf = surface(nb_wall)?;
    let (c_own, s_support, s_own) = end_ball(q.seed, own_wall_rail, &own_surf)?;
    let nb_seed = closest_point_on_curve(q.neighbour_rail, &q.seed);
    let (c_nb, _, s_nb) = end_ball(nb_seed, nb_wall_rail, &nb_surf)?;
    let faces = [(&s_surf, s_support), (&own_surf, s_own), (&nb_surf, s_nb)];
    let mut c = Point3::new(
        0.5 * (c_own.x + c_nb.x),
        0.5 * (c_own.y + c_nb.y),
        0.5 * (c_own.z + c_nb.z),
    );
    for _ in 0..60 {
        let mut rows = [[0.0; 3]; 3];
        let mut rhs = [0.0; 3];
        let mut worst = 0.0_f64;
        for (k, (g, s)) in faces.iter().enumerate() {
            let (f, n) = foot(g, c)?;
            let r = s * (c - f).dot(&n) - q.radius;
            worst = worst.max(r.abs());
            rows[k] = [s * n.x, s * n.y, s * n.z];
            rhs[k] = -r;
        }
        if tr {
            eprintln!(
                "[anchor] reflex mitre ball c=({:.4},{:.4},{:.4}) residual {worst:.3e}",
                c.x, c.y, c.z
            );
        }
        if worst <= q.gate * 1e-3 {
            let (m, _) = foot(&s_surf, c)?;
            return (m.distance_to(&q.seed) <= 4.0 * q.radius).then_some(m);
        }
        let det = |m: [[f64; 3]; 3]| {
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        };
        let d = det(rows);
        if d.abs() < 1e-12 {
            return None;
        }
        let mut step = [0.0; 3];
        for (i, si) in step.iter_mut().enumerate() {
            let mut m = rows;
            for k in 0..3 {
                m[k][i] = rhs[k];
            }
            *si = det(m) / d;
        }
        c = Point3::new(c.x + step[0], c.y + step[1], c.z + step[2]);
    }
    None
}

/// PLANAR-BEVEL anchor (chamfers): the strip's true end boundary is
/// where the BEVEL PLANE `(origin, normal)` crosses the perpendicular
/// edge `p0 → p1`. On orthogonal corners this equals the perp-walk
/// position exactly; on skewed corners the perp-walk anchor sits OFF
/// the bevel plane, the end transition gets merged into the planar
/// strip face anyway, and the face no longer contains its own boundary
/// (measured: a tilted boolean-cut bar's mitre corner 0.46 off the
/// recognized plane — the mesher and OCC's sewing both cracked along
/// both rails). Capped at 8× the chamfer distance so a near-parallel
/// edge falls back to perp-walk instead of a grotesque mitre. Returns
/// the crossing's edge parameter.
pub(super) fn bevel_plane_crossing(
    (po, pn): (Point3, popkern_core::geometry::Vec3),
    p0: Point3,
    p1: Point3,
    forward: bool,
    edge_len: f64,
    radius: f64,
) -> Option<f64> {
    let denom = pn.dot(&(p1 - p0));
    if denom.abs() <= 1e-12 {
        return None;
    }
    let s = pn.dot(&(po - p0)) / denom;
    let from_vertex = if forward {
        s * edge_len
    } else {
        (1.0 - s) * edge_len
    };
    ((0.001..0.999).contains(&s) && from_vertex > 1e-9 && from_vertex <= radius * 8.0).then_some(s)
}
