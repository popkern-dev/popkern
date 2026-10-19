//! Curve adoption helpers for boolean assembly (extracted from
//! `assemble.rs` under the file-size ratchet).

/// Subrange `stored` to the window spanning `p_start`→`p_end` (by
/// nearest-param projection) and HEAL the ends exactly onto the node
/// positions via degree-3 interpolation through samples of the true
/// curve, as dense as the curve needs. The chord fallback this replaces sat 0.17 off both adjacent
/// torus surfaces; a raw subrange window left the junction-slop
/// (measured 1.8e-3) end gap as a shell-gap validation failure.
pub(crate) fn subrange_healed_curve(
    stored: &popkern_core::geometry::CurveOnEdge,
    p_start: &popkern_core::geometry::Point3,
    p_end: &popkern_core::geometry::Point3,
) -> Option<popkern_core::geometry::CurveOnEdge> {
    use popkern_core::geometry::curve::CurveEval;
    let (t_a, d_a) = crate::merge::graph_split::nearest_param_on_curve(stored, p_start);
    let (t_b, d_b) = crate::merge::graph_split::nearest_param_on_curve(stored, p_end);
    let slop = 1.0e-2; // junction-slop class (measured 1.8e-3)
    if d_a >= slop || d_b >= slop || (t_a - t_b).abs() <= 1.0e-12 {
        return None;
    }
    let exact = 1.0e-8;
    if d_a < exact && d_b < exact {
        let mut piece = stored.clone();
        piece.t_start = t_a;
        piece.t_end = t_b;
        return Some(piece);
    }
    // Sample the true curve over the window, pin exact ends, and
    // interpolate; interior deviation stays at the slop scale near
    // the ends only. Samples double from 17 until the interpolant
    // follows the curve at every span midpoint clear of the two end
    // spans, to 1e-7 of the window's length: 17 fixed samples left a
    // 0.7-turn M12 thread flank x bore piece, about 22 long, 2e-4 off
    // both surfaces, past the 1e-4 precision bound (PK-898).
    let length = {
        let mut acc = 0.0;
        let mut prev = stored.geometry.point_at(t_a);
        for k in 1..=64 {
            let q = stored
                .geometry
                .point_at(t_a + (t_b - t_a) * k as f64 / 64.0);
            acc += prev.distance_to(&q);
            prev = q;
        }
        acc
    };
    let target = 1.0e-7 * length;
    let mut n = 16usize;
    loop {
        let at = |f: f64| t_a + (t_b - t_a) * f;
        let mut points = Vec::with_capacity(n + 1);
        let mut params = Vec::with_capacity(n + 1);
        for k in 0..=n {
            let f = k as f64 / n as f64;
            points.push(stored.geometry.point_at(at(f)));
            params.push(f);
        }
        points[0] = *p_start;
        points[n] = *p_end;
        let fit = popkern_core::geometry::fitting::interpolate_curve(&points, &params, 3).ok()?;
        let worst = (1..n.saturating_sub(1))
            .map(|k| {
                let f = (k as f64 + 0.5) / n as f64;
                fit.point_at(f)
                    .distance_to(&stored.geometry.point_at(at(f)))
            })
            .fold(0.0_f64, f64::max);
        if worst <= target || n >= 256 {
            return Some(popkern_core::geometry::CurveOnEdge {
                geometry: popkern_core::geometry::curve::CurveGeometry::Nurbs(fit),
                t_start: 0.0,
                t_end: 1.0,
            });
        }
        n *= 2;
    }
}

/// A CLOSED edge's direction is its canonical half-edge (PK-478): the
/// endpoints coincide, so no endpoint swap can express it. When the
/// curve belongs to the TWIN of the graph half-edge the edge was first
/// created from, the twin's BRep half-edge must be canonical. Left on
/// the first-needed side, a seam-fin union's housing claimed the top
/// rim's curve forward: the cap served inward, the wall's loop wound
/// 4π, and a later fillet turned the seam fin inside out (PK-730).
pub(super) fn canonical_curve_owner(
    brep: &mut popkern_core::topology::BRep,
    edge: popkern_core::topology::EdgeHandle,
    twin_half_edge: popkern_core::topology::HalfEdgeHandle,
    twin_owns_closed_curve: bool,
) {
    if twin_owns_closed_curve && let Some(e) = brep.edges.get_mut(edge) {
        e.half_edge = twin_half_edge;
    }
}
