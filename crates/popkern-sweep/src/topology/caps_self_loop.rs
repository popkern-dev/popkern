//! Half-edge direction along cap curves and the self-loop (full
//! circle) cap pcurve — split from `caps.rs` under the 500-line cap.

use popkern_math::PortableMath;
use std::f64::consts::FRAC_1_SQRT_2;

use popkern_core::constructors::builder::BRepBuilder;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{BRep, HalfEdgeHandle};

/// Whether half-edge `he` traverses its cap `curve` in increasing
/// parameter.
///
/// Cap curves are authored along the edge's CANONICAL half-edge
/// (`edge.half_edge`, the one `add_edge` returns first — the side
/// face's profile-direction half), which is the PK-303 convention.
/// On a SELF-LOOP edge (a full-circle profile: one vertex at both
/// curve ends) that convention is the only evidence: an
/// origin-vs-endpoint distance test ties and always read "forward",
/// so the start cap — whose half-edge is the twin — took the end
/// cap's winding and its normal pointed into the solid (PK-565). Open
/// edges keep the endpoint test, which measures the geometry
/// directly and holds for any producer.
pub(super) fn half_edge_runs_curve_forward(
    brep: &BRep,
    he: HalfEdgeHandle,
    curve: &NurbsCurve,
) -> bool {
    let half_edge = &brep.half_edges[he];
    let twin_origin = brep.half_edges[half_edge.twin].origin;
    if half_edge.origin == twin_origin {
        return brep.edges[half_edge.edge].half_edge == he;
    }
    let (t0, t1) = curve.parameter_range();
    let origin = brep.vertices[half_edge.origin].position;
    origin.distance_to(&curve.point_at(t0)) <= origin.distance_to(&curve.point_at(t1))
}

/// Author the pcurve of a self-loop cap half-edge: an exact rational
/// circle in the cap's UV plane that STARTS at the loop vertex and
/// winds in the half-edge's traversal sense.
///
/// The sense is measured, not assumed: the cap curve is sampled in
/// half-edge order (see [`half_edge_runs_curve_forward`]) and the
/// signed UV area decides CCW vs CW. An outer loop about an outward
/// normal in a right-handed frame measures CCW; a multi-wire cap's
/// inner (hole) loop measures CW. A pcurve starting anywhere but the
/// vertex disagrees with the loop's vertex chain at its ends.
pub(super) fn attach_self_loop_pcurve(
    builder: &mut BRepBuilder,
    he: HalfEdgeHandle,
    curve: &NurbsCurve,
    surface: &SurfaceGeometry,
) {
    const SAMPLES: usize = 16;
    let brep = builder.brep();
    let vertex = brep.vertices[brep.half_edges[he].origin].position;
    let forward = half_edge_runs_curve_forward(brep, he, curve);
    let (t0, t1) = curve.parameter_range();
    let mut ring: Vec<(f64, f64)> = Vec::with_capacity(SAMPLES);
    for k in 0..SAMPLES {
        let f = k as f64 / SAMPLES as f64;
        let f = if forward { f } else { 1.0 - f };
        match project_point_to_uv(surface, &curve.point_at(t0 + (t1 - t0) * f)) {
            Some(uv) => ring.push(uv),
            None => return,
        }
    }
    let Some((u_vertex, v_vertex)) = project_point_to_uv(surface, &vertex) else {
        return;
    };
    // Centre: the ring centroid. Samples are uniform in the curve
    // parameter, which for the quarter-arc rational circles the
    // profiles author (and for analytic arcs) places them
    // symmetrically about the centre.
    let n = ring.len() as f64;
    let cu = ring.iter().map(|p| p.0).sum::<f64>() / n;
    let cv = ring.iter().map(|p| p.1).sum::<f64>() / n;
    let (du, dv) = (u_vertex - cu, v_vertex - cv);
    let radius = (du * du + dv * dv).sqrt();
    if radius <= 1e-15 {
        return;
    }
    let signed_area: f64 = (0..ring.len())
        .map(|i| {
            let (a, b) = (ring[i], ring[(i + 1) % ring.len()]);
            a.0 * b.1 - b.0 * a.1
        })
        .sum();
    let pcurve = cap_circle_pcurve(radius, (cu, cv), dv.patan2(du), signed_area >= 0.0);
    builder.set_half_edge_pcurve(he, pcurve);
}

/// Degree-2 rational NURBS full circle in a cap's UV plane, centred
/// at `centre`, starting at polar angle `start_angle`, winding CCW
/// when `ccw` and CW otherwise. The quarter-arc control net is
/// rotated rigidly (an affine map of the control points reproduces
/// the rotated rational curve exactly); the CW circle is the CCW one
/// reversed (the knot vector is symmetric, so it maps onto itself).
pub(super) fn cap_circle_pcurve(
    radius: f64,
    centre: (f64, f64),
    start_angle: f64,
    ccw: bool,
) -> NurbsCurve {
    let (s, c) = start_angle.psin_cos();
    let w = FRAC_1_SQRT_2;
    let unit = [
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 0.0),
        (-1.0, -1.0),
        (0.0, -1.0),
        (1.0, -1.0),
        (1.0, 0.0),
    ];
    let mut control_points: Vec<Point3> = unit
        .iter()
        .map(|&(x, y)| {
            Point3::new(
                centre.0 + radius * (c * x - s * y),
                centre.1 + radius * (s * x + c * y),
                0.0,
            )
        })
        .collect();
    let mut weights = vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
    if !ccw {
        control_points.reverse();
        weights.reverse();
    }
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points,
        weights,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn circle_pcurve_starts_at_its_angle_and_winds_as_asked() {
        let (r, centre, a0) = (0.7, (0.3, -0.2), 2.1_f64);
        for ccw in [true, false] {
            let c = cap_circle_pcurve(r, centre, a0, ccw);
            let p0 = c.point_at(0.0);
            assert!((p0.x - (centre.0 + r * a0.cos())).abs() < 1e-12);
            assert!((p0.y - (centre.1 + r * a0.sin())).abs() < 1e-12);
            let p1 = c.point_at(1.0);
            assert!(p0.distance_to(&p1) < 1e-12, "closed");
            let mut area = 0.0;
            let pts: Vec<Point3> = (0..64).map(|k| c.point_at(k as f64 / 64.0)).collect();
            for i in 0..pts.len() {
                let (a, b) = (pts[i], pts[(i + 1) % pts.len()]);
                area += a.x * b.y - b.x * a.y;
                let rr = ((a.x - centre.0).powi(2) + (a.y - centre.1).powi(2)).sqrt();
                assert!((rr - r).abs() < 1e-12, "exact circle");
            }
            assert_eq!(area > 0.0, ccw);
        }
    }
}
