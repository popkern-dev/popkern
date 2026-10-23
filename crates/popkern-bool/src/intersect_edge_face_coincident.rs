//! Edges lying in a face's surface (PK-750), split from
//! `intersect_edge_face.rs` under the file-size budget.

use popkern_core::geometry::CurveEval;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;

/// Samples along the edge's own window.
const SAMPLES: usize = 9;

/// Whether the edge curve over `t_range` lies on `surface`: every sample
/// (ends included) within LINEAR tolerance of its closest surface point.
/// Coincidence is the model tolerance's decision; the 10 × linear
/// honesty band is too loose for it (measured: at a coarse 1e-2
/// tolerance every unit-cube edge stays within 0.1 of an r = 0.8 sphere
/// it crosses, and skipping them dropped every crossing).
///
/// Such an edge (a rim circle in a coplanar cap) has distance
/// identically zero to the surface: any root the solver returns is an
/// arbitrary point of a coincidence, not a crossing, and the pair is
/// skipped. Measured (PK-750): a rod's half-rim in a flush box's bottom
/// face hit at a point 4.2e-5 off the true junction, the spurious node
/// chopped the rim record into zero-length pieces, and the union
/// declined on curveless edges. Coincident contact is the overlap
/// records' business.
pub(crate) fn edge_lies_on_surface(
    curve: &CurveGeometry,
    t_range: (f64, f64),
    surface: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> bool {
    let band = tolerance.linear();
    (0..SAMPLES).all(|k| {
        let t = t_range.0 + (t_range.1 - t_range.0) * k as f64 / (SAMPLES - 1) as f64;
        let p = curve.point_at(t);
        surface
            .closest_point(&p)
            .is_some_and(|(u, v)| surface.point_at(u, v).distance_to(&p) <= band)
    })
}

/// Whether the hit at `t` is the contact at the nearer endpoint of the
/// window: the curve stays on `surface` (within the honesty band) all
/// the way from the hit to that endpoint. A curve TANGENT to the
/// surface at its own endpoint has a double root there, and a double
/// root's parameter is ill-conditioned: the solver's hit slides along
/// the curve by ~sqrt(solver residual) and passes a parameter endpoint
/// test (measured, PK-750: a rod's half-rim touching a flush box's
/// tangent face at its own vertex hit 1e-5 rad inside the window; the
/// node 4.2e-5 off the vertex chopped the rim into zero-length pieces).
pub(crate) fn touches_at_endpoint(
    curve: &CurveGeometry,
    (t_lo, t_hi): (f64, f64),
    t: f64,
    surface: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> bool {
    let end = if t - t_lo <= t_hi - t { t_lo } else { t_hi };
    edge_lies_on_surface(curve, (t, end), surface, tolerance)
}

/// Positions whose parameter box bounds `face` for the edge-face trim
/// test. On a PLANE the loop's vertices alone under-cover any curved
/// boundary: a seat disc with a notch has a few vertices and a rim
/// bulging far beyond them, and the hit of a pocket's seam line inside
/// the disc was rejected as off the face (measured, PK-751: the record
/// ring then crossed the seam unsplit and the pocket wall's chart
/// overwrapped). There, each conic edge adds its exact extremes along
/// the plane's OWN u and v axes over its window, and a NURBS edge its
/// (in-plane) control points: the box is the face's, in any rotation. A
/// 3D bounding box projected into a rotated plane's frame is far larger
/// than the face and re-admits the phantom hits the test exists to
/// reject (measured: a scene rotated 30° served 1.6e-3 too much
/// volume). Other charts keep their vertex positions (their u is
/// periodic and their trim test is tuned around it).
pub(crate) fn trim_box_positions(
    brep: &popkern_core::topology::BRep,
    face: popkern_core::topology::FaceHandle,
) -> Vec<popkern_core::geometry::Point3> {
    let mut points = brep.face_outer_loop_positions(face);
    let Some(f) = brep.faces.get(face) else {
        return points;
    };
    let Some(SurfaceGeometry::Plane { u_axis, v_axis, .. }) =
        f.surface.as_ref().map(|s| &s.geometry)
    else {
        return points;
    };
    let Some(l) = brep.loops.get(f.outer_loop) else {
        return points;
    };
    let start = l.half_edge;
    let mut he = start;
    for _ in 0..=brep.half_edges.len() {
        let Some(h) = brep.half_edges.get(he) else {
            break;
        };
        if let Some(c) = brep.edges.get(h.edge).and_then(|e| e.curve.as_ref()) {
            in_plane_extremes(c, (*u_axis, *v_axis), &mut points);
        }
        he = h.next;
        if he == start {
            break;
        }
    }
    points
}

/// The points of `c` extreme along each of `axes` over its window (exact
/// for circles, arcs and ellipses), or its control points (NURBS).
fn in_plane_extremes(
    c: &popkern_core::geometry::curve::CurveOnEdge,
    axes: (popkern_core::geometry::Vec3, popkern_core::geometry::Vec3),
    points: &mut Vec<popkern_core::geometry::Point3>,
) {
    use popkern_core::geometry::PortableMath as _;
    match &c.geometry {
        CurveGeometry::Circle { center, .. }
        | CurveGeometry::Arc { center, .. }
        | CurveGeometry::Ellipse { center, .. } => {
            // p(t) = centre + U cos t + V sin t, read back through the
            // curve's own evaluator.
            let u = c.geometry.point_at(0.0) - *center;
            let v = c.geometry.point_at(std::f64::consts::FRAC_PI_2) - *center;
            let (lo, hi) = (c.t_start.min(c.t_end), c.t_start.max(c.t_end));
            let tau = std::f64::consts::TAU;
            let full = !(lo.is_finite() && hi.is_finite()) || hi - lo >= tau;
            if !full {
                points.push(c.geometry.point_at(lo));
                points.push(c.geometry.point_at(hi));
            }
            for axis in [axes.0, axes.1] {
                let peak = v.dot(&axis).patan2(u.dot(&axis));
                for base in [peak, peak + std::f64::consts::PI] {
                    let t = if full {
                        base
                    } else {
                        base + ((lo - base) / tau).ceil() * tau
                    };
                    if full || t <= hi {
                        points.push(*center + u * t.pcos() + v * t.psin());
                    }
                }
            }
        }
        CurveGeometry::Nurbs(n) => points.extend(n.control_points.iter().copied()),
        _ => crate::intersect::conservative_curve_extent(c, points),
    }
}

#[cfg(test)]
#[path = "intersect_edge_face_coincident_tests.rs"]
mod tests;

/// The trim box of a face whose loops visit fewer than three distinct
/// vertices (rim-bounded walls and caps), from its loops' pcurves:
/// their vertex box collapses to the seam vertices, and with no box
/// every hit on the face's UNBOUNDED surface was accepted (measured,
/// PK-921: a blind radial drill's wall split a torus's equator where
/// its infinite cylinder met it, 0.5 behind the drill's axial range;
/// the result carried two phantom degree-2 vertices). The pcurve
/// window is sampled, so it is padded by a tenth of its extent per
/// axis (a full circle sampled 9 times sits within 8% of its radius of
/// the sample box). `None` unless every loop half-edge has a pcurve:
/// the vertex fallback is the collapsed box this replaces.
pub(crate) fn pcurve_trim_box(
    brep: &popkern_core::topology::BRep,
    face: popkern_core::topology::FaceHandle,
) -> Option<((f64, f64), (f64, f64))> {
    let f = brep.faces.get(face)?;
    for l in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let start = brep.loops.get(l)?.half_edge;
        let mut he = start;
        for _ in 0..=brep.half_edges.len() {
            let h = brep.half_edges.get(he)?;
            h.pcurve.as_ref()?;
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    let ((u0, u1), (v0, v1)) = brep.face_trim_uv_window(face)?;
    let (pu, pv) = (0.1 * (u1 - u0), 0.1 * (v1 - v0));
    Some(((u0 - pu, u1 + pu), (v0 - pv, v1 + pv)))
}

/// The face's vertex trim box (exact over `positions`) with every axis on which it collapses
/// (all the face's vertices at one coordinate) widened to the face's
/// padded pcurve window there. Measured (PK-922): the merged half
/// torus carries three distinct vertices, all on its outer equator, so
/// its vertex box spanned v = 0 only and rejected every pierce of a
/// radial drill's seam line through the tube walls; the drill wall's
/// seam stayed unsplit at the record nodes and its chart walk declined.
/// Faces whose vertices spread on both axes keep the vertex box.
pub(crate) fn vertex_trim_box(
    brep: &popkern_core::topology::BRep,
    face: popkern_core::topology::FaceHandle,
    positions: &[popkern_core::geometry::Point3],
    surface: &SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
    let vbox =
        popkern_core::geometry::uv_projection::exact_bounded_parameter_range(positions, surface);
    let collapsed = |(lo, hi): (f64, f64)| (hi - lo).abs() <= 1e-9 * (1.0 + lo.abs().max(hi.abs()));
    if !collapsed(vbox.0) && !collapsed(vbox.1) {
        return vbox;
    }
    let Some(window) = pcurve_trim_box(brep, face) else {
        return vbox;
    };
    (
        if collapsed(vbox.0) { window.0 } else { vbox.0 },
        if collapsed(vbox.1) { window.1 } else { vbox.1 },
    )
}
