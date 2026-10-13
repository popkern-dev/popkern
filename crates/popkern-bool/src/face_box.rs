//! Face boxes for the boolean's face-pair prefilter: every box must
//! CONTAIN its face (an under-covering box silently drops the pair's
//! intersections) and should be trim-scale (an over-covering one pairs
//! the face with everything near its untrimmed surface).

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_math::PortableMath;

/// Sample a face to compute its AABB.
///
/// For faces with finite parameter ranges (e.g. NURBS, sphere, torus), samples
/// a 5×5 grid over the surface. For faces with infinite ranges (e.g. planes),
/// falls back to collecting vertex positions from the face's outer loop.
///
/// Returns `None` if the face has no surface geometry or no valid vertices.
pub(crate) fn face_aabb(brep: &BRep, face_handle: FaceHandle) -> Option<Aabb> {
    let surface_on_face = brep.faces[face_handle].surface.as_ref()?;
    let surf = &surface_on_face.geometry;
    let ((u_min, u_max), (v_min, v_max)) = surf.parameter_range();

    let mut points = Vec::new();

    // Boundary CURVE samples, not just vertices: a planar face
    // bounded by a closed rim (cylinder cap) has ONE vertex, and a
    // vertex-only AABB collapses to the seam point — every face
    // pair with such a cap was pruned before intersection ran
    // (measured: pockets into swept-cylinder caps silently
    // no-opped). Curve samples span the true extent; for straight
    // polygonal loops they add interpolants on the chords, which
    // never widen the box.
    points.extend(brep.face_outer_loop_sample_positions(face_handle, 8));

    // CONSERVATIVE curved-edge extents. Finite sampling of a curved
    // boundary UNDER-approximates its extent — 8 samples of a
    // rational circle (whose parameter is not angle-uniform) missed
    // the 180° extreme by 2.2 mm on an r = 36 bore, the box stopped
    // at x = −58.8 instead of −61, the breached wall at x = −60 fell
    // outside it, and the BVH never paired the faces: the whole
    // plane×cylinder curve set went missing and the subtract failed
    // LowValence. An AABB used as a pruning gate must CONTAIN the
    // face; over-covering only costs candidate pairs. Analytic
    // curves get exact per-axis amplitude boxes (full-sweep
    // superset for arcs), NURBS curves their control-point hull
    // (positive weights ⇒ curve ⊂ hull), parabola/hyperbola their
    // endpoint + interior derivative-zero extremes.
    let mut extend_curve_extent =
        |c: &popkern_core::geometry::curve::CurveOnEdge| conservative_curve_extent(c, &mut points);
    let face = &brep.faces[face_handle];
    let mut loops = vec![face.outer_loop];
    loops.extend(face.inner_loops.iter().copied());
    for lh in loops {
        let Some(lp) = brep.loops.get(lh) else {
            continue;
        };
        let start = lp.half_edge;
        let mut he = start;
        for _ in 0..brep.half_edges.len() + 1 {
            let Some(h) = brep.half_edges.get(he) else {
                break;
            };
            if let Some(edge) = brep.edges.get(h.edge)
                && let Some(c) = edge.curve.as_ref()
            {
                extend_curve_extent(c);
            }
            he = h.next;
            if he == start {
                break;
            }
        }
    }

    // A NURBS face trimmed to part of its chart (a thread runout corner
    // cut from a multi-wrap band) takes the control hull of the
    // sub-patch over its padded trim window: rigorous (positive
    // weights) and trim-scale, where the 5x5 whole-surface grid boxed
    // the band (PK-493: runout corners 3x their trim).
    if let popkern_core::geometry::SurfaceGeometry::Nurbs(nurbs) = surf
        && let Some(hull) =
            trimmed_nurbs_hull(brep, face_handle, nurbs, (u_min, u_max), (v_min, v_max))
    {
        points.extend(hull);
        return Aabb::from_points(&points);
    }

    // If the surface has a finite parameter range, also sample a 5×5 grid for
    // curved surfaces where vertices may not capture the full extent.
    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
        const N: usize = 5;
        for i in 0..N {
            let u = u_min + (u_max - u_min) * (i as f64) / ((N - 1) as f64);
            for j in 0..N {
                let v = v_min + (v_max - v_min) * (j as f64) / ((N - 1) as f64);
                points.push(surf.point_at(u, v));
            }
        }
    }

    Aabb::from_points(&points)
}

/// The sub-patch hull of a NURBS face whose trim window (padded by
/// [`TRIM_WINDOW_PAD`] of its span: the window is pcurve-sampled) is a
/// proper part of the chart; `None` when it covers the chart or the
/// hull is unavailable.
fn trimmed_nurbs_hull(
    brep: &BRep,
    face: FaceHandle,
    nurbs: &popkern_core::geometry::NurbsSurface,
    u_domain: (f64, f64),
    v_domain: (f64, f64),
) -> Option<Vec<popkern_core::geometry::Point3>> {
    let ((u0, u1), (v0, v1)) = brep.face_trim_uv_window(face)?;
    let pad = |lo: f64, hi: f64, d: (f64, f64)| {
        let p = TRIM_WINDOW_PAD * (hi - lo).max(TRIM_WINDOW_PAD * (d.1 - d.0));
        ((lo - p).max(d.0), (hi + p).min(d.1))
    };
    let (u, v) = (pad(u0, u1, u_domain), pad(v0, v1, v_domain));
    let covers = |w: (f64, f64), d: (f64, f64)| w.0 <= d.0 && w.1 >= d.1;
    if covers(u, u_domain) && covers(v, v_domain) {
        return None;
    }
    nurbs.sub_patch_hull_points(u, v)
}

/// Pad of a sampled trim window, as a fraction of its span (a pcurve
/// bulges between its samples by far less).
const TRIM_WINDOW_PAD: f64 = 0.05;

/// Push CONSERVATIVE extent points for a curve into `points` — a
/// superset of the curve's true axis-aligned extent. Finite sampling
/// UNDER-approximates curved extents (a rational circle's parameter
/// is not angle-uniform; a closed rim's two vertices coincide at the
/// seam), and AABBs used as pruning gates must CONTAIN the geometry.
/// Analytic curves contribute exact per-axis amplitude boxes
/// (full-sweep superset for arcs), NURBS their control-point hull
/// (positive weights ⇒ curve ⊂ hull), parabola/hyperbola their
/// interior derivative-zero extremes.
pub(crate) fn conservative_curve_extent(
    c: &popkern_core::geometry::curve::CurveOnEdge,
    points: &mut Vec<popkern_core::geometry::Point3>,
) {
    use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
    match &c.geometry {
        CurveGeometry::Line { .. } => {} // endpoints suffice
        CurveGeometry::Circle { center, .. }
        | CurveGeometry::Arc { center, .. }
        | CurveGeometry::Ellipse { center, .. } => {
            // p(t) = centre + U cos t + V sin t for all three (U, V
            // read back through the curve's own evaluator). The extent
            // over the EDGE's range, not the full sweep (PK-731): a
            // plane nearly parallel to a hole's axis cuts it in an
            // ellipse a thousand times longer than the arc the edge
            // uses, and the full-sweep box made the hole's face box
            // span the whole housing — every later tool paired with
            // it in the SSI and every ray tested it.
            let u = c.geometry.point_at(0.0) - *center;
            let v = c.geometry.point_at(std::f64::consts::FRAC_PI_2) - *center;
            conic_arc_extent(*center, u, v, c.t_start, c.t_end, points);
        }
        CurveGeometry::Parabola {
            linear, quadratic, ..
        } => {
            // Interior extreme along axis w at linear·w + 2t·quadratic·w = 0.
            for (l, q) in [
                (linear.x, quadratic.x),
                (linear.y, quadratic.y),
                (linear.z, quadratic.z),
            ] {
                if q.abs() > 1e-30 {
                    let t = -l / (2.0 * q);
                    if t > c.t_start.min(c.t_end) && t < c.t_start.max(c.t_end) {
                        points.push(c.geometry.point_at(t));
                    }
                }
            }
        }
        CurveGeometry::Hyperbola {
            semi_major,
            semi_minor,
            ..
        } => {
            // Extreme along w where sinh(t)·(sM·w) + cosh(t)·(sm·w) = 0.
            for (a, b) in [
                (semi_major.x, semi_minor.x),
                (semi_major.y, semi_minor.y),
                (semi_major.z, semi_minor.z),
            ] {
                if a.abs() > 1e-30 && (b / a).abs() < 1.0 {
                    let t = (-b / a).patanh();
                    if t > c.t_start.min(c.t_end) && t < c.t_start.max(c.t_end) {
                        points.push(c.geometry.point_at(t));
                    }
                }
            }
        }
        // The hull of the EDGE's own span: the whole curve's hull boxed a
        // thread runout edge as its entire helix (PK-493).
        CurveGeometry::Nurbs(n) => {
            let sub = c.geometry.to_nurbs_in_range((c.t_start, c.t_end));
            if sub.weights.iter().all(|w| *w > 0.0) {
                points.extend(sub.control_points.iter().copied());
            } else {
                points.extend(n.control_points.iter().copied());
            }
        }
    }
}

/// Push the extreme points of `centre + u cos t + v sin t` over
/// `t ∈ [t0, t1]` (either order): both ends, and on each axis the
/// angles where that coordinate peaks (`atan2(v_w, u_w)` and its
/// antipode) when they fall inside the range. Their box CONTAINS the
/// arc; a range of a full turn or more yields the full-sweep extremes.
pub(crate) fn conic_arc_extent(
    center: popkern_core::geometry::Point3,
    u: popkern_core::geometry::Vec3,
    v: popkern_core::geometry::Vec3,
    t0: f64,
    t1: f64,
    points: &mut Vec<popkern_core::geometry::Point3>,
) {
    use popkern_core::geometry::PortableMath as _;
    let at = |t: f64| center + u * t.pcos() + v * t.psin();
    let (lo, hi) = (t0.min(t1), t0.max(t1));
    let tau = std::f64::consts::TAU;
    if !(lo.is_finite() && hi.is_finite()) || hi - lo >= tau {
        for (uw, vw) in [(u.x, v.x), (u.y, v.y), (u.z, v.z)] {
            let peak = vw.patan2(uw);
            points.push(at(peak));
            points.push(at(peak + std::f64::consts::PI));
        }
        return;
    }
    points.push(at(lo));
    points.push(at(hi));
    for (uw, vw) in [(u.x, v.x), (u.y, v.y), (u.z, v.z)] {
        let peak = vw.patan2(uw);
        for base in [peak, peak + std::f64::consts::PI] {
            // The representative of `base` at or above `lo`.
            let t = base + ((lo - base) / tau).ceil() * tau;
            if t <= hi {
                points.push(at(t));
            }
        }
    }
}
