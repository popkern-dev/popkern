//! Tangent curves traced from the face-pair pass's contacts (W3c,
//! tangency design §4.8): a contact classified `TangentCurve` is traced
//! along the curve (`contact_trace::trace_tangent_curve`), the trace
//! clipped to both faces, and each run that is a line or a circle within
//! the linear tolerance becomes a contact curve, exact. A run neither is
//! left to the general glue's detection until a fit with a bounded
//! deviation lands (the NURBS tangent-curve step).

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::intersection::contact_classify::{ContactScale, SurfaceContactKind};
use popkern_core::geometry::intersection::contact_discovery::SurfaceContact;
use popkern_core::geometry::intersection::contact_trace::trace_tangent_curve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{PortableMath, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

/// One face of a task: its surface, window, body and handle.
pub(super) type FaceSide<'a> = (
    &'a SurfaceGeometry,
    ((f64, f64), (f64, f64)),
    &'a BRep,
    FaceHandle,
);

/// The contact curves of one task: each tangent-curve contact within
/// `reach` traced once (a contact on a curve already traced is skipped),
/// clipped to both faces, its line or circle runs kept.
pub(super) fn contact_curves(
    contacts: &[SurfaceContact],
    (sa, wa, a, fa): FaceSide<'_>,
    (sb, wb, b, fb): FaceSide<'_>,
    scale: ContactScale,
    reach: f64,
) -> Vec<CurveOnEdge> {
    let mut traced: Vec<Vec<Point3>> = Vec::new();
    let mut out = Vec::new();
    for c in contacts {
        if c.gap.abs() > reach
            || !matches!(
                c.classification.kind,
                SurfaceContactKind::TangentCurve { .. }
            )
            || traced.iter().any(|t| {
                t.iter()
                    .any(|p| p.distance_to(&c.point) <= 2.0 * step_of(scale))
            })
        {
            continue;
        }
        let Some(trace) = trace_tangent_curve((sa, wa), (sb, wb), c, scale) else {
            continue;
        };
        let on = |p: &Point3| {
            let on_face = |brep: &BRep, f: FaceHandle| {
                brep.face_interior_holds(f, p, reach) || brep.point_near_face_boundary(f, p, reach)
            };
            on_face(a, fa) && on_face(b, fb)
        };
        let whole = trace.closed && trace.points.iter().all(on);
        if whole {
            out.extend(circle_through(&trace.points, true, reach));
        } else {
            for run in runs(&trace.points, &on) {
                out.extend(
                    line_through(&run, reach).or_else(|| circle_through(&run, false, reach)),
                );
            }
        }
        traced.push(trace.points);
    }
    out
}

/// The trace's step (`trace_tangent_curve`'s), for the skip test.
fn step_of(scale: ContactScale) -> f64 {
    (scale.extent / 64.0).max(10.0 * scale.linear)
}

/// The maximal runs of `points` (two or more) for which `on` holds.
fn runs(points: &[Point3], on: &impl Fn(&Point3) -> bool) -> Vec<Vec<Point3>> {
    let mut out: Vec<Vec<Point3>> = Vec::new();
    let mut cur: Vec<Point3> = Vec::new();
    for p in points {
        if on(p) {
            cur.push(*p);
        } else if !cur.is_empty() {
            out.push(std::mem::take(&mut cur));
        }
    }
    out.push(cur);
    out.retain(|r| r.len() >= 2);
    out
}

/// The run as a segment between its ends, when every point lies within
/// `reach` of it.
fn line_through(run: &[Point3], reach: f64) -> Option<CurveOnEdge> {
    let (p0, p1) = (*run.first()?, *run.last()?);
    let d = (p1 - p0).normalized()?;
    run.iter()
        .all(|p| {
            let w = *p - p0;
            (w - d * w.dot(&d)).length() <= reach
        })
        .then(|| CurveOnEdge::line_segment(p0, p1, reach))
        .flatten()
}

/// The run as an arc of the circle through its first, middle and last
/// points (the whole circle when `closed`), when every point lies within
/// `reach` of it.
fn circle_through(run: &[Point3], closed: bool, reach: f64) -> Option<CurveOnEdge> {
    if run.len() < 3 {
        return None;
    }
    let (a, b, c) = (run[0], run[run.len() / 2], run[run.len() - 1]);
    let (ab, ac) = (b - a, c - a);
    let normal = ab.cross(&ac);
    let n2 = normal.dot(&normal);
    if n2 <= 1e-24 {
        return None;
    }
    // Circumcentre of a, b, c.
    let to_centre =
        (ac.cross(&normal) * ab.dot(&ab) + normal.cross(&ab) * ac.dot(&ac)) * (0.5 / n2);
    let centre = a + to_centre;
    let radius = to_centre.length();
    let n = normal.normalized()?;
    let fits = run.iter().all(|p| {
        let w = *p - centre;
        let off_plane = w.dot(&n);
        let in_plane = (w - n * off_plane).length();
        off_plane.abs() <= reach && (in_plane - radius).abs() <= reach
    });
    if !fits {
        return None;
    }
    if closed {
        return Some(CurveOnEdge {
            geometry: CurveGeometry::Circle {
                center: centre,
                normal: n,
                radius,
            },
            t_start: 0.0,
            t_end: std::f64::consts::TAU,
        });
    }
    let angle = |p: &Point3| -> f64 { circle_angle(centre, n, p) };
    let (t0, t_mid, t1) = (angle(&a), angle(&b), angle(&c));
    // The arc from a through b to c, its parameters increasing.
    let span = |from: f64, to: f64| {
        popkern_core::geometry::chart::ops::forward_offset(to, from, std::f64::consts::TAU)
    };
    let (start, end) = if span(t0, t_mid) <= span(t0, t1) {
        (t0, t0 + span(t0, t1))
    } else {
        (t1, t1 + span(t1, t0))
    };
    Some(CurveOnEdge {
        geometry: CurveGeometry::Circle {
            center: centre,
            normal: n,
            radius,
        },
        t_start: start,
        t_end: end,
    })
}

/// The circle parameter of `p` on the circle about `centre` with normal
/// `n`, read as `CurveGeometry::Circle` reads it.
fn circle_angle(centre: Point3, n: Vec3, p: &Point3) -> f64 {
    let probe = CurveGeometry::Circle {
        center: centre,
        normal: n,
        radius: 1.0,
    };
    use popkern_core::geometry::curve::CurveEval;
    let x = probe.point_at(0.0) - centre;
    let y = probe.point_at(std::f64::consts::FRAC_PI_2) - centre;
    let w = *p - centre;
    w.dot(&y).patan2(w.dot(&x))
}
