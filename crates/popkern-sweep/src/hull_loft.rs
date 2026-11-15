//! The loft between two convex sections of different side counts in
//! parallel planes: their convex hull (PK-1197).
//!
//! Two sections that split into the same number of spans loft span by
//! span; two that do not (a square and a circle, a square and an
//! octagon) went to one skinned surface fitted through both, which
//! paired their points by chord length: a square-to-round read 15379
//! against the transition's 35069, a square-to-octagon 26409 against
//! 33760. The transition a sheet-metal shop makes, and the one the
//! sections determine, is their convex hull: every wall a straight
//! ruling joining the points of the two sections whose outward normals
//! point the same way (the support points in one direction). A side of
//! one section facing a corner of the other makes a flat triangle, two
//! parallel sides a flat trapezoid, a corner facing an arc a piece of an
//! oblique cone with its apex at the corner, and two arcs facing each
//! other a ruled patch between them (a cone, or a cylinder at equal
//! radii); the walls meet tangentially across their shared rulings.
//!
//! Construction is exact. Each section is read as its spans, each a
//! line or a circular arc recognised within the linear tolerance, and
//! must be convex and planar, the planes parallel; each feature (a
//! line, an arc, a corner between spans) covers its range of outward
//! normal angles in the first section's circle frame. The two ranges
//! merge into the cycle of rulings, one at each side of every event
//! angle, and between consecutive rulings lies one wall, its surface a
//! plane or the rational ruled surface between the arc pieces (both
//! written over the same angles, so their control points pair). The
//! body is assembled face by face, its pcurves completed, faced out of
//! its volume and validated; any section outside this class keeps the
//! skinned loft.

use popkern_core::constructors::BRepBuilder;
use popkern_core::geometry::analytic_curves::circle_basis;
use popkern_core::geometry::circle_recognition::{CircleFit, arc_between, recognize_circular};
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, PortableMath, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};

use crate::types::LoftSection;

/// Angles within this of one another are one event.
const ANGLE_EPS: f64 = 1e-9;
/// Samples per span for the planarity and convexity tests.
const SPAN_SAMPLES: usize = 16;

/// A section's feature over its normal angles, in its plane's frame.
#[derive(Debug, Clone, Copy)]
enum Feature {
    /// A straight side from `p` to `q` (counter-clockwise), normal `at`.
    Line { p: Point3, q: Point3, at: f64 },
    /// An arc about `centre` of `radius` from normal angle `from` to `to`.
    Arc {
        centre: Point3,
        radius: f64,
        from: f64,
        to: f64,
    },
    /// A corner at `v` turning the normal from `from` to `to`.
    Corner { v: Point3, from: f64, to: f64 },
}

/// One section read as features, counter-clockwise about the frame
/// normal, angles increasing over one turn.
struct Section {
    features: Vec<Feature>,
    frame: (Vec3, Vec3, Vec3),
}

impl Section {
    /// The point in direction `angle`, just before it (`after` false) or
    /// just after it: a side facing `angle` exactly answers its start or
    /// its end.
    fn support(&self, angle: f64, after: bool) -> Point3 {
        let wrap = |a: f64| a - std::f64::consts::TAU * (a / std::f64::consts::TAU).floor();
        let a = wrap(angle);
        let within = |from: f64, to: f64| {
            let (f, t) = (wrap(from), wrap(from) + (to - from));
            (a >= f - ANGLE_EPS && a <= t + ANGLE_EPS)
                || (a + std::f64::consts::TAU >= f - ANGLE_EPS
                    && a + std::f64::consts::TAU <= t + ANGLE_EPS)
        };
        for f in &self.features {
            if let Feature::Line { p, q, at } = *f
                && (wrap(at) - a)
                    .abs()
                    .min(std::f64::consts::TAU - (wrap(at) - a).abs())
                    <= ANGLE_EPS
            {
                return if after { q } else { p };
            }
        }
        for f in &self.features {
            match *f {
                Feature::Arc {
                    centre,
                    radius,
                    from,
                    to,
                } if within(from, to) => {
                    let (e1, e2, _) = self.frame;
                    return centre + (e1 * angle.pcos() + e2 * angle.psin()) * radius;
                }
                Feature::Corner { v, from, to } if within(from, to) => return v,
                _ => {}
            }
        }
        // Angles are covered by construction; the nearest corner stands.
        self.features
            .iter()
            .find_map(|f| match f {
                Feature::Corner { v, .. } => Some(*v),
                _ => None,
            })
            .unwrap_or(Point3::new(0.0, 0.0, 0.0))
    }

    /// The arc covering the open interval round `mid`, if any.
    fn arc_at(&self, mid: f64) -> Option<(Point3, f64)> {
        let wrap = |a: f64| a - std::f64::consts::TAU * (a / std::f64::consts::TAU).floor();
        self.features.iter().find_map(|f| match *f {
            Feature::Arc {
                centre,
                radius,
                from,
                to,
            } => {
                let (f0, m) = (wrap(from), wrap(mid));
                let span = to - from;
                let inside = (m > f0 && m < f0 + span) || (m + std::f64::consts::TAU < f0 + span);
                inside.then_some((centre, radius))
            }
            _ => None,
        })
    }

    /// Every angle at which a feature starts or ends.
    fn events(&self) -> Vec<f64> {
        let mut out = Vec::new();
        for f in &self.features {
            match *f {
                Feature::Line { at, .. } => out.push(at),
                Feature::Arc { from, to, .. } | Feature::Corner { from, to, .. } => {
                    out.push(from);
                    out.push(to);
                }
            }
        }
        out
    }
}

/// The hull loft of `sections`, when they are two convex planar
/// sections of lines and arcs in parallel planes whose span structures
/// differ; `None` keeps the skinned loft.
pub(crate) fn hull_loft(sections: &[LoftSection], tolerance: Tolerance) -> Option<BRep> {
    let [a, b] = sections else {
        return None;
    };
    if a.position.is_some() || b.position.is_some() {
        return None;
    }
    let lin = tolerance.linear();
    let normal = plane_normal(&a.curve, lin)?;
    let nb = plane_normal(&b.curve, lin)?;
    if normal.cross(&nb).length() > 1e-9 {
        return None;
    }
    let (e1, e2) = circle_basis(normal);
    let frame = (e1, e2, normal);
    let sa = read_section(&a.curve, frame, lin)?;
    let sb = read_section(&b.curve, frame, lin)?;
    let lift = (b.curve.point_at(b.curve.parameter_range().0)
        - a.curve.point_at(a.curve.parameter_range().0))
    .dot(&normal);
    if lift.abs() <= lin {
        return None;
    }
    let mut events: Vec<f64> = sa
        .events()
        .into_iter()
        .chain(sb.events())
        .map(|e| e - std::f64::consts::TAU * (e / std::f64::consts::TAU).floor())
        .collect();
    events.sort_by(f64::total_cmp);
    events.dedup_by(|x, y| (*x - *y).abs() <= ANGLE_EPS);
    if events.len() > 1
        && (events[0] + std::f64::consts::TAU - events[events.len() - 1]).abs() <= ANGLE_EPS
    {
        events.pop();
    }
    // The cycle of rulings: one before and one after each event.
    let mut rulings: Vec<(Point3, Point3, f64)> = Vec::new();
    for &e in &events {
        for after in [false, true] {
            let r = (sa.support(e, after), sb.support(e, after), e);
            let same = rulings
                .last()
                .is_some_and(|l| l.0.distance_to(&r.0) <= lin && l.1.distance_to(&r.1) <= lin);
            if !same {
                rulings.push(r);
            }
        }
    }
    if rulings.len() > 1 {
        let (f, l) = (rulings[0], rulings[rulings.len() - 1]);
        if f.0.distance_to(&l.0) <= lin && f.1.distance_to(&l.1) <= lin {
            rulings.pop();
        }
    }
    if rulings.len() < 3 {
        return None;
    }
    assemble(&rulings, &sa, &sb, normal, tolerance)
}

/// The unit normal of a planar closed curve (Newell), `None` when the
/// curve leaves its plane.
fn plane_normal(c: &NurbsCurve, lin: f64) -> Option<Vec3> {
    let (t0, t1) = c.parameter_range();
    let n = 256;
    let pts: Vec<Point3> = (0..n)
        .map(|k| c.point_at(t0 + (t1 - t0) * k as f64 / n as f64))
        .collect();
    let mut acc = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..n {
        let (p, q) = (pts[i], pts[(i + 1) % n]);
        acc = acc
            + Vec3::new(
                (p.y - q.y) * (p.z + q.z),
                (p.z - q.z) * (p.x + q.x),
                (p.x - q.x) * (p.y + q.y),
            );
    }
    let normal = acc.normalized()?;
    let o = pts[0];
    pts.iter()
        .all(|p| (*p - o).dot(&normal).abs() <= lin)
        .then_some(normal)
}

/// `curve` as features counter-clockwise about `frame.2`; `None` when a
/// span is neither a line nor a circular arc or the section is not
/// convex.
fn read_section(curve: &NurbsCurve, frame: (Vec3, Vec3, Vec3), lin: f64) -> Option<Section> {
    let (e1, e2, n) = frame;
    let mut spans = crate::profile::detect_and_split(curve).ok()?;
    let ccw = {
        let (t0, t1) = curve.parameter_range();
        let m = 256;
        let pts: Vec<Point3> = (0..m)
            .map(|k| curve.point_at(t0 + (t1 - t0) * k as f64 / m as f64))
            .collect();
        (0..m)
            .map(|i| {
                (pts[i] - pts[0])
                    .cross(&(pts[(i + 1) % m] - pts[0]))
                    .dot(&n)
            })
            .sum::<f64>()
            > 0.0
    };
    if !ccw {
        spans.reverse();
        for s in &mut spans {
            *s = reversed(s);
        }
    }
    let angle = |v: Vec3| v.dot(&e2).patan2(v.dot(&e1));
    // Each span's start and end normal, and its feature.
    let mut parts: Vec<(f64, f64, Feature)> = Vec::new();
    for s in &spans {
        let (t0, t1) = s.parameter_range();
        let (p, q) = (s.point_at(t0), s.point_at(t1));
        let samples: Vec<Point3> = (0..=SPAN_SAMPLES)
            .map(|k| s.point_at(t0 + (t1 - t0) * k as f64 / SPAN_SAMPLES as f64))
            .collect();
        // A closed span (a whole circle) has no chord and is no line.
        let chord = (q - p).normalized();
        let straight = chord.is_some_and(|chord| {
            samples.iter().all(|x| {
                let d = *x - p;
                (d - chord * d.dot(&chord)).length() <= lin
            })
        });
        if let (true, Some(chord)) = (straight, chord) {
            let out = chord.cross(&n);
            let at = angle(out);
            parts.push((at, at, Feature::Line { p, q, at }));
            continue;
        }
        let fit: CircleFit = recognize_circular(s, (t0, t1), lin)?;
        let from = angle(p - fit.center);
        let mut to = angle(q - fit.center);
        // Counter-clockwise about the frame normal: the arc's own
        // normal must agree (a convex section's arcs bulge outward).
        if fit.normal.dot(&n) <= 0.0 {
            return None;
        }
        while to <= from + ANGLE_EPS {
            to += std::f64::consts::TAU;
        }
        parts.push((
            from,
            to,
            Feature::Arc {
                centre: fit.center,
                radius: fit.radius,
                from,
                to,
            },
        ));
    }
    // Corners between spans: the normal turns forward by less than a
    // half turn at each (convexity), the turns summing to one turn.
    let mut features = Vec::new();
    let mut total = 0.0;
    for i in 0..parts.len() {
        let (_, end, f) = parts[i];
        let (start_next, _, _) = parts[(i + 1) % parts.len()];
        features.push(f);
        total += match f {
            Feature::Arc { from, to, .. } => to - from,
            _ => 0.0,
        };
        let mut turn = start_next - end;
        turn -= std::f64::consts::TAU * (turn / std::f64::consts::TAU).floor();
        if turn > std::f64::consts::TAU - ANGLE_EPS {
            turn = 0.0;
        }
        if turn >= std::f64::consts::PI {
            return None;
        }
        if turn > ANGLE_EPS {
            let v = match f {
                Feature::Line { q, .. } => q,
                Feature::Arc {
                    centre, radius, to, ..
                } => centre + (e1 * to.pcos() + e2 * to.psin()) * radius,
                Feature::Corner { v, .. } => v,
            };
            features.push(Feature::Corner {
                v,
                from: end,
                to: end + turn,
            });
            total += turn;
        }
    }
    ((total - std::f64::consts::TAU).abs() <= 1e-6).then_some(Section { features, frame })
}

/// `c` traversed the other way.
fn reversed(c: &NurbsCurve) -> NurbsCurve {
    let (lo, hi) = (c.knots[0], c.knots[c.knots.len() - 1]);
    NurbsCurve {
        degree: c.degree,
        knots: c.knots.iter().rev().map(|k| lo + hi - k).collect(),
        control_points: c.control_points.iter().rev().copied().collect(),
        weights: c.weights.iter().rev().copied().collect(),
    }
}

#[path = "hull_loft_assemble.rs"]
mod hull_loft_assemble;
use hull_loft_assemble::assemble;
