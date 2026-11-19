//! The chart descriptor (design §3): a surface's periods with their kind
//! and seam continuity, and its singular rows with their kind (chart or
//! geometric), measured from the surface, cached by content.

use popkern_math::PortableMath;
use std::cell::RefCell;
use std::rc::Rc;

use popkern_math::hash::FastHashMap as HashMap;

use crate::geometry::Point3;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

/// How a chart direction closes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PeriodKind {
    /// An analytic surface's periodicity (cylinder, cone, sphere, torus).
    Analytic,
    /// A NURBS stored in periodic form along the direction.
    PeriodicForm,
    /// A clamped NURBS whose end rows coincide: closed, not periodic.
    GeometricClosure,
}

/// The order of continuity across a seam.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Continuity {
    /// Position only: the seam is a crease of the chart.
    G0,
    /// Tangent planes agree across the seam.
    G1,
    /// Smooth (analytic, or a periodic form with simple seam knots).
    Smooth,
}

/// A closed direction of the chart.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChartPeriod {
    /// The period.
    pub period: f64,
    /// Where the chart's domain starts along the direction (the cut).
    pub origin: f64,
    /// How the direction closes.
    pub kind: PeriodKind,
    /// Continuity across the cut.
    pub continuity: Continuity,
}

/// The parameter direction a singular row runs along.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    /// The row is a line `u = at`.
    U,
    /// The row is a line `v = at`.
    V,
}

/// What kind of point a singular row maps to (design §2).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SingularKind {
    /// The surface is smooth there; only the chart degenerates (a
    /// sphere's pole). Removed by re-charting.
    Chart,
    /// The surface is singular there (a cone's apex: its normals sweep a
    /// cone about the point). A topological event, never smoothed.
    Apex,
}

/// A chart row that maps to one point.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SingularRow {
    /// The row is the line `axis = at`.
    pub axis: Axis,
    /// The row's parameter.
    pub at: f64,
    /// The point it maps to.
    pub point: Point3,
    /// Chart or geometric singularity, as measured.
    pub kind: SingularKind,
    /// The distortion radius (design §4): within this distance of the
    /// point the chart's speed across the row exceeds ten times its speed
    /// along it, and interpolation runs in the pole projection.
    pub distortion: f64,
}

/// A surface's chart artifacts.
#[derive(Debug, Clone, PartialEq)]
pub struct ChartDescriptor {
    /// The u direction's period, when the chart closes along u.
    pub u: Option<ChartPeriod>,
    /// The v direction's period, when the chart closes along v.
    pub v: Option<ChartPeriod>,
    /// The chart's singular rows.
    pub singular_rows: Vec<SingularRow>,
    /// The chart's domain (infinite along an unbounded direction).
    pub domain: ((f64, f64), (f64, f64)),
    /// The linear tolerance it was measured under: the reach of "at an
    /// artifact" for consumers asking in space.
    pub linear: f64,
}

fn size(s: &SurfaceGeometry, ((u0, u1), (v0, v1)): ((f64, f64), (f64, f64))) -> f64 {
    let clamp = |a: f64, b: f64| {
        if a.is_finite() && b.is_finite() {
            (a, b)
        } else {
            (-1.0, 1.0)
        }
    };
    let ((u0, u1), (v0, v1)) = (clamp(u0, u1), clamp(v0, v1));
    let pts = [
        s.point_at(u0, v0),
        s.point_at(u1, v1),
        s.point_at(0.5 * (u0 + u1), 0.5 * (v0 + v1)),
        s.point_at(u0, v1),
    ];
    pts.iter()
        .flat_map(|a| pts.iter().map(move |b| a.distance_to(b)))
        .fold(0.0, f64::max)
        .max(1e-300)
}

fn period_kind(s: &SurfaceGeometry, axis: Axis) -> PeriodKind {
    // A NURBS first: its `periodicity` reports its periodic form, which is
    // not an analytic period.
    match s {
        SurfaceGeometry::Nurbs(n) => {
            let form = match axis {
                Axis::U => n.u_period(),
                Axis::V => n.v_period(),
            };
            if form.is_some() {
                PeriodKind::PeriodicForm
            } else {
                PeriodKind::GeometricClosure
            }
        }
        _ => PeriodKind::Analytic,
    }
}

/// Continuity across the cut of a closed direction, measured: the unit
/// normals on either side at several stations along the cut.
fn seam_continuity(
    s: &SurfaceGeometry,
    axis: Axis,
    kind: PeriodKind,
    domain: ((f64, f64), (f64, f64)),
    angular: f64,
) -> Continuity {
    if kind == PeriodKind::Analytic {
        return Continuity::Smooth;
    }
    let ((u0, u1), (v0, v1)) = domain;
    let side = |at: f64, other: f64| -> Option<crate::geometry::Vec3> {
        let (u, v) = match axis {
            Axis::U => (at, other),
            Axis::V => (other, at),
        };
        let (su, sv) = s.first_partials(u, v);
        su.cross(&sv).normalized()
    };
    let (lo, hi, a, b) = match axis {
        Axis::U => (u0, u1, v0, v1),
        Axis::V => (v0, v1, u0, u1),
    };
    let mut worst = 0.0f64;
    for k in 1..8 {
        let t = a + (b - a) * k as f64 / 8.0;
        if let (Some(n0), Some(n1)) = (side(lo, t), side(hi, t)) {
            worst = worst.max(n0.dot(&n1).clamp(-1.0, 1.0).pacos());
        }
    }
    if worst > angular {
        Continuity::G0
    } else if kind == PeriodKind::PeriodicForm {
        Continuity::Smooth
    } else {
        Continuity::G1
    }
}

/// The singular rows: domain-boundary rows whose samples coincide,
/// classified by the spread of the unit normals on a ring about the point.
fn singular_rows(
    s: &SurfaceGeometry,
    domain: ((f64, f64), (f64, f64)),
    size: f64,
    tol: &Tolerance,
) -> Vec<SingularRow> {
    let ((u0, u1), (v0, v1)) = domain;
    let mut out = Vec::new();
    let finite = |a: f64, b: f64| a.is_finite() && b.is_finite();
    let mut candidates = Vec::new();
    if finite(u0, u1) && finite(v0, v1) {
        for at in [v0, v1] {
            candidates.push((Axis::V, at, (u0, u1)));
        }
        for at in [u0, u1] {
            candidates.push((Axis::U, at, (v0, v1)));
        }
    } else if finite(v0, v1) {
        // A cone's v is unbounded along the axis but its apex row is v = 0
        // when the chart is the cone's: probed below as a candidate row.
    }
    if let SurfaceGeometry::Cone { .. } = s {
        candidates.push((
            Axis::V,
            0.0,
            if finite(u0, u1) {
                (u0, u1)
            } else {
                (0.0, std::f64::consts::TAU)
            },
        ));
    }
    // A spindle torus (tube radius at least the major radius) collapses
    // its rows where the distance from the axis, R + r cos v, vanishes:
    // its two pinch points, apexes the tube passes through.
    if let SurfaceGeometry::Torus {
        major_radius,
        minor_radius,
        ..
    } = s
        && *minor_radius > 0.0
        && *minor_radius >= *major_radius
    {
        let v = (-major_radius / minor_radius).clamp(-1.0, 1.0).pacos();
        for at in [v, std::f64::consts::TAU - v] {
            candidates.push((Axis::V, at, (0.0, std::f64::consts::TAU)));
        }
    }
    for (axis, at, (a, b)) in candidates {
        let eval = |t: f64| match axis {
            Axis::V => s.point_at(t, at),
            Axis::U => s.point_at(at, t),
        };
        let p0 = eval(a + 0.11 * (b - a));
        let collapsed =
            (1..6).all(|k| eval(a + (b - a) * k as f64 / 6.0).distance_to(&p0) <= tol.linear());
        if !collapsed
            || out
                .iter()
                .any(|r: &SingularRow| r.point.distance_to(&p0) <= tol.linear())
        {
            continue;
        }
        out.push(SingularRow {
            axis,
            at,
            point: p0,
            kind: classify(s, axis, at, (a, b), domain, size, tol),
            distortion: 0.0,
        });
    }
    out
}

/// Chart or apex: the unit normals a small step off the row, around it.
fn classify(
    s: &SurfaceGeometry,
    axis: Axis,
    at: f64,
    (a, b): (f64, f64),
    domain: ((f64, f64), (f64, f64)),
    size: f64,
    tol: &Tolerance,
) -> SingularKind {
    let ((u0, u1), (v0, v1)) = domain;
    let (lo, hi) = match axis {
        Axis::V => (v0, v1),
        Axis::U => (u0, u1),
    };
    let span = if lo.is_finite() && hi.is_finite() {
        hi - lo
    } else {
        size
    };
    let inward = if (at - lo).abs() <= (at - hi).abs() {
        1.0
    } else {
        -1.0
    };
    let step = inward * 1e-4 * span;
    let normals: Vec<crate::geometry::Vec3> = (0..12)
        .filter_map(|k| {
            let t = a + (b - a) * (k as f64 + 0.5) / 12.0;
            let (u, v) = match axis {
                Axis::V => (t, at + step),
                Axis::U => (at + step, t),
            };
            let (su, sv) = s.first_partials(u, v);
            su.cross(&sv).normalized()
        })
        .collect();
    let Some(first) = normals.first() else {
        return SingularKind::Apex;
    };
    let spread = normals
        .iter()
        .map(|n| n.dot(first).clamp(-1.0, 1.0).pacos())
        .fold(0.0, f64::max);
    // A smooth surface's normals a small step off the point agree to the
    // step's own curvature turn; an apex's sweep its half-angle.
    if spread <= (10.0 * tol.angular()).max(1e-2) {
        SingularKind::Chart
    } else {
        SingularKind::Apex
    }
}

type Cache = RefCell<HashMap<(u64, u64, u64, u64), Rc<ChartDescriptor>>>;

thread_local! {
    static CACHE: Cache = RefCell::new(HashMap::default());
}

impl ChartDescriptor {
    /// The descriptor of `surface` under `tolerance`, measured once per
    /// surface content (NURBS are cached by their content fingerprint).
    pub fn of(surface: &SurfaceGeometry, tolerance: &Tolerance) -> Rc<ChartDescriptor> {
        let key = match surface {
            SurfaceGeometry::Nurbs(n) => {
                let (w, a, b) = n.content_fingerprint();
                Some((w, a, b, tolerance.linear().to_bits()))
            }
            _ => analytic_key(surface).map(|(w, a, b)| (w, a, b, tolerance.linear().to_bits())),
        };
        if let Some(k) = key
            && let Some(hit) = CACHE.with(|c| c.borrow().get(&k).cloned())
        {
            return hit;
        }
        let built = Rc::new(Self::measure(surface, tolerance));
        if let Some(k) = key {
            CACHE.with(|c| {
                let mut map = c.borrow_mut();
                if map.len() >= 256 {
                    map.clear();
                }
                map.insert(k, built.clone());
            });
        }
        built
    }

    fn measure(s: &SurfaceGeometry, tol: &Tolerance) -> ChartDescriptor {
        let domain = s.parameter_range();
        let size = size(s, domain);
        let (pu, pv) = super::periods::marching_periods(s);
        let period = |p: Option<f64>, axis: Axis| {
            p.map(|period| {
                let kind = period_kind(s, axis);
                let origin = match axis {
                    Axis::U => domain.0.0,
                    Axis::V => domain.1.0,
                };
                ChartPeriod {
                    period,
                    origin,
                    kind,
                    continuity: seam_continuity(s, axis, kind, domain, tol.angular()),
                }
            })
        };
        let mut d = ChartDescriptor {
            u: period(pu, Axis::U),
            v: period(pv, Axis::V),
            singular_rows: singular_rows(s, domain, size, tol),
            domain,
            linear: tol.linear(),
        };
        let radii: Vec<f64> = d
            .singular_rows
            .iter()
            .map(|row| super::ops::distortion_radius(&d, s, row))
            .collect();
        for (row, r) in d.singular_rows.iter_mut().zip(radii) {
            row.distortion = r;
        }
        d
    }

    /// Every closed direction's period (analytic, periodic form and a
    /// clamped NURBS's geometric closure): the chart's periods for
    /// chart-space reasoning. Equal to the former `geometric_chart_periods`
    /// and `marching_periods`.
    pub fn periods_of(surface: &SurfaceGeometry) -> (Option<f64>, Option<f64>) {
        super::periods::marching_periods(surface)
    }

    /// The periods along which evaluation itself is periodic (analytic,
    /// periodic form; not a geometric closure, which evaluates beyond its
    /// domain by extrapolation). Equal to the former `parameter_periodicity`.
    pub fn evaluation_periods_of<S: SurfaceEval + ?Sized>(
        surface: &S,
    ) -> (Option<f64>, Option<f64>) {
        surface.parameter_periodicity()
    }

    /// Whether the chart closes along u and v (any kind): the former
    /// `marching_periodicity`.
    pub fn closed_of(surface: &SurfaceGeometry) -> (bool, bool) {
        let (u, v) = Self::periods_of(surface);
        (u.is_some(), v.is_some())
    }

    /// A NURBS's periodic-form periods (evaluation-periodic directions).
    pub fn nurbs_evaluation_periods(
        n: &crate::geometry::nurbs_surface::NurbsSurface,
    ) -> (Option<f64>, Option<f64>) {
        (n.u_period(), n.v_period())
    }

    /// Whether a NURBS is in periodic form along either direction.
    pub fn is_periodic_form(n: &crate::geometry::nurbs_surface::NurbsSurface) -> bool {
        let (pu, pv) = Self::nurbs_evaluation_periods(n);
        pu.or(pv).is_some()
    }

    /// Whether the chart has a singular row under `linear` (a NURBS row of
    /// control points within it counted collapsed): the former
    /// `has_singular_row`, by the surface type and the net.
    pub fn has_singular_row(surface: &SurfaceGeometry, linear: f64) -> bool {
        super::singularity::has_singular_row(surface, linear)
    }

    /// The parameter wrap the surface's own evaluation applies, read once.
    pub fn evaluation_wrap_of<S: SurfaceEval + ?Sized>(surface: &S) -> super::wrap::EvaluationWrap {
        surface.evaluation_wrap()
    }

    /// The period along `axis`, when the chart closes along it.
    pub fn period(&self, axis: Axis) -> Option<ChartPeriod> {
        match axis {
            Axis::U => self.u,
            Axis::V => self.v,
        }
    }
}

/// A content key for an analytic surface: its kind and every defining
/// value's bits, folded two ways (so a descriptor is measured once per
/// distinct surface, as a NURBS's is by its fingerprint). `None` for the
/// patch kinds, which are measured each time.
fn analytic_key(s: &SurfaceGeometry) -> Option<(u64, u64, u64)> {
    let p = |q: &Point3| [q.x, q.y, q.z];
    let v = |q: &crate::geometry::Vec3| [q.x, q.y, q.z];
    let (kind, values): (u64, Vec<f64>) = match s {
        SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        } => (1, [p(origin), v(normal), v(u_axis), v(v_axis)].concat()),
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => (2, [p(origin).as_slice(), &v(axis), &[*radius]].concat()),
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => (3, [p(apex).as_slice(), &v(axis), &[*half_angle]].concat()),
        SurfaceGeometry::Sphere {
            center,
            radius,
            axis,
        } => (4, [p(center).as_slice(), &v(axis), &[*radius]].concat()),
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => (
            5,
            [
                p(center).as_slice(),
                &v(axis),
                &[*major_radius, *minor_radius],
            ]
            .concat(),
        ),
        _ => return None,
    };
    let (mut a, mut b) = (0xcbf2_9ce4_8422_2325_u64, 0x8422_2325_cbf2_9ce4_u64);
    for x in values {
        let bits = x.to_bits();
        a = (a ^ bits).wrapping_mul(0x0000_0100_0000_01b3);
        b = b.rotate_left(17) ^ bits.wrapping_mul(0x9e37_79b9_7f4a_7c15);
    }
    Some((kind | 0x8000_0000_0000_0000, a, b))
}
