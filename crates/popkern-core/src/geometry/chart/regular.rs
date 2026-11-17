//! Regular charts (design §4): around any point, a chart with no cut and
//! no singular row within a radius.

use popkern_math::PortableMath;

use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};

use super::descriptor::{Axis, ChartDescriptor, SingularKind, SingularRow};

/// Why a regular chart could not be used.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ChartError {
    /// The pole projection's inverse did not converge.
    #[error("pole chart inverse did not converge at ({x:e}, {y:e}): residual {residual:e}")]
    InverseDiverged {
        /// Local chart coordinates.
        x: f64,
        /// Local chart coordinates.
        y: f64,
        /// The final residual in the tangent plane.
        residual: f64,
    },
}

/// A chart free of artifacts around a point.
#[derive(Debug, Clone, PartialEq)]
pub enum RegularChart {
    /// The face chart itself: no artifact within the radius.
    Face,
    /// The face chart re-centred at `centre` along its closed directions:
    /// local coordinates are the parameters' representatives within half a
    /// period of the centre, so the cut lies half a period away.
    Shifted {
        /// The chart's centre.
        centre: (f64, f64),
        /// The periods along u and v (when closed).
        periods: (Option<f64>, Option<f64>),
        /// The domain (to evaluate in range).
        domain: ((f64, f64), (f64, f64)),
    },
    /// Near a chart singularity: coordinates are the orthogonal projection
    /// onto the tangent plane at the singular point.
    PoleProjection {
        /// The singular point.
        centre: Point3,
        /// An orthonormal tangent frame there.
        e1: Vec3,
        /// An orthonormal tangent frame there.
        e2: Vec3,
        /// The singular row.
        row: SingularRow,
    },
    /// A geometric singularity lies within the radius: a topological
    /// event the consumer must treat as one.
    Singular(SingularKind, Point3),
}

fn frame(n: Vec3) -> (Vec3, Vec3) {
    let a = if n.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let e1 = n.cross(&a).normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    (e1, n.cross(&e1))
}

/// The normal at a singular row's point, read a small step off the row.
fn row_normal(s: &SurfaceGeometry, row: &SingularRow, d: &ChartDescriptor) -> Option<Vec3> {
    let ((u0, u1), (v0, v1)) = d.domain;
    let (lo, hi, a, b) = match row.axis {
        Axis::V => (v0, v1, u0, u1),
        Axis::U => (u0, u1, v0, v1),
    };
    let span = if lo.is_finite() && hi.is_finite() {
        hi - lo
    } else {
        1.0
    };
    let inward = if (row.at - lo).abs() <= (row.at - hi).abs() {
        1.0
    } else {
        -1.0
    };
    let t = if a.is_finite() && b.is_finite() {
        0.5 * (a + b)
    } else {
        0.0
    };
    let (u, v) = match row.axis {
        Axis::V => (t, row.at + inward * 1e-7 * span),
        Axis::U => (row.at + inward * 1e-7 * span, t),
    };
    let (su, sv) = s.first_partials(u, v);
    let n = su.cross(&sv).normalized()?;
    // Point it away from the row's neighbourhood's centroid side: a
    // chart singularity's two sides agree, so either sign is consistent.
    Some(n)
}

/// The regular chart for a point at parameters `at` (design §4): a pole
/// projection within `radius` of a chart singularity, `Singular` within it
/// of a geometric one, `Shifted` when a cut lies within it, else `Face`.
pub fn regular_chart(
    d: &ChartDescriptor,
    s: &SurfaceGeometry,
    at: (f64, f64),
    radius: f64,
) -> RegularChart {
    let p = s.point_at(at.0, at.1);
    for row in &d.singular_rows {
        if row.point.distance_to(&p) <= radius {
            return match row.kind {
                SingularKind::Chart => match row_normal(s, row, d) {
                    Some(n) => {
                        let (e1, e2) = frame(n);
                        RegularChart::PoleProjection {
                            centre: row.point,
                            e1,
                            e2,
                            row: *row,
                        }
                    }
                    None => RegularChart::Singular(row.kind, row.point),
                },
                kind => RegularChart::Singular(kind, row.point),
            };
        }
    }
    let near_cut = |period: Option<super::descriptor::ChartPeriod>, value: f64, axis: Axis| {
        period.is_some_and(|c| {
            let ((u0, u1), (v0, v1)) = d.domain;
            let (lo, hi) = match axis {
                Axis::U => (u0, u1),
                Axis::V => (v0, v1),
            };
            // The cut's distance in space from the point.
            let (a, b) = match axis {
                Axis::U => (s.point_at(lo, at.1), s.point_at(hi, at.1)),
                Axis::V => (s.point_at(at.0, lo), s.point_at(at.0, hi)),
            };
            let _ = (value, c);
            a.distance_to(&p).min(b.distance_to(&p)) <= radius
        })
    };
    if near_cut(d.u, at.0, Axis::U) || near_cut(d.v, at.1, Axis::V) {
        return RegularChart::Shifted {
            centre: at,
            periods: (d.u.map(|c| c.period), d.v.map(|c| c.period)),
            domain: d.domain,
        };
    }
    RegularChart::Face
}

/// `t`'s representative within half a period of `centre`.
fn near(t: f64, centre: f64, period: Option<f64>) -> f64 {
    match period {
        Some(p) => centre + (t - centre + 0.5 * p).rem_euclid(p) - 0.5 * p,
        None => t,
    }
}

/// `t`'s representative in `[lo, hi]`.
fn in_domain(t: f64, (lo, hi): (f64, f64), period: Option<f64>) -> f64 {
    match period {
        Some(p) if lo.is_finite() && !(lo..=hi).contains(&t) => lo + (t - lo).rem_euclid(p),
        _ => t,
    }
}

impl RegularChart {
    /// Local coordinates of the face parameters `uv`.
    pub fn from_surface(&self, s: &SurfaceGeometry, uv: (f64, f64)) -> (f64, f64) {
        match self {
            RegularChart::Face | RegularChart::Singular(..) => uv,
            RegularChart::Shifted {
                centre, periods, ..
            } => (
                near(uv.0, centre.0, periods.0),
                near(uv.1, centre.1, periods.1),
            ),
            RegularChart::PoleProjection { centre, e1, e2, .. } => {
                let p = s.point_at(uv.0, uv.1);
                let d = p - *centre;
                (d.dot(e1), d.dot(e2))
            }
        }
    }

    /// Face parameters (in the domain) of local coordinates `xy`.
    pub fn to_surface(
        &self,
        s: &SurfaceGeometry,
        xy: (f64, f64),
    ) -> Result<(f64, f64), ChartError> {
        match self {
            RegularChart::Face | RegularChart::Singular(..) => Ok(xy),
            RegularChart::Shifted {
                periods, domain, ..
            } => Ok((
                in_domain(xy.0, domain.0, periods.0),
                in_domain(xy.1, domain.1, periods.1),
            )),
            RegularChart::PoleProjection { centre, e1, e2, .. } => {
                // Newton on the tangent-plane projection, seeded by the
                // surface's own closest point to the planar point.
                let target = *centre + *e1 * xy.0 + *e2 * xy.1;
                let (mut u, mut v) =
                    s.closest_point(&target)
                        .ok_or(ChartError::InverseDiverged {
                            x: xy.0,
                            y: xy.1,
                            residual: f64::INFINITY,
                        })?;
                let mut residual = f64::INFINITY;
                for _ in 0..30 {
                    let d = s.point_at(u, v) - *centre;
                    let (rx, ry) = (d.dot(e1) - xy.0, d.dot(e2) - xy.1);
                    residual = rx.phypot(ry);
                    let scale = (xy.0.phypot(xy.1))
                        .max(centre.x.abs().max(centre.y.abs()).max(centre.z.abs()))
                        .max(1.0);
                    if residual <= 1e-13 * scale {
                        return Ok((u, v));
                    }
                    let (su, sv) = s.first_partials(u, v);
                    let (a, b, c, e) = (su.dot(e1), sv.dot(e1), su.dot(e2), sv.dot(e2));
                    let det = a * e - b * c;
                    if det.abs() < 1e-300 {
                        break;
                    }
                    u -= (e * rx - b * ry) / det;
                    v -= (a * ry - c * rx) / det;
                }
                Err(ChartError::InverseDiverged {
                    x: xy.0,
                    y: xy.1,
                    residual,
                })
            }
        }
    }
}
