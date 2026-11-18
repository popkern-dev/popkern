//! Angular facet measures for the `max_facet_angle` bound (see
//! [`crate::TessellationOptions::max_facet_angle`]).
//!
//! Two measures, both derivative-free so they behave identically on
//! every curve and surface kind:
//!
//! - **Curve turn** of an interval from its end and midpoint samples:
//!   twice the angle between the two half-chords. On a circular arc
//!   this is exactly the arc's subtended angle (inscribed-angle
//!   theorem), independent of where the midpoint sits; on a general
//!   curve it is the second-order tangent-turn estimate.
//! - **Normal spread** of a chart cell: the largest angle between the
//!   unit normals at its four corners and centre — the facet-to-facet
//!   crease a flat triangle over that cell would carry.
//! - **Edge normal turn** ([`EdgeNormalProbe`]): how far the INCIDENT
//!   FACES' normals turn along an edge interval, read through the
//!   half-edge pcurves. This — not the edge's own tangent turn — is
//!   what a boundary facet must resolve: a fillet's end arc turns the
//!   fillet normal 90°, while a boolean intersection curve can carry
//!   a sub-chord-height hook whose tangent turns 40° over smooth
//!   surfaces (measured: threaded_stud's thread/shank seam, 2.8e-3
//!   off the flank at its runout). Chasing that tangent turn sampled
//!   INTO the hook, the off-surface gate declined the flank, and its
//!   50-turn fallback grid meshed 468k triangles with 2870 open
//!   edges.
//! - **Chord-foot error** ([`chord_foot_error`]): the angle between the
//!   smooth-shaded normal at a chord's midpoint and the surface normal
//!   at the surface point the midpoint actually lies over. Vertex
//!   normals cannot see a chord that crosses iso-lines: along a
//!   helix on a coil tube the normal barely turns at the tube's top,
//!   so a long seam chord passes the spread bound while its midpoint
//!   sits on the tube 0.016 in `u` off the vertices' iso-line, where
//!   the true normal is 5.2° away (PK-587). Radial sagitta (a chord
//!   across a cylinder or sphere) leaves the foot normal equal to the
//!   shaded one, so the measure is zero there by construction.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, Edge};
use popkern_math::PortableMath;

/// Angle (radians) between two vectors; zero when either is null.
fn angle_between(a: Vec3, b: Vec3) -> f64 {
    let (la, lb) = (a.length(), b.length());
    if la < 1e-300 || lb < 1e-300 {
        return 0.0;
    }
    (a.dot(&b) / (la * lb)).clamp(-1.0, 1.0).pacos()
}

/// Tangent turn (radians) across the curve interval sampled at
/// `pa`, `pm` (parametric midpoint), `pb`.
pub(crate) fn curve_turn(pa: Point3, pm: Point3, pb: Point3) -> f64 {
    let d1 = Vec3::new(pm.x - pa.x, pm.y - pa.y, pm.z - pa.z);
    let d2 = Vec3::new(pb.x - pm.x, pb.y - pm.y, pb.z - pm.z);
    2.0 * angle_between(d1, d2)
}

/// Largest pairwise angle (radians) between the surface normals at
/// the corners and centre of the chart cell `[u0,u1] × [v0,v1]`,
/// taking one-sided limit normals at chart singularities (poles).
pub(crate) fn normal_spread(
    (surface, chart): (
        &SurfaceGeometry,
        &popkern_core::geometry::chart::ChartDescriptor,
    ),
    u: (f64, f64),
    v: (f64, f64),
) -> f64 {
    let samples = [
        (u.0, v.0),
        (u.1, v.0),
        (u.0, v.1),
        (u.1, v.1),
        (0.5 * (u.0 + u.1), 0.5 * (v.0 + v.1)),
    ];
    let normals: Vec<Vec3> = samples
        .iter()
        .filter_map(|&(a, b)| crate::vertex_normal::limit_normal_at(chart, surface, (a, b)))
        .collect();
    let mut worst = 0.0_f64;
    for i in 0..normals.len() {
        for j in i + 1..normals.len() {
            worst = worst.max(angle_between(normals[i], normals[j]));
        }
    }
    worst
}

/// Share of the angular facet bound granted to [`chord_foot_error`]:
/// the displayed-normal error the zebra facet-accuracy gate accepts is
/// half the bound (PK-578), and a chord's foot error is that error at
/// its worst interior point.
pub(crate) const FOOT_SHARE: f64 = 0.5;

/// Gauss–Newton iterations for a chord midpoint's foot point. One step
/// from the parametric midpoint already lands within the chord's
/// second-order sagitta; the second absorbs chart distortion, leaving
/// at most 1.3e-4 rad of truncation on an exact rational cylinder, far
/// below any foot share. A step below `FOOT_STEP_TOL` of the chord
/// ends the iteration early.
const FOOT_ITERATIONS: usize = 2;
const FOOT_STEP_TOL: f64 = 1e-9;

/// Chord-foot normal error (radians) of the chord `pa`–`pb`, whose ends
/// carry unit surface normals `na`, `nb` and whose parametric midpoint
/// is `uv_mid`. The error is the angle between the shaded normal at the
/// chord midpoint (`na + nb`, normalised) and the surface normal at the
/// midpoint's foot point, which is the closest surface point, found by
/// Gauss–Newton from `uv_mid`. Sign-agnostic, since facing is the
/// caller's.
///
/// Planes, cylinders and spheres answer zero without evaluation: their
/// foot is the radial projection, whose normal is the shaded one
/// exactly. The refine pass asks this of every edge, so the analytic
/// cases matter (PK-589).
///
/// Returns `None`, so the chord does not vote, when the foot is not
/// local: it lies farther than half the chord from the midpoint
/// (across a chart singularity or seam), the chart is singular there,
/// or the shaded normal cancels.
pub(crate) fn chord_foot_error(
    (surface, chart): (
        &SurfaceGeometry,
        &popkern_core::geometry::chart::ChartDescriptor,
    ),
    uv_mid: (f64, f64),
    (pa, pb): (Point3, Point3),
    (na, nb): (Vec3, Vec3),
) -> Option<f64> {
    if matches!(
        surface,
        SurfaceGeometry::Plane { .. }
            | SurfaceGeometry::Cylinder { .. }
            | SurfaceGeometry::Sphere { .. }
    ) {
        return Some(0.0);
    }
    let q = Point3::new(
        0.5 * (pa.x + pb.x),
        0.5 * (pa.y + pb.y),
        0.5 * (pa.z + pb.z),
    );
    let len = pa.distance_to(&pb);
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    let (per_u, per_v) =
        popkern_core::geometry::chart::ChartDescriptor::evaluation_periods_of(surface);
    let clamp = |x: f64, lo: f64, hi: f64, periodic: bool| {
        if periodic || !(lo.is_finite() && hi.is_finite()) {
            x
        } else {
            x.clamp(lo, hi)
        }
    };
    let mut uv = uv_mid;
    for _ in 0..FOOT_ITERATIONS {
        let s = surface.point_at(uv.0, uv.1);
        let r = Vec3::new(s.x - q.x, s.y - q.y, s.z - q.z);
        let (su, sv) = surface.first_partials(uv.0, uv.1);
        let (a, b, c) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
        let det = a * c - b * b;
        if det.is_nan() || det <= 1e-24 * (a * c).max(1e-300) {
            return None;
        }
        let (g0, g1) = (r.dot(&su), r.dot(&sv));
        let (du, dv) = ((c * g0 - b * g1) / det, (a * g1 - b * g0) / det);
        uv = (
            clamp(uv.0 - du, u0, u1, per_u.is_some()),
            clamp(uv.1 - dv, v0, v1, per_v.is_some()),
        );
        if !(uv.0.is_finite() && uv.1.is_finite()) {
            return None;
        }
        if (du * du * a + dv * dv * c).sqrt() <= FOOT_STEP_TOL * len {
            break;
        }
    }
    if surface.point_at(uv.0, uv.1).distance_to(&q) > 0.5 * len {
        return None;
    }
    // The normal from the chart service (the row's at a chart
    // singularity, none at an apex or within its neighbourhood, where no
    // facet bound can be met: the old partial-ratio gate's reach, without
    // which edges beside a vanishing fillet strip's tip split without end).
    let at = popkern_core::geometry::chart::ops::evaluate(chart, surface, uv);
    if crate::vertex_normal::singular_without_unique_normal(chart, &at) {
        return None;
    }
    let cross = crate::vertex_normal::limit_normal(chart, surface, uv, at)?;
    let shaded = Vec3::new(na.x + nb.x, na.y + nb.y, na.z + nb.z);
    if shaded.length() < 1e-12 {
        return None;
    }
    let angle = angle_between(cross, shaded);
    Some(angle.min(std::f64::consts::PI - angle))
}

/// Square-lattice station count for a boundary-first chart under the
/// angular facet bound: a lattice cell's normal spread scales ~1/n,
/// so double `n` (up to `max_subdivisions`) until the chart-centre
/// cell meets the bound. Without a bound `n` is returned unchanged.
/// This only seeds the interior density — the conforming refinement
/// pass (`angular_refine`) enforces the bound per triangle.
pub(crate) fn densify_lattice(
    (surface, chart): (
        &SurfaceGeometry,
        &popkern_core::geometry::chart::ChartDescriptor,
    ),
    ((u0, u1), (v0, v1)): ((f64, f64), (f64, f64)),
    mut n: usize,
    options: &crate::TessellationOptions,
) -> usize {
    let Some(max_spread) = options.angular_bound() else {
        return n;
    };
    let (um, vm) = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
    while n < options.max_subdivisions {
        let (hu, hv) = (0.5 * (u1 - u0) / n as f64, 0.5 * (v1 - v0) / n as f64);
        if normal_spread((surface, chart), (um - hu, um + hu), (vm - hv, vm + hv)) <= max_spread {
            break;
        }
        n *= 2;
    }
    n
}

/// One incident face's view of an edge: its surface and the pcurve
/// oriented along the EDGE CURVE's parameter (fraction 0 → curve
/// start).
struct FaceTrace<'a> {
    surface: &'a SurfaceGeometry,
    chart: std::rc::Rc<popkern_core::geometry::chart::ChartDescriptor>,
    pcurve: &'a Pcurve,
    reversed: bool,
}

/// Incident-face normal-turn measure for one edge (see the module
/// docs). Faces without a pcurve are skipped; an edge with no usable
/// face trace answers `None` and callers fall back to [`curve_turn`].
pub(crate) struct EdgeNormalProbe<'a> {
    traces: Vec<FaceTrace<'a>>,
}

impl<'a> EdgeNormalProbe<'a> {
    /// Build the probe for `edge`, orienting each incident half-edge's
    /// pcurve against the edge curve by comparing their images at the
    /// quarter fraction (endpoints coincide on closed edges, so they
    /// cannot decide direction).
    pub(crate) fn new(brep: &'a BRep, edge: &'a Edge) -> Self {
        let mut traces = Vec::new();
        let Some(curve) = edge.curve.as_ref() else {
            return Self { traces };
        };
        let curve_at = |f: f64| {
            curve
                .geometry
                .point_at(curve.t_start + (curve.t_end - curve.t_start) * f)
        };
        let q = curve_at(0.25);
        let mut half_edges = vec![edge.half_edge];
        if let Some(mate) = brep
            .half_edges
            .get(edge.half_edge)
            .and_then(|he| he.mate())
            .filter(|m| *m != edge.half_edge)
        {
            half_edges.push(mate);
        }
        for heh in half_edges {
            let Some(he) = brep.half_edges.get(heh) else {
                continue;
            };
            let (Some(pcurve), Some(surface)) = (
                he.pcurve.as_ref(),
                brep.loops
                    .get(he.loop_)
                    .and_then(|l| brep.faces.get(l.face))
                    .and_then(|f| f.surface.as_ref())
                    .map(|sof| &sof.geometry),
            ) else {
                continue;
            };
            let (t0, t1) = pcurve.parameter_range();
            if !(t0.is_finite() && t1.is_finite()) {
                continue;
            }
            let image = |t: f64| {
                let uv = pcurve.point_at(t);
                surface.point_at(uv.x, uv.y).distance_to(&q)
            };
            let reversed = image(t1 - 0.25 * (t1 - t0)) < image(t0 + 0.25 * (t1 - t0));
            traces.push(FaceTrace {
                chart: popkern_core::geometry::chart::ChartDescriptor::of(surface, &brep.tolerance),
                surface,
                pcurve,
                reversed,
            });
        }
        Self { traces }
    }

    /// Largest [`chord_foot_error`] (radians) any incident face shows
    /// for the chord between curve fractions `fa` and `fb` (images
    /// `pa`, `pb`), or `None` when no face votes.
    pub(crate) fn foot_error(&self, fa: f64, fb: f64, pa: Point3, pb: Point3) -> Option<f64> {
        let mut worst: Option<f64> = None;
        for trace in &self.traces {
            let (t0, t1) = trace.pcurve.parameter_range();
            let uv_at = |f: f64| {
                let t = if trace.reversed {
                    t1 - (t1 - t0) * f
                } else {
                    t0 + (t1 - t0) * f
                };
                let uv = trace.pcurve.point_at(t);
                (uv.x, uv.y)
            };
            let (ua, ub, um) = (uv_at(fa), uv_at(fb), uv_at(0.5 * (fa + fb)));
            let normal = |uv: (f64, f64)| {
                crate::vertex_normal::limit_normal_at(&trace.chart, trace.surface, uv)
            };
            let (Some(na), Some(nb)) = (normal(ua), normal(ub)) else {
                continue;
            };
            if let Some(e) = chord_foot_error((trace.surface, &trace.chart), um, (pa, pb), (na, nb))
            {
                worst = Some(worst.map_or(e, |w| w.max(e)));
            }
        }
        worst
    }

    /// Whether the step between curve fractions `fa` and `fb`, on any
    /// incident face's periodic chart (read as the nearest period
    /// image), exceeds a quarter period. The face rings unwrap each
    /// station step as the nearest image, which is faithful only under
    /// half a period; near a pole the longitude sweeps half a turn over
    /// a short, straight stretch of the edge that no 3D criterion
    /// splits (PK-1059: a hole 0.005 from a sphere's poles meshed the
    /// hole instead of the face, its volume read 3.34 for 1.40).
    pub(crate) fn wide_uv_step(&self, fa: f64, fb: f64) -> bool {
        self.traces.iter().any(|trace| {
            let (t0, t1) = trace.pcurve.parameter_range();
            let uv_at = |f: f64| {
                let t = if trace.reversed {
                    t1 - (t1 - t0) * f
                } else {
                    t0 + (t1 - t0) * f
                };
                trace.pcurve.point_at(t)
            };
            let (a, b) = (uv_at(fa), uv_at(fb));
            let (period_u, period_v) =
                popkern_core::geometry::chart::ChartDescriptor::evaluation_periods_of(
                    trace.surface,
                );
            [(b.x - a.x, period_u), (b.y - a.y, period_v)]
                .iter()
                .any(|&(d, p)| {
                    p.is_some_and(|p| p > 0.0 && (d - (d / p).round() * p).abs() > 0.25 * p)
                })
        })
    }

    /// Largest normal turn (radians) any incident face shows across
    /// the curve fractions `fractions` (pairwise over the samples), or
    /// `None` when no face trace is available.
    pub(crate) fn turn(&self, fractions: &[f64]) -> Option<f64> {
        if self.traces.is_empty() {
            return None;
        }
        let mut worst = 0.0_f64;
        for trace in &self.traces {
            let (t0, t1) = trace.pcurve.parameter_range();
            let normals: Vec<Vec3> = fractions
                .iter()
                .map(|&f| {
                    let t = if trace.reversed {
                        t1 - (t1 - t0) * f
                    } else {
                        t0 + (t1 - t0) * f
                    };
                    let uv = trace.pcurve.point_at(t);
                    // LIMIT normal: an edge ending at a pole meets it
                    // where the raw ∂u×∂v vanishes; dropping that
                    // sample measured the interval without its pole
                    // end (a 21.6° last seam segment read as 10.8°).
                    crate::vertex_normal::limit_normal_at(&trace.chart, trace.surface, (uv.x, uv.y))
                })
                .collect::<Option<Vec<Vec3>>>()
                .unwrap_or_default();
            for i in 0..normals.len() {
                for j in i + 1..normals.len() {
                    worst = worst.max(angle_between(normals[i], normals[j]));
                }
            }
        }
        Some(worst)
    }
}

#[cfg(test)]
#[path = "facet_angle_tests.rs"]
mod tests;
