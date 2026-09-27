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

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, Edge};

/// Angle (radians) between two vectors; zero when either is null.
fn angle_between(a: Vec3, b: Vec3) -> f64 {
    let (la, lb) = (a.length(), b.length());
    if la < 1e-300 || lb < 1e-300 {
        return 0.0;
    }
    (a.dot(&b) / (la * lb)).clamp(-1.0, 1.0).acos()
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
pub(crate) fn normal_spread(surface: &SurfaceGeometry, u: (f64, f64), v: (f64, f64)) -> f64 {
    let samples = [
        (u.0, v.0),
        (u.1, v.0),
        (u.0, v.1),
        (u.1, v.1),
        (0.5 * (u.0 + u.1), 0.5 * (v.0 + v.1)),
    ];
    let normals: Vec<Vec3> = samples
        .iter()
        .filter_map(|&(a, b)| crate::vertex_normal::limit_normal(surface, a, b))
        .collect();
    let mut worst = 0.0_f64;
    for i in 0..normals.len() {
        for j in i + 1..normals.len() {
            worst = worst.max(angle_between(normals[i], normals[j]));
        }
    }
    worst
}

/// Square-lattice station count for a boundary-first chart under the
/// angular facet bound: a lattice cell's normal spread scales ~1/n,
/// so double `n` (up to `max_subdivisions`) until the chart-centre
/// cell meets the bound. Without a bound `n` is returned unchanged.
/// This only seeds the interior density — the conforming refinement
/// pass (`angular_refine`) enforces the bound per triangle.
pub(crate) fn densify_lattice(
    surface: &SurfaceGeometry,
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
        if normal_spread(surface, (um - hu, um + hu), (vm - hv, vm + hv)) <= max_spread {
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
        if let Some(he) = brep.half_edges.get(edge.half_edge)
            && he.twin != edge.half_edge
        {
            half_edges.push(he.twin);
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
                surface,
                pcurve,
                reversed,
            });
        }
        Self { traces }
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
                    crate::vertex_normal::limit_normal(trace.surface, uv.x, uv.y)
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
mod tests {
    use super::*;

    #[test]
    fn curve_turn_of_a_circular_arc_is_its_subtended_angle() {
        let r = 0.04;
        for &(a0, a1) in &[(0.0, 0.5 * std::f64::consts::PI), (0.3, 0.5), (1.0, 1.01)] {
            let at = |t: f64| Point3::new(r * f64::cos(t), r * f64::sin(t), 0.0);
            let turn = curve_turn(at(a0), at(0.5 * (a0 + a1)), at(a1));
            assert!((turn - (a1 - a0)).abs() < 1e-9, "{turn} vs {}", a1 - a0);
        }
    }

    #[test]
    fn curve_turn_of_a_line_is_zero() {
        let turn = curve_turn(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(2.0, 2.0, 2.0),
        );
        assert!(turn.abs() < 1e-7);
    }

    #[test]
    fn normal_spread_of_a_cylinder_cell_is_its_angular_width() {
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 0.1,
        };
        let spread = normal_spread(&cyl, (0.0, 0.4), (0.0, 5.0));
        assert!((spread - 0.4).abs() < 1e-9, "{spread}");
    }

    #[test]
    fn normal_spread_of_a_plane_is_zero() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        assert!(normal_spread(&plane, (-3.0, 3.0), (-3.0, 3.0)) < 1e-12);
    }
}
