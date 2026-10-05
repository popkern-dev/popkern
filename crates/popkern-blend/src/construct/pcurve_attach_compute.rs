//! Pcurve computation core for [`super::pcurve_attach`]: analytic
//! candidates, anchored projection walks with period unwrapping, and
//! the LSQ fit — split out for the 500-line structural gate. The
//! attach/wiring passes stay in `pcurve_attach.rs`.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_math::PortableMath;

use popkern_core::topology::{BRep, HalfEdgeHandle};

use super::pcurve_attach::{
    PCURVE_FIT_CONTROL_POINTS, PCURVE_FIT_DEGREE, PCURVE_SAMPLE_COUNT, project_to_uv_clamped,
};

/// Compute a pcurve for a single half-edge by projecting its 3D geometry onto a surface.
///
/// Samples the half-edge's underlying edge curve (or linear interpolation
/// between origin and destination vertices) at `PCURVE_SAMPLE_COUNT` points,
/// projects each into the surface's UV space, and fits a NURBS curve.
/// Analytic inversion for cylinder charts: generic Newton projection
/// wanders across period branches and can converge to the far-side
/// critical point (measured: 1.84 pcurve-image displacement on a
/// strip-end ruling lying ON the cylinder to machine zero). Points
/// at or near the surface invert exactly via atan2 in the cylinder's
/// own basis — the same cure the boundary-first mesher and the
/// projection walk use.
pub(super) fn analytic_cylinder_uv(surface: &SurfaceGeometry, p: &Point3) -> Option<(f64, f64)> {
    let SurfaceGeometry::Cylinder { origin, axis, .. } = surface else {
        return None;
    };
    let (perp1, perp2) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(*axis);
    let d = *p - *origin;
    Some((d.dot(&perp2).patan2(d.dot(&perp1)), d.dot(axis)))
}

pub(super) fn compute_halfedge_pcurve(
    brep: &BRep,
    he: HalfEdgeHandle,
    surface: &SurfaceGeometry,
) -> Option<NurbsCurve> {
    let he_data = brep.half_edges.get(he)?;
    let origin = brep.vertices[he_data.origin].position;

    // Find the destination vertex (next half-edge's origin).
    let next_data = brep.half_edges.get(he_data.next)?;
    let dest = brep.vertices[next_data.origin].position;

    // Try to get the 3D edge curve.
    let edge_data = brep.edges.get(he_data.edge)?;

    // Check half-edge direction vs edge direction.
    let edge_he = edge_data.half_edge;
    let forward = edge_he == he;

    // CONTINUATION projection with period unwrapping (PK-211): the
    // per-sample independent global projections branch-mix on
    // periodic charts — see compute_pcurve_on_surface.
    let n = PCURVE_SAMPLE_COUNT;
    let sample_at = |i: usize| -> Point3 {
        let t_frac = i as f64 / (n - 1) as f64;
        match &edge_data.curve {
            Some(curve) => {
                let t = if forward {
                    curve.t_start + t_frac * (curve.t_end - curve.t_start)
                } else {
                    curve.t_end - t_frac * (curve.t_end - curve.t_start)
                };
                curve.geometry.point_at(t)
            }
            None => Point3::new(
                origin.x + t_frac * (dest.x - origin.x),
                origin.y + t_frac * (dest.y - origin.y),
                origin.z + t_frac * (dest.z - origin.z),
            ),
        }
    };
    let (period_u, period_v) = surface.parameter_periodicity();
    let mut prev: Option<(f64, f64)> = None;
    let mut uv_points: Vec<Point3> = Vec::with_capacity(n);
    for i in 0..n {
        let pt = sample_at(i);
        let (mut u, mut v) = match analytic_cylinder_uv(surface, &pt) {
            Some(uv) => uv,
            None => match prev {
                Some((pu, pv)) => {
                    match popkern_core::geometry::uv_projection::project_point_to_uv_near(
                        surface, &pt, pu, pv,
                    ) {
                        Some(uv) => uv,
                        None => project_to_uv_clamped(surface, &pt),
                    }
                }
                None => project_to_uv_clamped(surface, &pt),
            },
        };
        if let Some((pu, pv)) = prev {
            if let Some(per) = period_u {
                while u - pu > per * 0.5 {
                    u -= per;
                }
                while u - pu < -per * 0.5 {
                    u += per;
                }
            }
            if let Some(per) = period_v {
                while v - pv > per * 0.5 {
                    v -= per;
                }
                while v - pv < -per * 0.5 {
                    v += per;
                }
            }
        }
        prev = Some((u, v));
        uv_points.push(Point3::new(u, v, 0.0));
    }

    if uv_points.len() < 2 {
        return None;
    }

    fit_pcurve(&uv_points)
}

/// Compute a pcurve by projecting a 3D NURBS curve onto a surface.
///
/// Samples the curve at `PCURVE_SAMPLE_COUNT` points, projects each onto the
/// surface's UV space, and fits a NURBS curve through the UV points.
pub(super) fn compute_pcurve_on_surface(
    curve_3d: &NurbsCurve,
    surface: &SurfaceGeometry,
) -> Option<NurbsCurve> {
    let (t_start, t_end) = curve_3d.parameter_range();
    let n = PCURVE_SAMPLE_COUNT;

    // CONTINUATION projection with period unwrapping (PK-211): each
    // sample independently globally projected mixes period branches
    // on a periodic chart (Newton lands wherever its internal seed
    // falls), and the chord-length cubic fit through branch-mixed
    // samples is garbage — measured: a strip-end ruling's pcurve
    // image 1.84 model units off its own line on offset_bracket,
    // which every later repair kept because the re-walk reproduced
    // the same mixing. Seed each sample from the previous result and
    // unwrap toward it — the same discipline as
    // pcurve_projection_walk and the boundary-first mesher.
    let (period_u, period_v) = surface.periodicity();
    let mut prev: Option<(f64, f64)> = None;
    let mut uv_points: Vec<Point3> = Vec::with_capacity(n);
    for i in 0..n {
        let t = t_start + (t_end - t_start) * (i as f64 / (n - 1) as f64);
        let pt = curve_3d.point_at(t);
        let (mut u, mut v) = match analytic_cylinder_uv(surface, &pt) {
            Some(uv) => uv,
            None => match prev {
                Some((pu, pv)) => {
                    match popkern_core::geometry::uv_projection::project_point_to_uv_near(
                        surface, &pt, pu, pv,
                    ) {
                        Some(uv) => uv,
                        None => project_to_uv_clamped(surface, &pt),
                    }
                }
                None => project_to_uv_clamped(surface, &pt),
            },
        };
        if let Some((pu, pv)) = prev {
            if let Some(per) = period_u {
                while u - pu > per * 0.5 {
                    u -= per;
                }
                while u - pu < -per * 0.5 {
                    u += per;
                }
            }
            if let Some(per) = period_v {
                while v - pv > per * 0.5 {
                    v -= per;
                }
                while v - pv < -per * 0.5 {
                    v += per;
                }
            }
        }
        prev = Some((u, v));
        uv_points.push(Point3::new(u, v, 0.0));
    }

    if uv_points.len() < 2 {
        return None;
    }

    fit_pcurve(&uv_points)
}

/// Fit a degree-3 NURBS curve through UV sample points.
///
/// Uses chord-length parameterization and [`popkern_core::geometry::fitting::fit_curve`].
/// Returns `None` if fitting fails.
pub(super) fn fit_pcurve(uv_points: &[Point3]) -> Option<NurbsCurve> {
    let n = uv_points.len();
    if n < 2 {
        return None;
    }

    // Chord-length parameterization.
    let mut params = vec![0.0_f64; n];
    for i in 1..n {
        params[i] = params[i - 1] + uv_points[i].distance_to(&uv_points[i - 1]);
    }
    let total = params[n - 1];
    if total < 1e-15 {
        // Degenerate: all points coincide. Return a constant pcurve.
        return Some(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![uv_points[0], uv_points[0]],
            weights: vec![1.0, 1.0],
        });
    }
    for p in &mut params {
        *p /= total;
    }

    let degree = PCURVE_FIT_DEGREE.min(n - 1);
    let num_cp = PCURVE_FIT_CONTROL_POINTS.min(n).max(degree + 1);

    popkern_core::geometry::fitting::fit_curve(uv_points, &params, degree, num_cp).ok()
}

/// Reverse a NURBS curve (swap traversal direction).
pub(super) fn reverse_nurbs(curve: &NurbsCurve) -> NurbsCurve {
    let (t_start, t_end) = curve.parameter_range();
    let mut reversed_knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|&k| t_start + t_end - k)
        .collect();
    reversed_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let mut reversed_cp = curve.control_points.clone();
    reversed_cp.reverse();

    let mut reversed_weights = curve.weights.clone();
    reversed_weights.reverse();

    NurbsCurve {
        degree: curve.degree,
        knots: reversed_knots,
        control_points: reversed_cp,
        weights: reversed_weights,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Vec3;

    #[test]
    fn compute_pcurve_on_plane() {
        let surface = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let curve_3d = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let result = compute_pcurve_on_surface(&curve_3d, &surface);
        assert!(result.is_some(), "should compute pcurve on plane");
        let pcurve = result.unwrap();
        // All control points should have z=0.
        for cp in &pcurve.control_points {
            assert!(cp.z.abs() < 1e-10, "pcurve z should be 0, got {}", cp.z);
        }
    }

    #[test]
    fn fit_pcurve_degenerate_points() {
        let pts = vec![Point3::new(0.5, 0.5, 0.0); 5];
        let result = fit_pcurve(&pts);
        assert!(
            result.is_some(),
            "degenerate points should produce a pcurve"
        );
    }

    #[test]
    fn fit_pcurve_too_few_points() {
        let result = fit_pcurve(&[]);
        assert!(result.is_none(), "empty points should return None");
        let result = fit_pcurve(&[Point3::new(0.0, 0.0, 0.0)]);
        assert!(result.is_none(), "single point should return None");
    }

    #[test]
    fn reverse_nurbs_roundtrip() {
        let curve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let reversed = reverse_nurbs(&curve);
        let (t0, t1) = reversed.parameter_range();
        let start = reversed.point_at(t0);
        let end = reversed.point_at(t1);
        assert!(start.distance_to(&Point3::new(1.0, 1.0, 0.0)) < 1e-10);
        assert!(end.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 1e-10);
    }
}
