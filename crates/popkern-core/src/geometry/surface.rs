//! Surface geometry types, evaluation trait, and topology binding.

use crate::geometry::gregory::GregoryPatch;
use crate::geometry::intersection::surface_subdivision::Aabb;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

/// Evaluation interface for parametric surfaces.
pub trait SurfaceEval {
    /// Point on the surface at parameters (u, v).
    fn point_at(&self, u: f64, v: f64) -> Point3;
    /// Partial derivative with respect to u.
    fn partial_u(&self, u: f64, v: f64) -> Vec3;
    /// Partial derivative with respect to v.
    fn partial_v(&self, u: f64, v: f64) -> Vec3;
    /// Surface normal at (u, v). May not be unit length for degenerate points.
    fn normal_at(&self, u: f64, v: f64) -> Vec3;
    /// Valid parameter range ((u_min, u_max), (v_min, v_max)).
    fn parameter_range(&self) -> ((f64, f64), (f64, f64));

    /// Second partial derivative d²S/du² at (u, v).
    fn second_partial_uu(&self, u: f64, v: f64) -> Vec3 {
        let h = 1e-7;
        let fwd = self.partial_u(u + h, v);
        let bwd = self.partial_u(u - h, v);
        (fwd - bwd) * (1.0 / (2.0 * h))
    }

    /// Mixed second partial derivative d²S/dudv at (u, v).
    fn second_partial_uv(&self, u: f64, v: f64) -> Vec3 {
        let h = 1e-7;
        let fwd = self.partial_u(u, v + h);
        let bwd = self.partial_u(u, v - h);
        (fwd - bwd) * (1.0 / (2.0 * h))
    }

    /// Second partial derivative d²S/dv² at (u, v).
    fn second_partial_vv(&self, u: f64, v: f64) -> Vec3 {
        let h = 1e-7;
        let fwd = self.partial_v(u, v + h);
        let bwd = self.partial_v(u, v - h);
        (fwd - bwd) * (1.0 / (2.0 * h))
    }

    /// Find the (u, v) parameters of the closest point on the surface to `point`.
    ///
    /// Uses Gauss-Newton projection. Returns `None` if the projection fails
    /// to converge (e.g. the surface is degenerate or the point is equidistant
    /// from multiple locations).
    fn closest_point(&self, point: &Point3) -> Option<(f64, f64)>;

    /// Compute an axis-aligned bounding box over the given parameter ranges.
    ///
    /// Samples the surface on a grid and returns the tight AABB of those samples.
    fn bounding_box(&self, u_range: (f64, f64), v_range: (f64, f64)) -> Aabb;
}

/// Gauss-Newton closest-point projection for any `SurfaceEval`.
///
/// Multi-start from 5 candidates (center + 4 quadrant midpoints) within the
/// parameter range. Returns the (u, v) that minimizes distance to `point`,
/// or `None` if all starts fail to converge.
pub fn closest_point_gauss_newton(surface: &dyn SurfaceEval, point: &Point3) -> Option<(f64, f64)> {
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    // Clamp infinite ranges for the search
    let u_lo = u0.max(-100.0);
    let u_hi = u1.min(100.0);
    let v_lo = v0.max(-100.0);
    let v_hi = v1.min(100.0);

    // 5x5 grid of starts plus domain corners for periodic surfaces.
    // The 5x5 grid covers the interior; corners catch periodic boundary
    // cases where the closest point is near u=0 or v=0 on a torus/sphere.
    let n_grid = 5;
    let mut starts = Vec::with_capacity(n_grid * n_grid + 4);
    for i in 0..n_grid {
        let u = u_lo + (u_hi - u_lo) * (i as f64 + 0.5) / (n_grid as f64);
        for j in 0..n_grid {
            let v = v_lo + (v_hi - v_lo) * (j as f64 + 0.5) / (n_grid as f64);
            starts.push((u, v));
        }
    }
    // Add domain corners (important for periodic surfaces)
    starts.push((u_lo, v_lo));
    starts.push((u_lo, v_hi));
    starts.push((u_hi, v_lo));
    starts.push((u_hi, v_hi));

    let max_iter = 20;
    let step_tol = 1e-12;
    let max_dist = 1e6;

    let mut best: Option<(f64, f64, f64)> = None; // (u, v, dist)

    for (u_start, v_start) in starts {
        let mut u = u_start;
        let mut v = v_start;

        for _ in 0..max_iter {
            let p = surface.point_at(u, v);
            let diff = p - *point;
            let dist = diff.length();
            if dist > max_dist {
                break;
            }

            let su = surface.partial_u(u, v);
            let sv = surface.partial_v(u, v);

            // 2x2 normal equation: J^T J [du, dv]^T = -J^T r
            let a00 = su.dot(&su);
            let a01 = su.dot(&sv);
            let a11 = sv.dot(&sv);
            let b0 = diff.dot(&su);
            let b1 = diff.dot(&sv);

            let det = a00 * a11 - a01 * a01;
            if det.abs() < 1e-30 {
                break;
            }

            let du = -(a11 * b0 - a01 * b1) / det;
            let dv = -(a00 * b1 - a01 * b0) / det;

            u += du;
            v += dv;

            if du * du + dv * dv < step_tol {
                let final_p = surface.point_at(u, v);
                let final_dist = final_p.distance_to(point);
                match &best {
                    Some((_, _, bd)) if final_dist >= *bd => {}
                    _ => best = Some((u, v, final_dist)),
                }
                break;
            }
        }

        // Even if we didn't converge, record the best distance found
        let p = surface.point_at(u, v);
        let dist = p.distance_to(point);
        match &best {
            Some((_, _, bd)) if dist >= *bd => {}
            _ => best = Some((u, v, dist)),
        }
    }

    best.map(|(u, v, _)| (u, v))
}

/// Compute AABB by sampling the surface on a grid.
pub fn bounding_box_by_sampling(
    surface: &dyn SurfaceEval,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> Aabb {
    let n = 8;
    let mut points = Vec::with_capacity((n + 1) * (n + 1));

    for i in 0..=n {
        let u = u_range.0 + (u_range.1 - u_range.0) * (i as f64) / (n as f64);
        for j in 0..=n {
            let v = v_range.0 + (v_range.1 - v_range.0) * (j as f64) / (n as f64);
            points.push(surface.point_at(u, v));
        }
    }

    Aabb::from_points(&points).unwrap_or(Aabb {
        min: points
            .first()
            .copied()
            .unwrap_or(Point3::new(0.0, 0.0, 0.0)),
        max: points
            .first()
            .copied()
            .unwrap_or(Point3::new(0.0, 0.0, 0.0)),
    })
}

/// Geometric surface types.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SurfaceGeometry {
    /// An infinite plane.
    Plane {
        /// A point on the plane.
        origin: Point3,
        /// Unit normal of the plane.
        normal: Vec3,
        /// Unit vector along the u parameter direction.
        u_axis: Vec3,
        /// Unit vector along the v parameter direction.
        v_axis: Vec3,
    },
    /// A cylinder (infinite along axis).
    Cylinder {
        /// Point on the axis (origin of the parameterization).
        origin: Point3,
        /// Unit vector along the cylinder axis.
        axis: Vec3,
        /// Radius of the cylinder.
        radius: f64,
    },
    /// A cone.
    Cone {
        /// Apex (tip) of the cone.
        apex: Point3,
        /// Unit vector along the cone axis (pointing away from apex into the body).
        axis: Vec3,
        /// Half-angle between the axis and the surface generator, in radians.
        half_angle: f64,
    },
    /// A sphere.
    Sphere {
        /// Center of the sphere.
        center: Point3,
        /// Radius of the sphere.
        radius: f64,
    },
    /// A torus.
    Torus {
        /// Center of the torus.
        center: Point3,
        /// Unit vector along the torus symmetry axis.
        axis: Vec3,
        /// Distance from the torus center to the center of the tube.
        major_radius: f64,
        /// Radius of the tube.
        minor_radius: f64,
    },
    /// A NURBS surface.
    Nurbs(NurbsSurface),
    /// A rational Gregory patch (N=3 triangular or N=4 rectangular).
    ///
    /// Used for §5.3 corner patches where a tensor-product NURBS
    /// can't represent the twist-compatibility-resolving rational
    /// inner control-point blend. See [`GregoryPatch`] for the
    /// Chiyokura-Kimura rectangular form and the Longhi triangular
    /// form.
    Gregory(GregoryPatch),
}

/// Binding between a face and its geometric surface.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SurfaceOnFace {
    /// The geometric surface.
    pub geometry: SurfaceGeometry,
    /// True if surface normal agrees with face normal.
    pub orientation: bool,
}

// ---------------------------------------------------------------------------
// Analytic-to-NURBS conversion
// ---------------------------------------------------------------------------

/// Build a uniform clamped knot vector for `n_pts` control points and `degree`.
///
/// The knot vector has `degree` repeated knots at each end and uniform interior
/// knots spanning `[range_start, range_end]`.
fn uniform_clamped_knots(
    n_pts: usize,
    degree: usize,
    range_start: f64,
    range_end: f64,
) -> Vec<f64> {
    // n_pts control points, degree `degree` => n_pts + degree + 1 knots
    // interior spans: n_pts - degree
    let n_interior = n_pts - degree; // number of distinct interior spans (>= 1)
    let mut knots = Vec::with_capacity(n_pts + degree + 1);
    // `degree + 1` repeated start knots
    for _ in 0..=degree {
        knots.push(range_start);
    }
    // interior knots
    for i in 1..n_interior {
        let t = range_start + (range_end - range_start) * (i as f64) / (n_interior as f64);
        knots.push(t);
    }
    // `degree + 1` repeated end knots
    for _ in 0..=degree {
        knots.push(range_end);
    }
    knots
}

/// Sample the surface on an n_u × n_v grid and use the sampled points as NURBS
/// control points (degree-3 interpolating approximation).
fn sample_to_nurbs(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> NurbsSurface {
    use crate::geometry::surface::SurfaceEval;

    let n_u: usize = 9;
    let n_v: usize = 9;
    let degree: usize = 3;

    // Build parameter sample positions (uniform over range)
    let u_params: Vec<f64> = (0..n_u)
        .map(|i| u_range.0 + (u_range.1 - u_range.0) * (i as f64) / ((n_u - 1) as f64))
        .collect();
    let v_params: Vec<f64> = (0..n_v)
        .map(|j| v_range.0 + (v_range.1 - v_range.0) * (j as f64) / ((n_v - 1) as f64))
        .collect();

    // Evaluate grid: control_points[i_u * n_v + i_v]
    let control_points: Vec<Point3> = u_params
        .iter()
        .flat_map(|&u| v_params.iter().map(move |&v| surface.point_at(u, v)))
        .collect();

    let weights = vec![1.0_f64; n_u * n_v];

    let knots_u = uniform_clamped_knots(n_u, degree, u_range.0, u_range.1);
    let knots_v = uniform_clamped_knots(n_v, degree, v_range.0, v_range.1);

    NurbsSurface {
        degree_u: degree,
        degree_v: degree,
        knots_u,
        knots_v,
        control_points,
        num_u: n_u,
        num_v: n_v,
        weights,
    }
}

impl SurfaceGeometry {
    /// Returns whether the surface is periodic in the u and v directions,
    /// as a tuple `(periodic_u, periodic_v)`.
    ///
    /// Periodicity means the surface wraps around in that parameter direction.
    /// For example, a cylinder wraps in u (the angular direction) but not in v
    /// (the axial direction). A torus wraps in both u and v.
    pub fn is_periodic(&self) -> (bool, bool) {
        match self {
            SurfaceGeometry::Plane { .. } => (false, false),
            SurfaceGeometry::Cylinder { .. } => (true, false),
            SurfaceGeometry::Cone { .. } => (true, false),
            SurfaceGeometry::Sphere { .. } => (true, false),
            SurfaceGeometry::Torus { .. } => (true, true),
            SurfaceGeometry::Nurbs(_) => (false, false),
            // Gregory patches are bounded (unit square / simplex)
            // with no parameter wrap-around in either direction.
            SurfaceGeometry::Gregory(_) => (false, false),
        }
    }

    /// Returns the periodicity of each parameter direction.
    ///
    /// `Some(period)` if the surface wraps with that period, `None` if not
    /// periodic. For example, a cylinder has `(Some(TAU), None)` because u
    /// is the angular direction with period 2pi while v is the axial
    /// direction with no periodicity.
    pub fn periodicity(&self) -> (Option<f64>, Option<f64>) {
        use std::f64::consts::TAU;
        match self {
            SurfaceGeometry::Plane { .. } => (None, None),
            SurfaceGeometry::Cylinder { .. } => (Some(TAU), None),
            SurfaceGeometry::Cone { .. } => (Some(TAU), None),
            SurfaceGeometry::Sphere { .. } => (Some(TAU), None),
            SurfaceGeometry::Torus { .. } => (Some(TAU), Some(TAU)),
            SurfaceGeometry::Nurbs(_) => (None, None),
            SurfaceGeometry::Gregory(_) => (None, None),
        }
    }

    /// Convert to a [`NurbsSurface`] covering the given parameter range.
    pub fn to_nurbs_in_range(&self, u_range: (f64, f64), v_range: (f64, f64)) -> NurbsSurface {
        match self {
            SurfaceGeometry::Plane {
                origin,
                u_axis,
                v_axis,
                ..
            } => {
                // Degree (1,1) bilinear patch with 4 control points at corners.
                // Control points stored row-major (u varies within each row, i.e. index = i_u*2 + i_v):
                //   (i_u=0,i_v=0)=p00, (i_u=0,i_v=1)=p01, (i_u=1,i_v=0)=p10, (i_u=1,i_v=1)=p11
                let p00 = *origin + *u_axis * u_range.0 + *v_axis * v_range.0;
                let p01 = *origin + *u_axis * u_range.0 + *v_axis * v_range.1;
                let p10 = *origin + *u_axis * u_range.1 + *v_axis * v_range.0;
                let p11 = *origin + *u_axis * u_range.1 + *v_axis * v_range.1;
                NurbsSurface {
                    degree_u: 1,
                    degree_v: 1,
                    knots_u: vec![u_range.0, u_range.0, u_range.1, u_range.1],
                    knots_v: vec![v_range.0, v_range.0, v_range.1, v_range.1],
                    control_points: vec![p00, p01, p10, p11],
                    num_u: 2,
                    num_v: 2,
                    weights: vec![1.0, 1.0, 1.0, 1.0],
                }
            }
            SurfaceGeometry::Nurbs(nurbs) => nurbs.clone(),
            // Rectangular rational Gregory has an exact bidegree-(7,7)
            // rational-Bezier representation obtained by clearing the
            // four corner denominators. Used when the requested range
            // is the natural `[0,1]²` parameter domain; otherwise fall
            // through to the sample-fit path.
            //
            // Triangular Gregory is fundamentally non-tensor-product
            // and still routes through the sample-fit path.
            SurfaceGeometry::Gregory(gp)
                if matches!(
                    gp,
                    crate::geometry::gregory::GregoryPatch::Rectangular(_)
                        | crate::geometry::gregory::GregoryPatch::RectangularDeg5(_)
                ) && u_range.0 == 0.0
                    && u_range.1 == 1.0
                    && v_range.0 == 0.0
                    && v_range.1 == 1.0 =>
            {
                match gp {
                    crate::geometry::gregory::GregoryPatch::Rectangular(r) => r.to_nurbs_exact(),
                    crate::geometry::gregory::GregoryPatch::RectangularDeg5(r) => {
                        r.to_nurbs_exact()
                    }
                    _ => unreachable!(),
                }
            }
            SurfaceGeometry::Gregory(_) => sample_to_nurbs(self, u_range, v_range),
            _ => sample_to_nurbs(self, u_range, v_range),
        }
    }

    /// Convert to a [`NurbsSurface`] using the default parameter range from
    /// [`SurfaceEval::parameter_range()`].
    ///
    /// For infinite surfaces such as planes and cylinders, a finite range must be
    /// specified manually via [`to_nurbs_in_range`](Self::to_nurbs_in_range).
    pub fn to_nurbs(&self) -> NurbsSurface {
        use crate::geometry::surface::SurfaceEval;
        let ((u0, u1), (v0, v1)) = self.parameter_range();
        self.to_nurbs_in_range((u0, u1), (v0, v1))
    }
}

#[cfg(test)]
#[path = "surface_tests.rs"]
mod tests;
