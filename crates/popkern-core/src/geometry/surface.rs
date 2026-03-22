//! Surface geometry types, evaluation trait, and topology binding.

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
mod tests {
    use super::*;
    use crate::geometry::surface::SurfaceEval;

    #[test]
    fn plane_to_nurbs_preserves_shape() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let nurbs = plane.to_nurbs_in_range((-1.0, 1.0), (-1.0, 1.0));
        for u in [-1.0_f64, -0.5, 0.0, 0.5, 0.99] {
            for v in [-1.0_f64, -0.5, 0.0, 0.5, 0.99] {
                let p1 = plane.point_at(u, v);
                let p2 = nurbs.point_at(u, v);
                assert!(
                    p1.distance_to(&p2) < 1e-8,
                    "Plane-to-NURBS mismatch at u={u}, v={v}: p1={p1:?}, p2={p2:?}"
                );
            }
        }
    }

    #[test]
    fn nurbs_to_nurbs_is_clone() {
        let nurbs = NurbsSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
            num_u: 2,
            num_v: 2,
            weights: vec![1.0, 1.0, 1.0, 1.0],
        };
        let geom = SurfaceGeometry::Nurbs(nurbs.clone());
        let promoted = geom.to_nurbs_in_range((0.0, 1.0), (0.0, 1.0));
        assert_eq!(promoted.degree_u, nurbs.degree_u);
    }

    #[test]
    fn sphere_to_nurbs_approximate() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let nurbs = sphere.to_nurbs_in_range((0.5, 2.5), (0.5, 2.5));
        for u in [0.6_f64, 1.0, 1.5, 2.0, 2.4] {
            for v in [0.6_f64, 1.0, 1.5, 2.0, 2.4] {
                let p = nurbs.point_at(u, v);
                let r = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
                assert!(
                    (r - 1.0).abs() < 0.1,
                    "Point not on unit sphere: r={r} at u={u}, v={v}"
                );
            }
        }
    }

    #[test]
    fn sphere_to_nurbs_default_range() {
        // Tests the to_nurbs() method (lines 195-198) which uses parameter_range()
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let nurbs = sphere.to_nurbs();
        // Sphere parameter range is [0, TAU] x [0, PI]
        let p = nurbs.point_at(std::f64::consts::FRAC_PI_2, std::f64::consts::FRAC_PI_2);
        let r = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
        assert!((r - 1.0).abs() < 0.15, "Point not on unit sphere: r={r}");
    }

    #[test]
    fn cylinder_to_nurbs_in_range() {
        // Cylinder's default range may be infinite, so use explicit range
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let nurbs = cyl.to_nurbs_in_range((0.0, std::f64::consts::TAU), (-1.0, 1.0));
        let p = nurbs.point_at(std::f64::consts::FRAC_PI_2, 0.0);
        let r = (p.x * p.x + p.y * p.y).sqrt();
        assert!((r - 1.0).abs() < 0.15, "Point not on cylinder: r={r}");
    }
}
