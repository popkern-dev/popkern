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

    /// Periodicity of each parameter, as `(u_period, v_period)`.
    /// `Some(p)` means the parameter is `p`-periodic — evaluating
    /// the surface at `u` and `u + k·p` gives the same point for
    /// any integer `k`. `None` means non-periodic. Default is
    /// non-periodic for both; analytic surfaces with rotational
    /// symmetry (cylinder, cone, sphere, torus) override this.
    fn parameter_periodicity(&self) -> (Option<f64>, Option<f64>) {
        (None, None)
    }

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

// `closest_point_gauss_newton` and `bounding_box_by_sampling` —
// numerical helpers reusable by any `SurfaceEval` — live in
// [`super::surface_helpers`].  `crop_nurbs_surface_to_range` and
// `sample_to_nurbs` are private helpers used by
// `SurfaceGeometry::to_nurbs_in_range` below.
pub use super::surface_helpers::{bounding_box_by_sampling, closest_point_gauss_newton};
use super::surface_helpers::{crop_nurbs_surface_to_range, sample_to_nurbs};

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
            // Exact rational-Bezier conversions for all Gregory
            // variants. The natural `[0,1]²` case calls the
            // type-specific exact converter. Sub-ranges go
            // exact → split 4× to extract the requested rectangle;
            // this preserves the rational form exactly, unlike the
            // `sample_to_nurbs` fallback which is degree-3 lossy.
            // §5.3 item #9 — Gregory sub-range exactness.
            SurfaceGeometry::Gregory(gp) => {
                let full: NurbsSurface = match gp {
                    crate::geometry::gregory::GregoryPatch::Rectangular(r) => r.to_nurbs_exact(),
                    crate::geometry::gregory::GregoryPatch::RectangularDeg5(r) => {
                        r.to_nurbs_exact()
                    }
                    crate::geometry::gregory::GregoryPatch::RectangularBiquinticFull(r) => {
                        r.to_nurbs_exact()
                    }
                    crate::geometry::gregory::GregoryPatch::RectangularBiseptic(r) => {
                        r.to_nurbs_exact()
                    }
                    crate::geometry::gregory::GregoryPatch::RectangularTwoStripe(r) => {
                        r.to_nurbs_exact()
                    }
                    crate::geometry::gregory::GregoryPatch::RectangularBiNonic(r) => {
                        r.to_nurbs_exact()
                    }
                    crate::geometry::gregory::GregoryPatch::Triangular(t) => t.to_nurbs_exact(),
                };
                if u_range.0 == 0.0 && u_range.1 == 1.0 && v_range.0 == 0.0 && v_range.1 == 1.0 {
                    return full;
                }
                crop_nurbs_surface_to_range(full, u_range, v_range)
            }
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

    /// Approximate this surface as a NURBS over the caller's numeric UV range.
    ///
    /// Unlike [`Self::to_nurbs_in_range`], this preserves the requested
    /// numeric parameter coordinates by sampling `self.point_at(u, v)` over
    /// that range. Repair code uses this only behind validation gates when a
    /// face's trims already depend on a wider bounded parameter domain.
    pub fn approximate_to_nurbs_in_range(
        &self,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> NurbsSurface {
        sample_to_nurbs(self, u_range, v_range)
    }
}

#[cfg(test)]
#[path = "surface_tests.rs"]
mod tests;
