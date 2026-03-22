//! NURBS surface evaluation via tensor-product de Boor algorithm.

use crate::error::CoreError;
use crate::geometry::nurbs_ops::{self, Homogeneous, de_boor_1d};
use crate::geometry::surface::SurfaceEval;
use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

/// A Non-Uniform Rational B-Spline surface.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NurbsSurface {
    /// Polynomial degree in the u direction.
    pub degree_u: usize,
    /// Polynomial degree in the v direction.
    pub degree_v: usize,
    /// Knot vector in the u direction.
    pub knots_u: Vec<f64>,
    /// Knot vector in the v direction.
    pub knots_v: Vec<f64>,
    /// Control point grid, stored row-major (u varies within each row).
    /// Length = num_u * num_v.
    pub control_points: Vec<Point3>,
    /// Number of control points in the u direction.
    pub num_u: usize,
    /// Number of control points in the v direction.
    pub num_v: usize,
    /// Weights, one per control point.
    pub weights: Vec<f64>,
}

impl NurbsSurface {
    /// Access control point at grid position (i_u, i_v) in the flat array.
    fn cp(&self, i_u: usize, i_v: usize) -> &Point3 {
        &self.control_points[i_u * self.num_v + i_v]
    }

    /// Access weight at grid position (i_u, i_v) in the flat array.
    fn wt(&self, i_u: usize, i_v: usize) -> f64 {
        self.weights[i_u * self.num_v + i_v]
    }

    /// Insert a single knot in the u direction. Surface shape is mathematically unchanged.
    pub fn insert_knot_u(&self, u: f64) -> Result<NurbsSurface, CoreError> {
        let ((u_min, u_max), _) = self.parameter_range();
        if u < u_min || u > u_max {
            return Err(CoreError::PreconditionFailed(format!(
                "u={} outside range [{}, {}]",
                u, u_min, u_max
            )));
        }

        // Process each column (fixed v_idx, varying i_u). Storage is row-major:
        // flat index = i_u * num_v + v_idx
        let mut new_knots_u = Vec::new();
        // After insertion, each column gains one point. We'll collect all results
        // and reassemble into a flat row-major grid.
        let new_num_u = self.num_u + 1;
        let new_num_v = self.num_v;
        let mut new_pts = vec![Point3::new(0.0, 0.0, 0.0); new_num_u * new_num_v];
        let mut new_wts = vec![0.0f64; new_num_u * new_num_v];

        for v_idx in 0..self.num_v {
            let col_pts: Vec<Point3> = (0..self.num_u)
                .map(|i_u| self.control_points[i_u * self.num_v + v_idx])
                .collect();
            let col_wts: Vec<f64> = (0..self.num_u)
                .map(|i_u| self.weights[i_u * self.num_v + v_idx])
                .collect();
            let (nk, np, nw) = nurbs_ops::knot_insertion::insert_knot_1d(
                self.degree_u,
                &self.knots_u,
                &col_pts,
                &col_wts,
                u,
            );
            if v_idx == 0 {
                new_knots_u = nk;
            }
            for (i_u, (pt, wt)) in np.into_iter().zip(nw.into_iter()).enumerate() {
                new_pts[i_u * new_num_v + v_idx] = pt;
                new_wts[i_u * new_num_v + v_idx] = wt;
            }
        }

        Ok(NurbsSurface {
            degree_u: self.degree_u,
            degree_v: self.degree_v,
            knots_u: new_knots_u,
            knots_v: self.knots_v.clone(),
            control_points: new_pts,
            num_u: new_num_u,
            num_v: new_num_v,
            weights: new_wts,
        })
    }

    /// Insert a single knot in the v direction. Surface shape is mathematically unchanged.
    pub fn insert_knot_v(&self, v: f64) -> Result<NurbsSurface, CoreError> {
        let (_, (v_min, v_max)) = self.parameter_range();
        if v < v_min || v > v_max {
            return Err(CoreError::PreconditionFailed(format!(
                "v={} outside range [{}, {}]",
                v, v_min, v_max
            )));
        }

        // Process each row (fixed i_u, varying v). Storage is row-major:
        // row i_u lives at flat indices [i_u*num_v .. (i_u+1)*num_v].
        let new_num_u = self.num_u;
        let new_num_v = self.num_v + 1;
        let mut new_pts = Vec::with_capacity(new_num_u * new_num_v);
        let mut new_wts = Vec::with_capacity(new_num_u * new_num_v);
        let mut new_knots_v = Vec::new();

        for i_u in 0..self.num_u {
            let row_pts: Vec<Point3> = (0..self.num_v)
                .map(|i_v| self.control_points[i_u * self.num_v + i_v])
                .collect();
            let row_wts: Vec<f64> = (0..self.num_v)
                .map(|i_v| self.weights[i_u * self.num_v + i_v])
                .collect();
            let (nk, np, nw) = nurbs_ops::knot_insertion::insert_knot_1d(
                self.degree_v,
                &self.knots_v,
                &row_pts,
                &row_wts,
                v,
            );
            if i_u == 0 {
                new_knots_v = nk;
            }
            new_pts.extend(np);
            new_wts.extend(nw);
        }

        Ok(NurbsSurface {
            degree_u: self.degree_u,
            degree_v: self.degree_v,
            knots_u: self.knots_u.clone(),
            knots_v: new_knots_v,
            control_points: new_pts,
            num_u: new_num_u,
            num_v: new_num_v,
            weights: new_wts,
        })
    }

    /// Insert multiple sorted knots in the u direction. Surface shape is mathematically unchanged.
    pub fn insert_knots_u(&self, knots_to_insert: &[f64]) -> Result<NurbsSurface, CoreError> {
        let ((u_min, u_max), _) = self.parameter_range();
        for &u in knots_to_insert {
            if u < u_min || u > u_max {
                return Err(CoreError::PreconditionFailed(format!(
                    "u={} outside range [{}, {}]",
                    u, u_min, u_max
                )));
            }
        }
        let mut surface = self.clone();
        for &u in knots_to_insert {
            surface = surface.insert_knot_u(u)?;
        }
        Ok(surface)
    }

    /// Insert multiple sorted knots in the v direction. Surface shape is mathematically unchanged.
    pub fn insert_knots_v(&self, knots_to_insert: &[f64]) -> Result<NurbsSurface, CoreError> {
        let (_, (v_min, v_max)) = self.parameter_range();
        for &v in knots_to_insert {
            if v < v_min || v > v_max {
                return Err(CoreError::PreconditionFailed(format!(
                    "v={} outside range [{}, {}]",
                    v, v_min, v_max
                )));
            }
        }
        let mut surface = self.clone();
        for &v in knots_to_insert {
            surface = surface.insert_knot_v(v)?;
        }
        Ok(surface)
    }

    /// Evaluate the surface at (u, v) using tensor-product de Boor in homogeneous coordinates.
    fn eval_homogeneous(&self, u: f64, v: f64) -> Homogeneous {
        // Step 1: For each u-row, evaluate the 1D B-spline in v direction.
        let mut intermediate: Vec<Homogeneous> = Vec::with_capacity(self.num_u);
        for i_u in 0..self.num_u {
            let row: Vec<Homogeneous> = (0..self.num_v)
                .map(|i_v| Homogeneous::from_point_weight(self.cp(i_u, i_v), self.wt(i_u, i_v)))
                .collect();
            intermediate.push(de_boor_1d(self.degree_v, &self.knots_v, &row, v));
        }

        // Step 2: Evaluate the 1D B-spline in u direction on the intermediate points.
        de_boor_1d(self.degree_u, &self.knots_u, &intermediate, u)
    }
}

impl SurfaceEval for NurbsSurface {
    fn point_at(&self, u: f64, v: f64) -> Point3 {
        self.eval_homogeneous(u, v).to_point3()
    }

    fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        let h = 1e-7;
        let p_fwd = self.point_at(u + h, v);
        let p_bwd = self.point_at(u - h, v);
        let inv_2h = 1.0 / (2.0 * h);
        Vec3::new(
            (p_fwd.x - p_bwd.x) * inv_2h,
            (p_fwd.y - p_bwd.y) * inv_2h,
            (p_fwd.z - p_bwd.z) * inv_2h,
        )
    }

    fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        let h = 1e-7;
        let p_fwd = self.point_at(u, v + h);
        let p_bwd = self.point_at(u, v - h);
        let inv_2h = 1.0 / (2.0 * h);
        Vec3::new(
            (p_fwd.x - p_bwd.x) * inv_2h,
            (p_fwd.y - p_bwd.y) * inv_2h,
            (p_fwd.z - p_bwd.z) * inv_2h,
        )
    }

    fn normal_at(&self, u: f64, v: f64) -> Vec3 {
        let du = self.partial_u(u, v);
        let dv = self.partial_v(u, v);
        du.cross(&dv)
    }

    fn parameter_range(&self) -> ((f64, f64), (f64, f64)) {
        let p = self.degree_u;
        let n_u = self.num_u - 1;
        let q = self.degree_v;
        let n_v = self.num_v - 1;
        (
            (self.knots_u[p], self.knots_u[n_u + 1]),
            (self.knots_v[q], self.knots_v[n_v + 1]),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::{Point3, SurfaceEval};

    /// Helper to build a bilinear patch for testing.
    fn bilinear_patch() -> NurbsSurface {
        NurbsSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                // Row u=0: (0,0,0), (0,1,0)
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                // Row u=1: (1,0,0), (1,1,0)
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
            num_u: 2,
            num_v: 2,
            weights: vec![1.0, 1.0, 1.0, 1.0],
        }
    }

    #[test]
    fn bilinear_patch_center() {
        let surface = bilinear_patch();
        let p = surface.point_at(0.5, 0.5);
        assert!((p.x - 0.5).abs() < 1e-10);
        assert!((p.y - 0.5).abs() < 1e-10);
        assert!((p.z).abs() < 1e-10);
    }

    #[test]
    fn bilinear_patch_corners() {
        let surface = bilinear_patch();
        let p00 = surface.point_at(0.0, 0.0);
        assert!((p00.x).abs() < 1e-10);
        assert!((p00.y).abs() < 1e-10);
        // Use slightly less than 1.0 to stay in valid knot span
        let p11 = surface.point_at(1.0 - 1e-12, 1.0 - 1e-12);
        assert!((p11.x - 1.0).abs() < 1e-6);
        assert!((p11.y - 1.0).abs() < 1e-6);
    }

    #[test]
    fn bilinear_patch_normal_is_z() {
        let surface = bilinear_patch();
        let n = surface.normal_at(0.5, 0.5);
        // Normal should point in +z or -z direction for a flat XY patch
        let nlen = n.length();
        assert!(nlen > 1e-10, "Normal should be non-zero");
        let nz = n.z.abs() / nlen;
        assert!((nz - 1.0).abs() < 1e-6, "Normal should be ±z, got {:?}", n);
    }

    #[test]
    fn parameter_range() {
        let surface = bilinear_patch();
        let ((u0, u1), (v0, v1)) = surface.parameter_range();
        assert_eq!(u0, 0.0);
        assert_eq!(u1, 1.0);
        assert_eq!(v0, 0.0);
        assert_eq!(v1, 1.0);
    }

    #[test]
    fn bilinear_patch_partial_u() {
        let surface = bilinear_patch();
        let du = surface.partial_u(0.5, 0.5);
        // For a bilinear patch, d/du should be (1, 0, 0)
        assert!((du.x - 1.0).abs() < 1e-4);
        assert!((du.y).abs() < 1e-4);
    }

    #[test]
    fn bilinear_patch_partial_v() {
        let surface = bilinear_patch();
        let dv = surface.partial_v(0.5, 0.5);
        // For a bilinear patch, d/dv should be (0, 1, 0)
        assert!((dv.x).abs() < 1e-4);
        assert!((dv.y - 1.0).abs() < 1e-4);
    }

    #[test]
    fn bilinear_patch_at_end() {
        let surface = bilinear_patch();
        // Evaluate at exactly t=1.0 (clamped end)
        let p = surface.point_at(1.0, 1.0);
        assert!((p.x - 1.0).abs() < 1e-6);
        assert!((p.y - 1.0).abs() < 1e-6);
    }

    #[test]
    fn bilinear_patch_at_before_start() {
        let surface = bilinear_patch();
        // Evaluate before start — clamped to degree span, extrapolates
        // Just verify it doesn't panic and returns finite values
        let p = surface.point_at(-0.1, -0.1);
        assert!(p.x.is_finite());
        assert!(p.y.is_finite());
    }

    #[test]
    fn surface_insert_knot_u_preserves_shape() {
        let surface = bilinear_patch();
        let new_surface = surface.insert_knot_u(0.5).unwrap();
        // num_u grows by 1
        assert_eq!(new_surface.num_u, surface.num_u + 1);
        assert_eq!(new_surface.num_v, surface.num_v);
        // Shape preserved at multiple sample points
        for &u in &[0.0, 0.25, 0.5, 0.75, 0.99] {
            for &v in &[0.0, 0.25, 0.5, 0.75, 0.99] {
                let p1 = surface.point_at(u, v);
                let p2 = new_surface.point_at(u, v);
                assert!(
                    (p1.x - p2.x).abs() < 1e-10,
                    "x mismatch at u={u},v={v}: {p1:?} vs {p2:?}"
                );
                assert!(
                    (p1.y - p2.y).abs() < 1e-10,
                    "y mismatch at u={u},v={v}: {p1:?} vs {p2:?}"
                );
            }
        }
    }

    #[test]
    fn surface_insert_knot_v_preserves_shape() {
        let surface = bilinear_patch();
        let new_surface = surface.insert_knot_v(0.5).unwrap();
        assert_eq!(new_surface.num_u, surface.num_u);
        assert_eq!(new_surface.num_v, surface.num_v + 1);
        for &u in &[0.0, 0.25, 0.5, 0.75, 0.99] {
            for &v in &[0.0, 0.25, 0.5, 0.75, 0.99] {
                let p1 = surface.point_at(u, v);
                let p2 = new_surface.point_at(u, v);
                assert!(
                    (p1.x - p2.x).abs() < 1e-10,
                    "x mismatch at u={u},v={v}: {p1:?} vs {p2:?}"
                );
                assert!(
                    (p1.y - p2.y).abs() < 1e-10,
                    "y mismatch at u={u},v={v}: {p1:?} vs {p2:?}"
                );
            }
        }
    }

    #[test]
    fn surface_insert_knot_u_out_of_range() {
        let surface = bilinear_patch();
        assert!(surface.insert_knot_u(2.0).is_err());
        assert!(surface.insert_knot_u(-0.5).is_err());
    }

    #[test]
    fn surface_insert_knot_v_out_of_range() {
        let surface = bilinear_patch();
        assert!(surface.insert_knot_v(2.0).is_err());
        assert!(surface.insert_knot_v(-0.5).is_err());
    }

    #[test]
    fn surface_insert_knots_u_preserves_shape() {
        let surface = bilinear_patch();
        let new_surface = surface.insert_knots_u(&[0.25, 0.75]).unwrap();
        assert_eq!(new_surface.num_u, surface.num_u + 2);
        let p1 = surface.point_at(0.5, 0.5);
        let p2 = new_surface.point_at(0.5, 0.5);
        assert!((p1.x - p2.x).abs() < 1e-10);
        assert!((p1.y - p2.y).abs() < 1e-10);
    }

    #[test]
    fn surface_insert_knots_v_preserves_shape() {
        let surface = bilinear_patch();
        let new_surface = surface.insert_knots_v(&[0.25, 0.75]).unwrap();
        assert_eq!(new_surface.num_v, surface.num_v + 2);
        let p1 = surface.point_at(0.5, 0.5);
        let p2 = new_surface.point_at(0.5, 0.5);
        assert!((p1.x - p2.x).abs() < 1e-10);
        assert!((p1.y - p2.y).abs() < 1e-10);
    }
}
