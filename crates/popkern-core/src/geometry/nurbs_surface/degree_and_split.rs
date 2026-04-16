//! Degree elevation and splitting methods for NurbsSurface.

use super::NurbsSurface;
use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::nurbs_ops;

/// A split segment: (left_pts, left_wts, right_pts, right_wts).
type SplitSegment = (Vec<Point3>, Vec<f64>, Vec<Point3>, Vec<f64>);

impl NurbsSurface {
    /// Elevate degree by 1 in the u direction. Shape unchanged.
    pub fn elevate_degree_u(&self) -> NurbsSurface {
        // Apply 1D elevation per column (fixed v_idx, varying i_u)
        let new_num_v = self.num_v;
        let mut new_knots_u = Vec::new();
        let mut new_num_u = 0;
        let mut col_results: Vec<(Vec<Point3>, Vec<f64>)> = Vec::with_capacity(self.num_v);

        for v_idx in 0..self.num_v {
            let col_pts: Vec<Point3> = (0..self.num_u)
                .map(|i_u| self.control_points[i_u * self.num_v + v_idx])
                .collect();
            let col_wts: Vec<f64> = (0..self.num_u)
                .map(|i_u| self.weights[i_u * self.num_v + v_idx])
                .collect();
            let (new_deg, nk, np, nw) = nurbs_ops::degree_elevation::elevate_degree_1d(
                self.degree_u,
                &self.knots_u,
                &col_pts,
                &col_wts,
            );
            if v_idx == 0 {
                new_knots_u = nk;
                new_num_u = np.len();
                let _ = new_deg; // degree_u + 1
            }
            col_results.push((np, nw));
        }

        let mut new_pts = vec![Point3::new(0.0, 0.0, 0.0); new_num_u * new_num_v];
        let mut new_wts = vec![0.0f64; new_num_u * new_num_v];
        for (v_idx, (col_pts, col_wts)) in col_results.into_iter().enumerate() {
            for (i_u, (pt, wt)) in col_pts.into_iter().zip(col_wts).enumerate() {
                new_pts[i_u * new_num_v + v_idx] = pt;
                new_wts[i_u * new_num_v + v_idx] = wt;
            }
        }

        NurbsSurface {
            degree_u: self.degree_u + 1,
            degree_v: self.degree_v,
            knots_u: new_knots_u,
            knots_v: self.knots_v.clone(),
            control_points: new_pts,
            num_u: new_num_u,
            num_v: new_num_v,
            weights: new_wts,
        }
    }

    /// Elevate degree by 1 in the v direction. Shape unchanged.
    pub fn elevate_degree_v(&self) -> NurbsSurface {
        // Apply 1D elevation per row (fixed i_u, varying v)
        let new_num_u = self.num_u;
        let mut new_knots_v = Vec::new();
        let mut new_num_v = 0;
        let mut new_pts_flat: Vec<Point3> = Vec::new();
        let mut new_wts_flat: Vec<f64> = Vec::new();

        for i_u in 0..self.num_u {
            let row_pts: Vec<Point3> = (0..self.num_v)
                .map(|i_v| self.control_points[i_u * self.num_v + i_v])
                .collect();
            let row_wts: Vec<f64> = (0..self.num_v)
                .map(|i_v| self.weights[i_u * self.num_v + i_v])
                .collect();
            let (new_deg, nk, np, nw) = nurbs_ops::degree_elevation::elevate_degree_1d(
                self.degree_v,
                &self.knots_v,
                &row_pts,
                &row_wts,
            );
            if i_u == 0 {
                new_knots_v = nk;
                new_num_v = np.len();
                let _ = new_deg; // degree_v + 1
            }
            new_pts_flat.extend(np);
            new_wts_flat.extend(nw);
        }

        NurbsSurface {
            degree_u: self.degree_u,
            degree_v: self.degree_v + 1,
            knots_u: self.knots_u.clone(),
            knots_v: new_knots_v,
            control_points: new_pts_flat,
            num_u: new_num_u,
            num_v: new_num_v,
            weights: new_wts_flat,
        }
    }

    /// Elevate degree in u to `target`. Returns `Err` if `target <= current degree_u`.
    pub fn elevate_degree_u_to(&self, target: usize) -> Result<NurbsSurface, CoreError> {
        if target <= self.degree_u {
            return Err(CoreError::PreconditionFailed(format!(
                "target degree_u {} must be > current degree_u {}",
                target, self.degree_u
            )));
        }
        let mut surface = self.clone();
        for _ in 0..(target - self.degree_u) {
            surface = surface.elevate_degree_u();
        }
        Ok(surface)
    }

    /// Elevate degree in v to `target`. Returns `Err` if `target <= current degree_v`.
    pub fn elevate_degree_v_to(&self, target: usize) -> Result<NurbsSurface, CoreError> {
        if target <= self.degree_v {
            return Err(CoreError::PreconditionFailed(format!(
                "target degree_v {} must be > current degree_v {}",
                target, self.degree_v
            )));
        }
        let mut surface = self.clone();
        for _ in 0..(target - self.degree_v) {
            surface = surface.elevate_degree_v();
        }
        Ok(surface)
    }

    /// Split the surface at parameter `u` into two independent surfaces.
    ///
    /// Returns `(left, right)` where `left` covers `[u_start, u]` and `right` covers `[u, u_end]`.
    /// The parameter `u` must be strictly inside the u parameter range.
    pub fn split_at_u(&self, u: f64) -> Result<(NurbsSurface, NurbsSurface), CoreError> {
        // For each column (fixed v_idx, varying i_u), split the column in the u direction.
        // Collect left and right column results, then reassemble into two surfaces.
        let mut left_knots_u: Option<Vec<f64>> = None;
        let mut right_knots_u: Option<Vec<f64>> = None;
        let mut left_num_u = 0usize;
        let mut right_num_u = 0usize;

        // col_splits[v_idx] = (left_col_pts, left_col_wts, right_col_pts, right_col_wts)
        let mut col_splits: Vec<SplitSegment> = Vec::with_capacity(self.num_v);

        for v_idx in 0..self.num_v {
            let col_pts: Vec<Point3> = (0..self.num_u)
                .map(|i_u| self.control_points[i_u * self.num_v + v_idx])
                .collect();
            let col_wts: Vec<f64> = (0..self.num_u)
                .map(|i_u| self.weights[i_u * self.num_v + v_idx])
                .collect();
            let (lk, lp, lw, rk, rp, rw) = nurbs_ops::splitting::split_1d(
                self.degree_u,
                &self.knots_u,
                &col_pts,
                &col_wts,
                u,
            )?;
            if v_idx == 0 {
                left_num_u = lp.len();
                right_num_u = rp.len();
                left_knots_u = Some(lk);
                right_knots_u = Some(rk);
            }
            col_splits.push((lp, lw, rp, rw));
        }

        let left_knots_u = left_knots_u.unwrap_or_default();
        let right_knots_u = right_knots_u.unwrap_or_default();

        // Reassemble left surface (row-major: flat index = i_u * num_v + v_idx)
        let mut left_pts = vec![Point3::new(0.0, 0.0, 0.0); left_num_u * self.num_v];
        let mut left_wts = vec![0.0f64; left_num_u * self.num_v];
        for (v_idx, (lp, lw, _, _)) in col_splits.iter().enumerate() {
            for (i_u, (&pt, &wt)) in lp.iter().zip(lw.iter()).enumerate() {
                left_pts[i_u * self.num_v + v_idx] = pt;
                left_wts[i_u * self.num_v + v_idx] = wt;
            }
        }

        // Reassemble right surface
        let mut right_pts = vec![Point3::new(0.0, 0.0, 0.0); right_num_u * self.num_v];
        let mut right_wts = vec![0.0f64; right_num_u * self.num_v];
        for (v_idx, (_, _, rp, rw)) in col_splits.iter().enumerate() {
            for (i_u, (&pt, &wt)) in rp.iter().zip(rw.iter()).enumerate() {
                right_pts[i_u * self.num_v + v_idx] = pt;
                right_wts[i_u * self.num_v + v_idx] = wt;
            }
        }

        Ok((
            NurbsSurface {
                degree_u: self.degree_u,
                degree_v: self.degree_v,
                knots_u: left_knots_u,
                knots_v: self.knots_v.clone(),
                control_points: left_pts,
                num_u: left_num_u,
                num_v: self.num_v,
                weights: left_wts,
            },
            NurbsSurface {
                degree_u: self.degree_u,
                degree_v: self.degree_v,
                knots_u: right_knots_u,
                knots_v: self.knots_v.clone(),
                control_points: right_pts,
                num_u: right_num_u,
                num_v: self.num_v,
                weights: right_wts,
            },
        ))
    }

    /// Split the surface at parameter `v` into two independent surfaces.
    ///
    /// Returns `(bottom, top)` where `bottom` covers `[v_start, v]` and `top` covers `[v, v_end]`.
    /// The parameter `v` must be strictly inside the v parameter range.
    pub fn split_at_v(&self, v: f64) -> Result<(NurbsSurface, NurbsSurface), CoreError> {
        // For each row (fixed i_u, varying v), split the row in the v direction.
        // Collect bottom and top row results, then reassemble into two surfaces.
        let mut bottom_knots_v: Option<Vec<f64>> = None;
        let mut top_knots_v: Option<Vec<f64>> = None;
        let mut bottom_num_v = 0usize;
        let mut top_num_v = 0usize;

        // row_splits[i_u] = (bottom_row_pts, bottom_row_wts, top_row_pts, top_row_wts)
        let mut row_splits: Vec<SplitSegment> = Vec::with_capacity(self.num_u);

        for i_u in 0..self.num_u {
            let row_pts: Vec<Point3> = (0..self.num_v)
                .map(|i_v| self.control_points[i_u * self.num_v + i_v])
                .collect();
            let row_wts: Vec<f64> = (0..self.num_v)
                .map(|i_v| self.weights[i_u * self.num_v + i_v])
                .collect();
            let (lk, lp, lw, rk, rp, rw) = nurbs_ops::splitting::split_1d(
                self.degree_v,
                &self.knots_v,
                &row_pts,
                &row_wts,
                v,
            )?;
            if i_u == 0 {
                bottom_num_v = lp.len();
                top_num_v = rp.len();
                bottom_knots_v = Some(lk);
                top_knots_v = Some(rk);
            }
            row_splits.push((lp, lw, rp, rw));
        }

        let bottom_knots_v = bottom_knots_v.unwrap_or_default();
        let top_knots_v = top_knots_v.unwrap_or_default();

        // Reassemble bottom surface (row-major: each row i_u has bottom_num_v points)
        let mut bottom_pts = Vec::with_capacity(self.num_u * bottom_num_v);
        let mut bottom_wts = Vec::with_capacity(self.num_u * bottom_num_v);
        for (lp, lw, _, _) in &row_splits {
            bottom_pts.extend_from_slice(lp);
            bottom_wts.extend_from_slice(lw);
        }

        // Reassemble top surface
        let mut top_pts = Vec::with_capacity(self.num_u * top_num_v);
        let mut top_wts = Vec::with_capacity(self.num_u * top_num_v);
        for (_, _, rp, rw) in &row_splits {
            top_pts.extend_from_slice(rp);
            top_wts.extend_from_slice(rw);
        }

        Ok((
            NurbsSurface {
                degree_u: self.degree_u,
                degree_v: self.degree_v,
                knots_u: self.knots_u.clone(),
                knots_v: bottom_knots_v,
                control_points: bottom_pts,
                num_u: self.num_u,
                num_v: bottom_num_v,
                weights: bottom_wts,
            },
            NurbsSurface {
                degree_u: self.degree_u,
                degree_v: self.degree_v,
                knots_u: self.knots_u.clone(),
                knots_v: top_knots_v,
                control_points: top_pts,
                num_u: self.num_u,
                num_v: top_num_v,
                weights: top_wts,
            },
        ))
    }
}
