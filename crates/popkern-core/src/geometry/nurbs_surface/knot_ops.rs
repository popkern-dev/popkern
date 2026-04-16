//! Knot insertion and removal methods for NurbsSurface.

use super::NurbsSurface;
use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::nurbs_ops;
use crate::geometry::surface::SurfaceEval;

impl NurbsSurface {
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
            for (i_u, (pt, wt)) in np.into_iter().zip(nw).enumerate() {
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

    /// Attempt to remove the knot at `knot_index` in the u direction.
    ///
    /// The removal is applied independently to every v-column (fixed `iv`, all `iu` values).
    /// Returns `Some(new_surface)` only if **all** columns succeed within `tolerance`.
    /// Returns `None` if any column cannot be removed within the tolerance.
    pub fn remove_knot_u(&self, knot_index: usize, tolerance: f64) -> Option<NurbsSurface> {
        let mut new_cols: Vec<(Vec<Point3>, Vec<f64>)> = Vec::with_capacity(self.num_v);
        let mut new_knots_u: Option<Vec<f64>> = None;

        for iv in 0..self.num_v {
            let col_pts: Vec<Point3> = (0..self.num_u).map(|iu| *self.cp(iu, iv)).collect();
            let col_wts: Vec<f64> = (0..self.num_u).map(|iu| self.wt(iu, iv)).collect();

            let (nk, np, nw) = nurbs_ops::knot_removal::remove_knot_1d(
                self.degree_u,
                &self.knots_u,
                &col_pts,
                &col_wts,
                knot_index,
                tolerance,
            )?;

            new_knots_u.get_or_insert(nk);
            new_cols.push((np, nw));
        }

        let new_knots_u = new_knots_u?;
        let new_num_u = self.num_u - 1;

        // Re-assemble row-major grid from per-column results
        let total = new_num_u * self.num_v;
        let mut control_points = Vec::with_capacity(total);
        let mut weights = Vec::with_capacity(total);
        for iu in 0..new_num_u {
            for col in new_cols.iter().take(self.num_v) {
                control_points.push(col.0[iu]);
                weights.push(col.1[iu]);
            }
        }

        Some(NurbsSurface {
            degree_u: self.degree_u,
            degree_v: self.degree_v,
            knots_u: new_knots_u,
            knots_v: self.knots_v.clone(),
            control_points,
            num_u: new_num_u,
            num_v: self.num_v,
            weights,
        })
    }

    /// Attempt to remove the knot at `knot_index` in the v direction.
    ///
    /// The removal is applied independently to every u-row (fixed `iu`, all `iv` values).
    /// Returns `Some(new_surface)` only if **all** rows succeed within `tolerance`.
    /// Returns `None` if any row cannot be removed within the tolerance.
    pub fn remove_knot_v(&self, knot_index: usize, tolerance: f64) -> Option<NurbsSurface> {
        let mut new_rows: Vec<(Vec<Point3>, Vec<f64>)> = Vec::with_capacity(self.num_u);
        let mut new_knots_v: Option<Vec<f64>> = None;

        for iu in 0..self.num_u {
            let row_pts: Vec<Point3> = (0..self.num_v).map(|iv| *self.cp(iu, iv)).collect();
            let row_wts: Vec<f64> = (0..self.num_v).map(|iv| self.wt(iu, iv)).collect();

            let (nk, np, nw) = nurbs_ops::knot_removal::remove_knot_1d(
                self.degree_v,
                &self.knots_v,
                &row_pts,
                &row_wts,
                knot_index,
                tolerance,
            )?;

            new_knots_v.get_or_insert(nk);
            new_rows.push((np, nw));
        }

        let new_knots_v = new_knots_v?;
        let new_num_v = self.num_v - 1;

        // Re-assemble row-major grid from per-row results
        let total = self.num_u * new_num_v;
        let mut control_points = Vec::with_capacity(total);
        let mut weights = Vec::with_capacity(total);
        for row in new_rows.iter().take(self.num_u) {
            for iv in 0..new_num_v {
                control_points.push(row.0[iv]);
                weights.push(row.1[iv]);
            }
        }

        Some(NurbsSurface {
            degree_u: self.degree_u,
            degree_v: self.degree_v,
            knots_u: self.knots_u.clone(),
            knots_v: new_knots_v,
            control_points,
            num_u: self.num_u,
            num_v: new_num_v,
            weights,
        })
    }

    /// Remove the u-knot at `knot_index` up to `times` times.
    ///
    /// Returns `(new_surface, actual_count)`. Stops early if removal fails for any column.
    pub fn remove_knot_u_times(
        &self,
        knot_index: usize,
        times: usize,
        tolerance: f64,
    ) -> (NurbsSurface, usize) {
        let mut current = self.clone();
        let mut count = 0;
        for _ in 0..times {
            match current.remove_knot_u(knot_index, tolerance) {
                Some(next) => {
                    current = next;
                    count += 1;
                    if knot_index >= current.knots_u.len() {
                        break;
                    }
                }
                None => break,
            }
        }
        (current, count)
    }

    /// Remove the v-knot at `knot_index` up to `times` times.
    ///
    /// Returns `(new_surface, actual_count)`. Stops early if removal fails for any row.
    pub fn remove_knot_v_times(
        &self,
        knot_index: usize,
        times: usize,
        tolerance: f64,
    ) -> (NurbsSurface, usize) {
        let mut current = self.clone();
        let mut count = 0;
        for _ in 0..times {
            match current.remove_knot_v(knot_index, tolerance) {
                Some(next) => {
                    current = next;
                    count += 1;
                    if knot_index >= current.knots_v.len() {
                        break;
                    }
                }
                None => break,
            }
        }
        (current, count)
    }
}

#[cfg(test)]
#[path = "knot_ops_tests.rs"]
mod knot_ops_tests;
