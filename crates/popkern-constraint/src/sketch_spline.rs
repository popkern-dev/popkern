//! Spline creation, modification, and deletion methods for [`Sketch`],
//! plus the intersection-point convenience method.
//!
//! Split from `sketch_ops.rs` to stay within the 500-line file limit.

use crate::constraint::{Constraint, ConstraintHandle};
use crate::elements::{ElementHandle, PointHandle, Spline, SplineHandle};
use crate::error::ConstraintError;
use crate::Sketch;

impl Sketch {
    /// Add a B-spline defined by the given control points and degree.
    ///
    /// A clamped uniform knot vector is computed automatically.  All
    /// weights default to `1.0` (non-rational).  The control points must
    /// already exist in the sketch; they are referenced by handle, not
    /// owned.
    pub fn add_spline(&mut self, control_points: &[PointHandle], degree: usize) -> SplineHandle {
        if self.recording {
            self.push_undo_snapshot();
        }
        let n = control_points.len();
        let weights = vec![1.0; n];
        let knots = clamped_uniform_knots(n, degree);
        self.splines.insert(Spline {
            control_points: control_points.to_vec(),
            degree,
            knots,
            weights,
            periodic: false,
            construction: false,
        })
    }

    /// Override the rational weight for one control point.
    ///
    /// Returns [`ConstraintError::InvalidHandle`] if `spline` does not
    /// exist or `cp_index` is out of bounds.
    pub fn set_spline_weight(
        &mut self,
        spline: SplineHandle,
        cp_index: usize,
        weight: f64,
    ) -> Result<(), ConstraintError> {
        if self.recording {
            self.push_undo_snapshot();
        }
        let sp = self
            .splines
            .get_mut(spline)
            .ok_or(ConstraintError::InvalidHandle)?;
        if cp_index >= sp.weights.len() {
            return Err(ConstraintError::InvalidHandle);
        }
        sp.weights[cp_index] = weight;
        Ok(())
    }

    /// Override a single knot value.
    ///
    /// Returns [`ConstraintError::InvalidHandle`] if `spline` does not
    /// exist or `knot_index` is out of bounds.
    pub fn set_spline_knot(
        &mut self,
        spline: SplineHandle,
        knot_index: usize,
        value: f64,
    ) -> Result<(), ConstraintError> {
        if self.recording {
            self.push_undo_snapshot();
        }
        let sp = self
            .splines
            .get_mut(spline)
            .ok_or(ConstraintError::InvalidHandle)?;
        if knot_index >= sp.knots.len() {
            return Err(ConstraintError::InvalidHandle);
        }
        sp.knots[knot_index] = value;
        Ok(())
    }

    /// Return the ordered slice of control-point handles for a spline, or
    /// `None` if the handle is invalid.
    pub fn spline_control_points(&self, spline: SplineHandle) -> Option<&[PointHandle]> {
        self.splines
            .get(spline)
            .map(|sp| sp.control_points.as_slice())
    }

    /// Create a B-spline that passes through the given points.
    ///
    /// Internally computes control points via cubic interpolation, adds
    /// them as sketch points, creates the spline, and constrains each
    /// through-point to lie on the spline at its interpolation parameter.
    ///
    /// Returns `(spline_handle, through_point_handles)`.
    pub fn add_spline_through_points(
        &mut self,
        through_points: &[(f64, f64)],
    ) -> Result<(SplineHandle, Vec<PointHandle>), ConstraintError> {
        use crate::spline_interp::interpolate_cubic;

        let (cp_positions, knots, t_values) =
            interpolate_cubic(through_points).map_err(ConstraintError::SolveFailed)?;

        // Create point handles for through-points.
        let tp_handles: Vec<PointHandle> = through_points
            .iter()
            .map(|&(x, y)| self.add_point(x, y))
            .collect();

        // Create point handles for computed control points.
        let cp_handles: Vec<PointHandle> = cp_positions
            .iter()
            .map(|&(x, y)| self.add_point(x, y))
            .collect();

        let degree = if cp_handles.len() >= 4 {
            3
        } else {
            cp_handles.len() - 1
        };

        // Create the spline.
        let sh = self.splines.insert(Spline {
            control_points: cp_handles.clone(),
            degree,
            knots,
            weights: vec![1.0; cp_handles.len()],
            periodic: false,
            construction: false,
        });

        // Add PointOnSpline constraints and set aux params to t_i.
        for (i, &tph) in tp_handles.iter().enumerate() {
            let ch = self.add_constraint(Constraint::PointOnSpline(tph, sh));
            if let Some(&aux_start) = self.aux_param_map.get(&ch) {
                self.params[aux_start] = t_values[i];
            }
        }

        Ok((sh, tp_handles))
    }

    /// Create a closed (periodic) B-spline that passes through the given
    /// points.
    ///
    /// Similar to [`add_spline_through_points`](Self::add_spline_through_points)
    /// but the spline wraps smoothly from the last point back to the first.
    ///
    /// Returns `(spline_handle, through_point_handles)`.
    pub fn add_closed_spline_through_points(
        &mut self,
        through_points: &[(f64, f64)],
    ) -> Result<(SplineHandle, Vec<PointHandle>), ConstraintError> {
        use crate::spline_interp::interpolate_cubic_periodic;

        let (cp_positions, knots, t_values) =
            interpolate_cubic_periodic(through_points).map_err(ConstraintError::SolveFailed)?;

        let tp_handles: Vec<PointHandle> = through_points
            .iter()
            .map(|&(x, y)| self.add_point(x, y))
            .collect();

        let cp_handles: Vec<PointHandle> = cp_positions
            .iter()
            .map(|&(x, y)| self.add_point(x, y))
            .collect();

        let sh = self.splines.insert(Spline {
            control_points: cp_handles.clone(),
            degree: 3,
            knots,
            weights: vec![1.0; cp_handles.len()],
            periodic: true,
            construction: false,
        });

        for (i, &tph) in tp_handles.iter().enumerate() {
            let ch = self.add_constraint(Constraint::PointOnSpline(tph, sh));
            if let Some(&aux_start) = self.aux_param_map.get(&ch) {
                self.params[aux_start] = t_values[i];
            }
        }

        Ok((sh, tp_handles))
    }

    /// Remove a spline and all constraints that reference it.
    ///
    /// Control points are **not** removed — they are shared elements that
    /// may be referenced by other geometry.
    pub(crate) fn remove_spline(&mut self, sh: SplineHandle) {
        let to_remove: Vec<ConstraintHandle> = self
            .constraints
            .iter()
            .filter(|(_, e)| e.constraint.references_spline(sh))
            .map(|(h, _)| h)
            .collect();
        for ch in to_remove {
            self.constraints.remove(ch);
        }
        self.splines.remove(sh);
    }

    /// Constrain a point to lie at the intersection of two curves.
    ///
    /// This is a convenience method that adds two point-on-curve constraints
    /// to the same point — one for each curve. The solver naturally finds
    /// the intersection. Supported curve types: Line, Arc, Ellipse, Parabola,
    /// Hyperbola, Spline.
    ///
    /// Returns the two constraint handles created.
    pub fn constrain_intersection(
        &mut self,
        point: PointHandle,
        curve1: ElementHandle,
        curve2: ElementHandle,
    ) -> Result<[ConstraintHandle; 2], ConstraintError> {
        let c1 = self.add_point_on_curve(point, curve1)?;
        let c2 = self.add_point_on_curve(point, curve2)?;
        Ok([c1, c2])
    }

    /// Add a point-on-curve constraint, dispatching on the curve element type.
    fn add_point_on_curve(
        &mut self,
        point: PointHandle,
        curve: ElementHandle,
    ) -> Result<ConstraintHandle, ConstraintError> {
        let constraint = match curve {
            ElementHandle::Line(lh) => {
                if !self.lines.contains_key(lh) {
                    return Err(ConstraintError::InvalidHandle);
                }
                Constraint::PointOnLine(point, lh)
            }
            ElementHandle::Arc(ah) => {
                if !self.arcs.contains_key(ah) {
                    return Err(ConstraintError::InvalidHandle);
                }
                Constraint::PointOnArc(point, ah)
            }
            ElementHandle::Ellipse(eh) => {
                if !self.ellipses.contains_key(eh) {
                    return Err(ConstraintError::InvalidHandle);
                }
                Constraint::PointOnEllipse(point, eh)
            }
            ElementHandle::Parabola(ph) => {
                if !self.parabolas.contains_key(ph) {
                    return Err(ConstraintError::InvalidHandle);
                }
                Constraint::PointOnParabola(point, ph)
            }
            ElementHandle::Hyperbola(hh) => {
                if !self.hyperbolas.contains_key(hh) {
                    return Err(ConstraintError::InvalidHandle);
                }
                Constraint::PointOnHyperbola(point, hh)
            }
            ElementHandle::Spline(sh) => {
                if !self.splines.contains_key(sh) {
                    return Err(ConstraintError::InvalidHandle);
                }
                Constraint::PointOnSpline(point, sh)
            }
            ElementHandle::Point(_) => {
                return Err(ConstraintError::InvalidHandle);
            }
        };
        Ok(self.add_constraint(constraint))
    }
}

/// Compute a clamped uniform knot vector for a B-spline with `n_cp`
/// control points and the given `degree`.
///
/// The vector has `n_cp + degree + 1` entries: `degree + 1` zeros at the
/// start, `degree + 1` ones at the end, and interior knots evenly spaced
/// in (0, 1).
fn clamped_uniform_knots(n_cp: usize, degree: usize) -> Vec<f64> {
    let m = n_cp + degree + 1;
    let mut knots = Vec::with_capacity(m);
    let n_interior = m.saturating_sub(2 * (degree + 1));
    knots.extend(std::iter::repeat_n(0.0, degree + 1));
    for i in 1..=n_interior {
        knots.push(i as f64 / (n_interior + 1) as f64);
    }
    knots.extend(std::iter::repeat_n(1.0, degree + 1));
    knots
}
