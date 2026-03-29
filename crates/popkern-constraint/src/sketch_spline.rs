//! Spline creation, modification, and deletion methods for [`Sketch`].
//!
//! Split from `sketch_ops.rs` to stay within the 500-line file limit.

use crate::constraint::ConstraintHandle;
use crate::elements::{PointHandle, Spline, SplineHandle};
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
