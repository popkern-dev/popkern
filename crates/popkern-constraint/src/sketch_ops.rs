//! Additional [`Sketch`] operations: driven values, dimensional value
//! modification, interactive drag, and element removal with cascade.
//!
//! Split from `lib.rs` to stay within the 500-line file limit.

use crate::constraint::{Constraint, ConstraintEntry, ConstraintHandle};
use crate::constraint::{ConstraintKind, ConstraintStatus};
use crate::elements::{ElementHandle, LineHandle, PointHandle, SplineHandle};
use crate::error::ConstraintError;
use crate::{equations, DragState, Sketch};

impl Sketch {
    /// Return the current geometric value that a driven constraint measures.
    ///
    /// For dimensional constraints (Distance, Angle, PointLineDistance) this
    /// evaluates the live geometric quantity. For other constraint kinds the
    /// residual magnitude is returned.
    ///
    /// Returns `None` if the handle is invalid or the referenced elements
    /// are missing.
    pub fn driven_value(&self, handle: ConstraintHandle) -> Option<f64> {
        let entry = self.constraints.get(handle)?;
        match &entry.constraint {
            Constraint::Distance(p1, p2, _) => {
                let (x1, y1) = self.point_position(*p1)?;
                let (x2, y2) = self.point_position(*p2)?;
                Some(((x2 - x1).powi(2) + (y2 - y1).powi(2)).sqrt())
            }
            Constraint::Angle(l1, l2, _) => self.compute_angle_value(*l1, *l2),
            Constraint::PointLineDistance(p, l, _) => self.compute_point_line_distance(*p, *l),
            Constraint::OrdinateX(p, _) => Some(self.point_position(*p)?.0),
            Constraint::OrdinateY(p, _) => Some(self.point_position(*p)?.1),
            Constraint::OffsetArc(a1, a2, _) => {
                let (_, _, r1, _, _) = self.arc_params(*a1)?;
                let (_, _, r2, _, _) = self.arc_params(*a2)?;
                Some((r1 - r2).abs())
            }
            _ => {
                let output = equations::evaluate_constraint(self, Some(handle), entry, 0);
                let norm_sq: f64 = output.residuals.iter().map(|r| r * r).sum();
                Some(norm_sq.sqrt())
            }
        }
    }

    /// Compute the signed angle between two lines.
    fn compute_angle_value(&self, l1: LineHandle, l2: LineHandle) -> Option<f64> {
        let line1 = self.lines.get(l1)?;
        let line2 = self.lines.get(l2)?;
        let (x11, y11) = self.point_position(line1.start)?;
        let (x12, y12) = self.point_position(line1.end)?;
        let (x21, y21) = self.point_position(line2.start)?;
        let (x22, y22) = self.point_position(line2.end)?;
        let dx1 = x12 - x11;
        let dy1 = y12 - y11;
        let dx2 = x22 - x21;
        let dy2 = y22 - y21;
        let cross = dx1 * dy2 - dy1 * dx2;
        let dot = dx1 * dx2 + dy1 * dy2;
        Some(cross.atan2(dot))
    }

    /// Compute the perpendicular distance from a point to a line.
    fn compute_point_line_distance(&self, p: PointHandle, l: LineHandle) -> Option<f64> {
        let (px, py) = self.point_position(p)?;
        let line = self.lines.get(l)?;
        let (x1, y1) = self.point_position(line.start)?;
        let (x2, y2) = self.point_position(line.end)?;
        let dx = x2 - x1;
        let dy = y2 - y1;
        let len = (dx * dx + dy * dy).sqrt();
        if len < 1e-15 {
            return Some(((px - x1).powi(2) + (py - y1).powi(2)).sqrt());
        }
        let signed = ((px - x1) * dy - (py - y1) * dx) / len;
        Some(signed.abs())
    }

    /// Modify the target value of a dimensional constraint.
    ///
    /// Only Distance, PointLineDistance, and Angle constraints have a
    /// modifiable value parameter. Returns an error for non-dimensional
    /// constraint kinds.
    pub fn set_constraint_value(
        &mut self,
        handle: ConstraintHandle,
        value: f64,
    ) -> Result<(), ConstraintError> {
        // Validate before taking snapshot.
        let entry = self
            .constraints
            .get(handle)
            .ok_or(ConstraintError::InvalidHandle)?;
        match &entry.constraint {
            Constraint::Distance(..)
            | Constraint::PointLineDistance(..)
            | Constraint::Angle(..)
            | Constraint::OffsetLine(..)
            | Constraint::OffsetArc(..)
            | Constraint::OrdinateX(..)
            | Constraint::OrdinateY(..) => {}
            _ => return Err(ConstraintError::NotDimensional),
        }
        if self.recording {
            self.push_undo_snapshot();
        }
        let entry = self
            .constraints
            .get_mut(handle)
            .ok_or(ConstraintError::InvalidHandle)?;
        match &mut entry.constraint {
            Constraint::Distance(_, _, d) => *d = value,
            Constraint::PointLineDistance(_, _, d) => *d = value,
            Constraint::Angle(_, _, a) => *a = value,
            Constraint::OffsetLine(_, _, d) => *d = value,
            Constraint::OffsetArc(_, _, d) => *d = value,
            Constraint::OrdinateX(_, v) => *v = value,
            Constraint::OrdinateY(_, v) => *v = value,
            _ => return Err(ConstraintError::NotDimensional),
        }
        Ok(())
    }

    /// Begin an interactive drag on a point.
    ///
    /// Saves the current parameter snapshot so subsequent
    /// [`drag_to`](Self::drag_to) calls can solve relative to the
    /// original state.
    pub fn begin_drag(&mut self, point: PointHandle) -> Result<(), ConstraintError> {
        if !self.points.contains_key(point) {
            return Err(ConstraintError::InvalidHandle);
        }
        self.drag_context = Some(DragState {
            point,
            saved_params: self.params.clone(),
            temp_constraint: None,
        });
        Ok(())
    }

    /// Move the dragged point to `(x, y)` and re-solve the sketch.
    ///
    /// Sets the point parameters to the target position (warm start), adds
    /// a temporary Fixed constraint, solves, then removes the constraint.
    /// The solver finds the nearest valid configuration to the drag target.
    pub fn drag_to(&mut self, x: f64, y: f64) -> Result<crate::SolveResult, ConstraintError> {
        let drag = self
            .drag_context
            .as_ref()
            .ok_or(ConstraintError::NoDragActive)?;
        let point = drag.point;

        let p = self
            .points
            .get(point)
            .ok_or(ConstraintError::InvalidHandle)?;
        let [ix, iy] = p.param_indices();
        self.params[ix] = x;
        self.params[iy] = y;

        let temp = self.constraints.insert(ConstraintEntry {
            constraint: Constraint::Fixed(point),
            kind: ConstraintKind::Driving,
            status: ConstraintStatus::Unknown,
            fixed_values: vec![x, y],
        });
        if let Some(drag) = self.drag_context.as_mut() {
            drag.temp_constraint = Some(temp);
        }

        let result = self.solve();

        self.constraints.remove(temp);
        if let Some(drag) = self.drag_context.as_mut() {
            drag.temp_constraint = None;
        }

        result
    }

    /// End the interactive drag, committing the current geometry.
    ///
    /// Does not restore the parameter snapshot -- the drag result is kept.
    pub fn end_drag(&mut self) {
        self.drag_context = None;
    }

    /// Add an ellipse defined by centre `(cx, cy)`, semi-axes `semi_a`
    /// and `semi_b`, and `rotation` angle (radians).
    ///
    /// The ellipse owns five consecutive entries in the parameter vector.
    pub fn add_ellipse(
        &mut self,
        cx: f64,
        cy: f64,
        semi_a: f64,
        semi_b: f64,
        rotation: f64,
    ) -> crate::EllipseHandle {
        if self.recording {
            self.push_undo_snapshot();
        }
        let param_start = self.params.len();
        self.params.push(cx);
        self.params.push(cy);
        self.params.push(semi_a);
        self.params.push(semi_b);
        self.params.push(rotation);
        self.ellipses.insert(crate::Ellipse {
            param_start,
            construction: false,
        })
    }

    /// Add a parabola defined by vertex `(vx, vy)`, `focal_dist`, and
    /// `rotation` angle (radians).
    ///
    /// The parabola owns four consecutive entries in the parameter vector.
    pub fn add_parabola(
        &mut self,
        vx: f64,
        vy: f64,
        focal_dist: f64,
        rotation: f64,
    ) -> crate::ParabolaHandle {
        if self.recording {
            self.push_undo_snapshot();
        }
        let param_start = self.params.len();
        self.params.push(vx);
        self.params.push(vy);
        self.params.push(focal_dist);
        self.params.push(rotation);
        self.parabolas.insert(crate::Parabola {
            param_start,
            construction: false,
        })
    }

    /// Add a hyperbola defined by centre `(cx, cy)`, semi-axes `semi_a`
    /// and `semi_b`, and `rotation` angle (radians).
    ///
    /// The hyperbola owns five consecutive entries in the parameter vector.
    pub fn add_hyperbola(
        &mut self,
        cx: f64,
        cy: f64,
        semi_a: f64,
        semi_b: f64,
        rotation: f64,
    ) -> crate::HyperbolaHandle {
        if self.recording {
            self.push_undo_snapshot();
        }
        let param_start = self.params.len();
        self.params.push(cx);
        self.params.push(cy);
        self.params.push(semi_a);
        self.params.push(semi_b);
        self.params.push(rotation);
        self.hyperbolas.insert(crate::Hyperbola {
            param_start,
            construction: false,
        })
    }

    /// Return the current ellipse parameters `(cx, cy, semi_a, semi_b,
    /// rotation)`, or `None` if the handle is invalid.
    pub fn ellipse_params(&self, h: crate::EllipseHandle) -> Option<(f64, f64, f64, f64, f64)> {
        let e = self.ellipses.get(h)?;
        let [icx, icy, ia, ib, ir] = e.param_indices();
        Some((
            self.params[icx],
            self.params[icy],
            self.params[ia],
            self.params[ib],
            self.params[ir],
        ))
    }

    /// Return the current parabola parameters `(vx, vy, focal_dist,
    /// rotation)`, or `None` if the handle is invalid.
    pub fn parabola_params(&self, h: crate::ParabolaHandle) -> Option<(f64, f64, f64, f64)> {
        let p = self.parabolas.get(h)?;
        let [ivx, ivy, ifd, ir] = p.param_indices();
        Some((
            self.params[ivx],
            self.params[ivy],
            self.params[ifd],
            self.params[ir],
        ))
    }

    /// Return the current hyperbola parameters `(cx, cy, semi_a, semi_b,
    /// rotation)`, or `None` if the handle is invalid.
    pub fn hyperbola_params(&self, h: crate::HyperbolaHandle) -> Option<(f64, f64, f64, f64, f64)> {
        let hp = self.hyperbolas.get(h)?;
        let [icx, icy, ia, ib, ir] = hp.param_indices();
        Some((
            self.params[icx],
            self.params[icy],
            self.params[ia],
            self.params[ib],
            self.params[ir],
        ))
    }

    /// Remove a geometric element and cascade-delete dependent elements
    /// and constraints.
    ///
    /// - **Point:** removes all lines that use this point as an endpoint
    ///   (recursively cascading their constraints), all constraints that
    ///   reference the point, then the point itself.
    /// - **Line:** removes all constraints that reference the line, then
    ///   the line (endpoint points are **not** removed).
    /// - **Arc:** removes all constraints that reference the arc, then
    ///   the arc.
    pub fn remove_element(&mut self, element: ElementHandle) {
        if self.recording {
            self.push_undo_snapshot();
        }
        match element {
            ElementHandle::Point(ph) => self.remove_point(ph),
            ElementHandle::Line(lh) => self.remove_line(lh),
            ElementHandle::Arc(ah) => self.remove_arc(ah),
            ElementHandle::Ellipse(eh) => self.remove_ellipse(eh),
            ElementHandle::Parabola(ph) => self.remove_parabola(ph),
            ElementHandle::Hyperbola(hh) => self.remove_hyperbola(hh),
            ElementHandle::Spline(sh) => self.remove_spline(sh),
        }
    }

    /// Remove a point and cascade to dependent lines, splines, and
    /// constraints.
    ///
    /// Splines that list this point as a control point are removed first
    /// (a spline is invalid without all its control points).  Their
    /// constraints are cascade-deleted by `remove_spline`.
    fn remove_point(&mut self, ph: PointHandle) {
        let dependent_lines: Vec<LineHandle> = self
            .lines
            .iter()
            .filter(|(_, l)| l.start == ph || l.end == ph)
            .map(|(h, _)| h)
            .collect();
        for lh in dependent_lines {
            self.remove_line(lh);
        }
        let dependent_splines: Vec<SplineHandle> = self
            .splines
            .iter()
            .filter(|(_, sp)| sp.control_points.contains(&ph))
            .map(|(h, _)| h)
            .collect();
        for sh in dependent_splines {
            self.remove_spline(sh);
        }
        let to_remove: Vec<ConstraintHandle> = self
            .constraints
            .iter()
            .filter(|(_, e)| e.constraint.references_point(ph))
            .map(|(h, _)| h)
            .collect();
        for ch in to_remove {
            self.constraints.remove(ch);
        }
        self.points.remove(ph);
    }

    /// Remove a line and its dependent constraints (but not endpoint
    /// points).
    fn remove_line(&mut self, lh: LineHandle) {
        let to_remove: Vec<ConstraintHandle> = self
            .constraints
            .iter()
            .filter(|(_, e)| e.constraint.references_line(lh))
            .map(|(h, _)| h)
            .collect();
        for ch in to_remove {
            self.constraints.remove(ch);
        }
        self.lines.remove(lh);
    }

    /// Remove an arc and its dependent constraints.
    fn remove_arc(&mut self, ah: crate::ArcHandle) {
        let to_remove: Vec<ConstraintHandle> = self
            .constraints
            .iter()
            .filter(|(_, e)| e.constraint.references_arc(ah))
            .map(|(h, _)| h)
            .collect();
        for ch in to_remove {
            self.constraints.remove(ch);
        }
        self.arcs.remove(ah);
    }

    /// Remove an ellipse and its dependent constraints.
    fn remove_ellipse(&mut self, eh: crate::EllipseHandle) {
        let to_remove: Vec<ConstraintHandle> = self
            .constraints
            .iter()
            .filter(|(_, e)| e.constraint.references_ellipse(eh))
            .map(|(h, _)| h)
            .collect();
        for ch in to_remove {
            self.constraints.remove(ch);
        }
        self.ellipses.remove(eh);
    }

    /// Remove a parabola and its dependent constraints.
    fn remove_parabola(&mut self, ph: crate::ParabolaHandle) {
        let to_remove: Vec<ConstraintHandle> = self
            .constraints
            .iter()
            .filter(|(_, e)| e.constraint.references_parabola(ph))
            .map(|(h, _)| h)
            .collect();
        for ch in to_remove {
            self.constraints.remove(ch);
        }
        self.parabolas.remove(ph);
    }

    /// Remove a hyperbola and its dependent constraints.
    fn remove_hyperbola(&mut self, hh: crate::HyperbolaHandle) {
        let to_remove: Vec<ConstraintHandle> = self
            .constraints
            .iter()
            .filter(|(_, e)| e.constraint.references_hyperbola(hh))
            .map(|(h, _)| h)
            .collect();
        for ch in to_remove {
            self.constraints.remove(ch);
        }
        self.hyperbolas.remove(hh);
    }
}
