pub(crate) mod analysis;
pub mod constraint;
pub mod elements;
pub(crate) mod equations;
pub mod error;
pub(crate) mod jacobian;
mod solver;

#[cfg(test)]
#[path = "solver_tests.rs"]
mod solver_tests;

pub use constraint::{
    Constraint, ConstraintEntry, ConstraintHandle, ConstraintKind, ConstraintStatus,
};
pub use elements::{Arc, ArcHandle, ElementHandle, Line, LineHandle, Point, PointHandle};
pub use error::{ConstraintError, FreedomInfo, SolveResult, SolveStatus};
pub use solver::SolverOptions;

use slotmap::SlotMap;

/// Transient drag state used during interactive point dragging.
#[derive(Debug, Clone)]
pub(crate) struct DragState {
    /// The point being dragged.
    point: PointHandle,
    /// Parameter snapshot taken at the start of the drag (available for
    /// future undo support).
    #[allow(dead_code)]
    saved_params: Vec<f64>,
    /// Temporary Fixed constraint added during each `drag_to` step.
    temp_constraint: Option<ConstraintHandle>,
}

/// A 2D parametric sketch: a collection of geometric elements connected by
/// constraints.
pub struct Sketch {
    /// Flat parameter vector shared by all elements.
    pub(crate) params: Vec<f64>,
    /// Point elements, keyed by [`PointHandle`].
    pub(crate) points: SlotMap<PointHandle, Point>,
    /// Line elements, keyed by [`LineHandle`].
    pub(crate) lines: SlotMap<LineHandle, Line>,
    /// Arc elements, keyed by [`ArcHandle`].
    pub(crate) arcs: SlotMap<ArcHandle, Arc>,
    /// Constraints, keyed by [`ConstraintHandle`].
    pub(crate) constraints: SlotMap<ConstraintHandle, ConstraintEntry>,
    /// Active drag context (transient, not serialised).
    pub(crate) drag_context: Option<DragState>,
}

impl Sketch {
    /// Create a new, empty sketch.
    pub fn new() -> Self {
        Self {
            params: Vec::new(),
            points: SlotMap::with_key(),
            lines: SlotMap::with_key(),
            arcs: SlotMap::with_key(),
            constraints: SlotMap::with_key(),
            drag_context: None,
        }
    }

    /// Add a free point at `(x, y)`.
    ///
    /// The point owns two consecutive entries in the parameter vector.
    pub fn add_point(&mut self, x: f64, y: f64) -> PointHandle {
        let param_start = self.params.len();
        self.params.push(x);
        self.params.push(y);
        self.points.insert(Point {
            param_start,
            construction: false,
        })
    }

    /// Add a line segment between two existing points.
    ///
    /// Lines share the parameter storage of their endpoint points; no new
    /// parameters are allocated.
    pub fn add_line(&mut self, start: PointHandle, end: PointHandle) -> LineHandle {
        self.lines.insert(Line {
            start,
            end,
            construction: false,
        })
    }

    /// Add an arc defined by centre `(cx, cy)`, radius `r`, and angular
    /// span `[start_angle, end_angle]` (radians).
    ///
    /// The arc owns five consecutive entries in the parameter vector.
    pub fn add_arc(
        &mut self,
        cx: f64,
        cy: f64,
        r: f64,
        start_angle: f64,
        end_angle: f64,
    ) -> ArcHandle {
        let param_start = self.params.len();
        self.params.push(cx);
        self.params.push(cy);
        self.params.push(r);
        self.params.push(start_angle);
        self.params.push(end_angle);
        self.arcs.insert(Arc {
            param_start,
            construction: false,
        })
    }

    /// Set the construction flag for any element (construction geometry is
    /// used for reference but not exported as model edges).
    pub fn set_construction(&mut self, element: ElementHandle, value: bool) {
        match element {
            ElementHandle::Point(h) => {
                if let Some(p) = self.points.get_mut(h) {
                    p.construction = value;
                }
            }
            ElementHandle::Line(h) => {
                if let Some(l) = self.lines.get_mut(h) {
                    l.construction = value;
                }
            }
            ElementHandle::Arc(h) => {
                if let Some(a) = self.arcs.get_mut(h) {
                    a.construction = value;
                }
            }
        }
    }

    /// Add a driving constraint to the sketch.
    ///
    /// Driving constraints determine the geometry (e.g. a distance fixes a
    /// length).
    pub fn add_constraint(&mut self, constraint: Constraint) -> ConstraintHandle {
        let fixed_values = self.capture_fixed_values(&constraint);
        self.constraints.insert(ConstraintEntry {
            constraint,
            kind: ConstraintKind::Driving,
            status: ConstraintStatus::Unknown,
            fixed_values,
        })
    }

    /// Add a driven constraint to the sketch.
    ///
    /// Driven constraints are display-only annotations; they do not influence
    /// the solver.
    pub fn add_driven(&mut self, constraint: Constraint) -> ConstraintHandle {
        let fixed_values = self.capture_fixed_values(&constraint);
        self.constraints.insert(ConstraintEntry {
            constraint,
            kind: ConstraintKind::Driven,
            status: ConstraintStatus::Unknown,
            fixed_values,
        })
    }

    /// Remove a constraint from the sketch.
    ///
    /// Does nothing if the handle is no longer valid.
    pub fn remove_constraint(&mut self, handle: ConstraintHandle) {
        self.constraints.remove(handle);
    }

    /// Return the current `(x, y)` position of a point, or `None` if the
    /// handle is invalid.
    pub fn point_position(&self, handle: PointHandle) -> Option<(f64, f64)> {
        let p = self.points.get(handle)?;
        let [xi, yi] = p.param_indices();
        Some((self.params[xi], self.params[yi]))
    }

    /// Return the cached status of a constraint, or `None` if the handle is
    /// invalid.
    pub fn constraint_status(&self, handle: ConstraintHandle) -> Option<ConstraintStatus> {
        Some(self.constraints.get(handle)?.status)
    }

    /// Solve the constraint system using default solver options.
    ///
    /// Updates the parameter vector in place so that all driving constraints
    /// are satisfied (if a solution exists).
    pub fn solve(&mut self) -> Result<SolveResult, ConstraintError> {
        self.solve_with_options(&SolverOptions::default())
    }

    /// Solve the constraint system with custom solver options.
    ///
    /// Updates the parameter vector in place so that all driving constraints
    /// are satisfied (if a solution exists).
    pub fn solve_with_options(
        &mut self,
        options: &SolverOptions,
    ) -> Result<SolveResult, ConstraintError> {
        solver::solve(self, options)
    }

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
            Constraint::Angle(l1, l2, _) => {
                let line1 = self.lines.get(*l1)?;
                let line2 = self.lines.get(*l2)?;
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
            Constraint::PointLineDistance(p, l, _) => {
                let (px, py) = self.point_position(*p)?;
                let line = self.lines.get(*l)?;
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
            _ => {
                // For non-dimensional constraints, return residual magnitude.
                let output = equations::evaluate_constraint(self, entry, 0);
                let norm_sq: f64 = output.residuals.iter().map(|r| r * r).sum();
                Some(norm_sq.sqrt())
            }
        }
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
        let entry = self
            .constraints
            .get_mut(handle)
            .ok_or(ConstraintError::InvalidHandle)?;
        match &mut entry.constraint {
            Constraint::Distance(_, _, d) => {
                *d = value;
                Ok(())
            }
            Constraint::PointLineDistance(_, _, d) => {
                *d = value;
                Ok(())
            }
            Constraint::Angle(_, _, a) => {
                *a = value;
                Ok(())
            }
            _ => Err(ConstraintError::NotDimensional),
        }
    }

    /// Begin an interactive drag on a point.
    ///
    /// Saves the current parameter snapshot so subsequent [`drag_to`](Self::drag_to)
    /// calls can solve relative to the original state.
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
    pub fn drag_to(&mut self, x: f64, y: f64) -> Result<SolveResult, ConstraintError> {
        let drag = self
            .drag_context
            .as_ref()
            .ok_or(ConstraintError::NoDragActive)?;
        let point = drag.point;

        // Set the point to the target position (warm start).
        let p = self
            .points
            .get(point)
            .ok_or(ConstraintError::InvalidHandle)?;
        let [ix, iy] = p.param_indices();
        self.params[ix] = x;
        self.params[iy] = y;

        // Add a temporary Fixed constraint at (x, y).
        let temp = self.constraints.insert(ConstraintEntry {
            constraint: Constraint::Fixed(point),
            kind: ConstraintKind::Driving,
            status: ConstraintStatus::Unknown,
            fixed_values: vec![x, y],
        });
        if let Some(drag) = self.drag_context.as_mut() {
            drag.temp_constraint = Some(temp);
        }

        // Solve — for over-determined systems the QR least-squares path
        // finds the nearest valid configuration.
        let result = self.solve();

        // Always remove the temporary constraint, even on failure.
        self.constraints.remove(temp);
        if let Some(drag) = self.drag_context.as_mut() {
            drag.temp_constraint = None;
        }

        result
    }

    /// End the interactive drag, committing the current geometry.
    ///
    /// Does not restore the parameter snapshot — the drag result is kept.
    pub fn end_drag(&mut self) {
        self.drag_context = None;
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
        match element {
            ElementHandle::Point(ph) => self.remove_point(ph),
            ElementHandle::Line(lh) => self.remove_line(lh),
            ElementHandle::Arc(ah) => self.remove_arc(ah),
        }
    }

    /// Remove a point and cascade to dependent lines and constraints.
    fn remove_point(&mut self, ph: PointHandle) {
        // Find all lines that use this point as an endpoint.
        let dependent_lines: Vec<LineHandle> = self
            .lines
            .iter()
            .filter(|(_, l)| l.start == ph || l.end == ph)
            .map(|(h, _)| h)
            .collect();
        for lh in dependent_lines {
            self.remove_line(lh);
        }
        // Remove constraints that directly reference this point.
        let to_remove: Vec<ConstraintHandle> = self
            .constraints
            .iter()
            .filter(|(_, e)| e.constraint.references_point(ph))
            .map(|(h, _)| h)
            .collect();
        for ch in to_remove {
            self.constraints.remove(ch);
        }
        // Remove the point itself.
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
    fn remove_arc(&mut self, ah: ArcHandle) {
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

    /// Snapshot the current parameter values that a `Fixed` or `FixedLine`
    /// constraint needs to lock in place.
    ///
    /// Returns an empty vector for all other constraint kinds.
    pub fn capture_fixed_values(&self, constraint: &Constraint) -> Vec<f64> {
        match constraint {
            Constraint::Fixed(ph) => {
                if let Some(p) = self.points.get(*ph) {
                    let [xi, yi] = p.param_indices();
                    vec![self.params[xi], self.params[yi]]
                } else {
                    vec![]
                }
            }
            Constraint::FixedLine(lh) => {
                if let Some(l) = self.lines.get(*lh) {
                    let mut vals = Vec::with_capacity(4);
                    for ph in [l.start, l.end] {
                        if let Some(p) = self.points.get(ph) {
                            let [xi, yi] = p.param_indices();
                            vals.push(self.params[xi]);
                            vals.push(self.params[yi]);
                        }
                    }
                    vals
                } else {
                    vec![]
                }
            }
            _ => vec![],
        }
    }
}

impl Default for Sketch {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
#[path = "sketch_tests.rs"]
mod tests;
