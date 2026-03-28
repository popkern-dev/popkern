pub mod constraint;
pub mod elements;
pub mod error;

pub use constraint::{
    Constraint, ConstraintEntry, ConstraintHandle, ConstraintKind, ConstraintStatus,
};
pub use elements::{Arc, ArcHandle, ElementHandle, Line, LineHandle, Point, PointHandle};
pub use error::{ConstraintError, FreedomInfo, SolveResult, SolveStatus};

use slotmap::SlotMap;

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
mod tests {
    use super::*;

    #[test]
    fn add_point_and_read_position() {
        let mut sketch = Sketch::new();
        let h = sketch.add_point(3.0, 4.0);
        assert_eq!(sketch.point_position(h), Some((3.0, 4.0)));
    }

    #[test]
    fn add_line_shares_point_params() {
        let mut sketch = Sketch::new();
        let a = sketch.add_point(0.0, 0.0);
        let b = sketch.add_point(1.0, 0.0);
        let _l = sketch.add_line(a, b);

        // Mutate the underlying param vector directly to verify sharing.
        let idx = sketch.points[a].param_indices()[0];
        sketch.params[idx] = 5.0;

        assert_eq!(sketch.point_position(a), Some((5.0, 0.0)));
    }

    #[test]
    fn add_and_remove_constraint() {
        let mut sketch = Sketch::new();
        let a = sketch.add_point(0.0, 0.0);
        let b = sketch.add_point(1.0, 0.0);
        let ch = sketch.add_constraint(Constraint::Coincident(a, b));
        assert!(sketch.constraint_status(ch).is_some());

        sketch.remove_constraint(ch);
        assert!(sketch.constraint_status(ch).is_none());
    }
}
