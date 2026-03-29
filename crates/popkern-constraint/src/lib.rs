pub(crate) mod analysis;
pub mod constraint;
pub mod elements;
pub(crate) mod equations;
pub mod error;
pub(crate) mod jacobian;
pub mod pattern;
pub(crate) mod pattern_clone;
pub(crate) mod pattern_geo_helpers;
pub(crate) mod pattern_remap;
mod sketch_ops;
mod sketch_pattern;
mod sketch_pattern_geo;
mod sketch_spline;
mod solver;
pub(crate) mod spline_interp;
pub(crate) mod undo;

#[cfg(test)]
#[path = "solver_tests.rs"]
mod solver_tests;

pub use constraint::{
    Constraint, ConstraintEntry, ConstraintHandle, ConstraintKind, ConstraintStatus,
};
pub use elements::{
    Arc, ArcHandle, ElementHandle, Ellipse, EllipseHandle, Hyperbola, HyperbolaHandle, Line,
    LineHandle, Parabola, ParabolaHandle, Point, PointHandle, Spline, SplineHandle,
};
pub use error::{ConstraintError, FreedomInfo, SolveResult, SolveStatus};
pub use pattern::{PatternDef, PatternHandle, PatternInstance, PatternKind};
pub use solver::SolverOptions;

use serde::{Deserialize, Serialize};
use slotmap::SlotMap;
use undo::Snapshot;

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
#[derive(Serialize, Deserialize)]
pub struct Sketch {
    /// Flat parameter vector shared by all elements.
    pub(crate) params: Vec<f64>,
    /// Point elements, keyed by [`PointHandle`].
    pub(crate) points: SlotMap<PointHandle, Point>,
    /// Line elements, keyed by [`LineHandle`].
    pub(crate) lines: SlotMap<LineHandle, Line>,
    /// Arc elements, keyed by [`ArcHandle`].
    pub(crate) arcs: SlotMap<ArcHandle, Arc>,
    /// Ellipse elements, keyed by [`EllipseHandle`].
    pub(crate) ellipses: SlotMap<EllipseHandle, Ellipse>,
    /// Parabola elements, keyed by [`ParabolaHandle`].
    pub(crate) parabolas: SlotMap<ParabolaHandle, Parabola>,
    /// Hyperbola elements, keyed by [`HyperbolaHandle`].
    pub(crate) hyperbolas: SlotMap<HyperbolaHandle, Hyperbola>,
    /// Spline elements, keyed by [`SplineHandle`].
    pub(crate) splines: SlotMap<SplineHandle, Spline>,
    /// Constraints, keyed by [`ConstraintHandle`].
    pub(crate) constraints: SlotMap<ConstraintHandle, ConstraintEntry>,
    /// Pattern definitions, keyed by [`PatternHandle`].
    pub(crate) patterns: SlotMap<PatternHandle, PatternDef>,
    /// Active drag context (transient, not serialised).
    #[serde(skip)]
    pub(crate) drag_context: Option<DragState>,
    /// Undo snapshot log (session-local, not serialised).
    #[serde(skip)]
    pub(crate) undo_log: Vec<Snapshot>,
    /// Cursor into [`undo_log`](Self::undo_log) — points to the next
    /// insert position. Snapshots before the cursor can be undone;
    /// snapshots at or after the cursor can be redone.
    #[serde(skip)]
    pub(crate) undo_cursor: usize,
    /// Maps constraint handles to the start index of their auxiliary
    /// parameters in the parameter vector. Transient — rebuilt as needed.
    #[serde(skip)]
    pub(crate) aux_param_map: std::collections::HashMap<ConstraintHandle, usize>,
    /// When `false`, mutations do not record snapshots (used by undo/redo
    /// to avoid infinite recursion).
    #[serde(skip, default = "default_recording")]
    pub(crate) recording: bool,
}

/// Default for the `recording` field when deserialising (always `true`).
fn default_recording() -> bool {
    true
}

impl Sketch {
    /// Create a new, empty sketch.
    pub fn new() -> Self {
        Self {
            params: Vec::new(),
            points: SlotMap::with_key(),
            lines: SlotMap::with_key(),
            arcs: SlotMap::with_key(),
            ellipses: SlotMap::with_key(),
            parabolas: SlotMap::with_key(),
            hyperbolas: SlotMap::with_key(),
            splines: SlotMap::with_key(),
            constraints: SlotMap::with_key(),
            patterns: SlotMap::with_key(),
            drag_context: None,
            aux_param_map: std::collections::HashMap::new(),
            undo_log: Vec::new(),
            undo_cursor: 0,
            recording: true,
        }
    }

    /// Capture the current sketch state and push it onto the undo log,
    /// truncating any redo history beyond the cursor.
    pub(crate) fn push_undo_snapshot(&mut self) {
        self.undo_log.truncate(self.undo_cursor);
        self.undo_log.push(Snapshot {
            params: self.params.clone(),
            points: self.points.clone(),
            lines: self.lines.clone(),
            arcs: self.arcs.clone(),
            ellipses: self.ellipses.clone(),
            parabolas: self.parabolas.clone(),
            hyperbolas: self.hyperbolas.clone(),
            splines: self.splines.clone(),
            constraints: self.constraints.clone(),
            patterns: self.patterns.clone(),
        });
        self.undo_cursor += 1;
    }

    /// Add a free point at `(x, y)`.
    ///
    /// The point owns two consecutive entries in the parameter vector.
    pub fn add_point(&mut self, x: f64, y: f64) -> PointHandle {
        if self.recording {
            self.push_undo_snapshot();
        }
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
        if self.recording {
            self.push_undo_snapshot();
        }
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
        if self.recording {
            self.push_undo_snapshot();
        }
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
        if self.recording {
            self.push_undo_snapshot();
        }
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
            ElementHandle::Ellipse(h) => {
                if let Some(e) = self.ellipses.get_mut(h) {
                    e.construction = value;
                }
            }
            ElementHandle::Parabola(h) => {
                if let Some(p) = self.parabolas.get_mut(h) {
                    p.construction = value;
                }
            }
            ElementHandle::Hyperbola(h) => {
                if let Some(hp) = self.hyperbolas.get_mut(h) {
                    hp.construction = value;
                }
            }
            ElementHandle::Spline(h) => {
                if let Some(sp) = self.splines.get_mut(h) {
                    sp.construction = value;
                }
            }
        }
    }

    /// Add a driving constraint to the sketch.
    ///
    /// Driving constraints determine the geometry (e.g. a distance fixes a
    /// length).
    pub fn add_constraint(&mut self, constraint: Constraint) -> ConstraintHandle {
        if self.recording {
            self.push_undo_snapshot();
        }
        let fixed_values = self.capture_fixed_values(&constraint);
        let aux_count = constraint.aux_param_count();
        let handle = self.constraints.insert(ConstraintEntry {
            constraint,
            kind: ConstraintKind::Driving,
            status: ConstraintStatus::Unknown,
            fixed_values,
        });
        if aux_count > 0 {
            let aux_start = self.params.len();
            self.params.extend(std::iter::repeat_n(0.0, aux_count));
            self.aux_param_map.insert(handle, aux_start);
        }
        handle
    }

    /// Add a driven constraint to the sketch.
    ///
    /// Driven constraints are display-only annotations; they do not influence
    /// the solver.
    pub fn add_driven(&mut self, constraint: Constraint) -> ConstraintHandle {
        if self.recording {
            self.push_undo_snapshot();
        }
        let fixed_values = self.capture_fixed_values(&constraint);
        let aux_count = constraint.aux_param_count();
        let handle = self.constraints.insert(ConstraintEntry {
            constraint,
            kind: ConstraintKind::Driven,
            status: ConstraintStatus::Unknown,
            fixed_values,
        });
        if aux_count > 0 {
            let aux_start = self.params.len();
            self.params.extend(std::iter::repeat_n(0.0, aux_count));
            self.aux_param_map.insert(handle, aux_start);
        }
        handle
    }

    /// Remove a constraint from the sketch.
    ///
    /// Does nothing if the handle is no longer valid.
    pub fn remove_constraint(&mut self, handle: ConstraintHandle) {
        if self.recording {
            self.push_undo_snapshot();
        }
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

    /// Return the current arc parameters `(cx, cy, r, start_angle,
    /// end_angle)`, or `None` if the handle is invalid.
    pub fn arc_params(&self, handle: ArcHandle) -> Option<(f64, f64, f64, f64, f64)> {
        let a = self.arcs.get(handle)?;
        let [icx, icy, ir, isa, iea] = a.param_indices();
        Some((
            self.params[icx],
            self.params[icy],
            self.params[ir],
            self.params[isa],
            self.params[iea],
        ))
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

    /// Return `true` if there is at least one operation that can be undone.
    pub fn can_undo(&self) -> bool {
        self.undo_cursor > 0
    }

    /// Return `true` if there is at least one operation that can be redone.
    pub fn can_redo(&self) -> bool {
        self.undo_cursor < self.undo_log.len()
    }

    /// Undo the most recent operation, returning `false` if nothing to
    /// undo.
    ///
    /// Restores the sketch state to the snapshot captured before the last
    /// mutation. The current (post-mutation) state is preserved in the
    /// redo log so it can be re-applied.
    pub fn undo(&mut self) -> bool {
        if !self.can_undo() {
            return false;
        }
        // If we are at the tip (no redo entries beyond cursor), save the
        // current state so redo can get back here.
        if self.undo_cursor == self.undo_log.len() {
            self.undo_log.push(Snapshot {
                params: self.params.clone(),
                points: self.points.clone(),
                lines: self.lines.clone(),
                arcs: self.arcs.clone(),
                ellipses: self.ellipses.clone(),
                parabolas: self.parabolas.clone(),
                hyperbolas: self.hyperbolas.clone(),
                splines: self.splines.clone(),
                constraints: self.constraints.clone(),
                patterns: self.patterns.clone(),
            });
        }
        self.undo_cursor -= 1;
        let snapshot = self.undo_log[self.undo_cursor].clone();
        self.restore_snapshot(snapshot);
        true
    }

    /// Redo the most recently undone operation, returning `false` if
    /// nothing to redo.
    pub fn redo(&mut self) -> bool {
        if self.undo_cursor + 1 >= self.undo_log.len() {
            return false;
        }
        self.undo_cursor += 1;
        let snapshot = self.undo_log[self.undo_cursor].clone();
        self.restore_snapshot(snapshot);
        true
    }

    /// Replace the sketch's mutable state from a snapshot.
    fn restore_snapshot(&mut self, snapshot: Snapshot) {
        self.params = snapshot.params;
        self.points = snapshot.points;
        self.lines = snapshot.lines;
        self.arcs = snapshot.arcs;
        self.ellipses = snapshot.ellipses;
        self.parabolas = snapshot.parabolas;
        self.hyperbolas = snapshot.hyperbolas;
        self.splines = snapshot.splines;
        self.constraints = snapshot.constraints;
        self.patterns = snapshot.patterns;
    }

    /// Snapshot the current parameter values that a fixed constraint needs
    /// to lock in place.
    ///
    /// Returns an empty vector for all non-fixed constraint kinds.
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
            Constraint::FixedEllipse(eh) => {
                if let Some(e) = self.ellipses.get(*eh) {
                    let idx = e.param_indices();
                    idx.iter().map(|&i| self.params[i]).collect()
                } else {
                    vec![]
                }
            }
            Constraint::FixedParabola(ph) => {
                if let Some(p) = self.parabolas.get(*ph) {
                    let idx = p.param_indices();
                    idx.iter().map(|&i| self.params[i]).collect()
                } else {
                    vec![]
                }
            }
            Constraint::FixedHyperbola(hh) => {
                if let Some(h) = self.hyperbolas.get(*hh) {
                    let idx = h.param_indices();
                    idx.iter().map(|&i| self.params[i]).collect()
                } else {
                    vec![]
                }
            }
            Constraint::SplineFixedEndTangent(_s, _end_idx, tx, ty) => {
                vec![*tx, *ty]
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
