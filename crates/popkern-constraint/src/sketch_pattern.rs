//! Pattern creation and deletion operations for [`Sketch`].
//!
//! Provides `linear_pattern`, `delete_pattern`, and pattern validation
//! as `impl Sketch` methods.

use crate::Sketch;
use crate::constraint::Constraint;
use crate::elements::ElementHandle;
use crate::error::ConstraintError;
use crate::pattern::{PatternDef, PatternHandle, PatternInstance, PatternKind};
use crate::pattern_clone::{
    clone_constraints, clone_elements, complete_seed_set, find_seed_constraints,
};

impl Sketch {
    /// Create a linear pattern replicating `seed` elements `count` times
    /// along `direction` with the given `spacing`.
    ///
    /// Instance 0 is the seed itself; instances 1..count are cloned copies
    /// offset by `k * spacing * direction` for k = 1..count-1.
    pub fn linear_pattern(
        &mut self,
        seed: &[ElementHandle],
        count: usize,
        direction: (f64, f64),
        spacing: f64,
    ) -> Result<PatternHandle, ConstraintError> {
        if count < 2 {
            return Err(ConstraintError::InvalidPattern("count must be >= 2".into()));
        }
        self.validate_seed_not_in_pattern(seed)?;
        if self.recording {
            self.push_undo_snapshot();
        }

        let seed_set = complete_seed_set(self, seed);
        let seed_constraints = find_seed_constraints(self, &seed_set);

        // Normalize direction
        let len = (direction.0 * direction.0 + direction.1 * direction.1).sqrt();
        let (dx, dy) = if len > 1e-15 {
            (direction.0 / len, direction.1 / len)
        } else {
            (1.0, 0.0)
        };

        let mut instances = Vec::with_capacity(count - 1);
        let was_recording = self.recording;
        self.recording = false;

        for k in 1..count {
            let offset_x = k as f64 * spacing * dx;
            let offset_y = k as f64 * spacing * dy;

            let (cloned_elems, remap, phantom_points) =
                clone_elements(self, seed, offset_x, offset_y);
            let internal = clone_constraints(self, &seed_constraints, &remap);

            // Add RelativePosition constraints for each point pair
            let mut offset_constraints = Vec::new();
            for (&orig_ph, &new_ph) in &remap.points {
                let ch = self.add_constraint(Constraint::RelativePosition(
                    orig_ph, new_ph, offset_x, offset_y,
                ));
                offset_constraints.push(ch);
            }

            // Add PinToCenter for each arc/conic phantom
            let pin_constraints = add_pin_constraints(self, seed, &remap, &phantom_points);

            instances.push(PatternInstance {
                elements: cloned_elems,
                phantom_points,
                internal_constraints: internal,
                offset_constraints,
                pin_constraints,
                angle_constraints: Vec::new(),
                suppressed: false,
            });
        }

        self.recording = was_recording;

        let handle = self.patterns.insert(PatternDef {
            kind: PatternKind::Linear {
                direction: (dx, dy),
                spacing,
            },
            seed_elements: seed.to_vec(),
            seed_internal_constraints: seed_constraints,
            instances,
        });
        Ok(handle)
    }

    /// Delete a pattern, removing all cloned elements and constraints.
    ///
    /// The seed elements and their constraints are left intact.
    pub fn delete_pattern(&mut self, handle: PatternHandle) -> Result<(), ConstraintError> {
        let pattern = self
            .patterns
            .remove(handle)
            .ok_or(ConstraintError::InvalidHandle)?;

        if self.recording {
            self.push_undo_snapshot();
        }

        let was_recording = self.recording;
        self.recording = false;

        for instance in &pattern.instances {
            // Remove all constraints first (offset, pin, angle, internal)
            for &ch in instance
                .offset_constraints
                .iter()
                .chain(&instance.pin_constraints)
                .chain(&instance.angle_constraints)
                .chain(&instance.internal_constraints)
            {
                self.constraints.remove(ch);
            }
            // Remove phantom points
            for &ph in &instance.phantom_points {
                self.points.remove(ph);
            }
            // Remove cloned elements
            for elem in &instance.elements {
                match elem {
                    ElementHandle::Point(ph) => {
                        self.points.remove(*ph);
                    }
                    ElementHandle::Line(lh) => {
                        self.lines.remove(*lh);
                    }
                    ElementHandle::Arc(ah) => {
                        self.arcs.remove(*ah);
                    }
                    ElementHandle::Ellipse(eh) => {
                        self.ellipses.remove(*eh);
                    }
                    ElementHandle::Parabola(ph) => {
                        self.parabolas.remove(*ph);
                    }
                    ElementHandle::Hyperbola(hh) => {
                        self.hyperbolas.remove(*hh);
                    }
                    ElementHandle::Spline(sh) => {
                        self.splines.remove(*sh);
                    }
                }
            }
        }

        self.recording = was_recording;
        Ok(())
    }

    /// Validate that no seed element is already part of a pattern instance.
    pub(crate) fn validate_seed_not_in_pattern(
        &self,
        seed: &[ElementHandle],
    ) -> Result<(), ConstraintError> {
        for (_ph, pattern) in self.patterns.iter() {
            for instance in &pattern.instances {
                for elem in &instance.elements {
                    if seed.contains(elem) {
                        return Err(ConstraintError::InvalidPattern(
                            "seed element is already a pattern instance".into(),
                        ));
                    }
                }
            }
        }
        Ok(())
    }

    /// Find the pattern that owns a given element (seed or instance).
    ///
    /// Returns `None` if the element is not part of any pattern.
    pub(crate) fn find_pattern_for_element(&self, elem: &ElementHandle) -> Option<PatternHandle> {
        for (ph, pattern) in self.patterns.iter() {
            if pattern.seed_elements.contains(elem) {
                return Some(ph);
            }
            for instance in &pattern.instances {
                if instance.elements.contains(elem) {
                    return Some(ph);
                }
            }
        }
        None
    }
}

/// Add PinToCenter constraints for phantom points paired with arc/conic
/// cloned elements. Returns the list of new constraint handles.
pub(crate) fn add_pin_constraints(
    sketch: &mut Sketch,
    seed: &[ElementHandle],
    remap: &crate::pattern::HandleRemap,
    phantom_points: &[crate::PointHandle],
) -> Vec<crate::ConstraintHandle> {
    let mut pin_constraints = Vec::new();
    let mut phantom_idx = 0;
    for elem in seed {
        match elem {
            ElementHandle::Arc(ah) => {
                if let Some(&new_ah) = remap.arcs.get(ah)
                    && phantom_idx < phantom_points.len()
                {
                    let phantom = phantom_points[phantom_idx];
                    let ch = sketch.add_constraint(Constraint::PinToCenter(
                        phantom,
                        ElementHandle::Arc(new_ah),
                    ));
                    pin_constraints.push(ch);
                    phantom_idx += 1;
                }
            }
            ElementHandle::Ellipse(eh) => {
                if let Some(&new_eh) = remap.ellipses.get(eh)
                    && phantom_idx < phantom_points.len()
                {
                    let phantom = phantom_points[phantom_idx];
                    let ch = sketch.add_constraint(Constraint::PinToCenter(
                        phantom,
                        ElementHandle::Ellipse(new_eh),
                    ));
                    pin_constraints.push(ch);
                    phantom_idx += 1;
                }
            }
            ElementHandle::Parabola(ph_key) => {
                if let Some(&new_ph) = remap.parabolas.get(ph_key)
                    && phantom_idx < phantom_points.len()
                {
                    let phantom = phantom_points[phantom_idx];
                    let ch = sketch.add_constraint(Constraint::PinToCenter(
                        phantom,
                        ElementHandle::Parabola(new_ph),
                    ));
                    pin_constraints.push(ch);
                    phantom_idx += 1;
                }
            }
            ElementHandle::Hyperbola(hh) => {
                if let Some(&new_hh) = remap.hyperbolas.get(hh)
                    && phantom_idx < phantom_points.len()
                {
                    let phantom = phantom_points[phantom_idx];
                    let ch = sketch.add_constraint(Constraint::PinToCenter(
                        phantom,
                        ElementHandle::Hyperbola(new_hh),
                    ));
                    pin_constraints.push(ch);
                    phantom_idx += 1;
                }
            }
            _ => {}
        }
    }
    pin_constraints
}
