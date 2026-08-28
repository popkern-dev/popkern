//! Pattern suppression, update, and refresh operations for [`Sketch`].
//!
//! Provides instance suppression/unsuppression, spacing/count/angle
//! updates, and constraint refresh as `impl Sketch` methods.

use crate::Sketch;
use crate::constraint::{Constraint, ConstraintKind};
use crate::elements::ElementHandle;
use crate::error::ConstraintError;
use crate::pattern::{PatternHandle, PatternKind};
use crate::pattern_clone::{
    clone_constraints, clone_elements, complete_seed_set, find_seed_constraints,
};
use crate::sketch_pattern::add_pin_constraints;

impl Sketch {
    /// Suppress a pattern instance (set all its constraints to Driven).
    ///
    /// Suppressed instances still exist but do not influence the solver.
    pub fn suppress_instance(
        &mut self,
        pattern: PatternHandle,
        index: usize,
    ) -> Result<(), ConstraintError> {
        let pat = self
            .patterns
            .get_mut(pattern)
            .ok_or(ConstraintError::InvalidHandle)?;
        let instance = pat
            .instances
            .get_mut(index)
            .ok_or(ConstraintError::InvalidPattern(
                "instance index out of range".into(),
            ))?;
        if instance.suppressed {
            return Ok(());
        }
        instance.suppressed = true;

        let all_constraints: Vec<_> = instance
            .offset_constraints
            .iter()
            .chain(&instance.pin_constraints)
            .chain(&instance.angle_constraints)
            .chain(&instance.internal_constraints)
            .copied()
            .collect();

        for ch in all_constraints {
            if let Some(entry) = self.constraints.get_mut(ch) {
                entry.kind = ConstraintKind::Driven;
            }
        }
        Ok(())
    }

    /// Unsuppress a pattern instance (set all its constraints back to Driving).
    pub fn unsuppress_instance(
        &mut self,
        pattern: PatternHandle,
        index: usize,
    ) -> Result<(), ConstraintError> {
        let pat = self
            .patterns
            .get_mut(pattern)
            .ok_or(ConstraintError::InvalidHandle)?;
        let instance = pat
            .instances
            .get_mut(index)
            .ok_or(ConstraintError::InvalidPattern(
                "instance index out of range".into(),
            ))?;
        if !instance.suppressed {
            return Ok(());
        }
        instance.suppressed = false;

        let all_constraints: Vec<_> = instance
            .offset_constraints
            .iter()
            .chain(&instance.pin_constraints)
            .chain(&instance.angle_constraints)
            .chain(&instance.internal_constraints)
            .copied()
            .collect();

        for ch in all_constraints {
            if let Some(entry) = self.constraints.get_mut(ch) {
                entry.kind = ConstraintKind::Driving;
            }
        }
        Ok(())
    }

    /// Update the spacing of a linear pattern.
    ///
    /// Recomputes RelativePosition constraint offsets for all instances.
    pub fn update_pattern_spacing(
        &mut self,
        pattern: PatternHandle,
        spacing: f64,
    ) -> Result<(), ConstraintError> {
        let pat = self
            .patterns
            .get_mut(pattern)
            .ok_or(ConstraintError::InvalidHandle)?;
        let (dx, dy) = match &mut pat.kind {
            PatternKind::Linear {
                direction,
                spacing: s,
            } => {
                *s = spacing;
                *direction
            }
            _ => {
                return Err(ConstraintError::InvalidPattern(
                    "update_pattern_spacing only applies to linear patterns".into(),
                ));
            }
        };

        for (k_minus_1, instance) in pat.instances.iter().enumerate() {
            let k = (k_minus_1 + 1) as f64;
            let offset_x = k * spacing * dx;
            let offset_y = k * spacing * dy;
            for &ch in &instance.offset_constraints {
                if let Some(entry) = self.constraints.get_mut(ch)
                    && let Constraint::RelativePosition(_, _, ref mut cdx, ref mut cdy) =
                        entry.constraint
                {
                    *cdx = offset_x;
                    *cdy = offset_y;
                }
            }
        }
        Ok(())
    }

    /// Update the angle of a circular pattern.
    ///
    /// Recomputes RotatedPosition and AngleOffset constraint values.
    pub fn update_pattern_angle(
        &mut self,
        pattern: PatternHandle,
        angle: f64,
    ) -> Result<(), ConstraintError> {
        let pat = self
            .patterns
            .get_mut(pattern)
            .ok_or(ConstraintError::InvalidHandle)?;
        let count = pat.instances.len() + 1;
        match &mut pat.kind {
            PatternKind::Circular { angle: a, .. } => *a = angle,
            _ => {
                return Err(ConstraintError::InvalidPattern(
                    "update_pattern_angle only applies to circular patterns".into(),
                ));
            }
        };

        for (k_minus_1, instance) in pat.instances.iter().enumerate() {
            let angle_k = (k_minus_1 + 1) as f64 * angle / (count - 1) as f64;
            for &ch in &instance.offset_constraints {
                if let Some(entry) = self.constraints.get_mut(ch)
                    && let Constraint::RotatedPosition(_, _, _, ref mut t) = entry.constraint
                {
                    *t = angle_k;
                }
            }
            for &ch in &instance.angle_constraints {
                if let Some(entry) = self.constraints.get_mut(ch) {
                    match &mut entry.constraint {
                        Constraint::AngleOffsetArc(_, _, t) => *t = angle_k,
                        Constraint::AngleOffsetConic(_, _, t) => *t = angle_k,
                        _ => {}
                    }
                }
            }
        }
        Ok(())
    }

    /// Update the instance count of a linear pattern.
    ///
    /// Increasing adds new instances; decreasing removes from the end.
    pub fn update_pattern_count(
        &mut self,
        pattern: PatternHandle,
        count: usize,
    ) -> Result<(), ConstraintError> {
        if count < 2 {
            return Err(ConstraintError::InvalidPattern("count must be >= 2".into()));
        }
        let pat = self
            .patterns
            .get(pattern)
            .ok_or(ConstraintError::InvalidHandle)?;
        let current_count = pat.instances.len() + 1;

        if count == current_count {
            return Ok(());
        }

        if count < current_count {
            // Remove instances from the end
            let to_remove = current_count - count;
            for _ in 0..to_remove {
                let pat = self
                    .patterns
                    .get_mut(pattern)
                    .ok_or(ConstraintError::InvalidHandle)?;
                if let Some(instance) = pat.instances.pop() {
                    remove_instance_entities(self, &instance);
                }
            }
        } else {
            // Add new instances
            let pat = self
                .patterns
                .get(pattern)
                .ok_or(ConstraintError::InvalidHandle)?;
            let (dx, dy, spacing) = match &pat.kind {
                PatternKind::Linear { direction, spacing } => (direction.0, direction.1, *spacing),
                _ => {
                    return Err(ConstraintError::InvalidPattern(
                        "update_pattern_count currently only supports linear patterns".into(),
                    ));
                }
            };
            let seed = pat.seed_elements.clone();
            let seed_constraints = pat.seed_internal_constraints.clone();

            let was_recording = self.recording;
            self.recording = false;

            for k in current_count..count {
                let offset_x = k as f64 * spacing * dx;
                let offset_y = k as f64 * spacing * dy;

                let (cloned_elems, remap, phantom_points) =
                    clone_elements(self, &seed, offset_x, offset_y);
                let internal = clone_constraints(self, &seed_constraints, &remap);

                let mut offset_constraints = Vec::new();
                for (&orig_ph, &new_ph) in &remap.points {
                    let ch = self.add_constraint(Constraint::RelativePosition(
                        orig_ph, new_ph, offset_x, offset_y,
                    ));
                    offset_constraints.push(ch);
                }

                let pin_constraints = add_pin_constraints(self, &seed, &remap, &phantom_points);

                let instance = crate::pattern::PatternInstance {
                    elements: cloned_elems,
                    phantom_points,
                    internal_constraints: internal,
                    offset_constraints,
                    pin_constraints,
                    angle_constraints: Vec::new(),
                    suppressed: false,
                };
                let pat = self
                    .patterns
                    .get_mut(pattern)
                    .ok_or(ConstraintError::InvalidHandle)?;
                pat.instances.push(instance);
            }

            self.recording = was_recording;
        }
        Ok(())
    }

    /// Refresh a pattern by re-identifying seed constraints and cloning
    /// any new ones to all non-suppressed instances.
    pub fn refresh_pattern(&mut self, pattern: PatternHandle) -> Result<(), ConstraintError> {
        let pat = self
            .patterns
            .get(pattern)
            .ok_or(ConstraintError::InvalidHandle)?;
        let seed = pat.seed_elements.clone();
        let old_seed_constraints: std::collections::HashSet<_> =
            pat.seed_internal_constraints.iter().copied().collect();

        let seed_set = complete_seed_set(self, &seed);
        let current = find_seed_constraints(self, &seed_set);

        // Find new constraints not in the original seed set
        let new_constraints: Vec<_> = current
            .iter()
            .filter(|ch| !old_seed_constraints.contains(ch))
            .copied()
            .collect();

        if new_constraints.is_empty() {
            return Ok(());
        }

        // Update seed_internal_constraints
        let pat = self
            .patterns
            .get_mut(pattern)
            .ok_or(ConstraintError::InvalidHandle)?;
        pat.seed_internal_constraints = current;

        // Collect instance remap data
        let instance_data: Vec<(bool, Vec<ElementHandle>)> = pat
            .instances
            .iter()
            .map(|inst| (inst.suppressed, inst.elements.clone()))
            .collect();

        let was_recording = self.recording;
        self.recording = false;

        for (i, (suppressed, _elements)) in instance_data.iter().enumerate() {
            if *suppressed {
                continue;
            }

            // Build remap from seed elements to instance elements
            let pat = self
                .patterns
                .get(pattern)
                .ok_or(ConstraintError::InvalidHandle)?;
            let remap = build_remap_from_instance(self, &pat.seed_elements, &pat.instances[i]);

            let cloned = clone_constraints(self, &new_constraints, &remap);

            // Add to instance's internal_constraints
            let pat = self
                .patterns
                .get_mut(pattern)
                .ok_or(ConstraintError::InvalidHandle)?;
            pat.instances[i].internal_constraints.extend(cloned);
        }

        self.recording = was_recording;
        Ok(())
    }
}

/// Remove all entities belonging to a pattern instance.
fn remove_instance_entities(sketch: &mut Sketch, instance: &crate::pattern::PatternInstance) {
    for &ch in instance
        .offset_constraints
        .iter()
        .chain(&instance.pin_constraints)
        .chain(&instance.angle_constraints)
        .chain(&instance.internal_constraints)
    {
        sketch.constraints.remove(ch);
    }
    for &ph in &instance.phantom_points {
        sketch.points.remove(ph);
    }
    for elem in &instance.elements {
        match elem {
            ElementHandle::Point(ph) => {
                sketch.points.remove(*ph);
            }
            ElementHandle::Line(lh) => {
                sketch.lines.remove(*lh);
            }
            ElementHandle::Arc(ah) => {
                sketch.arcs.remove(*ah);
            }
            ElementHandle::Ellipse(eh) => {
                sketch.ellipses.remove(*eh);
            }
            ElementHandle::Parabola(ph) => {
                sketch.parabolas.remove(*ph);
            }
            ElementHandle::Hyperbola(hh) => {
                sketch.hyperbolas.remove(*hh);
            }
            ElementHandle::Spline(sh) => {
                sketch.splines.remove(*sh);
            }
        }
    }
}

/// Build a HandleRemap from seed elements to an instance's elements.
fn build_remap_from_instance(
    _sketch: &Sketch,
    seed: &[ElementHandle],
    instance: &crate::pattern::PatternInstance,
) -> crate::pattern::HandleRemap {
    let mut remap = crate::pattern::HandleRemap::default();
    for (seed_elem, inst_elem) in seed.iter().zip(instance.elements.iter()) {
        match (seed_elem, inst_elem) {
            (ElementHandle::Point(s), ElementHandle::Point(i)) => {
                remap.points.insert(*s, *i);
            }
            (ElementHandle::Line(s), ElementHandle::Line(i)) => {
                remap.lines.insert(*s, *i);
            }
            (ElementHandle::Arc(s), ElementHandle::Arc(i)) => {
                remap.arcs.insert(*s, *i);
            }
            (ElementHandle::Ellipse(s), ElementHandle::Ellipse(i)) => {
                remap.ellipses.insert(*s, *i);
            }
            (ElementHandle::Parabola(s), ElementHandle::Parabola(i)) => {
                remap.parabolas.insert(*s, *i);
            }
            (ElementHandle::Hyperbola(s), ElementHandle::Hyperbola(i)) => {
                remap.hyperbolas.insert(*s, *i);
            }
            (ElementHandle::Spline(s), ElementHandle::Spline(i)) => {
                remap.splines.insert(*s, *i);
            }
            _ => {}
        }
    }
    remap
}
