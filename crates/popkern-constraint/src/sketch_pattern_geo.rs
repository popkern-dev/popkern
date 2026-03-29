//! Circular and mirror pattern operations for [`Sketch`].
//!
//! Split from `sketch_pattern.rs` to stay within the 500-line file limit.

use crate::constraint::Constraint;
use crate::elements::{ElementHandle, PointHandle};
use crate::error::ConstraintError;
use crate::pattern::{PatternDef, PatternHandle, PatternInstance, PatternKind};
use crate::pattern_clone::{
    clone_constraints, clone_elements, complete_seed_set, find_seed_constraints,
};
use crate::pattern_geo_helpers::{
    reflect_cloned_centres, reflect_point, rotate_cloned_centres, rotate_point, set_point_params,
    sync_phantom_positions,
};
use crate::sketch_pattern::add_pin_constraints;
use crate::Sketch;

impl Sketch {
    /// Create a circular pattern replicating `seed` elements `count` times
    /// around `center` spanning `angle` radians.
    ///
    /// Instance k (k = 1..count-1) is rotated by `k * angle / (count - 1)`.
    pub fn circular_pattern(
        &mut self,
        seed: &[ElementHandle],
        count: usize,
        center: PointHandle,
        angle: f64,
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
        let (cx, cy) = self
            .point_position(center)
            .ok_or(ConstraintError::InvalidHandle)?;

        let mut instances = Vec::with_capacity(count - 1);
        let was_recording = self.recording;
        self.recording = false;

        for k in 1..count {
            let angle_k = k as f64 * angle / (count - 1) as f64;
            let cos_a = angle_k.cos();
            let sin_a = angle_k.sin();

            let (cloned_elems, remap, phantom_points) = clone_elements(self, seed, 0.0, 0.0);

            // Set rotated initial positions for cloned points
            for (&orig_ph, &new_ph) in &remap.points {
                if let Some((px, py)) = self.point_position(orig_ph) {
                    let (rx, ry) = rotate_point(px, py, cx, cy, cos_a, sin_a);
                    set_point_params(self, new_ph, rx, ry);
                }
            }
            rotate_cloned_centres(self, seed, &remap, cx, cy, cos_a, sin_a);
            sync_phantom_positions(self, seed, &remap, &phantom_points);

            let internal = clone_constraints(self, &seed_constraints, &remap);

            let mut offset_constraints = Vec::new();
            for (&orig_ph, &new_ph) in &remap.points {
                let ch = self.add_constraint(Constraint::RotatedPosition(
                    orig_ph, new_ph, center, angle_k,
                ));
                offset_constraints.push(ch);
            }

            let pin_constraints = add_pin_constraints(self, seed, &remap, &phantom_points);
            let angle_constraints = add_angle_offset_constraints(self, seed, &remap, angle_k);

            instances.push(PatternInstance {
                elements: cloned_elems,
                phantom_points,
                internal_constraints: internal,
                offset_constraints,
                pin_constraints,
                angle_constraints,
                suppressed: false,
            });
        }
        self.recording = was_recording;

        let handle = self.patterns.insert(PatternDef {
            kind: PatternKind::Circular { center, angle },
            seed_elements: seed.to_vec(),
            seed_internal_constraints: seed_constraints,
            instances,
        });
        Ok(handle)
    }

    /// Create a mirror pattern reflecting `seed` elements across `line`.
    ///
    /// Always creates exactly one instance (the reflection).
    pub fn mirror_pattern(
        &mut self,
        seed: &[ElementHandle],
        line: crate::LineHandle,
    ) -> Result<PatternHandle, ConstraintError> {
        self.validate_seed_not_in_pattern(seed)?;
        if self.recording {
            self.push_undo_snapshot();
        }

        let seed_set = complete_seed_set(self, seed);
        let seed_constraints = find_seed_constraints(self, &seed_set);

        let ml = self.lines.get(line).ok_or(ConstraintError::InvalidHandle)?;
        let (lx1, ly1) = self
            .point_position(ml.start)
            .ok_or(ConstraintError::InvalidHandle)?;
        let (lx2, ly2) = self
            .point_position(ml.end)
            .ok_or(ConstraintError::InvalidHandle)?;

        let was_recording = self.recording;
        self.recording = false;

        let (cloned_elems, remap, phantom_points) = clone_elements(self, seed, 0.0, 0.0);

        // Set reflected positions for cloned points
        for (&orig_ph, &new_ph) in &remap.points {
            if let Some((px, py)) = self.point_position(orig_ph) {
                let (rx, ry) = reflect_point(px, py, lx1, ly1, lx2, ly2);
                set_point_params(self, new_ph, rx, ry);
            }
        }
        reflect_cloned_centres(self, seed, &remap, lx1, ly1, lx2, ly2);
        sync_phantom_positions(self, seed, &remap, &phantom_points);

        let internal = clone_constraints(self, &seed_constraints, &remap);

        let mut offset_constraints = Vec::new();
        for (&orig_ph, &new_ph) in &remap.points {
            let ch = self.add_constraint(Constraint::ReflectedPosition(orig_ph, new_ph, line));
            offset_constraints.push(ch);
        }

        let pin_constraints = add_pin_constraints(self, seed, &remap, &phantom_points);
        let angle_constraints = add_mirror_angle_constraints(self, seed, &remap, line);

        self.recording = was_recording;

        let handle = self.patterns.insert(PatternDef {
            kind: PatternKind::Mirror { line },
            seed_elements: seed.to_vec(),
            seed_internal_constraints: seed_constraints,
            instances: vec![PatternInstance {
                elements: cloned_elems,
                phantom_points,
                internal_constraints: internal,
                offset_constraints,
                pin_constraints,
                angle_constraints,
                suppressed: false,
            }],
        });
        Ok(handle)
    }
}

// ── Constraint-adding helpers ─────────────────────────────────────────────

/// Add AngleOffsetArc/Conic constraints for circular patterns.
fn add_angle_offset_constraints(
    sketch: &mut Sketch,
    seed: &[ElementHandle],
    remap: &crate::pattern::HandleRemap,
    angle_k: f64,
) -> Vec<crate::ConstraintHandle> {
    let mut constraints = Vec::new();
    for elem in seed {
        match elem {
            ElementHandle::Arc(ah) => {
                if let Some(&new_ah) = remap.arcs.get(ah) {
                    constraints.push(
                        sketch.add_constraint(Constraint::AngleOffsetArc(*ah, new_ah, angle_k)),
                    );
                }
            }
            ElementHandle::Ellipse(eh) => {
                if let Some(&new_eh) = remap.ellipses.get(eh) {
                    constraints.push(sketch.add_constraint(Constraint::AngleOffsetConic(
                        ElementHandle::Ellipse(*eh),
                        ElementHandle::Ellipse(new_eh),
                        angle_k,
                    )));
                }
            }
            ElementHandle::Parabola(ph) => {
                if let Some(&new_ph) = remap.parabolas.get(ph) {
                    constraints.push(sketch.add_constraint(Constraint::AngleOffsetConic(
                        ElementHandle::Parabola(*ph),
                        ElementHandle::Parabola(new_ph),
                        angle_k,
                    )));
                }
            }
            ElementHandle::Hyperbola(hh) => {
                if let Some(&new_hh) = remap.hyperbolas.get(hh) {
                    constraints.push(sketch.add_constraint(Constraint::AngleOffsetConic(
                        ElementHandle::Hyperbola(*hh),
                        ElementHandle::Hyperbola(new_hh),
                        angle_k,
                    )));
                }
            }
            _ => {}
        }
    }
    constraints
}

/// Add MirrorAnglesArc/Conic constraints for mirror patterns.
fn add_mirror_angle_constraints(
    sketch: &mut Sketch,
    seed: &[ElementHandle],
    remap: &crate::pattern::HandleRemap,
    line: crate::LineHandle,
) -> Vec<crate::ConstraintHandle> {
    let mut constraints = Vec::new();
    for elem in seed {
        match elem {
            ElementHandle::Arc(ah) => {
                if let Some(&new_ah) = remap.arcs.get(ah) {
                    constraints.push(
                        sketch.add_constraint(Constraint::MirrorAnglesArc(*ah, new_ah, line)),
                    );
                }
            }
            ElementHandle::Ellipse(eh) => {
                if let Some(&new_eh) = remap.ellipses.get(eh) {
                    constraints.push(sketch.add_constraint(Constraint::MirrorAnglesConic(
                        ElementHandle::Ellipse(*eh),
                        ElementHandle::Ellipse(new_eh),
                        line,
                    )));
                }
            }
            ElementHandle::Parabola(ph) => {
                if let Some(&new_ph) = remap.parabolas.get(ph) {
                    constraints.push(sketch.add_constraint(Constraint::MirrorAnglesConic(
                        ElementHandle::Parabola(*ph),
                        ElementHandle::Parabola(new_ph),
                        line,
                    )));
                }
            }
            ElementHandle::Hyperbola(hh) => {
                if let Some(&new_hh) = remap.hyperbolas.get(hh) {
                    constraints.push(sketch.add_constraint(Constraint::MirrorAnglesConic(
                        ElementHandle::Hyperbola(*hh),
                        ElementHandle::Hyperbola(new_hh),
                        line,
                    )));
                }
            }
            _ => {}
        }
    }
    constraints
}
