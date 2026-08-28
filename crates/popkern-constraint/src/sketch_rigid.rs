//! [`Sketch`] methods for creating and deleting rigid sets.
//!
//! A rigid set groups elements into a rigid body with 3 DOF (tx, ty, theta)
//! instead of 2N individual point DOFs.

use crate::Sketch;
use crate::elements::ElementHandle;
use crate::error::ConstraintError;
use crate::rigid_set::{RigidSet, RigidSetHandle};

impl Sketch {
    /// Create a rigid set from the given elements.
    ///
    /// Collects all points referenced by the elements (line endpoints,
    /// spline control points, standalone points) and arc/conic centres.
    /// Validates that no element is already in a rigid set. Computes the
    /// centroid and reference offsets, allocates 3 rigid-body parameters
    /// (tx, ty, theta), and tags all elements.
    ///
    /// Returns the handle to the new rigid set.
    pub fn create_rigid_set(
        &mut self,
        elements: &[ElementHandle],
    ) -> Result<RigidSetHandle, ConstraintError> {
        if self.recording {
            self.push_undo_snapshot();
        }

        // Validate no element already in a rigid set.
        for elem in elements {
            if self.element_in_rigid_set(elem) {
                return Err(ConstraintError::InvalidPattern(
                    "element already belongs to a rigid set".to_owned(),
                ));
            }
        }

        // Collect points and center info.
        let mut points = Vec::new();
        let mut positions = Vec::new();
        let mut center_entries: Vec<(ElementHandle, (f64, f64))> = Vec::new();

        for elem in elements {
            match elem {
                ElementHandle::Point(ph) => {
                    if let Some((x, y)) = self.point_position(*ph)
                        && !points.contains(ph)
                    {
                        points.push(*ph);
                        positions.push((x, y));
                    }
                }
                ElementHandle::Line(lh) => {
                    if let Some(line) = self.lines.get(*lh) {
                        let start = line.start;
                        let end = line.end;
                        for ph in [start, end] {
                            if let Some((x, y)) = self.point_position(ph)
                                && !points.contains(&ph)
                            {
                                points.push(ph);
                                positions.push((x, y));
                            }
                        }
                    }
                }
                ElementHandle::Arc(ah) => {
                    if let Some(a) = self.arcs.get(*ah) {
                        let ps = a.param_start;
                        let cx = self.params[ps];
                        let cy = self.params[ps + 1];
                        center_entries.push((*elem, (cx, cy)));
                    }
                }
                ElementHandle::Ellipse(eh) => {
                    if let Some(e) = self.ellipses.get(*eh) {
                        let ps = e.param_start;
                        let cx = self.params[ps];
                        let cy = self.params[ps + 1];
                        center_entries.push((*elem, (cx, cy)));
                    }
                }
                ElementHandle::Parabola(ph) => {
                    if let Some(p) = self.parabolas.get(*ph) {
                        let ps = p.param_start;
                        let vx = self.params[ps];
                        let vy = self.params[ps + 1];
                        center_entries.push((*elem, (vx, vy)));
                    }
                }
                ElementHandle::Hyperbola(hh) => {
                    if let Some(h) = self.hyperbolas.get(*hh) {
                        let ps = h.param_start;
                        let cx = self.params[ps];
                        let cy = self.params[ps + 1];
                        center_entries.push((*elem, (cx, cy)));
                    }
                }
                ElementHandle::Spline(sh) => {
                    if let Some(sp) = self.splines.get(*sh) {
                        let cps: Vec<_> = sp.control_points.clone();
                        for ph in cps {
                            if let Some((x, y)) = self.point_position(ph)
                                && !points.contains(&ph)
                            {
                                points.push(ph);
                                positions.push((x, y));
                            }
                        }
                    }
                }
            }
        }

        // Compute centroid over all positions (points + centers).
        let mut all_positions = positions.clone();
        for (_, pos) in &center_entries {
            all_positions.push(*pos);
        }
        if all_positions.is_empty() {
            return Err(ConstraintError::InvalidPattern(
                "rigid set has no positions".to_owned(),
            ));
        }
        let n = all_positions.len() as f64;
        let cx: f64 = all_positions.iter().map(|(x, _)| x).sum::<f64>() / n;
        let cy: f64 = all_positions.iter().map(|(_, y)| y).sum::<f64>() / n;

        // Compute reference offsets for points.
        let reference_offsets: Vec<(f64, f64)> =
            positions.iter().map(|(x, y)| (x - cx, y - cy)).collect();

        // Compute center reference offsets.
        let center_ref_offsets: Vec<(f64, f64)> = center_entries
            .iter()
            .map(|(_, (x, y))| (x - cx, y - cy))
            .collect();
        let center_offsets: Vec<(ElementHandle, usize)> = center_entries
            .iter()
            .enumerate()
            .map(|(i, (elem, _))| (*elem, i))
            .collect();

        // Allocate 3 rigid-body parameters: tx=cx, ty=cy, theta=0.
        let param_start = self.params.len();
        self.params.push(cx);
        self.params.push(cy);
        self.params.push(0.0);

        let rsh = self.rigid_sets.insert(RigidSet {
            elements: elements.to_vec(),
            points: points.clone(),
            param_start,
            reference_offsets,
            center_offsets: center_offsets.clone(),
            center_ref_offsets,
            initial_centroid: (cx, cy),
        });

        // Tag elements.
        for &ph in &points {
            if let Some(pt) = self.points.get_mut(ph) {
                pt.rigid_set = Some(rsh);
            }
        }
        for (elem, _) in &center_offsets {
            self.tag_element_rigid_set(elem, Some(rsh));
        }

        Ok(rsh)
    }

    /// Delete a rigid set, untagging all elements.
    ///
    /// The rigid-body parameters remain in the parameter vector (indices
    /// are stable). Does nothing if the handle is invalid.
    pub fn delete_rigid_set(&mut self, handle: RigidSetHandle) {
        if self.recording {
            self.push_undo_snapshot();
        }
        let rs = match self.rigid_sets.remove(handle) {
            Some(rs) => rs,
            None => return,
        };
        for &ph in &rs.points {
            if let Some(pt) = self.points.get_mut(ph) {
                pt.rigid_set = None;
            }
        }
        for (elem, _) in &rs.center_offsets {
            self.tag_element_rigid_set(elem, None);
        }
    }

    /// Check if an element is already in a rigid set.
    fn element_in_rigid_set(&self, elem: &ElementHandle) -> bool {
        match elem {
            ElementHandle::Point(ph) => self.points.get(*ph).and_then(|p| p.rigid_set).is_some(),
            ElementHandle::Line(lh) => {
                if let Some(line) = self.lines.get(*lh) {
                    self.points
                        .get(line.start)
                        .and_then(|p| p.rigid_set)
                        .is_some()
                        || self
                            .points
                            .get(line.end)
                            .and_then(|p| p.rigid_set)
                            .is_some()
                } else {
                    false
                }
            }
            ElementHandle::Arc(ah) => self.arcs.get(*ah).and_then(|a| a.rigid_set).is_some(),
            ElementHandle::Ellipse(eh) => {
                self.ellipses.get(*eh).and_then(|e| e.rigid_set).is_some()
            }
            ElementHandle::Parabola(ph) => {
                self.parabolas.get(*ph).and_then(|p| p.rigid_set).is_some()
            }
            ElementHandle::Hyperbola(hh) => {
                self.hyperbolas.get(*hh).and_then(|h| h.rigid_set).is_some()
            }
            ElementHandle::Spline(sh) => {
                if let Some(sp) = self.splines.get(*sh) {
                    sp.control_points
                        .iter()
                        .any(|cp| self.points.get(*cp).and_then(|p| p.rigid_set).is_some())
                } else {
                    false
                }
            }
        }
    }

    /// Set or clear the rigid_set tag on an arc/conic element.
    fn tag_element_rigid_set(&mut self, elem: &ElementHandle, rsh: Option<RigidSetHandle>) {
        match elem {
            ElementHandle::Arc(ah) => {
                if let Some(a) = self.arcs.get_mut(*ah) {
                    a.rigid_set = rsh;
                }
            }
            ElementHandle::Ellipse(eh) => {
                if let Some(e) = self.ellipses.get_mut(*eh) {
                    e.rigid_set = rsh;
                }
            }
            ElementHandle::Parabola(ph) => {
                if let Some(p) = self.parabolas.get_mut(*ph) {
                    p.rigid_set = rsh;
                }
            }
            ElementHandle::Hyperbola(hh) => {
                if let Some(h) = self.hyperbolas.get_mut(*hh) {
                    h.rigid_set = rsh;
                }
            }
            _ => {}
        }
    }
}
