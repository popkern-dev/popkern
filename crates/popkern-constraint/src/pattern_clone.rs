//! Element and constraint cloning helpers for pattern creation.
//!
//! These functions duplicate seed geometry at offset positions and
//! remap constraint handles from seed space to instance space.

use crate::Sketch;
use crate::constraint::ConstraintHandle;
use crate::elements::{ElementHandle, PointHandle};
use crate::pattern::HandleRemap;
use crate::pattern_remap;
use std::collections::HashSet;

/// Compute the "complete seed set" -- seed elements plus all points
/// they reference (line endpoints, arc centres via phantom points,
/// spline control points).
///
/// Returns the union of `seed` with any `PointHandle`s referenced by
/// non-point seed elements.
pub(crate) fn complete_seed_set(sketch: &Sketch, seed: &[ElementHandle]) -> HashSet<ElementHandle> {
    let mut set: HashSet<ElementHandle> = seed.iter().copied().collect();
    for elem in seed {
        match elem {
            ElementHandle::Line(lh) => {
                if let Some(l) = sketch.lines.get(*lh) {
                    set.insert(ElementHandle::Point(l.start));
                    set.insert(ElementHandle::Point(l.end));
                }
            }
            ElementHandle::Spline(sh) => {
                if let Some(sp) = sketch.splines.get(*sh) {
                    for &cp in &sp.control_points {
                        set.insert(ElementHandle::Point(cp));
                    }
                }
            }
            // Points, arcs, ellipses, parabolas, hyperbolas own their
            // own params -- no extra point references to add.
            _ => {}
        }
    }
    set
}

/// Clone all seed elements at a translational offset, returning the
/// remap tables and any phantom points created for arc/conic centres.
///
/// The recording flag is temporarily disabled so that individual
/// element additions do not generate undo snapshots (the caller is
/// responsible for the snapshot).
pub(crate) fn clone_elements(
    sketch: &mut Sketch,
    seed: &[ElementHandle],
    offset_x: f64,
    offset_y: f64,
) -> (Vec<ElementHandle>, HandleRemap, Vec<PointHandle>) {
    let mut remap = HandleRemap::default();
    let mut phantom_points: Vec<PointHandle> = Vec::new();
    let mut cloned_elements: Vec<ElementHandle> = Vec::new();

    // Pass 1: collect unique points from seed (line endpoints, spline CPs,
    // explicit point seeds). Clone each once at offset.
    let mut point_order: Vec<PointHandle> = Vec::new();
    for elem in seed {
        match elem {
            ElementHandle::Point(ph) => {
                point_order.push(*ph);
            }
            ElementHandle::Line(lh) => {
                if let Some(l) = sketch.lines.get(*lh) {
                    point_order.push(l.start);
                    point_order.push(l.end);
                }
            }
            ElementHandle::Spline(sh) => {
                if let Some(sp) = sketch.splines.get(*sh) {
                    for &cp in &sp.control_points {
                        point_order.push(cp);
                    }
                }
            }
            _ => {}
        }
    }
    // Deduplicate while preserving insertion order
    let mut seen_points: HashSet<PointHandle> = HashSet::new();
    point_order.retain(|ph| seen_points.insert(*ph));

    // Clone each unique point
    let was_recording = sketch.recording;
    sketch.recording = false;
    for &ph in &point_order {
        if let Some((px, py)) = sketch.point_position(ph) {
            let new_ph = sketch.add_point(px + offset_x, py + offset_y);
            remap.points.insert(ph, new_ph);
        }
    }

    // Pass 2: clone each seed element
    for elem in seed {
        match elem {
            ElementHandle::Point(ph) => {
                // Already cloned in pass 1
                if let Some(&new_ph) = remap.points.get(ph) {
                    cloned_elements.push(ElementHandle::Point(new_ph));
                }
            }
            ElementHandle::Line(lh) => {
                if let Some(l) = sketch.lines.get(*lh) {
                    let new_start = remap.points[&l.start];
                    let new_end = remap.points[&l.end];
                    let new_lh = sketch.add_line(new_start, new_end);
                    remap.lines.insert(*lh, new_lh);
                    cloned_elements.push(ElementHandle::Line(new_lh));
                }
            }
            ElementHandle::Arc(ah) => {
                if let Some((cx, cy, r, sa, ea)) = sketch.arc_params(*ah) {
                    let phantom = sketch.add_point(cx + offset_x, cy + offset_y);
                    phantom_points.push(phantom);
                    let new_ah = sketch.add_arc(cx + offset_x, cy + offset_y, r, sa, ea);
                    remap.arcs.insert(*ah, new_ah);
                    cloned_elements.push(ElementHandle::Arc(new_ah));
                }
            }
            ElementHandle::Ellipse(eh) => {
                if let Some((cx, cy, sa, sb, rot)) = sketch.ellipse_params(*eh) {
                    let phantom = sketch.add_point(cx + offset_x, cy + offset_y);
                    phantom_points.push(phantom);
                    let new_eh = sketch.add_ellipse(cx + offset_x, cy + offset_y, sa, sb, rot);
                    remap.ellipses.insert(*eh, new_eh);
                    cloned_elements.push(ElementHandle::Ellipse(new_eh));
                }
            }
            ElementHandle::Parabola(ph) => {
                if let Some((vx, vy, fd, rot)) = sketch.parabola_params(*ph) {
                    let phantom = sketch.add_point(vx + offset_x, vy + offset_y);
                    phantom_points.push(phantom);
                    let new_ph = sketch.add_parabola(vx + offset_x, vy + offset_y, fd, rot);
                    remap.parabolas.insert(*ph, new_ph);
                    cloned_elements.push(ElementHandle::Parabola(new_ph));
                }
            }
            ElementHandle::Hyperbola(hh) => {
                if let Some((cx, cy, sa, sb, rot)) = sketch.hyperbola_params(*hh) {
                    let phantom = sketch.add_point(cx + offset_x, cy + offset_y);
                    phantom_points.push(phantom);
                    let new_hh = sketch.add_hyperbola(cx + offset_x, cy + offset_y, sa, sb, rot);
                    remap.hyperbolas.insert(*hh, new_hh);
                    cloned_elements.push(ElementHandle::Hyperbola(new_hh));
                }
            }
            ElementHandle::Spline(sh) => {
                if let Some(sp) = sketch.splines.get(*sh) {
                    let new_cps: Vec<PointHandle> = sp
                        .control_points
                        .iter()
                        .map(|cp| remap.points[cp])
                        .collect();
                    let degree = sp.degree;
                    let knots = sp.knots.clone();
                    let weights = sp.weights.clone();
                    let periodic = sp.periodic;
                    let new_sh = sketch.add_spline(&new_cps, degree);
                    if let Some(new_sp) = sketch.splines.get_mut(new_sh) {
                        new_sp.knots = knots;
                        new_sp.weights = weights;
                        new_sp.periodic = periodic;
                    }
                    remap.splines.insert(*sh, new_sh);
                    cloned_elements.push(ElementHandle::Spline(new_sh));
                }
            }
        }
    }
    sketch.recording = was_recording;

    (cloned_elements, remap, phantom_points)
}

/// Find all constraints whose referenced elements are entirely within
/// the given seed set (i.e. "internal" constraints).
pub(crate) fn find_seed_constraints(
    sketch: &Sketch,
    seed_set: &HashSet<ElementHandle>,
) -> Vec<ConstraintHandle> {
    let mut result = Vec::new();
    for (ch, entry) in sketch.constraints.iter() {
        if pattern_remap::constraint_elements_in_set(&entry.constraint, sketch, seed_set) {
            result.push(ch);
        }
    }
    result
}

/// Clone a set of constraints, remapping all handles from seed space
/// to instance space using the provided remap tables.
///
/// Returns the list of newly created constraint handles. Constraints
/// with handles that don't appear in the remap are skipped (this
/// handles the rare case where a constraint references an element not
/// in the seed set).
pub(crate) fn clone_constraints(
    sketch: &mut Sketch,
    constraint_handles: &[ConstraintHandle],
    remap: &HandleRemap,
) -> Vec<ConstraintHandle> {
    // Gather constraint data first to avoid borrow conflicts.
    let constraint_data: Vec<_> = constraint_handles
        .iter()
        .filter_map(|&ch| {
            sketch
                .constraints
                .get(ch)
                .map(|e| (e.constraint.clone(), e.fixed_values.clone()))
        })
        .collect();

    let was_recording = sketch.recording;
    sketch.recording = false;

    let mut result = Vec::new();
    for (constraint, _fixed_values) in constraint_data {
        if let Some(remapped) = pattern_remap::remap_constraint(&constraint, remap) {
            let new_fixed = sketch.capture_fixed_values(&remapped);
            let ch = sketch.add_constraint(remapped);
            if let Some(entry) = sketch.constraints.get_mut(ch) {
                entry.fixed_values = new_fixed;
            }
            result.push(ch);
        }
    }

    sketch.recording = was_recording;
    result
}
