//! Rigid-set data model and synchronisation.
//!
//! A rigid set groups sketch elements into a rigid body that moves with three
//! degrees of freedom: translation (tx, ty) and rotation (theta). Individual
//! point parameters within the set are derived from the rigid-body parameters
//! and per-point reference offsets.

use crate::Sketch;
use crate::elements::ElementHandle;
use crate::elements::PointHandle;
use popkern_math::PortableMath;
use serde::{Deserialize, Serialize};

slotmap::new_key_type! {
    /// Opaque handle to a [`RigidSet`] inside a [`Sketch`](crate::Sketch).
    pub struct RigidSetHandle;
}

/// A group of sketch elements that move together as a 2D rigid body.
///
/// The rigid body has three parameters stored in the sketch's parameter
/// vector at `param_start`: `tx`, `ty`, `theta`. Each member point's
/// position is derived as:
///
/// ```text
/// P_i = (tx + rx*cos(theta) - ry*sin(theta),
///        ty + rx*sin(theta) + ry*cos(theta))
/// ```
///
/// where `(rx, ry)` is the point's fixed offset from the initial centroid.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RigidSet {
    /// Elements belonging to this rigid set.
    pub elements: Vec<ElementHandle>,
    /// Points controlled by this rigid set.
    pub points: Vec<PointHandle>,
    /// Index of `tx` in the sketch parameter vector. `ty` is at
    /// `param_start + 1` and `theta` at `param_start + 2`.
    pub param_start: usize,
    /// Per-point offset from the initial centroid, indexed in parallel
    /// with [`points`](Self::points).
    pub reference_offsets: Vec<(f64, f64)>,
    /// Arc/conic centre parameters that also move rigidly. Each entry is
    /// `(element, ref_offset_index)` where `ref_offset_index` indexes into
    /// an auxiliary offset list stored in [`center_ref_offsets`](Self::center_ref_offsets).
    pub center_offsets: Vec<(ElementHandle, usize)>,
    /// Reference offsets for arc/conic centres, indexed by the second
    /// field of [`center_offsets`](Self::center_offsets).
    pub center_ref_offsets: Vec<(f64, f64)>,
    /// The centroid position at the time the rigid set was created.
    pub initial_centroid: (f64, f64),
}

/// Synchronise all rigid sets: write derived point positions from the
/// rigid-body parameters (tx, ty, theta).
///
/// This must be called after each Newton parameter update step, before the
/// next residual evaluation.
pub(crate) fn sync_rigid_sets(sketch: &mut Sketch) {
    // Collect keys first to avoid borrow-checker issues.
    let keys: Vec<RigidSetHandle> = sketch.rigid_sets.keys().collect();
    for key in keys {
        // Re-borrow each iteration; the set itself is read-only here.
        let rs = match sketch.rigid_sets.get(key) {
            Some(rs) => rs.clone(),
            None => continue,
        };

        let tx = sketch.params[rs.param_start];
        let ty = sketch.params[rs.param_start + 1];
        let theta = sketch.params[rs.param_start + 2];
        let (sin_t, cos_t) = theta.psin_cos();

        // Write point positions.
        for (i, &ph) in rs.points.iter().enumerate() {
            let (rx, ry) = rs.reference_offsets[i];
            let new_x = tx + rx * cos_t - ry * sin_t;
            let new_y = ty + rx * sin_t + ry * cos_t;
            if let Some(pt) = sketch.points.get(ph) {
                let [ix, iy] = pt.param_indices();
                sketch.params[ix] = new_x;
                sketch.params[iy] = new_y;
            }
        }

        // Write arc/conic centre positions.
        for (elem, ref_idx) in &rs.center_offsets {
            let (rx, ry) = rs.center_ref_offsets[*ref_idx];
            let new_x = tx + rx * cos_t - ry * sin_t;
            let new_y = ty + rx * sin_t + ry * cos_t;
            write_center_params(sketch, elem, new_x, new_y);
        }
    }
}

/// Write the centre (cx, cy) parameters of an arc or conic element.
fn write_center_params(sketch: &mut Sketch, elem: &ElementHandle, cx: f64, cy: f64) {
    match elem {
        ElementHandle::Arc(ah) => {
            if let Some(a) = sketch.arcs.get(*ah) {
                let ps = a.param_start;
                sketch.params[ps] = cx;
                sketch.params[ps + 1] = cy;
            }
        }
        ElementHandle::Ellipse(eh) => {
            if let Some(e) = sketch.ellipses.get(*eh) {
                let ps = e.param_start;
                sketch.params[ps] = cx;
                sketch.params[ps + 1] = cy;
            }
        }
        ElementHandle::Parabola(ph) => {
            if let Some(p) = sketch.parabolas.get(*ph) {
                let ps = p.param_start;
                sketch.params[ps] = cx;
                sketch.params[ps + 1] = cy;
            }
        }
        ElementHandle::Hyperbola(hh) => {
            if let Some(h) = sketch.hyperbolas.get(*hh) {
                let ps = h.param_start;
                sketch.params[ps] = cx;
                sketch.params[ps + 1] = cy;
            }
        }
        _ => {}
    }
}

/// Build a lookup table mapping individual parameter indices to their
/// rigid-body parameter info.
///
/// Returns a map from `param_index` to `(rigid_param_start, ref_x, ref_y)`.
/// This is used by the Jacobian remapping to redirect columns.
pub(crate) fn build_rigid_param_lookup(
    sketch: &Sketch,
) -> std::collections::HashMap<usize, (usize, f64, f64)> {
    let mut map = std::collections::HashMap::new();
    for rs in sketch.rigid_sets.values() {
        // Point params.
        for (i, &ph) in rs.points.iter().enumerate() {
            if let Some(pt) = sketch.points.get(ph) {
                let [ix, iy] = pt.param_indices();
                let (rx, ry) = rs.reference_offsets[i];
                map.insert(ix, (rs.param_start, rx, ry));
                map.insert(iy, (rs.param_start, rx, ry));
            }
        }
        // Center params.
        for (elem, ref_idx) in &rs.center_offsets {
            let (rx, ry) = rs.center_ref_offsets[*ref_idx];
            if let Some(ps) = center_param_start(sketch, elem) {
                map.insert(ps, (rs.param_start, rx, ry));
                map.insert(ps + 1, (rs.param_start, rx, ry));
            }
        }
    }
    map
}

/// Get the `param_start` of an element's centre parameters.
fn center_param_start(sketch: &Sketch, elem: &ElementHandle) -> Option<usize> {
    match elem {
        ElementHandle::Arc(ah) => sketch.arcs.get(*ah).map(|a| a.param_start),
        ElementHandle::Ellipse(eh) => sketch.ellipses.get(*eh).map(|e| e.param_start),
        ElementHandle::Parabola(ph) => sketch.parabolas.get(*ph).map(|p| p.param_start),
        ElementHandle::Hyperbola(hh) => sketch.hyperbolas.get(*hh).map(|h| h.param_start),
        _ => None,
    }
}

/// Check whether all element handles referenced by a constraint belong
/// to the same rigid set AND the constraint references at least two
/// distinct elements.
///
/// Single-element constraints (Fixed, OrdinateX, etc.) are never
/// considered intra-set even if their element is in a rigid set.
///
/// Returns `Some(handle)` if all elements are in the same set and there
/// are at least two, `None` otherwise.
pub(crate) fn all_in_same_rigid_set(
    sketch: &Sketch,
    constraint: &crate::constraint::Constraint,
) -> Option<RigidSetHandle> {
    let elems = constraint_element_handles(constraint);
    // Need at least 2 element references to be an "intra-set" constraint.
    if elems.len() < 2 {
        return None;
    }
    let first_set = element_rigid_set(sketch, &elems[0])?;
    for elem in &elems[1..] {
        let set = element_rigid_set(sketch, elem)?;
        if set != first_set {
            return None;
        }
    }
    Some(first_set)
}

/// Look up which rigid set an element belongs to (if any).
fn element_rigid_set(sketch: &Sketch, elem: &ElementHandle) -> Option<RigidSetHandle> {
    match elem {
        ElementHandle::Point(ph) => sketch.points.get(*ph)?.rigid_set,
        ElementHandle::Line(lh) => {
            let line = sketch.lines.get(*lh)?;
            let s = sketch.points.get(line.start)?.rigid_set?;
            let e = sketch.points.get(line.end)?.rigid_set?;
            if s == e { Some(s) } else { None }
        }
        ElementHandle::Arc(ah) => sketch.arcs.get(*ah)?.rigid_set,
        ElementHandle::Ellipse(eh) => sketch.ellipses.get(*eh)?.rigid_set,
        ElementHandle::Parabola(ph) => sketch.parabolas.get(*ph)?.rigid_set,
        ElementHandle::Hyperbola(hh) => sketch.hyperbolas.get(*hh)?.rigid_set,
        ElementHandle::Spline(sh) => {
            let sp = sketch.splines.get(*sh)?;
            let mut iter = sp.control_points.iter();
            let first = sketch.points.get(*iter.next()?)?.rigid_set?;
            for &cp in iter {
                if sketch.points.get(cp)?.rigid_set? != first {
                    return None;
                }
            }
            Some(first)
        }
    }
}

/// Extract the element handles referenced by a constraint.
///
/// Returns a flat list of all elements involved (points, lines, arcs,
/// conics, splines).
pub(crate) fn constraint_element_handles(c: &crate::constraint::Constraint) -> Vec<ElementHandle> {
    use crate::constraint::Constraint;
    match c {
        Constraint::Coincident(a, b) | Constraint::Distance(a, b, _) => {
            vec![ElementHandle::Point(*a), ElementHandle::Point(*b)]
        }
        Constraint::Fixed(p) | Constraint::OrdinateX(p, _) | Constraint::OrdinateY(p, _) => {
            vec![ElementHandle::Point(*p)]
        }
        Constraint::PointOnLine(p, l) | Constraint::PointLineDistance(p, l, _) => {
            vec![ElementHandle::Point(*p), ElementHandle::Line(*l)]
        }
        Constraint::PointOnArc(p, a) => {
            vec![ElementHandle::Point(*p), ElementHandle::Arc(*a)]
        }
        Constraint::PointOnEllipse(p, e) => {
            vec![ElementHandle::Point(*p), ElementHandle::Ellipse(*e)]
        }
        Constraint::PointOnParabola(p, pa) => {
            vec![ElementHandle::Point(*p), ElementHandle::Parabola(*pa)]
        }
        Constraint::PointOnHyperbola(p, h) => {
            vec![ElementHandle::Point(*p), ElementHandle::Hyperbola(*h)]
        }
        Constraint::Horizontal(l) | Constraint::Vertical(l) | Constraint::FixedLine(l) => {
            vec![ElementHandle::Line(*l)]
        }
        Constraint::Angle(l1, l2, _)
        | Constraint::Parallel(l1, l2)
        | Constraint::Perpendicular(l1, l2)
        | Constraint::EqualLength(l1, l2)
        | Constraint::Collinear(l1, l2)
        | Constraint::OffsetLine(l1, l2, _) => {
            vec![ElementHandle::Line(*l1), ElementHandle::Line(*l2)]
        }
        Constraint::Tangent(l, a) => {
            vec![ElementHandle::Line(*l), ElementHandle::Arc(*a)]
        }
        Constraint::ArcTangent(a1, a2)
        | Constraint::EqualRadius(a1, a2)
        | Constraint::Concentric(a1, a2)
        | Constraint::Coradial(a1, a2)
        | Constraint::OffsetArc(a1, a2, _) => {
            vec![ElementHandle::Arc(*a1), ElementHandle::Arc(*a2)]
        }
        Constraint::Symmetric(a, b, l) => {
            vec![
                ElementHandle::Point(*a),
                ElementHandle::Point(*b),
                ElementHandle::Line(*l),
            ]
        }
        Constraint::Midpoint(p, l) => {
            vec![ElementHandle::Point(*p), ElementHandle::Line(*l)]
        }
        Constraint::TangentLineEllipse(l, e) => {
            vec![ElementHandle::Line(*l), ElementHandle::Ellipse(*e)]
        }
        Constraint::TangentLineParabola(l, p) => {
            vec![ElementHandle::Line(*l), ElementHandle::Parabola(*p)]
        }
        Constraint::TangentLineHyperbola(l, h) => {
            vec![ElementHandle::Line(*l), ElementHandle::Hyperbola(*h)]
        }
        Constraint::TangentArcEllipse(a, e) => {
            vec![ElementHandle::Arc(*a), ElementHandle::Ellipse(*e)]
        }
        Constraint::TangentArcParabola(a, p) => {
            vec![ElementHandle::Arc(*a), ElementHandle::Parabola(*p)]
        }
        Constraint::TangentArcHyperbola(a, h) => {
            vec![ElementHandle::Arc(*a), ElementHandle::Hyperbola(*h)]
        }
        Constraint::ConcentricEllipse(e1, e2) | Constraint::EqualSemiAxis(e1, e2) => {
            vec![ElementHandle::Ellipse(*e1), ElementHandle::Ellipse(*e2)]
        }
        Constraint::ConcentricArcEllipse(a, e) => {
            vec![ElementHandle::Arc(*a), ElementHandle::Ellipse(*e)]
        }
        Constraint::FocalDistance(e, _) | Constraint::EllipseEccentricity(e, _) => {
            vec![ElementHandle::Ellipse(*e)]
        }
        Constraint::HyperbolaEccentricity(h, _) => {
            vec![ElementHandle::Hyperbola(*h)]
        }
        Constraint::FixedEllipse(e) => vec![ElementHandle::Ellipse(*e)],
        Constraint::FixedParabola(p) => vec![ElementHandle::Parabola(*p)],
        Constraint::FixedHyperbola(h) => vec![ElementHandle::Hyperbola(*h)],
        Constraint::PointOnSpline(p, s) => {
            vec![ElementHandle::Point(*p), ElementHandle::Spline(*s)]
        }
        Constraint::SplineEndTangent(s, _, l) => {
            vec![ElementHandle::Spline(*s), ElementHandle::Line(*l)]
        }
        Constraint::SplineEndTangentToArc(s, _, a) => {
            vec![ElementHandle::Spline(*s), ElementHandle::Arc(*a)]
        }
        Constraint::SplineFixedEndTangent(s, _, _, _)
        | Constraint::SplineFixedEndCurvature(s, _, _)
        | Constraint::SplineLength(s, _) => {
            vec![ElementHandle::Spline(*s)]
        }
        Constraint::SplineTangentToLine(s, l) => {
            vec![ElementHandle::Spline(*s), ElementHandle::Line(*l)]
        }
        Constraint::SplineTangentToArc(s, a) => {
            vec![ElementHandle::Spline(*s), ElementHandle::Arc(*a)]
        }
        Constraint::SplineTangentToEllipse(s, e) => {
            vec![ElementHandle::Spline(*s), ElementHandle::Ellipse(*e)]
        }
        Constraint::SplineTangentToParabola(s, p) => {
            vec![ElementHandle::Spline(*s), ElementHandle::Parabola(*p)]
        }
        Constraint::SplineTangentToHyperbola(s, h) => {
            vec![ElementHandle::Spline(*s), ElementHandle::Hyperbola(*h)]
        }
        Constraint::SplineG1(s1, _, s2, _) | Constraint::SplineG2(s1, _, s2, _) => {
            vec![ElementHandle::Spline(*s1), ElementHandle::Spline(*s2)]
        }
        Constraint::RelativePosition(a, b, _, _) => {
            vec![ElementHandle::Point(*a), ElementHandle::Point(*b)]
        }
        Constraint::RotatedPosition(a, b, c, _) => {
            vec![
                ElementHandle::Point(*a),
                ElementHandle::Point(*b),
                ElementHandle::Point(*c),
            ]
        }
        Constraint::ReflectedPosition(a, b, l) => {
            vec![
                ElementHandle::Point(*a),
                ElementHandle::Point(*b),
                ElementHandle::Line(*l),
            ]
        }
        Constraint::PinToCenter(p, e) => {
            vec![ElementHandle::Point(*p), *e]
        }
        Constraint::AngleOffsetArc(a1, a2, _) | Constraint::MirrorAnglesArc(a1, a2, _) => {
            vec![ElementHandle::Arc(*a1), ElementHandle::Arc(*a2)]
        }
        Constraint::AngleOffsetConic(e1, e2, _) | Constraint::MirrorAnglesConic(e1, e2, _) => {
            vec![*e1, *e2]
        }
    }
}

/// Collect the set of parameter indices that belong to rigid-set members.
///
/// Used by the DOF analysis to filter out individual point freedoms for
/// rigid-set members.
pub(crate) fn rigid_set_point_params(sketch: &Sketch) -> std::collections::HashSet<usize> {
    let mut set = std::collections::HashSet::new();
    for rs in sketch.rigid_sets.values() {
        for &ph in &rs.points {
            if let Some(pt) = sketch.points.get(ph) {
                let [ix, iy] = pt.param_indices();
                set.insert(ix);
                set.insert(iy);
            }
        }
        for (elem, _) in &rs.center_offsets {
            if let Some(ps) = center_param_start(sketch, elem) {
                set.insert(ps);
                set.insert(ps + 1);
            }
        }
    }
    set
}

#[cfg(test)]
#[path = "rigid_set_tests.rs"]
mod tests;
