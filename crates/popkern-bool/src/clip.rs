//! Curve clipping to face boundaries.
//!
//! Clips intersection curves to the polygonal boundaries of BRep faces.
//! Used by the intersection stage to trim curves that extend beyond face
//! boundaries due to the UV-box clipping limitations for non-orthogonal
//! plane parameterizations.
//!
//! Line (degree-1) clipping is in the [`clip_line`] submodule.

mod clip_line;

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::intersection::intersect_curves;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

pub(crate) use clip_line::{clip_curve_to_single_face, clip_line_to_face};

/// Two-pi constant for periodic seam handling.
const TWO_PI: f64 = std::f64::consts::TAU;

/// Returns true if any control point of a UV curve has a u-coordinate
/// near the periodic seam (u < pi/4 or u > 7*pi/4).
pub(crate) fn is_near_seam(curve: &NurbsCurve) -> bool {
    let quarter = std::f64::consts::FRAC_PI_4;
    let seven_quarter = 7.0 * quarter;
    curve
        .control_points
        .iter()
        .any(|p| p.x < quarter || p.x > seven_quarter)
}

/// Shift all control point x-coordinates by `shift`, wrapping to [0, 2*pi].
pub(crate) fn shift_curve_u(curve: &NurbsCurve, shift: f64) -> NurbsCurve {
    let shifted_points = curve
        .control_points
        .iter()
        .map(|p| {
            let mut u = p.x + shift;
            // Wrap to [0, 2*pi).
            u = u.rem_euclid(TWO_PI);
            Point3::new(u, p.y, p.z)
        })
        .collect();
    NurbsCurve {
        degree: curve.degree,
        knots: curve.knots.clone(),
        control_points: shifted_points,
        weights: curve.weights.clone(),
    }
}

/// Collect vertex positions from the outer loop of a face by walking half-edges.
///
/// Returns an empty `Vec` if the face's loop has no valid half-edges.
/// Uses a bounded iteration limit to prevent infinite loops on broken topology.
pub(crate) fn face_vertex_positions(brep: &BRep, face_handle: FaceHandle) -> Vec<Point3> {
    let face = &brep.faces[face_handle];
    let loop_handle = face.outer_loop;
    let loop_data = match brep.loops.get(loop_handle) {
        Some(l) => l,
        None => return Vec::new(),
    };
    let start_he = loop_data.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return Vec::new();
    }

    let mut positions = Vec::new();
    let mut current = start_he;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        let he = &brep.half_edges[current];
        positions.push(brep.vertices[he.origin].position);
        current = he.next;
        if current == start_he {
            break;
        }
    }
    positions
}

/// Extract the UV boundary of a face from its half-edge pcurves.
///
/// Returns a list of `(half_edge_handle, pcurve)` forming the outer loop
/// boundary in UV space. Returns `None` if any half-edge lacks a pcurve.
pub(crate) fn face_uv_boundary(
    brep: &BRep,
    face: FaceHandle,
) -> Option<Vec<(HalfEdgeHandle, NurbsCurve)>> {
    let face_data = brep.faces.get(face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return None;
    }
    let mut boundary = Vec::new();
    let mut cursor = start;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        let he = brep.half_edges.get(cursor)?;
        let pcurve = he.pcurve.as_ref()?.clone();
        boundary.push((cursor, pcurve));
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    Some(boundary)
}

/// A crossing between a UV intersection curve and a face's UV boundary.
pub(crate) struct UvCrossing {
    /// Parameter on the intersection curve at the crossing.
    pub t_curve: f64,
    /// The boundary half-edge that was crossed.
    pub boundary_he: popkern_core::topology::HalfEdgeHandle,
    /// Parameter on the boundary half-edge's pcurve at the crossing.
    pub t_boundary: f64,
}

/// Find where a UV intersection curve crosses the face's UV boundary.
///
/// Returns crossings sorted by `t_curve`. Returns `None` if the face has
/// no UV boundary (missing pcurves).
pub(crate) fn find_uv_crossings(
    brep: &BRep,
    face: FaceHandle,
    curve_uv: &NurbsCurve,
    tolerance: &Tolerance,
) -> Option<Vec<UvCrossing>> {
    let boundary = face_uv_boundary(brep, face)?;
    let mut crossings = Vec::new();

    let curve_range = curve_uv.parameter_range();
    // Skip if curve has zero-length parameter range.
    if (curve_range.1 - curve_range.0).abs() < tolerance.numerical_floor() {
        return Some(crossings);
    }

    // Detect if curves are near the periodic seam (u ≈ 0 or u ≈ 2π) and
    // need a u-shift to avoid the discontinuity. Only shift when curves are
    // in the periodic UV range [0, 2π] — NOT for small UV values from
    // non-periodic faces (e.g., plane faces with UV range [0, 1]).
    let in_periodic_range = |curve: &NurbsCurve| -> bool {
        curve
            .control_points
            .iter()
            .any(|p| p.x > std::f64::consts::PI)
    };
    let needs_shift = (is_near_seam(curve_uv) || boundary.iter().any(|(_, pc)| is_near_seam(pc)))
        && (in_periodic_range(curve_uv) || boundary.iter().any(|(_, pc)| in_periodic_range(pc)));
    let shift = std::f64::consts::PI;

    let working_curve = if needs_shift {
        shift_curve_u(curve_uv, shift)
    } else {
        curve_uv.clone()
    };

    for (he_handle, boundary_pcurve) in &boundary {
        let boundary_range = boundary_pcurve.parameter_range();
        if (boundary_range.1 - boundary_range.0).abs() < tolerance.numerical_floor() {
            continue;
        }

        // Skip seam edges: if this half-edge's twin is on the same face,
        // the edge is a parameterization seam, not a physical boundary.
        // Crossings with the seam don't divide the face into inside/outside.
        let twin_he = brep.half_edges[*he_handle].twin;
        if brep.half_edges.contains_key(twin_he)
            && brep.half_edges[twin_he].loop_ == brep.half_edges[*he_handle].loop_
        {
            continue;
        }

        let working_boundary = if needs_shift {
            shift_curve_u(boundary_pcurve, shift)
        } else {
            boundary_pcurve.clone()
        };

        let result = intersect_curves(
            &CurveGeometry::Nurbs(working_boundary),
            boundary_range,
            &CurveGeometry::Nurbs(working_curve.clone()),
            curve_range,
            tolerance,
        );

        if let Ok(result) = result {
            use popkern_core::geometry::intersection::CurveIntersection;
            for hit in &result.intersections {
                match hit {
                    CurveIntersection::Point { t1, t2, .. } => {
                        crossings.push(UvCrossing {
                            t_curve: *t2,
                            boundary_he: *he_handle,
                            t_boundary: *t1,
                        });
                    }
                    CurveIntersection::Overlap { t1_range, t2_range } => {
                        // For near-tangent overlaps, use the midpoint of the
                        // overlap range as a crossing point.
                        let t1_mid = (t1_range.0 + t1_range.1) * 0.5;
                        let t2_mid = (t2_range.0 + t2_range.1) * 0.5;
                        crossings.push(UvCrossing {
                            t_curve: t2_mid,
                            boundary_he: *he_handle,
                            t_boundary: t1_mid,
                        });
                    }
                }
            }
        }
    }

    // Sort by t_curve.
    crossings.sort_by(|a, b| {
        a.t_curve
            .partial_cmp(&b.t_curve)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    Some(crossings)
}

/// Clip a NURBS curve to the interior segment between two parameter values.
///
/// For degree-1 curves, evaluates endpoints and builds a new line segment.
/// For higher-degree curves, uses `NurbsCurve::split_at` to extract the
/// sub-curve between `t_enter` and `t_exit`.
pub(crate) fn clip_curve_to_crossings(
    curve: &NurbsCurve,
    t_enter: f64,
    t_exit: f64,
) -> Option<NurbsCurve> {
    let (t_min, t_max) = curve.parameter_range();
    let t0 = t_enter.max(t_min);
    let t1 = t_exit.min(t_max);
    if t1 - t0 < 1e-15 {
        return None;
    }
    if curve.degree == 1 {
        let p0 = curve.point_at(t0);
        let p1 = curve.point_at(t1);
        return Some(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![p0, p1],
            weights: vec![1.0, 1.0],
        });
    }
    let eps = (t_max - t_min) * 1e-10;
    let right = if t0 - t_min > eps {
        curve.split_at(t0).ok()?.1
    } else {
        curve.clone()
    };
    let (_, r_max) = right.parameter_range();
    if t1 < r_max - eps {
        right.split_at(t1).ok().map(|(left, _)| left)
    } else {
        Some(right)
    }
}

/// Compute the interior parameter range from UV-space crossings for one face.
///
/// Returns `Some((t_enter, t_exit))` if crossings define a valid interior
/// segment. Returns `None` if the curve is outside the face or ambiguous.
/// For closed curves with no crossings, returns the full parameter range
/// (the curve is entirely inside the face).
pub(crate) fn uv_range_from_crossings(
    crossings: &[UvCrossing],
    curve: &NurbsCurve,
    is_closed: bool,
) -> Option<(f64, f64)> {
    let (t_min, t_max) = curve.parameter_range();
    match crossings.len() {
        0 if is_closed => Some((t_min, t_max)),
        0 => None,
        1 => None, // tangent — ambiguous
        _ => {
            // Use first and last crossings as the interior segment.
            let t_enter = crossings[0].t_curve;
            let t_exit = crossings[crossings.len() - 1].t_curve;
            if t_exit - t_enter < 1e-15 {
                // Near-zero range — crossings are at the same point.
                // For closed curves this means tangent touch (still inside).
                if is_closed {
                    return Some((t_min, t_max));
                }
                return None;
            }
            Some((t_enter, t_exit))
        }
    }
}

/// Result of UV-space range computation for one face.
enum UvClipResult {
    /// Face has no pcurves — cannot clip in UV space.
    NoPcurves,
    /// Curve is outside the face (or ambiguous single crossing).
    Outside,
    /// Valid interior parameter range.
    Range(f64, f64),
}

/// Compute the UV clip result for a single face.
fn uv_clip_one_face(
    brep: &BRep,
    face: FaceHandle,
    curve_uv: &NurbsCurve,
    is_closed: bool,
) -> UvClipResult {
    match find_uv_crossings(brep, face, curve_uv, &brep.tolerance) {
        None => UvClipResult::NoPcurves,
        Some(crossings) => match uv_range_from_crossings(&crossings, curve_uv, is_closed) {
            Some((t_enter, t_exit)) => UvClipResult::Range(t_enter, t_exit),
            None => UvClipResult::Outside,
        },
    }
}

/// Clip a non-linear intersection curve using UV-space boundary crossings.
///
/// For each face, finds where the curve's UV representation crosses the face's
/// UV boundary, then intersects the two parameter ranges and extracts the
/// interior segment. Falls back to returning the curve unclipped if both
/// faces lack pcurves.
fn clip_curve_uv_space(
    brep_a: &BRep,
    face_a: FaceHandle,
    brep_b: &BRep,
    face_b: FaceHandle,
    curve: &CurveOnSurface,
) -> Option<CurveOnSurface> {
    let clip_a = uv_clip_one_face(brep_a, face_a, &curve.curve_uv1, curve.is_closed);
    let clip_b = uv_clip_one_face(brep_b, face_b, &curve.curve_uv2, curve.is_closed);

    // If either face says the curve is outside, return None.
    if matches!(clip_a, UvClipResult::Outside) || matches!(clip_b, UvClipResult::Outside) {
        return None;
    }

    // If both faces lack pcurves, return unclipped.
    if matches!(clip_a, UvClipResult::NoPcurves) && matches!(clip_b, UvClipResult::NoPcurves) {
        return Some(curve.clone());
    }

    let (t_min, t_max) = curve.curve_3d.parameter_range();

    // Default to full range if a face lacks pcurves.
    let (a_enter, a_exit) = match clip_a {
        UvClipResult::Range(enter, exit) => (enter, exit),
        _ => (t_min, t_max),
    };
    let (b_enter, b_exit) = match clip_b {
        UvClipResult::Range(enter, exit) => (enter, exit),
        _ => (t_min, t_max),
    };

    // Intersect the two parameter ranges.
    let t_enter = a_enter.max(b_enter);
    let t_exit = a_exit.min(b_exit);

    if t_exit - t_enter < 1e-15 {
        return None;
    }

    // For closed curves that span the full range, return as-is.
    if curve.is_closed && (t_enter - t_min).abs() < 1e-12 && (t_exit - t_max).abs() < 1e-12 {
        return Some(curve.clone());
    }

    // Clip all three curves (3D, UV1, UV2) to the parameter range.
    let clipped_3d = clip_curve_to_crossings(&curve.curve_3d, t_enter, t_exit)?;
    let clipped_uv1 = clip_curve_to_crossings(&curve.curve_uv1, t_enter, t_exit)?;
    let clipped_uv2 = clip_curve_to_crossings(&curve.curve_uv2, t_enter, t_exit)?;

    Some(CurveOnSurface {
        curve_3d: clipped_3d,
        curve_uv1: clipped_uv1,
        curve_uv2: clipped_uv2,
        is_closed: false,
    })
}

/// Clip an intersection curve to both face boundaries and produce a trimmed
/// `CurveOnSurface`. Returns `None` if the curve does not pass through both faces.
pub(crate) fn clip_curve_to_faces(
    brep_a: &BRep,
    face_a: FaceHandle,
    brep_b: &BRep,
    face_b: FaceHandle,
    curve: &CurveOnSurface,
) -> Option<CurveOnSurface> {
    // For non-linear curves, use UV-space clipping.
    if curve.curve_3d.degree != 1 {
        return clip_curve_uv_space(brep_a, face_a, brep_b, face_b, curve);
    }

    let range_a = clip_line_to_face(brep_a, face_a, &curve.curve_3d)?;
    let range_b = clip_line_to_face(brep_b, face_b, &curve.curve_3d)?;

    // Intersect the two ranges.
    let t_start = range_a.0.max(range_b.0);
    let t_end = range_a.1.min(range_b.1);

    if t_start >= t_end {
        return None;
    }

    // Build a new degree-1 NURBS for the clipped segment.
    let p_start = curve.curve_3d.point_at(t_start);
    let p_end = curve.curve_3d.point_at(t_end);

    let clipped_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p_start, p_end],
        weights: vec![1.0, 1.0],
    };

    // Rebuild UV curves by evaluating the originals at the clipped parameters.
    let uv1_start = curve.curve_uv1.point_at(t_start);
    let uv1_end = curve.curve_uv1.point_at(t_end);
    let clipped_uv1 = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![uv1_start, uv1_end],
        weights: vec![1.0, 1.0],
    };

    let uv2_start = curve.curve_uv2.point_at(t_start);
    let uv2_end = curve.curve_uv2.point_at(t_end);
    let clipped_uv2 = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![uv2_start, uv2_end],
        weights: vec![1.0, 1.0],
    };

    Some(CurveOnSurface {
        curve_3d: clipped_3d,
        curve_uv1: clipped_uv1,
        curve_uv2: clipped_uv2,
        is_closed: false,
    })
}

#[cfg(test)]
#[path = "clip_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "clip_uv_tests.rs"]
mod uv_tests;
