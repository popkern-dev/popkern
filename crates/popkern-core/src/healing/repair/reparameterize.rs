//! Reparameterization repair: construct a new UV mapping for NURBS surfaces
//! that improves parameter-space uniformity while preserving 3D shape.
//!
//! Uses arc-length reparameterization: resample the surface at approximately
//! equal 3D intervals and refit as a new NURBS surface.

use crate::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{NurbsSurface, Point3};
use crate::healing::curvature::fundamental_forms::first_fundamental_form;
use crate::topology::{BRep, FaceHandle};

/// Result of reparameterizing a surface.
#[derive(Debug, Clone)]
pub struct ReparamResult {
    /// The face that was reparameterized.
    pub face: FaceHandle,
    /// Area element ratio before reparameterization.
    pub distortion_before: f64,
    /// Area element ratio after reparameterization.
    pub distortion_after: f64,
    /// Maximum 3D deviation between old and new surface.
    pub max_deviation: f64,
    /// Whether reparameterization improved uniformity.
    pub improved: bool,
}

/// Number of samples for arc-length table computation.
const ARC_LENGTH_SAMPLES: usize = 100;

/// Deviation verification grid size (per dimension).
const DEVIATION_GRID: usize = 20;

/// Threshold below which an area element is considered degenerate.
const AREA_ELEMENT_EPSILON: f64 = 1e-20;

/// Measure the parameterization distortion of a surface.
///
/// Returns the area element ratio: max(area_element) / min(area_element)
/// across a sampled grid. A ratio of 1.0 means perfectly uniform.
/// Returns `None` if the surface has degenerate points (area_element near 0
/// everywhere).
pub fn measure_distortion(
    surface: &dyn SurfaceEval,
    u_range: (f64, f64),
    v_range: (f64, f64),
    grid_size: usize,
) -> Option<f64> {
    if grid_size < 2 {
        return None;
    }

    let u_span = u_range.1 - u_range.0;
    let v_span = v_range.1 - v_range.0;
    if u_span <= 0.0 || v_span <= 0.0 {
        return None;
    }

    let mut min_ae = f64::INFINITY;
    let mut max_ae = 0.0_f64;

    for i in 0..grid_size {
        let u = u_range.0 + u_span * (i as f64) / ((grid_size - 1) as f64);
        for j in 0..grid_size {
            let v = v_range.0 + v_span * (j as f64) / ((grid_size - 1) as f64);
            let ff = first_fundamental_form(surface, u, v);
            let ae = ff.area_element();

            // Skip near-zero values (singularities, e.g. poles of a sphere).
            if ae < AREA_ELEMENT_EPSILON {
                continue;
            }

            if ae < min_ae {
                min_ae = ae;
            }
            if ae > max_ae {
                max_ae = ae;
            }
        }
    }

    // If all samples were degenerate, return None.
    if min_ae == f64::INFINITY || max_ae < AREA_ELEMENT_EPSILON {
        return None;
    }

    Some(max_ae / min_ae)
}

/// Compute arc-length parameterization along one parameter direction.
///
/// Returns a vector of (cumulative_arc_length, parameter) pairs for
/// interpolation. The first entry is always (0.0, varying_range.0).
fn arc_length_table(
    surface: &dyn SurfaceEval,
    fixed_param: f64,
    varying_range: (f64, f64),
    is_u_varying: bool,
    num_samples: usize,
) -> Vec<(f64, f64)> {
    let mut table = Vec::with_capacity(num_samples);
    let span = varying_range.1 - varying_range.0;

    let param_at =
        |i: usize| -> f64 { varying_range.0 + span * (i as f64) / ((num_samples - 1) as f64) };

    let point_at = |t: f64| -> Point3 {
        if is_u_varying {
            surface.point_at(t, fixed_param)
        } else {
            surface.point_at(fixed_param, t)
        }
    };

    let mut cumulative = 0.0;
    let mut prev_point = point_at(param_at(0));
    table.push((0.0, param_at(0)));

    for i in 1..num_samples {
        let t = param_at(i);
        let p = point_at(t);
        cumulative += prev_point.distance_to(&p);
        table.push((cumulative, t));
        prev_point = p;
    }

    table
}

/// Look up the parameter value for a target arc length.
///
/// Binary searches the arc-length table and linearly interpolates
/// between the two bracketing entries.
fn arc_length_to_param(table: &[(f64, f64)], target_s: f64) -> f64 {
    if table.is_empty() {
        return 0.0;
    }

    // Clamp to range.
    let total = table.last().map_or(0.0, |e| e.0);
    if target_s <= 0.0 {
        return table[0].1;
    }
    if target_s >= total {
        return table.last().map_or(0.0, |e| e.1);
    }

    // Binary search for the interval containing target_s.
    let mut lo = 0;
    let mut hi = table.len() - 1;
    while lo + 1 < hi {
        let mid = (lo + hi) / 2;
        if table[mid].0 <= target_s {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    // Linearly interpolate between table[lo] and table[hi].
    let s0 = table[lo].0;
    let s1 = table[hi].0;
    let t0 = table[lo].1;
    let t1 = table[hi].1;

    let ds = s1 - s0;
    if ds.abs() < 1e-30 {
        return t0;
    }

    t0 + (t1 - t0) * (target_s - s0) / ds
}

/// Build a uniform clamped knot vector for the `[0, 1]` domain.
///
/// Produces `num_points + degree + 1` knots with `degree + 1` repeated
/// values at each end and uniform interior spacing.
fn uniform_clamped_knots(num_points: usize, degree: usize) -> Vec<f64> {
    let n_interior = num_points.saturating_sub(degree);
    let total_knots = num_points + degree + 1;
    let mut knots = Vec::with_capacity(total_knots);

    knots.extend(std::iter::repeat_n(0.0, degree + 1));

    if n_interior > 1 {
        for i in 1..n_interior {
            knots.push(i as f64 / n_interior as f64);
        }
    }

    knots.extend(std::iter::repeat_n(1.0, degree + 1));

    knots
}

/// Reparameterize a NURBS surface for improved UV uniformity.
///
/// Resamples the surface at approximately equal 3D arc-length intervals
/// and constructs a new NURBS surface through the resampled points.
/// The new surface has the same 3D shape (within `max_deviation`) but
/// a more uniform UV mapping.
///
/// Returns `None` if the distortion is already below `target_ratio`,
/// or if reparameterization would exceed `max_deviation`.
pub fn reparameterize_surface(
    surface: &NurbsSurface,
    u_range: (f64, f64),
    v_range: (f64, f64),
    target_ratio: f64,
    max_deviation: f64,
) -> Option<(NurbsSurface, ReparamResult)> {
    // Measure current distortion.
    let distortion_before = measure_distortion(surface, u_range, v_range, 10)?;

    if distortion_before <= target_ratio {
        return None;
    }

    // Compute arc-length parameterization in u-direction.
    // Fix v at the midpoint of the range.
    let v_mid = (v_range.0 + v_range.1) * 0.5;
    let u_table = arc_length_table(surface, v_mid, u_range, true, ARC_LENGTH_SAMPLES);
    let u_total = u_table.last().map_or(0.0, |e| e.0);

    // Compute arc-length parameterization in v-direction.
    // Fix u at the midpoint of the range.
    let u_mid = (u_range.0 + u_range.1) * 0.5;
    let v_table = arc_length_table(surface, u_mid, v_range, false, ARC_LENGTH_SAMPLES);
    let v_total = v_table.last().map_or(0.0, |e| e.0);

    if u_total < 1e-30 || v_total < 1e-30 {
        return None;
    }

    // Resample on arc-length-uniform grid.
    let n_u = surface.num_u;
    let n_v = surface.num_v;

    // Guard against degenerate surfaces with fewer than 2 control points.
    if n_u < 2 || n_v < 2 {
        return None;
    }

    let mut resampled = Vec::with_capacity(n_u * n_v);
    for i in 0..n_u {
        let s_u = u_total * (i as f64) / ((n_u - 1) as f64);
        let u_old = arc_length_to_param(&u_table, s_u);
        for j in 0..n_v {
            let s_v = v_total * (j as f64) / ((n_v - 1) as f64);
            let v_old = arc_length_to_param(&v_table, s_v);
            resampled.push(surface.point_at(u_old, v_old));
        }
    }

    // Construct new NURBS surface.
    let degree_u = surface.degree_u;
    let degree_v = surface.degree_v;

    // Ensure we have enough control points for the requested degree.
    if n_u <= degree_u || n_v <= degree_v {
        return None;
    }

    let knots_u = uniform_clamped_knots(n_u, degree_u);
    let knots_v = uniform_clamped_knots(n_v, degree_v);
    let weights = vec![1.0; n_u * n_v];

    let new_surface = NurbsSurface {
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        control_points: resampled,
        num_u: n_u,
        num_v: n_v,
        weights,
    };

    // Verify deviation: sample both surfaces at DEVIATION_GRID x DEVIATION_GRID.
    let new_range = new_surface.parameter_range();
    let mut worst_deviation = 0.0_f64;

    for i in 0..DEVIATION_GRID {
        let u_new = new_range.0.0
            + (new_range.0.1 - new_range.0.0) * (i as f64) / ((DEVIATION_GRID - 1) as f64);
        for j in 0..DEVIATION_GRID {
            let v_new = new_range.1.0
                + (new_range.1.1 - new_range.1.0) * (j as f64) / ((DEVIATION_GRID - 1) as f64);
            let p_new = new_surface.point_at(u_new, v_new);

            // Find closest point on the old surface.
            let dist = match surface.closest_point(&p_new) {
                Some((cu, cv)) => {
                    let p_old = surface.point_at(cu, cv);
                    p_new.distance_to(&p_old)
                }
                None => {
                    // Fallback: sample old surface at a corresponding
                    // parameter position and use direct distance.
                    let u_frac = if (new_range.0.1 - new_range.0.0).abs() > 1e-30 {
                        (u_new - new_range.0.0) / (new_range.0.1 - new_range.0.0)
                    } else {
                        0.0
                    };
                    let v_frac = if (new_range.1.1 - new_range.1.0).abs() > 1e-30 {
                        (v_new - new_range.1.0) / (new_range.1.1 - new_range.1.0)
                    } else {
                        0.0
                    };
                    let u_old = u_range.0 + u_frac * (u_range.1 - u_range.0);
                    let v_old = v_range.0 + v_frac * (v_range.1 - v_range.0);
                    let p_old = surface.point_at(u_old, v_old);
                    p_new.distance_to(&p_old)
                }
            };

            if dist > worst_deviation {
                worst_deviation = dist;
            }
        }
    }

    if worst_deviation > max_deviation {
        return None;
    }

    // Measure new distortion.
    let distortion_after =
        measure_distortion(&new_surface, new_range.0, new_range.1, 10).unwrap_or(distortion_before);

    let improved = distortion_after < distortion_before;

    // Use a dummy face handle for the surface-level result.
    let result = ReparamResult {
        face: FaceHandle::default(),
        distortion_before,
        distortion_after,
        max_deviation: worst_deviation,
        improved,
    };

    Some((new_surface, result))
}

/// Reparameterize a face's NURBS surface in place.
///
/// Measures distortion, reparameterizes if above `target_ratio`,
/// and replaces the face's surface geometry if improved.
/// Pcurves on the face's half-edges are invalidated (set to `None`)
/// because they reference the old UV space.
pub fn reparameterize_face(
    brep: &mut BRep,
    face: FaceHandle,
    target_ratio: f64,
    max_deviation: f64,
) -> Option<ReparamResult> {
    // Get face surface; must be NURBS.
    let face_data = brep.faces.get(face)?;
    let sof = face_data.surface.as_ref()?;

    let nurbs = match &sof.geometry {
        SurfaceGeometry::Nurbs(n) => n.clone(),
        _ => return None,
    };
    let orientation = sof.orientation;

    // Get bounded parameter range from face vertices.
    let positions = brep.face_outer_loop_positions(face);
    if positions.is_empty() {
        return None;
    }
    let (u_range, v_range) =
        crate::geometry::uv_projection::tight_bounded_parameter_range(&positions, &sof.geometry);

    // Attempt reparameterization.
    let (new_nurbs, mut result) =
        reparameterize_surface(&nurbs, u_range, v_range, target_ratio, max_deviation)?;

    result.face = face;

    if !result.improved {
        return Some(result);
    }

    // Replace the face's surface geometry.
    if let Some(f) = brep.faces.get_mut(face) {
        f.surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Nurbs(new_nurbs),
            orientation,
        });
    }

    // Invalidate pcurves on all half-edges of this face: they reference
    // the old UV space and are no longer valid.
    clear_face_pcurves(brep, face);

    Some(result)
}

/// Clear pcurves on all half-edges belonging to all loops of a face.
fn clear_face_pcurves(brep: &mut BRep, face: FaceHandle) {
    let Some(face_data) = brep.faces.get(face) else {
        return;
    };

    // Collect all loop handles (outer + inner).
    let mut loop_handles = vec![face_data.outer_loop];
    loop_handles.extend_from_slice(&face_data.inner_loops);

    for lh in loop_handles {
        let Some(loop_) = brep.loops.get(lh) else {
            continue;
        };
        let start_he = loop_.half_edge;
        if !brep.half_edges.contains_key(start_he) {
            continue;
        }

        // Walk the loop and collect half-edge handles.
        let mut he_handles = Vec::new();
        let mut current = start_he;
        loop {
            he_handles.push(current);
            let next = brep.half_edges[current].next;
            current = next;
            if current == start_he {
                break;
            }
        }

        // Clear pcurves.
        for heh in he_handles {
            if let Some(he) = brep.half_edges.get_mut(heh) {
                he.pcurve = None;
            }
        }
    }
}

/// Reparameterize all distorted NURBS faces in a BRep.
///
/// Iterates every face and calls [`reparameterize_face`]. Returns a list
/// of results for faces that were examined (both improved and unimproved).
pub fn reparameterize_all(
    brep: &mut BRep,
    target_ratio: f64,
    max_deviation: f64,
) -> Vec<ReparamResult> {
    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
    let mut results = Vec::new();
    for fh in face_handles {
        if let Some(result) = reparameterize_face(brep, fh, target_ratio, max_deviation) {
            results.push(result);
        }
    }
    results
}

#[cfg(test)]
#[path = "reparameterize_tests.rs"]
mod tests;
