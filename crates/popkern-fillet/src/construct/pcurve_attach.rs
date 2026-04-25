//! Attach pcurves to half-edges produced by the fillet pipeline.
//!
//! After topology construction (face splits + merge), half-edges on the fillet
//! face lack pcurves relative to the fillet surface, and the trim edge
//! half-edges on the kept faces may have only approximate pcurves. This module
//! computes and wires the correct pcurves for both sides.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::helpers::project_to_uv;

/// Project a 3D point to UV and clamp to the surface's valid parameter range.
///
/// Infinite parameter bounds (e.g. planes) are left unclamped.
fn project_to_uv_clamped(surface: &SurfaceGeometry, pt: &Point3) -> (f64, f64) {
    let (u, v) = project_to_uv(surface, pt);
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    let u = if u_min.is_finite() && u_max.is_finite() {
        u.clamp(u_min, u_max)
    } else {
        u
    };
    let v = if v_min.is_finite() && v_max.is_finite() {
        v.clamp(v_min, v_max)
    } else {
        v
    };
    (u, v)
}

/// Number of sample points for projecting a 3D curve onto a surface UV space.
const PCURVE_SAMPLE_COUNT: usize = 20;

/// Degree for fitted pcurves.
const PCURVE_FIT_DEGREE: usize = 3;

/// Number of control points for fitted pcurves.
const PCURVE_FIT_CONTROL_POINTS: usize = 8;

/// Parameters for attaching pcurves to fillet face half-edges.
pub(crate) struct FilletPcurveParams<'a> {
    /// Handle of the fillet face.
    pub fillet_face: FaceHandle,
    /// The fillet surface geometry.
    pub fillet_surface: &'a SurfaceGeometry,
    /// Trim edges created by splitting face 1.
    pub trim_edges_1: &'a [EdgeHandle],
    /// Trim edges created by splitting face 2.
    pub trim_edges_2: &'a [EdgeHandle],
    /// EdgeGeometry pcurve on face 1's surface.
    pub eg_pcurve_1: &'a NurbsCurve,
    /// EdgeGeometry pcurve on face 2's surface.
    pub eg_pcurve_2: &'a NurbsCurve,
    /// 3D trim curve on face 1.
    pub trim_3d_1: &'a NurbsCurve,
    /// 3D trim curve on face 2.
    pub trim_3d_2: &'a NurbsCurve,
}

/// Attach pcurves to the fillet face and its neighboring kept faces.
///
/// For each trim edge (created by [`split_face`]):
/// - The half-edge on the **kept face** side receives the EdgeGeometry pcurve
///   (`eg_pcurve`), which is in the kept face's surface UV space.
/// - The half-edge on the **fillet face** side receives a pcurve computed by
///   sampling the 3D trim curve and projecting onto the fillet surface.
///
/// For other half-edges on the fillet face (inherited boundary edges), pcurves
/// are computed by projecting the 3D edge geometry onto the fillet surface.
///
/// Errors during pcurve computation are silently ignored (the half-edge keeps
/// whatever pcurve it had, or `None`).
pub(crate) fn attach_fillet_pcurves(brep: &mut BRep, params: &FilletPcurveParams<'_>) {
    // Wire EdgeGeometry pcurves onto the kept-face side of trim edges,
    // and fillet-surface pcurves onto the fillet-face side.
    wire_trim_edge_pcurves(
        brep,
        params.fillet_face,
        params.fillet_surface,
        params.trim_edges_1,
        params.eg_pcurve_1,
        params.trim_3d_1,
    );
    wire_trim_edge_pcurves(
        brep,
        params.fillet_face,
        params.fillet_surface,
        params.trim_edges_2,
        params.eg_pcurve_2,
        params.trim_3d_2,
    );

    // Compute pcurves for all remaining half-edges on the fillet face
    // that still lack a pcurve.
    attach_remaining_fillet_pcurves(brep, params.fillet_face, params.fillet_surface);
}

/// Wire pcurves for a set of trim edges.
///
/// For each trim edge, determines which half-edge is on the fillet face and
/// which is on the kept face. The kept-face half-edge receives `eg_pcurve`,
/// and the fillet-face half-edge receives a fitted pcurve on the fillet surface.
fn wire_trim_edge_pcurves(
    brep: &mut BRep,
    fillet_face: FaceHandle,
    fillet_surface: &SurfaceGeometry,
    trim_edges: &[EdgeHandle],
    eg_pcurve: &NurbsCurve,
    trim_3d: &NurbsCurve,
) {
    for &edge in trim_edges {
        let Some(edge_data) = brep.edges.get(edge) else {
            continue;
        };
        let he_fwd = edge_data.half_edge;
        let Some(he_fwd_data) = brep.half_edges.get(he_fwd) else {
            continue;
        };
        let he_rev = he_fwd_data.twin;

        let (fillet_he, kept_he) = classify_trim_halfedges(brep, fillet_face, he_fwd, he_rev);

        // Wire the EdgeGeometry pcurve onto the kept-face half-edge.
        if let Some(kh) = kept_he {
            brep.half_edges[kh].pcurve = Some(eg_pcurve.clone().into());
        }

        // Compute and wire a pcurve on the fillet surface for the fillet-face half-edge.
        if let Some(fh) = fillet_he
            && let Some(pcurve) = compute_pcurve_on_surface(trim_3d, fillet_surface)
        {
            // Orient the pcurve to match the half-edge traversal direction.
            let pcurve = orient_pcurve_to_halfedge(brep, fh, trim_3d, &pcurve);
            brep.half_edges[fh].pcurve = Some(pcurve.into());
        }
    }
}

/// Determine which half-edge of a trim edge is on the fillet face.
///
/// Returns `(fillet_side, kept_side)` where either may be `None` if the
/// half-edge's loop/face cannot be resolved.
fn classify_trim_halfedges(
    brep: &BRep,
    fillet_face: FaceHandle,
    he_a: HalfEdgeHandle,
    he_b: HalfEdgeHandle,
) -> (Option<HalfEdgeHandle>, Option<HalfEdgeHandle>) {
    let face_of = |he: HalfEdgeHandle| -> Option<FaceHandle> {
        let he_data = brep.half_edges.get(he)?;
        let loop_data = brep.loops.get(he_data.loop_)?;
        Some(loop_data.face)
    };

    let face_a = face_of(he_a);
    let face_b = face_of(he_b);

    match (face_a, face_b) {
        (Some(fa), Some(_)) if fa == fillet_face => (Some(he_a), Some(he_b)),
        (Some(_), Some(fb)) if fb == fillet_face => (Some(he_b), Some(he_a)),
        _ => (None, None),
    }
}

/// Orient a pcurve to match the direction of a half-edge.
///
/// If the half-edge origin is closer to the trim curve end than to its start,
/// the pcurve is reversed so it follows the half-edge direction.
fn orient_pcurve_to_halfedge(
    brep: &BRep,
    he: HalfEdgeHandle,
    trim_3d: &NurbsCurve,
    pcurve: &NurbsCurve,
) -> NurbsCurve {
    let Some(he_data) = brep.half_edges.get(he) else {
        return pcurve.clone();
    };
    let origin_pos = brep.vertices[he_data.origin].position;
    let (t_start, t_end) = trim_3d.parameter_range();
    let curve_start = trim_3d.point_at(t_start);
    let curve_end = trim_3d.point_at(t_end);

    let d_to_start = origin_pos.distance_to(&curve_start);
    let d_to_end = origin_pos.distance_to(&curve_end);

    if d_to_end < d_to_start {
        reverse_nurbs(pcurve)
    } else {
        pcurve.clone()
    }
}

/// Compute pcurves for half-edges on the fillet face that do not yet have one.
///
/// Walks all loops on the fillet face (outer loop and inner loops). For each
/// half-edge without a pcurve, samples the 3D edge curve, projects onto the
/// fillet surface, and fits a NURBS pcurve.
fn attach_remaining_fillet_pcurves(
    brep: &mut BRep,
    fillet_face: FaceHandle,
    fillet_surface: &SurfaceGeometry,
) {
    let Some(face_data) = brep.faces.get(fillet_face) else {
        return;
    };
    let outer_loop = face_data.outer_loop;
    let inner_loops = face_data.inner_loops.clone();

    // Collect half-edge handles that need pcurves from a single loop.
    let collect_from_loop =
        |brep: &BRep, loop_handle: popkern_core::topology::LoopHandle| -> Vec<HalfEdgeHandle> {
            let Some(loop_data) = brep.loops.get(loop_handle) else {
                return Vec::new();
            };
            let start_he = loop_data.half_edge;
            if !brep.half_edges.contains_key(start_he) {
                return Vec::new();
            }
            let mut handles = Vec::new();
            let mut cursor = start_he;
            let limit = brep.half_edges.len() + 1;
            for _ in 0..limit {
                if brep.half_edges[cursor].pcurve.is_none() {
                    handles.push(cursor);
                }
                cursor = brep.half_edges[cursor].next;
                if cursor == start_he {
                    break;
                }
            }
            handles
        };

    // Collect from the outer loop.
    let mut he_handles = collect_from_loop(brep, outer_loop);

    // Also process inner loops.
    for &inner_loop_handle in &inner_loops {
        he_handles.extend(collect_from_loop(brep, inner_loop_handle));
    }

    // Compute pcurves for the collected half-edges.
    for heh in he_handles {
        if let Some(pcurve) = compute_halfedge_pcurve(brep, heh, fillet_surface) {
            brep.half_edges[heh].pcurve = Some(pcurve.into());
        }
    }
}

/// Compute a pcurve for a single half-edge by projecting its 3D geometry onto a surface.
///
/// Samples the half-edge's underlying edge curve (or linear interpolation
/// between origin and destination vertices) at `PCURVE_SAMPLE_COUNT` points,
/// projects each into the surface's UV space, and fits a NURBS curve.
fn compute_halfedge_pcurve(
    brep: &BRep,
    he: HalfEdgeHandle,
    surface: &SurfaceGeometry,
) -> Option<NurbsCurve> {
    let he_data = brep.half_edges.get(he)?;
    let origin = brep.vertices[he_data.origin].position;

    // Find the destination vertex (next half-edge's origin).
    let next_data = brep.half_edges.get(he_data.next)?;
    let dest = brep.vertices[next_data.origin].position;

    // Try to get the 3D edge curve.
    let edge_data = brep.edges.get(he_data.edge)?;

    // Check half-edge direction vs edge direction.
    let edge_he = edge_data.half_edge;
    let forward = edge_he == he;

    let uv_points: Vec<Point3> = match &edge_data.curve {
        Some(curve) => {
            let n = PCURVE_SAMPLE_COUNT;
            (0..n)
                .map(|i| {
                    let t_frac = i as f64 / (n - 1) as f64;
                    let t = if forward {
                        curve.t_start + t_frac * (curve.t_end - curve.t_start)
                    } else {
                        curve.t_end - t_frac * (curve.t_end - curve.t_start)
                    };
                    let pt = curve.geometry.point_at(t);
                    let (u, v) = project_to_uv_clamped(surface, &pt);
                    Point3::new(u, v, 0.0)
                })
                .collect()
        }
        None => {
            // Linear edge: sample between origin and destination.
            let n = PCURVE_SAMPLE_COUNT;
            (0..n)
                .map(|i| {
                    let t = i as f64 / (n - 1) as f64;
                    let pt = Point3::new(
                        origin.x + t * (dest.x - origin.x),
                        origin.y + t * (dest.y - origin.y),
                        origin.z + t * (dest.z - origin.z),
                    );
                    let (u, v) = project_to_uv_clamped(surface, &pt);
                    Point3::new(u, v, 0.0)
                })
                .collect()
        }
    };

    if uv_points.len() < 2 {
        return None;
    }

    fit_pcurve(&uv_points)
}

/// Compute a pcurve by projecting a 3D NURBS curve onto a surface.
///
/// Samples the curve at `PCURVE_SAMPLE_COUNT` points, projects each onto the
/// surface's UV space, and fits a NURBS curve through the UV points.
fn compute_pcurve_on_surface(
    curve_3d: &NurbsCurve,
    surface: &SurfaceGeometry,
) -> Option<NurbsCurve> {
    let (t_start, t_end) = curve_3d.parameter_range();
    let n = PCURVE_SAMPLE_COUNT;

    let uv_points: Vec<Point3> = (0..n)
        .map(|i| {
            let t = t_start + (t_end - t_start) * (i as f64 / (n - 1) as f64);
            let pt = curve_3d.point_at(t);
            let (u, v) = project_to_uv_clamped(surface, &pt);
            Point3::new(u, v, 0.0)
        })
        .collect();

    if uv_points.len() < 2 {
        return None;
    }

    fit_pcurve(&uv_points)
}

/// Fit a degree-3 NURBS curve through UV sample points.
///
/// Uses chord-length parameterization and [`popkern_core::geometry::fitting::fit_curve`].
/// Returns `None` if fitting fails.
fn fit_pcurve(uv_points: &[Point3]) -> Option<NurbsCurve> {
    let n = uv_points.len();
    if n < 2 {
        return None;
    }

    // Chord-length parameterization.
    let mut params = vec![0.0_f64; n];
    for i in 1..n {
        params[i] = params[i - 1] + uv_points[i].distance_to(&uv_points[i - 1]);
    }
    let total = params[n - 1];
    if total < 1e-15 {
        // Degenerate: all points coincide. Return a constant pcurve.
        return Some(NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![uv_points[0], uv_points[0]],
            weights: vec![1.0, 1.0],
        });
    }
    for p in &mut params {
        *p /= total;
    }

    let degree = PCURVE_FIT_DEGREE.min(n - 1);
    let num_cp = PCURVE_FIT_CONTROL_POINTS.min(n).max(degree + 1);

    popkern_core::geometry::fitting::fit_curve(uv_points, &params, degree, num_cp).ok()
}

/// Reverse a NURBS curve (swap traversal direction).
fn reverse_nurbs(curve: &NurbsCurve) -> NurbsCurve {
    let (t_start, t_end) = curve.parameter_range();
    let mut reversed_knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|&k| t_start + t_end - k)
        .collect();
    reversed_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let mut reversed_cp = curve.control_points.clone();
    reversed_cp.reverse();

    let mut reversed_weights = curve.weights.clone();
    reversed_weights.reverse();

    NurbsCurve {
        degree: curve.degree,
        knots: reversed_knots,
        control_points: reversed_cp,
        weights: reversed_weights,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Vec3;

    #[test]
    fn compute_pcurve_on_plane() {
        let surface = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let curve_3d = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let result = compute_pcurve_on_surface(&curve_3d, &surface);
        assert!(result.is_some(), "should compute pcurve on plane");
        let pcurve = result.unwrap();
        // All control points should have z=0.
        for cp in &pcurve.control_points {
            assert!(cp.z.abs() < 1e-10, "pcurve z should be 0, got {}", cp.z);
        }
    }

    #[test]
    fn fit_pcurve_degenerate_points() {
        let pts = vec![Point3::new(0.5, 0.5, 0.0); 5];
        let result = fit_pcurve(&pts);
        assert!(
            result.is_some(),
            "degenerate points should produce a pcurve"
        );
    }

    #[test]
    fn fit_pcurve_too_few_points() {
        let result = fit_pcurve(&[]);
        assert!(result.is_none(), "empty points should return None");
        let result = fit_pcurve(&[Point3::new(0.0, 0.0, 0.0)]);
        assert!(result.is_none(), "single point should return None");
    }

    #[test]
    fn reverse_nurbs_roundtrip() {
        let curve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let reversed = reverse_nurbs(&curve);
        let (t0, t1) = reversed.parameter_range();
        let start = reversed.point_at(t0);
        let end = reversed.point_at(t1);
        assert!(start.distance_to(&Point3::new(1.0, 1.0, 0.0)) < 1e-10);
        assert!(end.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 1e-10);
    }
}
