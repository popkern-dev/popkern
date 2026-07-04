//! UV-space splitting path for face splitting.
//!
//! Split out of `split_helpers.rs` to keep that file under the 500-line
//! structural cap. Handles both the "curve crosses boundary" case and
//! the "closed curve entirely inside face" case.

use popkern_core::{
    geometry::Point3,
    geometry::curve::CurveEval,
    geometry::nurbs_curve::NurbsCurve,
    tolerance::Tolerance,
    topology::{BRep, FaceHandle},
};

use crate::clip::{clip_curve_to_crossings, find_uv_crossings};
use crate::error::BoolDiagnostics;
use crate::intersect::IntersectionRecord;
use crate::split_inner_loop::point_inside_face_uv;

use super::FragmentMap;
use super::helpers::{ensure_crossing_vertices, split_err, walk_face_loop};

/// Try UV-space splitting for a curved face: find where the intersection
/// curve crosses the face boundary in UV space, clip to interior, and split.
///
/// Returns `true` if the face was split, `false` if it should be skipped.
pub(super) fn try_uv_split(
    brep: &mut BRep,
    face_handle: FaceHandle,
    curve_3d: &NurbsCurve,
    curve_uv: &NurbsCurve,
    tolerance: &Tolerance,
    fragments: &mut FragmentMap,
    records: &[IntersectionRecord],
) -> Result<bool, BoolDiagnostics> {
    // Check if the original curve is closed (start ≈ end).
    let (t0_orig, t1_orig) = curve_3d.parameter_range();
    let is_closed =
        tolerance.points_coincident(&curve_3d.point_at(t0_orig), &curve_3d.point_at(t1_orig));

    let crossings = match find_uv_crossings(brep, face_handle, curve_uv, tolerance) {
        Some(c) => c,
        None => {
            if std::env::var_os("PK_BOOL_TRACE").is_some() {
                eprintln!("  [bool] face {face_handle:?}: no uv crossings result");
            }
            return Ok(false);
        }
    };
    if std::env::var_os("PK_BOOL_TRACE").is_some() {
        eprintln!(
            "  [bool] face {face_handle:?}: crossings={} is_closed={is_closed}",
            crossings.len()
        );
    }

    // Sub-problem A: curve crosses boundary (>= 2 crossings).
    // Clip to the interior segment and split the face.
    if crossings.len() >= 2 {
        let t_enter = crossings[0].t_curve;
        let t_exit = crossings[crossings.len() - 1].t_curve;
        let clipped_3d = match clip_curve_to_crossings(curve_3d, t_enter, t_exit) {
            Some(c) => c,
            None => return Ok(false),
        };
        let clipped_uv = match clip_curve_to_crossings(curve_uv, t_enter, t_exit) {
            Some(c) => c,
            None => return Ok(false),
        };
        let (ts, te) = clipped_3d.parameter_range();
        let start_pt = clipped_3d.point_at(ts);
        let end_pt = clipped_3d.point_at(te);
        // Use UV crossing edge identity for all face types. With orthonormal
        // plane UV axes (from plane_from_points), the pcurve parameter maps
        // exactly to the edge parameter for both planar and curved faces.
        let endpoints_ok = ensure_crossing_vertices(
            brep,
            face_handle,
            &[&crossings[0], &crossings[crossings.len() - 1]],
            curve_3d,
            tolerance,
        );
        match endpoints_ok {
            Ok(true) => {
                // Rebuild the clipped 3D curve endpoints from the actual vertex
                // positions on the face boundary (which may differ slightly from
                // curve_3d.point_at(t) due to NURBS evaluation precision).
                let pairs = walk_face_loop(brep, face_handle)
                    .map_err(|e| split_err(face_handle, e.to_string(), records))?;
                let snap_to_vertex = |pt: Point3| -> Point3 {
                    pairs
                        .iter()
                        .map(|&(_, vh)| brep.vertices[vh].position)
                        .min_by(|a, b| {
                            a.distance_to(&pt)
                                .partial_cmp(&b.distance_to(&pt))
                                .unwrap_or(std::cmp::Ordering::Equal)
                        })
                        .unwrap_or(pt)
                };
                let snapped_start = snap_to_vertex(start_pt);
                let snapped_end = snap_to_vertex(end_pt);
                let snapped_3d = NurbsCurve {
                    degree: 1,
                    knots: vec![0.0, 0.0, 1.0, 1.0],
                    control_points: vec![snapped_start, snapped_end],
                    weights: vec![1.0, 1.0],
                };
                match brep.split_face(face_handle, &clipped_uv, &snapped_3d, tolerance) {
                    Ok(result) => {
                        fragments.origin.insert(result.face_a, face_handle);
                        fragments.origin.insert(result.face_b, face_handle);
                        for &edge in &result.new_edges {
                            fragments.intersection_edges.insert(edge);
                        }
                        return Ok(true);
                    }
                    Err(e) => return Err(split_err(face_handle, e.to_string(), records)),
                }
            }
            Ok(false) if is_closed => {
                // Endpoints not on boundary but curve is closed — the UV
                // crossings are spurious (near-tangent). Fall through to
                // the inner-loop path below.
            }
            Ok(false) => return Ok(false),
            Err(e) => {
                return Err(split_err(
                    face_handle,
                    format!("pre-split edge failed: {e}"),
                    records,
                ));
            }
        }
    }

    // Sub-problem B: closed curve entirely inside face (0 real crossings).
    // Also reached when crossings were spurious (closed curve that fell through above).
    if is_closed {
        // Use the 3D midpoint projected to the face's UV space, rather than
        // evaluating the UV curve directly. For non-orthogonal UV axes, the
        // UV curve's midpoint may fall outside the face's UV polygon even when
        // the 3D point is inside the face.
        let mid_3d = curve_3d.point_at((t0_orig + t1_orig) / 2.0);
        let mid_uv = if let Some(sof) = brep.faces[face_handle].surface.as_ref() {
            popkern_core::geometry::uv_projection::project_point_to_uv(&sof.geometry, &mid_3d)
                .map(|(u, v)| Point3::new(u, v, 0.0))
                .unwrap_or_else(|| curve_uv.point_at((t0_orig + t1_orig) / 2.0))
        } else {
            curve_uv.point_at((t0_orig + t1_orig) / 2.0)
        };
        if point_inside_face_uv(brep, face_handle, mid_uv, tolerance) {
            let n = (curve_3d.control_points.len() * 2).max(12);
            match brep.add_inner_loop(face_handle, curve_3d, curve_uv, n, tolerance) {
                Ok(result) => {
                    for &edge in &result.edges {
                        fragments.intersection_edges.insert(edge);
                    }
                    return Ok(true);
                }
                Err(e) => {
                    let msg = format!("add_inner_loop failed: {e}");
                    return Err(split_err(face_handle, msg, records));
                }
            }
        }
    }

    Ok(false)
}
