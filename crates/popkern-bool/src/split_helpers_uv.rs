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

    let mut crossings = match find_uv_crossings(brep, face_handle, curve_uv, tolerance) {
        Some(c) => c,
        None => {
            if std::env::var_os("PK_BOOL_TRACE").is_some() {
                eprintln!("  [bool] face {face_handle:?}: no uv crossings result");
            }
            return Ok(false);
        }
    };
    // The crossing list accumulates per boundary edge in loop order:
    // adjacent edges each report a shared corner, near-parallel pairs
    // report repeatedly, and nothing guarantees monotonicity in the
    // curve parameter. Taking `first`/`last` from that soup as
    // enter/exit clipped garbage chords (measured: a flush cutter's
    // boundary-collinear flank trim collected 36 hits and the split
    // corrupted the face for every later record). Sort by curve
    // parameter and collapse duplicates.
    // Scoped to OPEN curves: a closed record's enter/exit pair may
    // legitimately be the complement of the sorted pair (the interior
    // span wraps through the closure), and the discovery order
    // encodes that (measured: sorting the drill-circle crossings of a
    // cube through-hole flipped the clipped arc and reopened the
    // hole).
    if !is_closed {
        crossings.sort_by(|a, b| {
            a.t_curve
                .partial_cmp(&b.t_curve)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let param_eps = (t1_orig - t0_orig).abs() * 1e-9 + 1e-15;
        crossings.dedup_by(|a, b| (a.t_curve - b.t_curve).abs() < param_eps);
    }
    if std::env::var_os("PK_BOOL_TRACE").is_some() {
        eprintln!(
            "  [bool] face {face_handle:?}: crossings={} is_closed={is_closed}",
            crossings.len()
        );
    }

    // A trim whose crossings are confined to its ENDPOINTS and whose
    // midpoint lies ON the face boundary runs ALONG the rim — it
    // partitions nothing. This covers both a coplanar cutter's record
    // duplicating the target's own edge (first pass) and a record
    // re-visited after its own split created the edge it now lies on
    // (retry pass). Splitting instead corrupted the face's junction
    // network for every later record.
    if !crossings.is_empty() {
        let interior = crossings.iter().any(|c| {
            let span = (t1_orig - t0_orig).abs().max(1e-30);
            let f = (c.t_curve - t0_orig.min(t1_orig)) / span;
            f > 1e-6 && f < 1.0 - 1e-6
        });
        if !interior
            && uv_point_on_face_boundary(
                brep,
                face_handle,
                &curve_uv.point_at(0.5 * (t0_orig + t1_orig)),
                tolerance,
            )
        {
            if std::env::var_os("PK_SPLIT_TRACE").is_some() {
                eprintln!("  [split1] {face_handle:?} rim-collinear trim — consumed as no-op");
            }
            return Ok(true);
        }
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

/// True when the UV point lies within tolerance of the face's outer
/// boundary, sampled through the boundary pcurves (16-station chords
/// per edge — exact for linear pcurves).
fn uv_point_on_face_boundary(
    brep: &BRep,
    face: FaceHandle,
    p: &popkern_core::geometry::Point3,
    tolerance: &Tolerance,
) -> bool {
    let Some(boundary) = crate::clip::face_uv_boundary(brep, face) else {
        return false;
    };
    let band = tolerance.linear() * 10.0;
    for (_, pc) in &boundary {
        let (t0, t1) = (pc.t_start, pc.t_end);
        let mut prev = pc.point_at(t0);
        for i in 1..=16 {
            let t = t0 + (t1 - t0) * i as f64 / 16.0;
            let q = pc.point_at(t);
            let ab = popkern_core::geometry::Vec3::new(q.x - prev.x, q.y - prev.y, 0.0);
            let l2 = ab.dot(&ab);
            let d = if l2 < 1e-24 {
                ((p.x - prev.x).powi(2) + (p.y - prev.y).powi(2)).sqrt()
            } else {
                let ap = popkern_core::geometry::Vec3::new(p.x - prev.x, p.y - prev.y, 0.0);
                let tt = (ap.dot(&ab) / l2).clamp(0.0, 1.0);
                ((p.x - (prev.x + ab.x * tt)).powi(2) + (p.y - (prev.y + ab.y * tt)).powi(2)).sqrt()
            };
            if d < band {
                return true;
            }
            prev = q;
        }
    }
    false
}
