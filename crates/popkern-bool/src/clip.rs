//! Curve clipping to face boundaries.
//!
//! Clips intersection curves to the polygonal boundaries of BRep faces.
//! Used by the intersection stage to trim curves that extend beyond face
//! boundaries due to the UV-box clipping limitations for non-orthogonal
//! plane parameterizations.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

/// Collect vertex positions from the outer loop of a face by walking half-edges.
///
/// Returns an empty `Vec` if the face's loop has no valid half-edges.
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
    loop {
        let he = &brep.half_edges[current];
        positions.push(brep.vertices[he.origin].position);
        current = he.next;
        if current == start_he {
            break;
        }
    }
    positions
}

/// Clip a degree-1 NURBS line to the interior of a face's boundary polygon.
///
/// Uses Cyrus-Beck line clipping against each edge of the face's outer loop.
/// Returns `Some((t_enter, t_exit))` if the line segment intersects the polygon,
/// or `None` if the line is entirely outside.
fn clip_line_to_face(
    brep: &BRep,
    face_handle: FaceHandle,
    curve: &NurbsCurve,
) -> Option<(f64, f64)> {
    let positions = face_vertex_positions(brep, face_handle);
    if positions.len() < 3 {
        return None;
    }

    let (t_min, t_max) = curve.parameter_range();
    let p0 = curve.point_at(t_min);
    let p1 = curve.point_at(t_max);
    let d = p1 - p0; // direction vector

    // Compute polygon centroid.
    let n = positions.len() as f64;
    let centroid = Point3::new(
        positions.iter().map(|p| p.x).sum::<f64>() / n,
        positions.iter().map(|p| p.y).sum::<f64>() / n,
        positions.iter().map(|p| p.z).sum::<f64>() / n,
    );

    // Compute the face normal from the first three vertices.
    let e1 = positions[1] - positions[0];
    let e2 = positions[2] - positions[0];
    let face_normal = e1.cross(&e2);

    // Use Cyrus-Beck in the [0, 1] parameter space of the line segment.
    let mut s_enter: f64 = 0.0;
    let mut s_exit: f64 = 1.0;

    for i in 0..positions.len() {
        let j = (i + 1) % positions.len();
        let edge_dir = positions[j] - positions[i];
        // Compute in-plane perpendicular to the edge.
        let mut inward = face_normal.cross(&edge_dir);
        // Ensure it points inward by checking against the centroid.
        let to_centroid = centroid - positions[i];
        if to_centroid.dot(&inward) < 0.0 {
            inward = Vec3::new(-inward.x, -inward.y, -inward.z);
        }

        let w = p0 - positions[i];
        let denom = d.dot(&inward);
        let numer = -(w.dot(&inward));

        if denom.abs() < 1e-15 {
            // Line is parallel to this edge.
            // numer = -dot(w, inward): positive means p0 is OUTSIDE.
            if numer > 1e-12 {
                return None; // outside this edge
            }
            continue;
        }

        let s = numer / denom;
        if denom > 0.0 {
            // Entering half-plane.
            if s > s_enter {
                s_enter = s;
            }
        } else {
            // Exiting half-plane.
            if s < s_exit {
                s_exit = s;
            }
        }

        if s_enter > s_exit + 1e-12 {
            return None;
        }
    }

    if s_enter >= s_exit - 1e-12 {
        return None;
    }

    // Convert from [0,1] line parameter to the curve's parameter domain.
    let t_enter = t_min + s_enter * (t_max - t_min);
    let t_exit = t_min + s_exit * (t_max - t_min);

    Some((t_enter, t_exit))
}

/// Clip a degree-1 intersection curve to a single face's boundary polygon.
///
/// Extends the underlying 3D line to the full extent of the face boundary,
/// producing a new `CurveOnSurface` whose endpoints lie on the face's edges.
/// Returns `None` if the curve direction does not intersect the face polygon.
///
/// This is used by the split stage to ensure each face gets a curve whose
/// endpoints land on that face's boundary, rather than on the other face's
/// boundary.
pub(crate) fn clip_curve_to_single_face(
    brep: &BRep,
    face_handle: FaceHandle,
    curve: &CurveOnSurface,
) -> Option<CurveOnSurface> {
    // Only handle degree-1 (line) curves for now.
    if curve.curve_3d.degree != 1 {
        return Some(curve.clone());
    }

    // We need to clip the 3D line to this face's polygon boundary.
    // To do this properly, we extend the line beyond its current endpoints
    // and clip to the face. We build a temporary extended curve for clipping.
    let (t_min, t_max) = curve.curve_3d.parameter_range();
    let p0 = curve.curve_3d.point_at(t_min);
    let p1 = curve.curve_3d.point_at(t_max);
    let dir = p1 - p0;
    let len = dir.length();
    if len < 1e-15 {
        return None;
    }

    // Extend the line by 10x in both directions to ensure we capture the
    // full intersection with the face boundary.
    let ext_p0 = p0 - Vec3::new(dir.x * 10.0, dir.y * 10.0, dir.z * 10.0);
    let ext_p1 = p1 + Vec3::new(dir.x * 10.0, dir.y * 10.0, dir.z * 10.0);

    let extended = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![ext_p0, ext_p1],
        weights: vec![1.0, 1.0],
    };

    // Verify the line actually lies on the face's plane. Cyrus-Beck only
    // tests against in-plane edge normals, so a line parallel to the face
    // but offset from it would falsely pass all half-plane tests.
    let face_positions = face_vertex_positions(brep, face_handle);
    if face_positions.len() >= 3 {
        let fe1 = face_positions[1] - face_positions[0];
        let fe2 = face_positions[2] - face_positions[0];
        let face_n = fe1.cross(&fe2);
        let n_len = face_n.length();
        if n_len > 1e-15 {
            let unit_n = Vec3::new(face_n.x / n_len, face_n.y / n_len, face_n.z / n_len);
            let w = p0 - face_positions[0];
            let dist = w.dot(&unit_n).abs();
            if dist > 1e-4 {
                return None; // line is not on the face's plane
            }
        }
    }

    let (s_enter, s_exit) = clip_line_to_face(brep, face_handle, &extended)?;

    // Convert clipped parameters back to 3D points.
    let clipped_start = extended.point_at(s_enter);
    let clipped_end = extended.point_at(s_exit);

    let clipped_3d = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![clipped_start, clipped_end],
        weights: vec![1.0, 1.0],
    };

    // Rebuild UV curves by projecting the clipped 3D endpoints.
    // Map from extended-curve parameter to original-curve parameter
    // to evaluate the UV curves.
    // Extended curve maps: ext_p0 corresponds to original t_min - 10*range,
    // ext_p1 corresponds to original t_max + 10*range.
    // s in [0,1] on extended → t on original = t_min + (s * 21 - 10) * range
    // But this can go outside [t_min, t_max], so we clamp and extrapolate UV linearly.
    // For UV curves, we use linear extrapolation from the original curve.
    let uv1_p0 = curve.curve_uv1.point_at(t_min);
    let uv1_p1 = curve.curve_uv1.point_at(t_max);
    let uv1_dir = uv1_p1 - uv1_p0;

    let uv2_p0 = curve.curve_uv2.point_at(t_min);
    let uv2_p1 = curve.curve_uv2.point_at(t_max);
    let uv2_dir = uv2_p1 - uv2_p0;

    // Map clipped 3D points back to parameter on the original line direction.
    let start_along = if len > 1e-15 {
        let v = clipped_start - p0;
        v.dot(&dir) / (len * len)
    } else {
        0.0
    };
    let end_along = if len > 1e-15 {
        let v = clipped_end - p0;
        v.dot(&dir) / (len * len)
    } else {
        1.0
    };

    let clipped_uv1_start = Point3::new(
        uv1_p0.x + uv1_dir.x * start_along,
        uv1_p0.y + uv1_dir.y * start_along,
        uv1_p0.z + uv1_dir.z * start_along,
    );
    let clipped_uv1_end = Point3::new(
        uv1_p0.x + uv1_dir.x * end_along,
        uv1_p0.y + uv1_dir.y * end_along,
        uv1_p0.z + uv1_dir.z * end_along,
    );
    let clipped_uv2_start = Point3::new(
        uv2_p0.x + uv2_dir.x * start_along,
        uv2_p0.y + uv2_dir.y * start_along,
        uv2_p0.z + uv2_dir.z * start_along,
    );
    let clipped_uv2_end = Point3::new(
        uv2_p0.x + uv2_dir.x * end_along,
        uv2_p0.y + uv2_dir.y * end_along,
        uv2_p0.z + uv2_dir.z * end_along,
    );

    let clipped_uv1 = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![clipped_uv1_start, clipped_uv1_end],
        weights: vec![1.0, 1.0],
    };
    let clipped_uv2 = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![clipped_uv2_start, clipped_uv2_end],
        weights: vec![1.0, 1.0],
    };

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
    // Only handle degree-1 (line) curves for now.
    if curve.curve_3d.degree != 1 {
        return Some(curve.clone());
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
