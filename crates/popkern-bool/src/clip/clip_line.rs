//! Line (degree-1) clipping helpers.
//!
//! Contains Cyrus-Beck polygon clipping for degree-1 NURBS lines and the
//! single-face line clipping used by the split stage.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use super::face_vertex_positions;

/// Outer-loop outline as an ordered polygon: every loop vertex plus
/// interior samples of each curved edge, in traversal order (each
/// half-edge's samples follow its own direction along the edge curve).
///
/// The VERTEX polygon is an under-approximation of a face bounded by
/// curved edges — it drops the sagitta bulge between junction
/// vertices. Cyrus-Beck against the chord polygon made the clip
/// window depend on WHICH junction vertices happened to exist: a
/// grind-plane trim across a coil end-cap disc is collinear with the
/// chord between its own two junction points, so with both junctions
/// already minted the trim lay ON the hull edge and clipped to an
/// empty window ([cb] s=[0.524,0.524]) — the split then failed and
/// the whole boolean declined, flipping on record ORDER (measured,
/// PK-460: coil_spring top grind). Sampling the true boundary curves
/// restores the bulge; the disc strictly contains its chords and the
/// clip is order-independent.
fn face_outline_sample_positions(
    brep: &BRep,
    face_handle: FaceHandle,
    samples_per_edge: usize,
) -> Vec<Point3> {
    let face = &brep.faces[face_handle];
    let Some(loop_data) = brep.loops.get(face.outer_loop) else {
        return Vec::new();
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
        let origin = brep.vertices[he.origin].position;
        positions.push(origin);
        if let Some(c) = brep.edges[he.edge].curve.as_ref() {
            // Straight edges carry no bulge; skip their samples so
            // planar fixtures keep byte-identical polygons.
            let curved = !matches!(
                c.geometry,
                popkern_core::geometry::curve::CurveGeometry::Line { .. }
            );
            if curved {
                // The half-edge may traverse the edge curve in either
                // direction: sample from whichever end sits at this
                // half-edge's origin.
                let forward = c.geometry.point_at(c.t_start).distance_to(&origin)
                    <= c.geometry.point_at(c.t_end).distance_to(&origin);
                let n = samples_per_edge.max(2);
                for i in 1..n {
                    let f = i as f64 / n as f64;
                    let t = if forward {
                        c.t_start + (c.t_end - c.t_start) * f
                    } else {
                        c.t_end + (c.t_start - c.t_end) * f
                    };
                    positions.push(c.geometry.point_at(t));
                }
            }
        }
        current = he.next;
        if current == start_he {
            break;
        }
    }
    positions
}

/// Clip a degree-1 NURBS line to the interior of a face's boundary polygon.
///
/// Uses Cyrus-Beck line clipping against each edge of the face's outer
/// loop VERTEX polygon (the long-tuned path). When that polygon clips
/// the line to a DEGENERATE window, the clip retries against the
/// curve-sampled outline: the vertex polygon under-covers a curved
/// boundary by the sagitta, so a trim collinear with the chord between
/// its own two junction vertices lies ON the hull edge and clips to
/// nothing — whether it happens depends on WHICH junction vertices
/// exist, i.e. record order (measured, PK-460: the coil end-cap disc's
/// grind trim clipped to s=[0.524,0.524] and the whole boolean
/// declined). The rescue is scoped to the degenerate outcome so every
/// healthy clip keeps its byte-identical tuned window.
///
/// Returns `Some((t_enter, t_exit))` if the line segment intersects the polygon,
/// or `None` if the line is entirely outside.
pub(crate) fn clip_line_to_face(
    brep: &BRep,
    face_handle: FaceHandle,
    curve: &NurbsCurve,
) -> Option<(f64, f64)> {
    let vertex_polygon = face_vertex_positions(brep, face_handle);
    match clip_line_to_polygon(&vertex_polygon, curve, face_handle) {
        Some((s_enter, s_exit)) if s_exit - s_enter < 1.0e-3 => {
            // Degenerate window — retry on the curve-sampled outline
            // (identical to the vertex polygon when every boundary
            // edge is straight, so only curved boundaries re-clip).
            let sampled = face_outline_sample_positions(brep, face_handle, 8);
            let rescued = if sampled.len() > vertex_polygon.len() {
                clip_line_to_polygon(&sampled, curve, face_handle)
            } else {
                None
            };
            let (s_enter, s_exit) = rescued.unwrap_or((s_enter, s_exit));
            let (t_min, t_max) = curve.parameter_range();
            Some((
                t_min + s_enter * (t_max - t_min),
                t_min + s_exit * (t_max - t_min),
            ))
        }
        Some((s_enter, s_exit)) => {
            let (t_min, t_max) = curve.parameter_range();
            Some((
                t_min + s_enter * (t_max - t_min),
                t_min + s_exit * (t_max - t_min),
            ))
        }
        None => None,
    }
}

/// Cyrus-Beck clip of the line against one polygon; returns the
/// `[0, 1]` line-parameter window.
fn clip_line_to_polygon(
    positions: &[Point3],
    curve: &NurbsCurve,
    face_handle: FaceHandle,
) -> Option<(f64, f64)> {
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

        // Parallelism and side tests on UNNORMALIZED products are
        // scale-dependent knife edges: |inward| carries the face
        // normal times the edge length, so a line lying exactly ON a
        // boundary edge produced denom ≈ 1e-12 — past the absolute
        // 1e-15 gate — and the division fabricated an enter parameter
        // of ~100, emptying the range (measured: a coplanar cutter's
        // top×flank record silently discarded, the wall never split
        // at the relief depth). Normalize to an angular test.
        let scale = d.length().max(1e-30) * inward.length().max(1e-30);
        if denom.abs() / scale < 1e-12 {
            // Line is parallel to this edge. Compare the SIGNED
            // IN-PLANE DISTANCE of the line from the edge — `numer`
            // is scaled by the unnormalized inward vector, so an
            // absolute threshold is a scale-dependent knife edge. A
            // line lying ON the edge (a coplanar cutter's face whose
            // intersection with the neighbour runs along this face's
            // own boundary) must clip to the edge's span via the
            // remaining edges, not vanish (measured: the flush
            // relief's cutter-top × flank-wall record was silently
            // discarded and the wall never split at the relief
            // depth).
            let inward_len = inward.length().max(1e-30);
            if numer / inward_len > 1e-9 {
                return None; // strictly outside this edge
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

    if popkern_core::trace::trace_enabled("clip") {
        eprintln!(
            "  [cb] {face_handle:?}: s=[{s_enter:.3},{s_exit:.3}] n_verts={}",
            positions.len()
        );
    }
    Some((s_enter, s_exit))
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
