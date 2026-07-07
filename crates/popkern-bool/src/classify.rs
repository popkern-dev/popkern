//! Stage 3: Classify face fragments as inside/outside/on-boundary.
//!
//! Uses multi-ray point-in-solid voting to determine the spatial relationship
//! of each face fragment to the opposite solid. Each ray is cast from a sample
//! point on the fragment, and transversal crossings with the solid's faces are
//! counted. Odd crossings indicate inside, even indicate outside.

use std::collections::{HashMap, HashSet};

use popkern_core::{
    geometry::{
        Point3, Vec3,
        curve::CurveGeometry,
        intersection::{
            SurfaceCurveIntersection, intersect_surface_curve, surface_subdivision::Aabb,
        },
        surface::{SurfaceEval, SurfaceGeometry},
    },
    tolerance::Tolerance,
    topology::{BRep, EdgeHandle, FaceHandle},
};

use crate::error::{BoolError, Classification, PointClassification};

// ---------------------------------------------------------------------------
// 1. Ray-AABB intersection (slab test)
// ---------------------------------------------------------------------------

/// Test whether a ray intersects an axis-aligned bounding box using the slab method.
///
/// Returns `true` if the ray `origin + t * direction` (for t >= 0) passes through
/// or touches `aabb`. Handles axis-aligned rays (direction component = 0) correctly.
fn ray_intersects_aabb(origin: &Point3, direction: &Vec3, aabb: &Aabb) -> bool {
    let mut t_min = f64::NEG_INFINITY;
    let mut t_max = f64::INFINITY;

    let axes = [
        (origin.x, direction.x, aabb.min.x, aabb.max.x),
        (origin.y, direction.y, aabb.min.y, aabb.max.y),
        (origin.z, direction.z, aabb.min.z, aabb.max.z),
    ];

    for (o, d, bmin, bmax) in axes {
        if d.abs() < 1e-15 {
            // Ray is parallel to this slab — miss if origin outside.
            // (1e-15: numerical floor, no tolerance in scope)
            if o < bmin || o > bmax {
                return false;
            }
        } else {
            let inv_d = 1.0 / d;
            let mut t1 = (bmin - o) * inv_d;
            let mut t2 = (bmax - o) * inv_d;
            if t1 > t2 {
                std::mem::swap(&mut t1, &mut t2);
            }
            t_min = t_min.max(t1);
            t_max = t_max.min(t2);
            if t_min > t_max {
                return false;
            }
        }
    }

    // Ray hits if the exit is at t >= 0.
    t_max >= 0.0_f64.max(t_min)
}

// ---------------------------------------------------------------------------
// 2. Face vertex bounding range helper
// ---------------------------------------------------------------------------

/// Compute a finite UV bounding range for a face by examining its vertex
/// positions projected onto the surface parameter space.
///
/// For infinite-range surfaces (planes), projects face loop vertices onto
/// the surface's u/v axes to determine a finite bounding rectangle with
/// generous padding. Returns the surface's own parameter range if it is
/// already finite.
fn face_parameter_range(brep: &BRep, face_handle: FaceHandle) -> ((f64, f64), (f64, f64)) {
    let face = &brep.faces[face_handle];
    let surface_on_face = match face.surface.as_ref() {
        Some(s) => s,
        None => return ((-10.0, 10.0), (-10.0, 10.0)),
    };
    let ((u_min, u_max), (v_min, v_max)) = surface_on_face.geometry.parameter_range();

    // If both ranges are finite, use them directly.
    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
        return ((u_min, u_max), (v_min, v_max));
    }

    // For infinite surfaces (planes), project boundary SAMPLE
    // positions (vertices + edge-curve samples) onto the surface's
    // u/v axes. Vertex-only projection collapses to a point for a
    // face bounded by a closed rim — a cylinder cap has ONE vertex —
    // so ray classification never counted cap crossings and every
    // point below such a body classified Inside via the lone wall
    // crossing (measured: (12, 0, −0.5) under a Ø70 cylinder).
    // Vertex positions bound the box for ordinary faces — for the
    // ray-crossing counter the box acts as the face-trim proxy, so a
    // WIDER box creates phantom crossings (measured: curve-sampling
    // every face broke disjoint-cube CSG). Only faces whose vertex
    // polygon is degenerate (closed rims: a cylinder cap has ONE
    // vertex) need curve samples — without them their box collapses
    // to the seam point and cap crossings never count (a point below
    // such a body classified Inside via the lone wall crossing).
    let mut positions = face_vertex_positions(brep, face_handle);
    if positions.len() < 3 {
        positions = brep.face_outer_loop_sample_positions(face_handle, 8);
    }
    if positions.is_empty() {
        return ((-100.0, 100.0), (-100.0, 100.0));
    }

    // Extract plane origin and axes for projection.
    let (plane_origin, u_axis, v_axis) = match &surface_on_face.geometry {
        SurfaceGeometry::Plane {
            origin,
            u_axis,
            v_axis,
            ..
        } => (*origin, *u_axis, *v_axis),
        _ => {
            // Non-plane infinite surface — use large default range.
            return ((-100.0, 100.0), (-100.0, 100.0));
        }
    };

    // Project each vertex onto the plane's u/v axes.
    let mut u_lo = f64::INFINITY;
    let mut u_hi = f64::NEG_INFINITY;
    let mut v_lo = f64::INFINITY;
    let mut v_hi = f64::NEG_INFINITY;

    for p in &positions {
        let d = *p - plane_origin;
        let u_val = d.dot(&u_axis);
        let v_val = d.dot(&v_axis);
        u_lo = u_lo.min(u_val);
        u_hi = u_hi.max(u_val);
        v_lo = v_lo.min(v_val);
        v_hi = v_hi.max(v_val);
    }

    // Add generous padding.
    let u_extent = (u_hi - u_lo).max(1.0);
    let v_extent = (v_hi - v_lo).max(1.0);
    u_lo -= u_extent;
    u_hi += u_extent;
    v_lo -= v_extent;
    v_hi += v_extent;

    // Override only the infinite components.
    let u_range = (
        if u_min.is_finite() { u_min } else { u_lo },
        if u_max.is_finite() { u_max } else { u_hi },
    );
    let v_range = (
        if v_min.is_finite() { v_min } else { v_lo },
        if v_max.is_finite() { v_max } else { v_hi },
    );

    (u_range, v_range)
}

/// Collect vertex positions from the outer loop of a face by walking half-edges.
///
/// Uses a bounded iteration limit to prevent infinite loops on broken topology.
fn face_vertex_positions(brep: &BRep, face_handle: FaceHandle) -> Vec<Point3> {
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

// ---------------------------------------------------------------------------
// 3. Single-ray crossing count
// ---------------------------------------------------------------------------

/// Cast a single ray and count transversal crossings with the faces of a solid.
///
/// Returns `Ok(count)` with the number of transversal face crossings, or
/// `Err(())` if any intersection is ambiguous (near-tangent or near an edge).
fn ray_crossing_count(
    origin: &Point3,
    direction: &Vec3,
    brep: &BRep,
    face_aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> Result<usize, ()> {
    let ray_curve = CurveGeometry::Line {
        origin: *origin,
        direction: *direction,
    };
    let t_range = (0.0, 1e6);
    let mut count = 0usize;

    for (face_handle, aabb) in face_aabbs {
        if !ray_intersects_aabb(origin, direction, aabb) {
            continue;
        }

        let face = &brep.faces[*face_handle];
        let surface_on_face = match face.surface.as_ref() {
            Some(s) => s,
            None => continue,
        };
        let surface = &surface_on_face.geometry;

        let ((u_min, u_max), (v_min, v_max)) = face_parameter_range(brep, *face_handle);

        let result = intersect_surface_curve(
            surface,
            (u_min, u_max),
            (v_min, v_max),
            &ray_curve,
            t_range,
            tolerance,
        )
        .map_err(|_| ())?;

        for isect in &result.intersections {
            match isect {
                SurfaceCurveIntersection::Point { t, u, v, .. } => {
                    // Only forward intersections (t > 0).
                    if *t < tolerance.linear() {
                        continue;
                    }

                    let normal = surface.normal_at(*u, *v);
                    let n_len = normal.length();
                    if n_len < tolerance.numerical_floor() {
                        // Degenerate normal — skip.
                        continue;
                    }
                    let n_unit = normal * (1.0 / n_len);
                    let d_len = direction.length();
                    if d_len < tolerance.numerical_floor() {
                        return Err(());
                    }
                    let d_unit = *direction * (1.0 / d_len);

                    let dot = n_unit.dot(&d_unit).abs();

                    // If the hit is nearly tangent, the ray is inconclusive.
                    if dot < tolerance.angular() {
                        return Err(());
                    }

                    // Transversal crossing.
                    count += 1;
                }
                SurfaceCurveIntersection::Overlap { .. } => {
                    // Ray lies on the surface — inconclusive.
                    return Err(());
                }
            }
        }
    }

    Ok(count)
}

// ---------------------------------------------------------------------------
// 4. Point-in-solid with multi-ray voting
// ---------------------------------------------------------------------------

/// Directions for the multi-ray voting scheme.
///
/// These are intentionally non-axis-aligned to reduce the chance of hitting
/// edges, vertices, or coincident surfaces.
const RAY_DIRECTIONS: [Vec3; 6] = [
    Vec3 {
        x: 1.0,
        y: 0.1,
        z: 0.05,
    },
    Vec3 {
        x: 0.05,
        y: 1.0,
        z: 0.1,
    },
    Vec3 {
        x: 0.1,
        y: 0.05,
        z: 1.0,
    },
    Vec3 {
        x: -0.7,
        y: 0.5,
        z: 0.3,
    },
    Vec3 {
        x: 0.3,
        y: -0.7,
        z: 0.5,
    },
    Vec3 {
        x: 0.5,
        y: 0.3,
        z: -0.7,
    },
];

/// Classify a 3D point as inside, outside, or on the boundary of a solid.
///
/// Uses multi-ray voting: casts several rays from `point` and counts
/// transversal crossings with the solid's faces. An odd count means inside,
/// even means outside. A majority vote from conclusive rays determines the
/// result.
///
/// # Errors
///
/// Returns [`BoolError::ClassificationFailed`] if no conclusive ray can be
/// found after exhausting all candidate directions.
pub fn point_in_solid(
    point: &Point3,
    solid: &BRep,
    face_aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> Result<PointClassification, BoolError> {
    let mut inside_votes = 0u32;
    let mut outside_votes = 0u32;

    for direction in &RAY_DIRECTIONS {
        match ray_crossing_count(point, direction, solid, face_aabbs, tolerance) {
            Ok(count) => {
                if count % 2 == 1 {
                    inside_votes += 1;
                } else {
                    outside_votes += 1;
                }

                // Early exit: if we have a clear majority from 3+ rays, stop.
                if inside_votes >= 2 || outside_votes >= 2 {
                    break;
                }
            }
            Err(()) => {
                // Inconclusive ray — try next direction.
                continue;
            }
        }
    }

    if inside_votes > outside_votes {
        Ok(PointClassification::Inside)
    } else if outside_votes > inside_votes {
        Ok(PointClassification::Outside)
    } else {
        // Tied or no conclusive rays at all — treat as on-boundary.
        // This is a heuristic; a more robust approach would add more rays.
        Ok(PointClassification::OnBoundary)
    }
}

// ---------------------------------------------------------------------------
// 5. Classify all fragments
// ---------------------------------------------------------------------------

/// Classify every face in `fragments_brep` relative to `other_brep`.
///
/// For each face, picks a sample point (surface midpoint or first vertex) and
/// calls [`point_in_solid`] against the other solid. The resulting
/// [`PointClassification`] is mapped to [`Classification`]:
///
/// - `Inside` -> `Classification::Inside`
/// - `Outside` -> `Classification::Outside`
/// - `OnBoundary` -> `Classification::OnBoundaryCooriented` (conservative default)
///
/// # Arguments
///
/// * `fragments_brep` — the BRep whose faces are being classified.
/// * `other_brep` — the solid to test containment against.
/// * `other_face_aabbs` — precomputed face AABBs for `other_brep`.
/// * `intersection_edges` — edges that lie on the intersection (used for
///   future propagation optimization; currently unused).
/// * `tolerance` — geometric tolerance for all distance and coincidence tests.
pub(crate) fn classify_fragments(
    fragments_brep: &BRep,
    other_brep: &BRep,
    other_face_aabbs: &[(FaceHandle, Aabb)],
    intersection_edges: &HashSet<EdgeHandle>,
    tolerance: &Tolerance,
) -> Result<HashMap<FaceHandle, Classification>, BoolError> {
    // Use topological propagation: classify one seed per connected component,
    // then propagate across non-intersection edges. Falls back to per-face
    // ray casting if propagation fails.
    match crate::classify_propagation::classify_with_propagation(
        fragments_brep,
        other_brep,
        other_face_aabbs,
        intersection_edges,
        tolerance,
    ) {
        Ok((result, _ray_count)) => Ok(result),
        Err(_) => {
            // Propagation failed — fall back to per-face classification
            let mut result = HashMap::new();
            for face_handle in fragments_brep.faces.keys() {
                let sample = pick_sample_point(fragments_brep, face_handle);
                let pc = point_in_solid(&sample, other_brep, other_face_aabbs, tolerance)?;
                let class = match pc {
                    PointClassification::Inside => Classification::Inside,
                    PointClassification::Outside => Classification::Outside,
                    PointClassification::OnBoundary => Classification::OnBoundaryCooriented,
                };
                result.insert(face_handle, class);
            }
            Ok(result)
        }
    }
}

/// Pick a representative sample point on a face for classification.
///
/// Tries UV-centroid sampling first (using pcurves to find a point guaranteed
/// inside the face boundary), then falls back to edge midpoint (for faces with
/// inner loops) or vertex centroid (for simple faces). Falls back to the surface
/// parameter midpoint for faces without vertices.
pub(crate) fn pick_sample_point(brep: &BRep, face_handle: FaceHandle) -> Point3 {
    let face = &brep.faces[face_handle];

    // Try UV-centroid sampling for faces with surface geometry and pcurves.
    if let Some(pt) = uv_centroid_sample(brep, face_handle) {
        return pt;
    }

    // For faces with inner loops, use the midpoint of an outer-loop edge.
    // The centroid can fall inside a hole, giving wrong classification.
    if !face.inner_loops.is_empty() {
        let loop_data = match brep.loops.get(face.outer_loop) {
            Some(l) => l,
            None => return Point3::new(0.0, 0.0, 0.0),
        };
        let start_he = loop_data.half_edge;
        if let Some(he) = brep.half_edges.get(start_he) {
            let origin = brep.vertices[he.origin].position;
            if let Some(next_he) = brep.half_edges.get(he.next) {
                let dest = brep.vertices[next_he.origin].position;
                return Point3::new(
                    (origin.x + dest.x) * 0.5,
                    (origin.y + dest.y) * 0.5,
                    (origin.z + dest.z) * 0.5,
                );
            }
        }
    }

    // Use vertex centroid — reliable for convex polygonal faces without holes.
    let positions = face_vertex_positions(brep, face_handle);
    if !positions.is_empty() {
        let n = positions.len() as f64;
        let sum = positions.iter().fold(Vec3::new(0.0, 0.0, 0.0), |acc, p| {
            Vec3::new(acc.x + p.x, acc.y + p.y, acc.z + p.z)
        });
        return Point3::new(sum.x / n, sum.y / n, sum.z / n);
    }

    // Fall back to surface parameter midpoint.
    let face = &brep.faces[face_handle];
    if let Some(surface_on_face) = face.surface.as_ref() {
        let ((u_min, u_max), (v_min, v_max)) = face_parameter_range(brep, face_handle);
        let u_mid = (u_min + u_max) * 0.5;
        let v_mid = (v_min + v_max) * 0.5;
        return surface_on_face.geometry.point_at(u_mid, v_mid);
    }

    // Last resort — origin.
    Point3::new(0.0, 0.0, 0.0)
}

#[path = "classify_uv_sample.rs"]
mod uv_sample;
use uv_sample::uv_centroid_sample;

#[cfg(test)]
#[path = "classify_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "classify_edge_tests.rs"]
mod edge_tests;
