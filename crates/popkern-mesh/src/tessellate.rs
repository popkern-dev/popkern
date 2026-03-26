//! Core tessellation: per-face meshing and top-level BRep tessellation.

use std::collections::HashMap;

use popkern_core::geometry::{SurfaceEval, SurfaceGeometry, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use crate::adaptive::adaptive_subdivide;
use crate::mesh::{FaceMesh, MeshVertex, TriMesh};
use crate::{MeshError, TessellationOptions};

/// Tessellate a single BRep face into a triangle mesh.
pub(crate) fn tessellate_face(
    brep: &BRep,
    face_handle: FaceHandle,
    options: &TessellationOptions,
) -> Result<FaceMesh, MeshError> {
    let face = brep
        .faces
        .get(face_handle)
        .ok_or_else(|| MeshError::TessellationFailed("invalid face handle".into()))?;

    let surface_on_face = face.surface.as_ref().ok_or(MeshError::NoSurfaceGeometry)?;

    let surface = &surface_on_face.geometry;
    let orientation = surface_on_face.orientation;

    // Get bounded parameter range (handles infinite planes).
    let (u_range, v_range) = bounded_parameter_range(brep, face_handle, surface);

    // Adaptive subdivision of UV space.
    let cells = adaptive_subdivide(
        surface,
        u_range,
        v_range,
        options.min_subdivisions,
        options.max_subdivisions,
        options.chord_height,
    );

    // Deduplicate vertices by UV coordinate using bit-exact keys.
    let mut vertex_map: HashMap<(u64, u64), usize> = HashMap::new();
    let mut vertices: Vec<MeshVertex> = Vec::new();
    let mut triangles: Vec<[usize; 3]> = Vec::new();

    for cell in &cells {
        // Four corners of the quad cell.
        let corners = [
            (cell.u_min, cell.v_min),
            (cell.u_max, cell.v_min),
            (cell.u_min, cell.v_max),
            (cell.u_max, cell.v_max),
        ];

        let mut indices = [0_usize; 4];
        for (i, &(u, v)) in corners.iter().enumerate() {
            let key = (u.to_bits(), v.to_bits());
            indices[i] = *vertex_map.entry(key).or_insert_with(|| {
                let pos = surface.point_at(u, v);
                let mut normal = surface.normal_at(u, v);
                // Normalize the normal vector.
                let len = normal.length();
                if len > 1e-12 {
                    normal = Vec3::new(normal.x / len, normal.y / len, normal.z / len);
                }
                // Flip normal if face orientation is reversed.
                if !orientation {
                    normal = -normal;
                }
                let idx = vertices.len();
                vertices.push(MeshVertex {
                    position: pos,
                    normal,
                    uv: (u, v),
                });
                idx
            });
        }

        // Split quad into 2 triangles along the shorter diagonal.
        // Indices: 0=bottom-left, 1=bottom-right, 2=top-left, 3=top-right
        let diag_02_13 = vertices[indices[0]]
            .position
            .distance_to(&vertices[indices[3]].position);
        let diag_01_23 = vertices[indices[1]]
            .position
            .distance_to(&vertices[indices[2]].position);

        if diag_02_13 <= diag_01_23 {
            // Split along 0-3 diagonal.
            triangles.push([indices[0], indices[1], indices[3]]);
            triangles.push([indices[0], indices[3], indices[2]]);
        } else {
            // Split along 1-2 diagonal.
            triangles.push([indices[0], indices[1], indices[2]]);
            triangles.push([indices[1], indices[3], indices[2]]);
        }
    }

    // Trim: discard tessellation triangles whose UV centroids fall outside
    // the face boundary. This handles self-loop caps (circular trim on
    // plane), NURBS body surfaces (trim at sweep endpoints), and any face
    // whose UV parameter range is larger than its actual boundary.
    // Skip trimming for degenerate boundaries (collinear seam edges on
    // spheres, single-point loops) where the polygon has negligible area.
    let boundary_polygon = face_boundary_uv_polygon(brep, face_handle, surface);
    let poly_area = polygon_area(&boundary_polygon);
    let uv_area = (u_range.1 - u_range.0) * (v_range.1 - v_range.0);
    if boundary_polygon.len() >= 3 && poly_area > uv_area * 0.01 {
        triangles.retain(|tri| {
            let c_u = (vertices[tri[0]].uv.0 + vertices[tri[1]].uv.0 + vertices[tri[2]].uv.0) / 3.0;
            let c_v = (vertices[tri[0]].uv.1 + vertices[tri[1]].uv.1 + vertices[tri[2]].uv.1) / 3.0;
            point_in_polygon(c_u, c_v, &boundary_polygon)
        });
    }

    Ok(FaceMesh {
        face: face_handle,
        vertices,
        triangles,
    })
}

/// Sample the face's outer loop boundary as a UV polygon.
///
/// Walks the outer loop's edges, projects curve sample points into UV
/// space, and returns the polygon vertices. Used for trimming
/// tessellation triangles that fall outside the face boundary.
fn face_boundary_uv_polygon(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
) -> Vec<(f64, f64)> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;

    let Some(face) = brep.faces.get(face_handle) else {
        return Vec::new();
    };
    let Some(loop_) = brep.loops.get(face.outer_loop) else {
        return Vec::new();
    };
    let start_he = loop_.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return Vec::new();
    }

    let mut uv_points = Vec::new();
    let mut current = start_he;
    let samples_per_edge = 32;

    loop {
        let he = &brep.half_edges[current];
        let edge = &brep.edges[he.edge];

        if let Some(ref curve) = edge.curve {
            // Determine traversal direction: compare this half-edge's origin
            // vertex with the curve's start point.  If they match, sample
            // forward (t_start → t_end); otherwise backward.
            let curve_start = curve.geometry.point_at(curve.t_start);
            let he_origin = brep.vertices[he.origin].position;
            let forward = curve_start.distance_to(&he_origin) < 1e-6;
            let (t_lo, t_hi) = if forward {
                (curve.t_start, curve.t_end)
            } else {
                (curve.t_end, curve.t_start)
            };

            for i in 0..samples_per_edge {
                let t = t_lo + (t_hi - t_lo) * (i as f64 / samples_per_edge as f64);
                let p = curve.geometry.point_at(t);
                if let Some(uv) = project_point_to_uv(surface, &p) {
                    uv_points.push(uv);
                }
            }
        } else {
            // No curve — just use the vertex position.
            let p = brep.vertices[he.origin].position;
            if let Some(uv) = project_point_to_uv(surface, &p) {
                uv_points.push(uv);
            }
        }

        current = he.next;
        if current == start_he {
            break;
        }
    }

    uv_points
}

/// Winding-number point-in-polygon test in 2D.
///
/// Returns true if the point `(px, py)` is inside the polygon defined
/// by `polygon` (a closed sequence of (x, y) vertices).
fn point_in_polygon(px: f64, py: f64, polygon: &[(f64, f64)]) -> bool {
    if polygon.len() < 3 {
        return true; // Degenerate polygon — don't trim.
    }
    let mut winding = 0i32;
    let n = polygon.len();
    for i in 0..n {
        let (x1, y1) = polygon[i];
        let (x2, y2) = polygon[(i + 1) % n];
        if y1 <= py {
            if y2 > py && cross_2d(x2 - x1, y2 - y1, px - x1, py - y1) > 0.0 {
                winding += 1;
            }
        } else if y2 <= py && cross_2d(x2 - x1, y2 - y1, px - x1, py - y1) < 0.0 {
            winding -= 1;
        }
    }
    winding != 0
}

/// Signed area of a 2D polygon via the shoelace formula.
fn polygon_area(polygon: &[(f64, f64)]) -> f64 {
    if polygon.len() < 3 {
        return 0.0;
    }
    let mut area = 0.0;
    let n = polygon.len();
    for i in 0..n {
        let (x1, y1) = polygon[i];
        let (x2, y2) = polygon[(i + 1) % n];
        area += x1 * y2 - x2 * y1;
    }
    (area / 2.0).abs()
}

/// 2D cross product of vectors (ax, ay) and (bx, by).
fn cross_2d(ax: f64, ay: f64, bx: f64, by: f64) -> f64 {
    ax * by - ay * bx
}

/// Compute a bounded parameter range for a face's surface, tight enough for
/// tessellation.
///
/// Uses edge curve sample points (not just vertex positions) to determine
/// the face's UV extent. This handles self-loop edges (e.g., circular caps)
/// where the vertex count is 1 but the curve covers a full circle.
/// Uses 5% padding to ensure the UV grid covers the face plus a small
/// margin. The trimming step then clips triangles outside the boundary.
fn bounded_parameter_range(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
    let positions = brep.face_boundary_sample_points(face_handle, 16);
    popkern_core::geometry::uv_projection::tight_bounded_parameter_range(&positions, surface)
}

/// Tessellate all faces of a BRep into a combined [`TriMesh`].
pub(crate) fn tessellate_brep(
    brep: &BRep,
    options: &TessellationOptions,
) -> Result<TriMesh, MeshError> {
    let mut face_meshes = Vec::new();
    for face_handle in brep.faces.keys() {
        if brep.faces[face_handle].surface.is_none() {
            continue;
        }
        face_meshes.push(tessellate_face(brep, face_handle, options)?);
    }
    Ok(TriMesh::from_face_meshes(face_meshes))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tessellate;
    use popkern_core::tolerance::Tolerance;

    #[test]
    fn tessellate_unit_cube_correct_counts() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let mesh = tessellate(&cube, &TessellationOptions::default()).unwrap();
        // With min_subdivisions=2, each planar face gets a 2x2 grid of quads,
        // each split into 2 triangles: 6 faces * 4 quads * 2 = 48 triangles.
        assert_eq!(mesh.triangles.len(), 48);
        assert_eq!(mesh.face_map.len(), 48);
        assert_eq!(mesh.face_meshes.len(), 6);
    }

    #[test]
    fn tessellate_unit_cube_min_subdiv_1() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let opts = TessellationOptions {
            min_subdivisions: 1,
            ..Default::default()
        };
        let mesh = tessellate(&cube, &opts).unwrap();
        // With min_subdivisions=1, each planar face gets 1 quad = 2 triangles.
        // 6 faces * 2 = 12 triangles.
        assert_eq!(mesh.triangles.len(), 12);
        assert_eq!(mesh.face_map.len(), 12);
        assert_eq!(mesh.face_meshes.len(), 6);
    }

    #[test]
    fn tessellate_sphere_adaptive() {
        let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
        let mesh = tessellate(&sphere, &TessellationOptions::default()).unwrap();
        assert!(mesh.triangles.len() > 12);
    }

    #[test]
    fn tessellate_cube_normals_axis_aligned() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let mesh = tessellate(&cube, &TessellationOptions::default()).unwrap();
        for v in &mesh.vertices {
            let abs_sum = v.normal.x.abs() + v.normal.y.abs() + v.normal.z.abs();
            assert!(
                (abs_sum - 1.0).abs() < 0.1,
                "Normal should be roughly unit axis-aligned, got ({}, {}, {})",
                v.normal.x,
                v.normal.y,
                v.normal.z
            );
        }
    }

    #[test]
    fn tessellate_sphere_vertices_on_surface() {
        let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
        let mesh = tessellate(&sphere, &TessellationOptions::default()).unwrap();
        for v in &mesh.vertices {
            let dist = (v.position.x.powi(2) + v.position.y.powi(2) + v.position.z.powi(2)).sqrt();
            assert!((dist - 1.0).abs() < 1e-6);
        }
    }

    #[test]
    fn tighter_chord_height_more_triangles() {
        let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
        let coarse = tessellate(
            &sphere,
            &TessellationOptions {
                chord_height: 0.1,
                ..Default::default()
            },
        )
        .unwrap();
        let fine = tessellate(
            &sphere,
            &TessellationOptions {
                chord_height: 0.001,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(fine.triangles.len() > coarse.triangles.len());
    }

    #[test]
    fn all_indices_in_bounds() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let mesh = tessellate(&cube, &TessellationOptions::default()).unwrap();
        for tri in &mesh.triangles {
            for &idx in tri {
                assert!(idx < mesh.vertices.len());
            }
        }
    }
}
