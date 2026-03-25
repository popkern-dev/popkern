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

    Ok(FaceMesh {
        face: face_handle,
        vertices,
        triangles,
    })
}

/// Compute a bounded parameter range for a face's surface, tight enough for
/// tessellation.
///
/// The shared `popkern_core::geometry::uv_projection::bounded_parameter_range`
/// uses generous padding (50%, min 1.0) for boolean intersection solvers. For
/// tessellation we need tight bounds — vertices are placed directly on the UV
/// grid, so excess padding produces mesh geometry outside the face boundary.
///
/// Uses 5 % padding (minimum 0.05) on projected vertex extents.
fn bounded_parameter_range(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
    use popkern_core::geometry::uv_projection::project_point_to_uv;

    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();

    // If both ranges are already finite, return them as-is.
    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
        return ((u_min, u_max), (v_min, v_max));
    }

    let positions = brep.face_outer_loop_positions(face_handle);
    if positions.is_empty() {
        return ((u_min, u_max), (v_min, v_max));
    }

    // Project each vertex to UV and accumulate the bounding box.
    let mut uv_lo = (f64::INFINITY, f64::INFINITY);
    let mut uv_hi = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut any_projected = false;

    for pos in &positions {
        if let Some((u, v)) = project_point_to_uv(surface, pos) {
            uv_lo.0 = uv_lo.0.min(u);
            uv_lo.1 = uv_lo.1.min(v);
            uv_hi.0 = uv_hi.0.max(u);
            uv_hi.1 = uv_hi.1.max(v);
            any_projected = true;
        }
    }

    if !any_projected {
        return ((u_min, u_max), (v_min, v_max));
    }

    // Tight padding for tessellation: 5 % of the span, minimum 0.05.
    let u_pad = ((uv_hi.0 - uv_lo.0) * 0.05).max(0.05);
    let v_pad = ((uv_hi.1 - uv_lo.1) * 0.05).max(0.05);

    // Only override dimensions that were infinite.
    let u_range = if u_min.is_finite() && u_max.is_finite() {
        (u_min, u_max)
    } else {
        (uv_lo.0 - u_pad, uv_hi.0 + u_pad)
    };
    let v_range = if v_min.is_finite() && v_max.is_finite() {
        (v_min, v_max)
    } else {
        (uv_lo.1 - v_pad, uv_hi.1 + v_pad)
    };

    (u_range, v_range)
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
