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

    // Trim: use pcurve-based adaptive polylines for precise boundary clipping.
    let trim_loops =
        crate::trim::build_trim_loops(brep, face_handle, surface, options.chord_height);
    if !trim_loops.is_empty() {
        let poly_area = crate::trim::trim_loop_area(&trim_loops);
        let uv_area = (u_range.1 - u_range.0) * (v_range.1 - v_range.0);
        // Skip trimming for degenerate boundaries (sphere seams, etc.)
        if poly_area > uv_area * 0.01 {
            // Only use the full clipping pipeline when ALL loops were built
            // from pcurves (reliable UV coordinates).  For loops built via
            // 3D→UV projection (unreliable for non-planar surfaces), fall
            // back to centroid-based discard to avoid generating bad geometry.
            let all_from_pcurves = trim_loops.iter().all(|l| l.from_pcurves);
            if all_from_pcurves {
                let (trimmed_verts, trimmed_tris) = crate::trim_clip::apply_trim(
                    &vertices,
                    &triangles,
                    &trim_loops,
                    surface,
                    orientation,
                );
                vertices = trimmed_verts;
                triangles = trimmed_tris;
            } else {
                // Centroid-based discard fallback for projection-based loops.
                // Only apply if the trim polygon area is reasonable (within
                // 2× the UV grid area). Larger areas indicate unreliable
                // UV projections that would create bad geometry.
                if poly_area < uv_area * 2.0 {
                    let outer = trim_loops.iter().find(|l| l.is_outer);
                    if let Some(outer_loop) = outer {
                        triangles.retain(|tri| {
                            let cu = (vertices[tri[0]].uv.0
                                + vertices[tri[1]].uv.0
                                + vertices[tri[2]].uv.0)
                                / 3.0;
                            let cv = (vertices[tri[0]].uv.1
                                + vertices[tri[1]].uv.1
                                + vertices[tri[2]].uv.1)
                                / 3.0;
                            crate::trim::winding_number(cu, cv, &outer_loop.points) != 0
                        });
                    }
                }
            }
        }
    }

    Ok(FaceMesh {
        face: face_handle,
        vertices,
        triangles,
    })
}

/// Compute a bounded parameter range for a face's surface for tessellation.
///
/// Uses the pcurve-based trim loop bounding box when pcurves are available,
/// giving exact UV bounds that match the trim boundary. Falls back to
/// sample-and-project with 5% padding for faces without pcurves, where
/// the padding ensures coverage for degenerate topologies (e.g., sphere
/// seam edges where all vertices project to u=0).
fn bounded_parameter_range(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
    // Compute UV range from 3D boundary sampling (always reliable).
    let positions = brep.face_boundary_sample_points(face_handle, 16);
    let fallback =
        popkern_core::geometry::uv_projection::tight_bounded_parameter_range(&positions, surface);

    // If pcurve-based bounds are available, intersect with the fallback
    // to prevent pcurve UV overflow (e.g., fillet trim curves that extend
    // beyond the face in UV space).
    if let Some(pc_range) = pcurve_parameter_range(brep, face_handle) {
        let u = (
            pc_range.0.0.max(fallback.0.0),
            pc_range.0.1.min(fallback.0.1),
        );
        let v = (
            pc_range.1.0.max(fallback.1.0),
            pc_range.1.1.min(fallback.1.1),
        );
        // Only use the intersection if it's non-degenerate.
        if u.0 < u.1 && v.0 < v.1 {
            return (u, v);
        }
    }

    fallback
}

/// Compute UV parameter range from the outer loop's pcurves.
///
/// Returns `Some((u_range, v_range))` if every half-edge in the outer loop
/// has a pcurve; `None` otherwise, signaling the caller to fall back to
/// 3D projection.
fn pcurve_parameter_range(
    brep: &BRep,
    face_handle: FaceHandle,
) -> Option<((f64, f64), (f64, f64))> {
    use popkern_core::geometry::curve::CurveEval;

    let face = brep.faces.get(face_handle)?;
    let loop_ = brep.loops.get(face.outer_loop)?;
    let start_he = loop_.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return None;
    }

    let mut u_min = f64::INFINITY;
    let mut u_max = f64::NEG_INFINITY;
    let mut v_min = f64::INFINITY;
    let mut v_max = f64::NEG_INFINITY;
    let samples = 16_usize;

    let mut current = start_he;
    loop {
        let he = &brep.half_edges[current];
        let pcurve = he.pcurve.as_ref()?; // Return None if any half-edge lacks a pcurve.

        let (t0, t1) = pcurve.parameter_range();
        for i in 0..=samples {
            let t = t0 + (t1 - t0) * (i as f64 / samples as f64);
            let pt = pcurve.point_at(t);
            u_min = u_min.min(pt.x);
            u_max = u_max.max(pt.x);
            v_min = v_min.min(pt.y);
            v_max = v_max.max(pt.y);
        }

        current = he.next;
        if current == start_he {
            break;
        }
    }

    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
        Some(((u_min, u_max), (v_min, v_max)))
    } else {
        None
    }
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
