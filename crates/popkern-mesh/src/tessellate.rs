//! Core tessellation: per-face meshing and top-level BRep tessellation.

use std::collections::HashMap;

use popkern_core::geometry::{SurfaceEval, SurfaceGeometry, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use crate::adaptive::adaptive_subdivide;
use crate::mesh::{FaceMesh, MeshVertex, TriMesh};
use crate::{MeshError, TessellationOptions};

/// Tessellate a single BRep face into a triangle mesh.
pub fn tessellate_face(
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
                if len > brep.tolerance.numerical_floor() {
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
            // All construction pipelines (constructors, sweep, boolean,
            // fillet, shell, draft) now generate pcurves via
            // recompute_all_pcurves. The full clipping pipeline runs
            // unconditionally. The 3D validation in apply_trim handles
            // any residual UV imprecision.
            let (trimmed_verts, trimmed_tris) = crate::trim_clip::apply_trim(
                &vertices,
                &triangles,
                &trim_loops,
                surface,
                orientation,
            );
            vertices = trimmed_verts;
            triangles = trimmed_tris;
        }
    }

    // When face orientation is reversed, flip triangle winding to match.
    // The adaptive subdivision always produces CCW winding from the surface
    // geometry's normal direction. If orientation=false, the effective face
    // normal is flipped, so triangles must be rewound to remain front-facing
    // under CCW back-face culling.
    if !orientation {
        for tri in &mut triangles {
            tri.swap(1, 2);
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
///
/// On periodic surfaces (cylinder, torus), if the raw UV range spans
/// more than half the period, the boundary likely crosses the seam.
/// In that case, shift all UV samples by half-period and recompute
/// min/max to find the shorter arc that doesn't straddle the seam.
fn bounded_parameter_range(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
    let positions = brep.face_boundary_sample_points(face_handle, 16);
    let mut result =
        popkern_core::geometry::uv_projection::tight_bounded_parameter_range(&positions, surface);

    // If pcurve-based bounds are available, intersect with the fallback.
    if let Some(pc_range) = pcurve_parameter_range(brep, face_handle) {
        let u = (pc_range.0.0.max(result.0.0), pc_range.0.1.min(result.0.1));
        let v = (pc_range.1.0.max(result.1.0), pc_range.1.1.min(result.1.1));
        if u.0 < u.1 && v.0 < v.1 {
            result = (u, v);
        }
    }

    result
}

/// Compute UV parameter range from the outer loop's pcurves.
///
/// Returns `Some((u_range, v_range))` if every half-edge in the outer loop
/// has a pcurve; `None` otherwise, signaling the caller to fall back to
/// 3D projection.
///
/// On periodic surfaces, UV values are unwrapped relative to the first
/// sample to prevent seam crossings from inflating the range to the
/// full period.
fn pcurve_parameter_range(
    brep: &BRep,
    face_handle: FaceHandle,
) -> Option<((f64, f64), (f64, f64))> {
    let face = brep.faces.get(face_handle)?;
    let surface = face.surface.as_ref()?;
    let (period_u, period_v) = surface.geometry.periodicity();

    let loop_ = brep.loops.get(face.outer_loop)?;
    let start_he = loop_.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return None;
    }

    let samples = 16_usize;
    let mut all_u = Vec::new();
    let mut all_v = Vec::new();

    let mut current = start_he;
    loop {
        let he = &brep.half_edges[current];
        let pcurve = he.pcurve.as_ref()?;
        let (t0, t1) = pcurve.parameter_range();
        for i in 0..=samples {
            let t = t0 + (t1 - t0) * (i as f64 / samples as f64);
            let pt = pcurve.point_at(t);
            all_u.push(pt.x);
            all_v.push(pt.y);
        }
        current = he.next;
        if current == start_he {
            break;
        }
    }

    if all_u.is_empty() {
        return None;
    }

    // Unwrap periodic values relative to the first sample.
    unwrap_periodic(&mut all_u, period_u);
    unwrap_periodic(&mut all_v, period_v);

    let u_min = all_u.iter().cloned().fold(f64::INFINITY, f64::min);
    let u_max = all_u.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let v_min = all_v.iter().cloned().fold(f64::INFINITY, f64::min);
    let v_max = all_v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
        Some(((u_min, u_max), (v_min, v_max)))
    } else {
        None
    }
}

/// Unwrap a sequence of values on a periodic domain so consecutive
/// values don't jump across the seam. Adjusts values in-place by
/// adding/subtracting the period to stay within half-period of the
/// previous value.
fn unwrap_periodic(values: &mut [f64], period: Option<f64>) {
    let Some(p) = period else { return };
    if p <= 0.0 || values.len() < 2 {
        return;
    }
    let half = p * 0.5;
    for i in 1..values.len() {
        let diff = values[i] - values[i - 1];
        if diff > half {
            values[i] -= p;
        } else if diff < -half {
            values[i] += p;
        }
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
        // With min_subdivisions=4, each planar face gets a 4x4 grid of quads,
        // each split into 2 triangles: 6 faces * 16 quads * 2 = 192 triangles.
        assert_eq!(mesh.triangles.len(), 192);
        assert_eq!(mesh.face_map.len(), 192);
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

#[cfg(test)]
#[path = "tessellate_edge_tests.rs"]
mod tessellate_edge_tests;
