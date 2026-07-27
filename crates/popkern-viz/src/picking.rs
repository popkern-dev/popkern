//! CPU-based face picking via ray-triangle intersection.
//!
//! Casts a ray from the camera through the click position and tests it
//! against every triangle in the mesh. The closest hit is mapped back to
//! its originating BRep [`FaceHandle`] via the mesh's face map.

use glam::Vec3;
use popkern_core::topology::FaceHandle;
use popkern_mesh::TriMesh;

use crate::camera::Camera;

/// Result of a pick operation.
#[derive(Debug, Clone)]
pub struct PickHit {
    /// The face that was hit.
    pub face: FaceHandle,
    /// The world-space position of the intersection.
    pub position: Vec3,
    /// Distance from the ray origin to the hit.
    pub distance: f32,
    /// Index of the triangle within the mesh that was hit.
    pub triangle_index: usize,
}

/// Encode a face index as an RGBA color using the lower 24 bits.
///
/// Used for GPU-based pick-buffer rendering (future use). Each face gets
/// a unique color that can be read back from a framebuffer pixel.
pub fn encode_face_id(id: u32) -> [f32; 4] {
    let r = ((id >> 16) & 0xFF) as f32 / 255.0;
    let g = ((id >> 8) & 0xFF) as f32 / 255.0;
    let b = (id & 0xFF) as f32 / 255.0;
    [r, g, b, 1.0]
}

/// Decode a framebuffer pixel back to a face index.
///
/// Returns `None` if the pixel has alpha == 0 (background).
pub fn decode_face_id(pixel: [u8; 4]) -> Option<u32> {
    if pixel[3] == 0 {
        return None;
    }
    Some(((pixel[0] as u32) << 16) | ((pixel[1] as u32) << 8) | pixel[2] as u32)
}

/// Cast a ray from the camera through a click position and find the closest
/// face intersection.
///
/// # Arguments
///
/// * `camera` - The camera defining the view and projection.
/// * `click_ndc` - Click position in normalized device coordinates, where
///   (0, 0) is the top-left and (1, 1) is the bottom-right of the viewport.
/// * `mesh` - The tessellated triangle mesh to test against.
///
/// Returns the closest [`PickHit`] or `None` if no triangle was hit.
pub fn pick_face(camera: &Camera, click_ndc: (f32, f32), mesh: &TriMesh) -> Option<PickHit> {
    let (origin, direction) = unproject_ray(camera, click_ndc);

    let mut closest: Option<PickHit> = None;

    for (tri_idx, tri) in mesh.triangles.iter().enumerate() {
        let v0 = mesh_vertex_pos(mesh, tri[0]);
        let v1 = mesh_vertex_pos(mesh, tri[1]);
        let v2 = mesh_vertex_pos(mesh, tri[2]);

        if let Some(t) = ray_triangle_intersect(origin, direction, v0, v1, v2)
            && t > 0.0
            && closest.as_ref().is_none_or(|prev| t < prev.distance)
        {
            let face = mesh.face_map.get(tri_idx).copied().unwrap_or_default();
            closest = Some(PickHit {
                face,
                position: origin + direction * t,
                distance: t,
                triangle_index: tri_idx,
            });
        }
    }

    closest
}

/// Extract a vertex position from the mesh as a [`glam::Vec3`].
fn mesh_vertex_pos(mesh: &TriMesh, idx: usize) -> Vec3 {
    let v = &mesh.vertices[idx];
    Vec3::new(
        v.position.x as f32,
        v.position.y as f32,
        v.position.z as f32,
    )
}

/// Unproject a normalized click position into a world-space ray.
///
/// Returns `(origin, direction)` where `direction` is normalized.
fn unproject_ray(camera: &Camera, click_ndc: (f32, f32)) -> (Vec3, Vec3) {
    // Convert from [0,1] viewport coords to [-1,1] clip space.
    // click_ndc: (0,0) = top-left, (1,1) = bottom-right
    let clip_x = click_ndc.0 * 2.0 - 1.0;
    let clip_y = 1.0 - click_ndc.1 * 2.0; // flip Y: top is +1 in clip space

    let inv_vp = camera.view_projection().inverse();

    // Near plane point in clip space → world
    let near_clip = glam::Vec4::new(clip_x, clip_y, 0.0, 1.0);
    let near_world = inv_vp * near_clip;
    let near_world = Vec3::new(
        near_world.x / near_world.w,
        near_world.y / near_world.w,
        near_world.z / near_world.w,
    );

    // Far plane point in clip space → world
    let far_clip = glam::Vec4::new(clip_x, clip_y, 1.0, 1.0);
    let far_world = inv_vp * far_clip;
    let far_world = Vec3::new(
        far_world.x / far_world.w,
        far_world.y / far_world.w,
        far_world.z / far_world.w,
    );

    let direction = (far_world - near_world).normalize();
    (near_world, direction)
}

/// Moller-Trumbore ray-triangle intersection.
///
/// Returns the distance `t` along the ray if the ray hits the triangle,
/// or `None` if it misses.
fn ray_triangle_intersect(origin: Vec3, dir: Vec3, v0: Vec3, v1: Vec3, v2: Vec3) -> Option<f32> {
    const EPSILON: f32 = 1e-8;

    let edge1 = v1 - v0;
    let edge2 = v2 - v0;
    let h = dir.cross(edge2);
    let a = edge1.dot(h);

    if a.abs() < EPSILON {
        return None; // Ray is parallel to the triangle
    }

    let f = 1.0 / a;
    let s = origin - v0;
    let u = f * s.dot(h);

    if !(0.0..=1.0).contains(&u) {
        return None;
    }

    let q = s.cross(edge1);
    let v = f * dir.dot(q);

    if v < 0.0 || u + v > 1.0 {
        return None;
    }

    let t = f * edge2.dot(q);
    if t > EPSILON { Some(t) } else { None }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_decode_roundtrip() {
        for id in [0, 1, 255, 256, 65535, 0xFF_FFFF] {
            let color = encode_face_id(id);
            let pixel = [
                (color[0] * 255.0).round() as u8,
                (color[1] * 255.0).round() as u8,
                (color[2] * 255.0).round() as u8,
                (color[3] * 255.0).round() as u8,
            ];
            let decoded = decode_face_id(pixel);
            assert_eq!(decoded, Some(id), "roundtrip failed for id={id}");
        }
    }

    #[test]
    fn decode_background_returns_none() {
        assert_eq!(decode_face_id([0, 0, 0, 0]), None);
    }

    #[test]
    fn encode_face_id_zero() {
        let c = encode_face_id(0);
        assert!((c[0]).abs() < f32::EPSILON);
        assert!((c[1]).abs() < f32::EPSILON);
        assert!((c[2]).abs() < f32::EPSILON);
        assert!((c[3] - 1.0).abs() < f32::EPSILON);
    }

    #[test]
    fn encode_face_id_max_24bit() {
        let c = encode_face_id(0xFF_FFFF);
        assert!((c[0] - 1.0).abs() < 0.01);
        assert!((c[1] - 1.0).abs() < 0.01);
        assert!((c[2] - 1.0).abs() < 0.01);
    }

    #[test]
    fn ray_hits_triangle() {
        let v0 = Vec3::new(-1.0, -1.0, 0.0);
        let v1 = Vec3::new(1.0, -1.0, 0.0);
        let v2 = Vec3::new(0.0, 1.0, 0.0);

        let origin = Vec3::new(0.0, 0.0, 5.0);
        let dir = Vec3::new(0.0, 0.0, -1.0);

        let t = ray_triangle_intersect(origin, dir, v0, v1, v2);
        assert!(t.is_some());
        assert!((t.unwrap() - 5.0).abs() < 1e-5);
    }

    #[test]
    fn ray_misses_triangle() {
        let v0 = Vec3::new(-1.0, -1.0, 0.0);
        let v1 = Vec3::new(1.0, -1.0, 0.0);
        let v2 = Vec3::new(0.0, 1.0, 0.0);

        // Ray pointing away
        let origin = Vec3::new(0.0, 0.0, 5.0);
        let dir = Vec3::new(0.0, 0.0, 1.0);

        let t = ray_triangle_intersect(origin, dir, v0, v1, v2);
        assert!(t.is_none());
    }

    #[test]
    fn ray_misses_outside_triangle() {
        let v0 = Vec3::new(-1.0, -1.0, 0.0);
        let v1 = Vec3::new(1.0, -1.0, 0.0);
        let v2 = Vec3::new(0.0, 1.0, 0.0);

        let origin = Vec3::new(5.0, 5.0, 5.0);
        let dir = Vec3::new(0.0, 0.0, -1.0);

        let t = ray_triangle_intersect(origin, dir, v0, v1, v2);
        assert!(t.is_none());
    }

    #[test]
    fn ray_parallel_to_triangle() {
        let v0 = Vec3::new(-1.0, -1.0, 0.0);
        let v1 = Vec3::new(1.0, -1.0, 0.0);
        let v2 = Vec3::new(0.0, 1.0, 0.0);

        let origin = Vec3::new(0.0, 0.0, 0.0);
        let dir = Vec3::new(1.0, 0.0, 0.0);

        let t = ray_triangle_intersect(origin, dir, v0, v1, v2);
        assert!(t.is_none());
    }

    #[test]
    fn unproject_center_points_forward() {
        let camera = Camera::default_iso();
        let (origin, direction) = unproject_ray(&camera, (0.5, 0.5));

        // The ray should point roughly from eye toward target
        let expected_dir = (camera.target - camera.eye).normalize();
        let dot = direction.dot(expected_dir);
        assert!(
            dot > 0.9,
            "center ray should point toward target, dot={dot}"
        );

        // Origin should be near the camera eye (on the near plane)
        let dist_to_eye = (origin - camera.eye).length();
        assert!(
            dist_to_eye < 1.0,
            "ray origin should be near camera eye, dist={dist_to_eye}"
        );
    }

    #[test]
    fn pick_face_no_mesh_triangles() {
        let camera = Camera::default_iso();
        let mesh = TriMesh::default();
        let hit = pick_face(&camera, (0.5, 0.5), &mesh);
        assert!(hit.is_none());
    }

    #[test]
    fn pick_face_hits_centered_triangle() {
        use popkern_core::geometry::{Point3, Vec3 as CoreVec3};
        use popkern_mesh::MeshVertex;

        let camera = Camera {
            eye: glam::Vec3::new(0.0, 0.0, 5.0),
            target: glam::Vec3::ZERO,
            up: glam::Vec3::Y,
            fov_y: std::f32::consts::FRAC_PI_4,
            aspect: 1.0,
            near: 0.1,
            far: 100.0,
        };

        let mesh = TriMesh {
            vertices: vec![
                MeshVertex {
                    position: Point3::new(-1.0, -1.0, 0.0),
                    normal: CoreVec3::new(0.0, 0.0, 1.0),
                    uv: (0.0, 0.0),
                },
                MeshVertex {
                    position: Point3::new(1.0, -1.0, 0.0),
                    normal: CoreVec3::new(0.0, 0.0, 1.0),
                    uv: (1.0, 0.0),
                },
                MeshVertex {
                    position: Point3::new(0.0, 1.0, 0.0),
                    normal: CoreVec3::new(0.0, 0.0, 1.0),
                    uv: (0.5, 1.0),
                },
            ],
            triangles: vec![[0, 1, 2]],
            face_map: vec![FaceHandle::default()],
            face_meshes: vec![],
        };

        // Click center of viewport
        let hit = pick_face(&camera, (0.5, 0.5), &mesh);
        assert!(hit.is_some(), "should hit centered triangle");
        let hit = hit.unwrap();
        assert_eq!(hit.triangle_index, 0);
        assert!(hit.distance > 0.0);
    }

    #[test]
    fn pick_face_misses_off_center() {
        use popkern_core::geometry::{Point3, Vec3 as CoreVec3};
        use popkern_mesh::MeshVertex;

        let camera = Camera {
            eye: glam::Vec3::new(0.0, 0.0, 5.0),
            target: glam::Vec3::ZERO,
            up: glam::Vec3::Y,
            fov_y: std::f32::consts::FRAC_PI_4,
            aspect: 1.0,
            near: 0.1,
            far: 100.0,
        };

        // Small triangle at origin
        let mesh = TriMesh {
            vertices: vec![
                MeshVertex {
                    position: Point3::new(-0.01, -0.01, 0.0),
                    normal: CoreVec3::new(0.0, 0.0, 1.0),
                    uv: (0.0, 0.0),
                },
                MeshVertex {
                    position: Point3::new(0.01, -0.01, 0.0),
                    normal: CoreVec3::new(0.0, 0.0, 1.0),
                    uv: (1.0, 0.0),
                },
                MeshVertex {
                    position: Point3::new(0.0, 0.01, 0.0),
                    normal: CoreVec3::new(0.0, 0.0, 1.0),
                    uv: (0.5, 1.0),
                },
            ],
            triangles: vec![[0, 1, 2]],
            face_map: vec![FaceHandle::default()],
            face_meshes: vec![],
        };

        // Click far corner — should miss the tiny triangle
        let hit = pick_face(&camera, (0.0, 0.0), &mesh);
        assert!(hit.is_none(), "should miss small triangle at corner click");
    }

    #[test]
    fn pick_face_returns_closest() {
        use popkern_core::geometry::{Point3, Vec3 as CoreVec3};
        use popkern_mesh::MeshVertex;

        let camera = Camera {
            eye: glam::Vec3::new(0.0, 0.0, 10.0),
            target: glam::Vec3::ZERO,
            up: glam::Vec3::Y,
            fov_y: std::f32::consts::FRAC_PI_4,
            aspect: 1.0,
            near: 0.1,
            far: 100.0,
        };

        let make_tri = |z: f64| -> Vec<MeshVertex> {
            vec![
                MeshVertex {
                    position: Point3::new(-2.0, -2.0, z),
                    normal: CoreVec3::new(0.0, 0.0, 1.0),
                    uv: (0.0, 0.0),
                },
                MeshVertex {
                    position: Point3::new(2.0, -2.0, z),
                    normal: CoreVec3::new(0.0, 0.0, 1.0),
                    uv: (1.0, 0.0),
                },
                MeshVertex {
                    position: Point3::new(0.0, 2.0, z),
                    normal: CoreVec3::new(0.0, 0.0, 1.0),
                    uv: (0.5, 1.0),
                },
            ]
        };

        let mut verts = make_tri(2.0); // closer
        verts.extend(make_tri(-2.0)); // farther

        let mesh = TriMesh {
            vertices: verts,
            triangles: vec![[0, 1, 2], [3, 4, 5]],
            face_map: vec![FaceHandle::default(), FaceHandle::default()],
            face_meshes: vec![],
            zipped_t_junctions: 0,
        };

        let hit = pick_face(&camera, (0.5, 0.5), &mesh);
        assert!(hit.is_some());
        // Should pick triangle 0 (at z=2, closer to camera at z=10)
        assert_eq!(hit.unwrap().triangle_index, 0);
    }
}
