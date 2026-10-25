//! Face sense, face area and solid orientation from the welded mesh.
//!
//! Measured at the pinned popkern (popkern ledger PK-772): planar faces
//! store an inward plane normal with `orientation = true`, while the
//! tessellation is outward. The outward sense of a plane is therefore
//! taken from its own triangles, which also checks that they lie on the
//! plane. Face areas come from the same triangles (PK-771: `face_area`
//! drops inner loops on planar faces), so they are bounded by the chord
//! tolerance like every other mesh quantity.

use popkern_analysis::{FaceId, FaceInfo, SolidId, Surface};

use crate::adapter::AdapterError;
use crate::geometry::{cross, dot, norm, sub};

/// Relative tolerance on |n·Σ n_tri| / 2 = area for a planar face.
pub const PLANE_AREA_TOL: f64 = 1.0e-9;

/// Sum of triangle cross products (twice the vector area) and total
/// area of each face in `faces`, from the solid's welded mesh.
fn face_measures(
    vertices: &[[f64; 3]],
    triangles: &[[u32; 3]],
    triangle_face: &[FaceId],
    faces: &[FaceInfo],
) -> Vec<([f64; 3], f64)> {
    let first = faces.first().map_or(0, |f| f.id.0 as usize);
    let mut sums = vec![([0.0; 3], 0.0); faces.len()];
    for (t, face) in triangles.iter().zip(triangle_face) {
        let [a, b, c] = t.map(|v| vertices[v as usize]);
        let x = cross(sub(b, a), sub(c, a));
        if let Some((vector, area)) = (face.0 as usize)
            .checked_sub(first)
            .and_then(|i| sums.get_mut(i))
        {
            for k in 0..3 {
                vector[k] += x[k];
            }
            *area += 0.5 * norm(x);
        }
    }
    sums
}

/// Signed volume enclosed by a closed mesh (positive when outward).
pub fn signed_volume(vertices: &[[f64; 3]], triangles: &[[u32; 3]]) -> f64 {
    let six: f64 = triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|v| vertices[v as usize]);
            dot(a, cross(b, c))
        })
        .sum();
    six / 6.0
}

/// Sets each face's area from the mesh and orients each plane outward.
/// `faces` are this solid's faces, with consecutive ids.
pub fn orient_faces(
    solid: SolidId,
    vertices: &[[f64; 3]],
    triangles: &[[u32; 3]],
    triangle_face: &[FaceId],
    faces: &mut [FaceInfo],
) -> Result<(), AdapterError> {
    if signed_volume(vertices, triangles) <= 0.0 {
        return Err(AdapterError::InsideOut { solid });
    }
    let measures = face_measures(vertices, triangles, triangle_face, faces);
    for (face, (vector, area)) in faces.iter_mut().zip(measures) {
        if area <= 0.0 {
            return Err(AdapterError::UntessellatedFace { face: face.id });
        }
        face.area = area;
        if let Surface::Plane { normal, offset } = &mut face.surface {
            let along = dot(*normal, vector);
            if ((along.abs() / 2.0) - area).abs() > PLANE_AREA_TOL * area {
                return Err(AdapterError::OffPlane {
                    face: face.id,
                    projected: along.abs() / 2.0,
                    area,
                });
            }
            if along < 0.0 {
                *normal = normal.map(|c| -c);
                *offset = -*offset;
            }
        }
    }
    Ok(())
}
