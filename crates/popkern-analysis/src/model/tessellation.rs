//! Faces whose tessellation does not match the face.
//!
//! popkern's tessellation of boss_plate's sphere corners covers most of
//! each sphere instead of its trimmed corner: 0.070 mm² of triangles
//! for a 0.010 mm² face, stitched in closed and consistently oriented,
//! so no edge census sees it (PK-787). The kernel's closed-form area
//! of a curved analytic face reveals it: taken over the face's
//! parameter box, it is exact for a rectangular domain and too large
//! otherwise (1372 against 31 mm² on threaded_stud's helically trimmed
//! cylinder), so it bounds the true area from above, and only a
//! tessellation larger than it is wrong. A face whose tessellation
//! folds against its surface (`FaceInfo::folded_area`, from the
//! kernel's exact normals) is mismatched too. Such faces make results
//! nearby unreliable, and wall readings ignore them. A tessellation
//! whose volume falls short of the kernel's (`volume_warning`) carries
//! chord error into mass and stiffness.

use super::{FaceInfo, SolidMesh, Surface};

/// A curved face's tessellated area may exceed the kernel's bound by
/// this factor before the face is mismatched. Chord error stays within
/// a few percent; a wrongly trimmed face is off by far more (7x on
/// boss_plate; no other part of the 30-part gallery exceeds 1.5).
pub const AREA_MISMATCH: f64 = 1.5;

/// A face is mismatched when more than this fraction of its
/// tessellated area is folded (none folds over the 30-part gallery).
pub const FOLD_FRACTION: f64 = 0.01;

fn triangle_area(mesh: &SolidMesh, t: [u32; 3]) -> f64 {
    let [a, b, c] = t.map(|i| mesh.vertices[i as usize]);
    let (u, v) = (
        [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
        [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
    );
    let n = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    0.5 * (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
}

/// Faces of `mesh` whose tessellated area exceeds the kernel's bound
/// (curved analytic faces) by more than [`AREA_MISMATCH`], or more
/// than [`FOLD_FRACTION`] of which is folded: (face id, tessellated
/// area, folded area), ascending by id.
pub fn mismatched_faces(mesh: &SolidMesh, faces: &[FaceInfo]) -> Vec<(u32, f64, f64)> {
    let mut tessellated: std::collections::BTreeMap<u32, f64> = std::collections::BTreeMap::new();
    for (t, f) in mesh.triangles.iter().zip(&mesh.triangle_face) {
        *tessellated.entry(f.0).or_default() += triangle_area(mesh, *t);
    }
    tessellated
        .into_iter()
        .filter_map(|(id, area)| {
            let face = faces.iter().find(|f| f.id.0 == id)?;
            let curved = matches!(
                face.surface,
                Surface::Cylinder { .. }
                    | Surface::Cone { .. }
                    | Surface::Sphere { .. }
                    | Surface::Torus { .. }
            );
            let folded = face.folded_area.unwrap_or(0.0);
            let beyond = face
                .exact_area
                .is_some_and(|exact| curved && exact > 0.0 && area > AREA_MISMATCH * exact);
            (beyond || folded > FOLD_FRACTION * area).then_some((id, area, folded))
        })
        .collect()
}

/// The `TessellationMismatch` warning for `mesh`, if any face is
/// mismatched, with a flag per triangle on such faces.
pub fn mismatch_warning(
    mesh: &SolidMesh,
    faces: &[FaceInfo],
) -> (Vec<bool>, Option<crate::report::Warning>) {
    let bad = mismatched_faces(mesh, faces);
    let flags = mesh
        .triangle_face
        .iter()
        .map(|f| bad.iter().any(|b| b.0 == f.0))
        .collect();
    let warning = (!bad.is_empty()).then(|| crate::report::Warning::TessellationMismatch {
        solid: mesh.id.0,
        faces: bad.iter().map(|b| b.0).collect(),
        tessellated_area: bad.iter().map(|b| b.1).sum(),
        folded_area: bad.iter().map(|b| b.2).sum(),
    });
    (flags, warning)
}

/// A solid's tessellated volume may differ from the kernel's volume by
/// this fraction before the solve warns. The moments match the
/// tessellation (X5), so mass and stiffness carry its chord error.
/// Measured per solid over the gallery at the display tessellation
/// (PH-094): centrifugal_impeller's thin curved blades 1.35% short,
/// ball_bearing's balls 0.78% over, piano_hinge's 2 mm pin 0.60% short,
/// every other solid within 0.18%. The kernel's volume of a curved solid is itself
/// a fine tessellation, within 0.6% of OpenCASCADE (PK-807): the
/// impeller is 1.44% short of OpenCASCADE's, so the warning holds, but
/// deficits near the threshold are only as good as that reference.
pub const VOLUME_DEFICIT: f64 = 5.0e-3;

/// The `TessellationVolume` warning for `mesh` when its enclosed volume
/// and the kernel's volume differ by more than [`VOLUME_DEFICIT`].
pub fn volume_warning(mesh: &SolidMesh) -> Option<crate::report::Warning> {
    let kernel = mesh.exact_volume?;
    let tessellated = mesh.enclosed_volume();
    (kernel > 0.0 && (tessellated - kernel).abs() > VOLUME_DEFICIT * kernel).then_some(
        crate::report::Warning::TessellationVolume {
            solid: mesh.id.0,
            tessellated,
            kernel,
        },
    )
}

#[cfg(test)]
#[path = "tessellation_tests.rs"]
mod tests;
