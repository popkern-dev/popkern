//! Bi-nonic surface accessors for the KP20 support-cell
//! reconciliation route: read a face's patch as a plain grid
//! (dual-twist patches read back with pair-midpoint corners —
//! boundary rows are plain controls either way), and write
//! plain or dual-twist surfaces onto faces.

use popkern_core::error::CoreError;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{BRep, FaceHandle};

pub(super) fn rectangular_binonic_on_face(
    brep: &BRep,
    face: FaceHandle,
) -> Result<popkern_core::geometry::gregory::RectangularGregoryBiNonic, CoreError> {
    match brep
        .faces
        .get(face)
        .and_then(|face_data| face_data.surface.as_ref())
        .map(|surface| &surface.geometry)
    {
        Some(SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonic(patch))) => {
            Ok((**patch).clone())
        }
        // Dual-twist patches read back as their plain grid (the
        // four corner slots hold the pair midpoints). Boundary
        // rows are plain controls, so trim measurement is exact;
        // a constructive re-author pass recomputes the pairs.
        Some(SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonicDualTwist(patch))) => {
            Ok(popkern_core::geometry::gregory::RectangularGregoryBiNonic { grid: patch.grid })
        }
        _ => Err(CoreError::PreconditionFailed(
            "kp20 support-cell boundary reconciliation requires rectangular bi-nonic surfaces"
                .into(),
        )),
    }
}

pub(super) fn set_rectangular_binonic_on_face(
    brep: &mut BRep,
    face: FaceHandle,
    patch: popkern_core::geometry::gregory::RectangularGregoryBiNonic,
) -> Result<(), CoreError> {
    let Some(face_data) = brep.faces.get_mut(face) else {
        return Err(CoreError::InvalidHandle(
            "kp20 support-cell boundary reconciliation face invalid".into(),
        ));
    };
    face_data.surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonic(Box::new(patch))),
        orientation: true,
    });
    Ok(())
}

/// Write a face's surface as a dual-twist bi-nonic: the plain
/// grid plus per-corner candidate pairs; corners without an
/// authored pair keep coincident candidates (plain behavior).
pub(super) fn set_rectangular_binonic_dual_twist_on_face(
    brep: &mut BRep,
    face: FaceHandle,
    patch: popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    corners: &[Option<(
        popkern_core::geometry::Point3,
        popkern_core::geometry::Point3,
    )>; 4],
) -> Result<(), CoreError> {
    let mut dual =
        popkern_core::geometry::gregory::RectangularGregoryBiNonicDualTwist::from_plain(&patch);
    for (slot, pair) in corners.iter().enumerate() {
        if let Some(pair) = pair {
            dual.inner_cps[slot] = *pair;
        }
    }
    let Some(face_data) = brep.faces.get_mut(face) else {
        return Err(CoreError::InvalidHandle(
            "kp20 support-cell boundary reconciliation face invalid".into(),
        ));
    };
    face_data.surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonicDualTwist(Box::new(
            dual,
        ))),
        orientation: true,
    });
    Ok(())
}
