//! Trim a sheet body to the portion lying on one side of a solid (§6.x).
//!
//! Used by [`crate::sheet::boolean_with_sheets`] to implement the two
//! operations that cut a sheet against a solid: `SheetIntersectSolid`
//! keeps the inside fragments, while `SolidUnionSheet` keeps the
//! outside fragments as an open-shell appendage.

use popkern_core::topology::BRep;

/// Which half of the sheet to keep when trimming against a solid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum KeepSide {
    /// Keep the portion of the sheet inside the solid.
    Inside,
    /// Keep the portion of the sheet outside the solid.
    Outside,
}

/// Trim a sheet body to the portion on the requested side of `solid`.
///
/// Walks every face of the (already-imprinted) sheet, samples a
/// centroid point, and tests it via point-in-solid against the
/// solid. Faces whose classification does not match `keep` are
/// removed from the sheet via slotmap delete; orphaned edges and
/// vertices are pruned afterwards.
///
/// The result is an open-shell BRep — possibly with multiple
/// disconnected components — containing only the surviving fragments.
/// This is the §6.x boundary-relative classification path: instead
/// of asking "what's the volume of the result", we ask "which
/// pieces of the sheet survived the cut".
pub(super) fn trim_sheet(
    mut sheet: BRep,
    solid: &BRep,
    tolerance: &popkern_core::tolerance::Tolerance,
    keep: KeepSide,
) -> BRep {
    use crate::clip::face_vertex_positions;
    use crate::error::PointClassification;
    use crate::point_in_solid;
    use popkern_core::topology::FaceHandle;
    use std::collections::HashSet;

    let solid_aabbs = crate::intersect::compute_face_aabbs(solid);

    // Classify every sheet face as inside / outside / on by sampling
    // its centroid and collect the ones to remove.
    let mut faces_to_remove: HashSet<FaceHandle> = HashSet::new();
    for fh in sheet.faces.keys() {
        let positions = face_vertex_positions(&sheet, fh);
        if positions.is_empty() {
            faces_to_remove.insert(fh);
            continue;
        }
        let mut cx = 0.0;
        let mut cy = 0.0;
        let mut cz = 0.0;
        for p in &positions {
            cx += p.x;
            cy += p.y;
            cz += p.z;
        }
        let n = positions.len() as f64;
        let centroid = popkern_core::geometry::Point3::new(cx / n, cy / n, cz / n);
        let classification = point_in_solid(&centroid, solid, &solid_aabbs, tolerance)
            .unwrap_or(PointClassification::Outside);
        let discard = match keep {
            KeepSide::Inside => matches!(classification, PointClassification::Outside),
            KeepSide::Outside => matches!(
                classification,
                PointClassification::Inside | PointClassification::OnBoundary
            ),
        };
        if discard {
            faces_to_remove.insert(fh);
        }
    }

    super::dedupe::drop_faces_and_prune(&mut sheet, &faces_to_remove);
    sheet
}
