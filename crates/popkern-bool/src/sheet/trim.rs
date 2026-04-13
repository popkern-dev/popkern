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

    let solid_aabbs = crate::intersect::compute_face_aabbs(solid);

    // Step 1: classify every sheet face as inside / outside / on by
    // sampling its centroid.
    let mut faces_to_remove: Vec<FaceHandle> = Vec::new();
    for fh in sheet.faces.keys() {
        let positions = face_vertex_positions(&sheet, fh);
        if positions.is_empty() {
            faces_to_remove.push(fh);
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
            faces_to_remove.push(fh);
        }
    }

    // Step 2: remove outside faces via slotmap delete; clean up the
    // orphaned half-edges / edges / vertices left behind.
    for fh in &faces_to_remove {
        if let Some(face) = sheet.faces.get(*fh) {
            let outer_loop = face.outer_loop;
            let inner_loops = face.inner_loops.clone();
            sheet.loops.remove(outer_loop);
            for il in inner_loops {
                sheet.loops.remove(il);
            }
            sheet.faces.remove(*fh);
        }
    }

    // Step 3: prune orphaned half-edges.
    let orphan_he: Vec<_> = sheet
        .half_edges
        .iter()
        .filter_map(|(heh, he)| {
            if sheet.loops.get(he.loop_).is_none() {
                Some(heh)
            } else {
                None
            }
        })
        .collect();
    for heh in orphan_he {
        sheet.half_edges.remove(heh);
    }

    // Step 4: prune orphaned edges.
    let orphan_edges: Vec<_> = sheet
        .edges
        .iter()
        .filter_map(|(eh, edge)| {
            if sheet.half_edges.get(edge.half_edge).is_none() {
                Some(eh)
            } else {
                None
            }
        })
        .collect();
    for eh in orphan_edges {
        sheet.edges.remove(eh);
    }

    // Step 5: prune orphaned vertices.
    let mut vertex_in_use = std::collections::HashSet::new();
    for (_, he) in &sheet.half_edges {
        vertex_in_use.insert(he.origin);
    }
    let orphan_vertices: Vec<_> = sheet
        .vertices
        .keys()
        .filter(|vh| !vertex_in_use.contains(vh))
        .collect();
    for vh in orphan_vertices {
        sheet.vertices.remove(vh);
    }

    sheet
}
