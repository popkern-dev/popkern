//! Thicken operation: convert a sheet body into a solid.
//!
//! Takes a sheet body (a BRep with open boundary edges) and thickens it
//! by offsetting in both directions and connecting the free boundaries
//! with side wall faces.

use crate::ShellWarning;
use crate::assemble::assemble_shell;
use crate::classify::classify_faces;
use crate::edge_recompute::recompute_edges;
use crate::error::ShellError;
use crate::inner_skin::{apply_offset_geometry, build_inner_skin};
use crate::skin_map::SkinMap;
use crate::surface_offset::compute_offset_surfaces;
use crate::types::{ShellDirection, ShellOptions, ShellResult};
use crate::vertex_recompute::recompute_vertices;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Thicken a sheet body into a solid.
///
/// Offsets the input sheet body by `thickness / 2` in both directions
/// (outward and inward relative to face normals), then connects the
/// free boundary edges with side wall faces.
///
/// The input BRep should be an open sheet body. All faces are kept;
/// side walls are generated at free boundary edges.
///
/// # Errors
///
/// Returns [`ShellError`] if the BRep has no faces or offset computation
/// fails.
pub fn thicken(
    brep: &BRep,
    thickness: f64,
    options: &ShellOptions,
    tolerance: Tolerance,
) -> Result<ShellResult, ShellError> {
    let half = thickness / 2.0;
    let stage_tol = Tolerance::default().with_linear(tolerance.linear / 4.0);

    // Offset outward to create the outer skin.
    let class_out = classify_faces(brep, half, &[], ShellDirection::Outward, &[])?;
    let offset_out = compute_offset_surfaces(brep, &class_out, options.allow_partial, &stage_tol)?;

    let (mut outer, outer_map) = build_inner_skin(brep, &[])?;
    recompute_vertices(
        &mut outer,
        &offset_out.surfaces,
        &outer_map,
        brep,
        &stage_tol,
    )?;
    recompute_edges(
        &mut outer,
        &offset_out.surfaces,
        &outer_map,
        brep,
        &stage_tol,
    )?;
    apply_offset_geometry(&mut outer, &offset_out.surfaces, &outer_map)?;

    // Restore orientation for outer skin (outward, not flipped).
    for (_, face) in &mut outer.faces {
        if let Some(sof) = &mut face.surface {
            sof.orientation = !sof.orientation;
        }
    }

    // Offset inward to create the inner skin.
    let class_in = classify_faces(brep, half, &[], ShellDirection::Inward, &[])?;
    let offset_in = compute_offset_surfaces(brep, &class_in, options.allow_partial, &stage_tol)?;

    let (mut inner, inner_map) = build_inner_skin(brep, &[])?;
    recompute_vertices(
        &mut inner,
        &offset_in.surfaces,
        &inner_map,
        brep,
        &stage_tol,
    )?;
    recompute_edges(
        &mut inner,
        &offset_in.surfaces,
        &inner_map,
        brep,
        &stage_tol,
    )?;
    apply_offset_geometry(&mut inner, &offset_in.surfaces, &inner_map)?;

    // Build a combined skin map from outer to inner.
    let combined_skin_map = build_combined_skin_map(&outer_map, &inner_map);

    // Assemble: outer + inner, no removed faces (side walls at free edges).
    let result = assemble_shell(&outer, &inner, &[], &combined_skin_map, tolerance)?;

    let mut warnings = offset_out.warnings;
    warnings.extend(offset_in.warnings);
    for face in offset_out
        .auto_removed
        .iter()
        .chain(&offset_in.auto_removed)
    {
        warnings.push(ShellWarning::FaceAutoRemoved {
            face: *face,
            reason: "offset failed".into(),
        });
    }

    Ok(ShellResult {
        brep: result,
        warnings,
    })
}

/// Build a combined skin map mapping outer skin handles to inner skin handles.
///
/// Both `outer_map` and `inner_map` map from the original BRep to their
/// respective skins. We need a map from outer skin to inner skin.
fn build_combined_skin_map(outer_map: &SkinMap, inner_map: &SkinMap) -> SkinMap {
    let reverse_outer_vertices = outer_map.reverse_vertex_map();
    let reverse_outer_edges = outer_map.reverse_edge_map();
    let reverse_outer_he = outer_map
        .half_edges
        .iter()
        .map(|(&k, &v)| (v, k))
        .collect::<std::collections::HashMap<_, _>>();
    let reverse_outer_faces = outer_map.reverse_face_map();

    let mut vertices = std::collections::HashMap::new();
    let mut edges = std::collections::HashMap::new();
    let mut half_edges = std::collections::HashMap::new();
    let mut faces = std::collections::HashMap::new();

    for (&outer_skin_vh, &orig_vh) in &reverse_outer_vertices {
        if let Some(&inner_skin_vh) = inner_map.vertices.get(&orig_vh) {
            vertices.insert(outer_skin_vh, inner_skin_vh);
        }
    }
    for (&outer_skin_eh, &orig_eh) in &reverse_outer_edges {
        if let Some(&inner_skin_eh) = inner_map.edges.get(&orig_eh) {
            edges.insert(outer_skin_eh, inner_skin_eh);
        }
    }
    for (&outer_skin_heh, &orig_heh) in &reverse_outer_he {
        if let Some(&inner_skin_heh) = inner_map.half_edges.get(&orig_heh) {
            half_edges.insert(outer_skin_heh, inner_skin_heh);
        }
    }
    for (&outer_skin_fh, &orig_fh) in &reverse_outer_faces {
        if let Some(&inner_skin_fh) = inner_map.faces.get(&orig_fh) {
            faces.insert(outer_skin_fh, inner_skin_fh);
        }
    }

    SkinMap {
        vertices,
        edges,
        half_edges,
        faces,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn thicken_cube_produces_result() {
        // A cube is a closed solid, not a sheet, but thicken should still
        // work (producing outer + inner skins with no side walls).
        let brep = BRep::unit_cube(tol()).unwrap();
        let result = thicken(&brep, 0.2, &ShellOptions::default(), tol()).unwrap();
        // Outer (6) + inner (6) = 12 faces.
        assert_eq!(result.brep.faces.len(), 12);
    }

    #[test]
    fn thicken_empty_brep_errors() {
        let brep = BRep::new(tol());
        let result = thicken(&brep, 0.1, &ShellOptions::default(), tol());
        assert!(matches!(result, Err(ShellError::NoFaces)));
    }
}
