//! Face classification for shell operations.

use std::collections::{HashMap, HashSet};

use crate::error::ShellError;
use crate::types::ShellDirection;
use popkern_core::topology::{BRep, FaceHandle};

/// Result of classifying faces for a shell operation.
#[derive(Debug, Clone)]
pub struct FaceClassification {
    /// Faces to be offset: (face_handle, signed_distance, inward_flag).
    pub offset_faces: Vec<(FaceHandle, f64, bool)>,
    /// Faces to be removed (open faces of the shell).
    pub removed_faces: Vec<FaceHandle>,
}

/// Classify faces for a shell operation.
///
/// Each face is either marked for removal or assigned an offset distance
/// and direction. The `face_thicknesses` map overrides `default_thickness`
/// for specific faces.
///
/// Directions are MATERIAL directions: `Inward` moves the skin into the
/// body's material (the hollow-shell void skin), `Outward` moves it out
/// of the material. popkern has no fixed normal-vs-material convention
/// (PK-253: a component's flag-folded normals may face either way — the
/// unit cube folds inward, the cylinder outward), so the mapping from a
/// material direction to a per-face geometric offset sign resolves the
/// component's [`popkern_core::validation::MaterialSense`] (wedge-integral
/// sign) and folds it with the per-face orientation flag. The historical
/// flag-only mapping served hollow bodies whose void skin sat OUTSIDE
/// the body.
///
/// `Midplane` halves the thickness and classifies the inner skin inward;
/// the caller owns moving the outer skin.
///
/// # Errors
///
/// Returns [`ShellError::NoFaces`] if the BRep has no faces, or
/// [`ShellError::MaterialSenseUnresolved`] when a component's material
/// side cannot be measured (offset direction would be meaningless).
pub fn classify_faces(
    brep: &BRep,
    default_thickness: f64,
    face_thicknesses: &[(FaceHandle, f64)],
    direction: ShellDirection,
    removed: &[FaceHandle],
) -> Result<FaceClassification, ShellError> {
    if brep.faces.is_empty() {
        return Err(ShellError::NoFaces);
    }

    let thickness_map: HashMap<FaceHandle, f64> = face_thicknesses.iter().copied().collect();
    let outward_multipliers = face_outward_multipliers(brep)?;

    let mut offset_faces = Vec::new();
    let mut removed_faces = Vec::new();

    for (fh, face) in &brep.faces {
        if removed.contains(&fh) {
            removed_faces.push(fh);
            continue;
        }

        let thickness = thickness_map.get(&fh).copied().unwrap_or(default_thickness);

        // Material direction requested by the caller.
        let (distance, base_inward) = match direction {
            ShellDirection::Inward => (thickness, true),
            ShellDirection::Outward => (thickness, false),
            ShellDirection::Midplane => (thickness / 2.0, true),
        };

        // Out-of-material coefficient along the STORED geometric normal:
        // outward = c_out * n, where c_out folds the component's material
        // sense with the face's orientation flag.
        let flag_sign = face
            .surface
            .as_ref()
            .map(|sof| if sof.orientation { 1.0 } else { -1.0 })
            .unwrap_or(1.0);
        let mult = outward_multipliers
            .get(&fh)
            .copied()
            .ok_or(ShellError::MaterialSenseUnresolved { face: fh })?;
        let c_out = mult * flag_sign;

        // The geometric `inward` flag consumed by the offset arms means
        // "offset along -n". Material-inward = -c_out·n, so the flag is
        // set exactly when the outward coefficient is positive.
        let inward = base_inward == (c_out > 0.0);

        offset_faces.push((fh, distance, inward));
    }

    Ok(FaceClassification {
        offset_faces,
        removed_faces,
    })
}

/// Resolve the out-of-material multiplier PER FACE via the
/// traversal-coherent normal (see [`crate::classify_sense`]).
///
/// The component-wide wedge multiplier assumed flag-folded normals
/// are uniform per component; served boolean outputs REFUTE that
/// (measured: the sphere-capped-rod union folds its cylinder wall and
/// bottom cap outward but its sphere cap center-ward — every face
/// winding-honest — so the component verdict moved the wall skin
/// OUTWARD under an "Inward" shell and served the mixed body
/// structurally clean at small thickness). The traversal-coherent
/// normal (outer loop CCW in UV, measured from the pcurve polygons)
/// is on ONE global side per component; the component wedge fixes
/// which. Faces whose winding vote is inconclusive fall back to the
/// legacy component verdict.
///
/// # Errors
///
/// [`ShellError::MaterialSenseUnresolved`] naming the first face
/// neither resolution can serve.
fn face_outward_multipliers(brep: &BRep) -> Result<HashMap<FaceHandle, f64>, ShellError> {
    let mut senses: HashMap<FaceHandle, crate::classify_sense::FaceSense> = HashMap::new();
    let mut face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
    face_handles.sort();
    for &fh in &face_handles {
        if let Some(s) = crate::classify_sense::face_sense(brep, fh) {
            senses.insert(fh, s);
        }
    }

    let mut multipliers: HashMap<FaceHandle, f64> = HashMap::new();
    let mut unresolved: Vec<FaceHandle> = Vec::new();
    let mut assigned: HashSet<FaceHandle> = HashSet::new();
    for &seed in &face_handles {
        if assigned.contains(&seed) {
            continue;
        }
        let component = shared_edge_component(brep, seed, &mut assigned);
        match crate::classify_sense::component_stored_normal_out_signs(brep, &component, &senses) {
            Some(signs) => {
                for &fh in &component {
                    match signs.get(&fh) {
                        // The sign applies to the STORED normal; the
                        // caller folds by the flag, so pre-divide it
                        // out (flag² = 1).
                        Some(&s) => {
                            let flag_sign = brep
                                .faces
                                .get(fh)
                                .and_then(|f| f.surface.as_ref())
                                .map(|sof| if sof.orientation { 1.0 } else { -1.0 })
                                .unwrap_or(1.0);
                            if std::env::var_os("PK_SENSE_TRACE").is_some() {
                                eprintln!("[sense-mult] {fh:?} sense s={s} flag={flag_sign}");
                            }
                            multipliers.insert(fh, s * flag_sign);
                        }
                        None => unresolved.push(fh),
                    }
                }
            }
            None => unresolved.extend(component.iter().copied()),
        }
    }

    // Fallback: the legacy component wedge verdict for faces the
    // winding vote could not measure (degenerate polygons, open
    // sheets whose wedge is the only signal).
    if !unresolved.is_empty() {
        let component = component_outward_multipliers(brep)?;
        for fh in unresolved {
            let mult = component
                .get(&fh)
                .copied()
                .ok_or(ShellError::MaterialSenseUnresolved { face: fh })?;
            if std::env::var_os("PK_SENSE_TRACE").is_some() {
                eprintln!("[sense-mult] {fh:?} LEGACY mult={mult}");
            }
            multipliers.insert(fh, mult);
        }
    }
    Ok(multipliers)
}

/// Faces reachable from `seed` over shared edges (marking them in
/// `assigned`), defensively walked.
fn shared_edge_component(
    brep: &BRep,
    seed: FaceHandle,
    assigned: &mut HashSet<FaceHandle>,
) -> Vec<FaceHandle> {
    let mut component = Vec::new();
    let mut queue = vec![seed];
    while let Some(fh) = queue.pop() {
        if !assigned.insert(fh) {
            continue;
        }
        component.push(fh);
        let Some(face) = brep.faces.get(fh) else {
            continue;
        };
        // ALL loops: a face bounded by a hole connects to its bore
        // wall through the INNER rim only (measured: the bracket's
        // outer-loop-only walk split {planes} from {bore, cone} and
        // the two components resolved OPPOSITE global sides).
        let loop_handles: Vec<_> = std::iter::once(face.outer_loop)
            .chain(face.inner_loops.iter().copied())
            .collect();
        for lh in loop_handles {
            let Some(loop_) = brep.loops.get(lh) else {
                continue;
            };
            let start = loop_.half_edge;
            let mut he = start;
            while let Some(he_data) = brep.half_edges.get(he) {
                if let Some(twin_data) = brep.half_edges.get(he_data.twin)
                    && let Some(partner_loop) = brep.loops.get(twin_data.loop_)
                {
                    let partner = partner_loop.face;
                    if !assigned.contains(&partner) {
                        queue.push(partner);
                    }
                }
                he = he_data.next;
                if he == start {
                    break;
                }
            }
        }
    }
    component.sort();
    component
}

/// Resolve the material-sense outward multiplier for every face, one
/// wedge-integral resolution per shared-edge-connected component
/// (fallback arm of [`face_outward_multipliers`]).
///
/// # Errors
///
/// [`ShellError::MaterialSenseUnresolved`] naming a seed face of the
/// first component whose sense cannot be measured.
fn component_outward_multipliers(brep: &BRep) -> Result<HashMap<FaceHandle, f64>, ShellError> {
    let mut multipliers = HashMap::new();
    let mut assigned: HashSet<FaceHandle> = HashSet::new();

    for (seed, _) in &brep.faces {
        if assigned.contains(&seed) {
            continue;
        }
        // Enumerate the component (faces reachable over shared edges) —
        // the same reachability component_material_sense resolves over.
        let mut component = Vec::new();
        let mut queue = vec![seed];
        while let Some(fh) = queue.pop() {
            if !assigned.insert(fh) {
                continue;
            }
            component.push(fh);
            let Some(face) = brep.faces.get(fh) else {
                continue;
            };
            let Some(loop_) = brep.loops.get(face.outer_loop) else {
                continue;
            };
            let start = loop_.half_edge;
            let mut he = start;
            // Free rims (open sheets) carry twin-less half-edges;
            // defensive access instead of indexing (PK-328: the
            // pre-flight query reaches classify on sheet bodies).
            while let Some(he_data) = brep.half_edges.get(he) {
                if let Some(twin_data) = brep.half_edges.get(he_data.twin)
                    && let Some(partner_loop) = brep.loops.get(twin_data.loop_)
                {
                    let partner = partner_loop.face;
                    if !assigned.contains(&partner) {
                        queue.push(partner);
                    }
                }
                he = he_data.next;
                if he == start {
                    break;
                }
            }
        }

        let sense = brep
            .component_material_sense(seed)
            .ok_or(ShellError::MaterialSenseUnresolved { face: seed })?;
        for fh in component {
            multipliers.insert(fh, sense.outward_multiplier());
        }
    }
    Ok(multipliers)
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    #[test]
    fn empty_brep_returns_no_faces() {
        let brep = BRep::new(tol());
        let result = classify_faces(&brep, 0.1, &[], ShellDirection::Inward, &[]);
        assert!(matches!(result, Err(ShellError::NoFaces)));
    }

    #[test]
    fn cube_all_faces_offset() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let result = classify_faces(&brep, 0.1, &[], ShellDirection::Inward, &[]).unwrap();
        assert_eq!(result.offset_faces.len(), 6);
        assert_eq!(result.removed_faces.len(), 0);
        // The unit cube stores plane normals pointing INTO the material
        // (component folds inward, orientation=true → c_out = −1), so
        // material-Inward maps to the along-+n geometric flag:
        // inward=false. The old flag-only pin (inward=true) encoded the
        // inverted semantics that served void-outside hollow bodies.
        for &(_, dist, inward) in &result.offset_faces {
            assert!((dist - 0.1).abs() < 1e-12);
            assert!(!inward, "material-Inward on the cube offsets along +n");
        }
    }

    #[test]
    fn cube_one_face_removed() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let first_face = brep.faces.keys().next().unwrap();
        let result =
            classify_faces(&brep, 0.1, &[], ShellDirection::Inward, &[first_face]).unwrap();
        assert_eq!(result.offset_faces.len(), 5);
        assert_eq!(result.removed_faces.len(), 1);
        assert_eq!(result.removed_faces[0], first_face);
    }

    #[test]
    fn midplane_halves_thickness() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let result = classify_faces(&brep, 0.2, &[], ShellDirection::Midplane, &[]).unwrap();
        for &(_, dist, inward) in &result.offset_faces {
            assert!((dist - 0.1).abs() < 1e-12);
            // Material-Inward on the cube = geometric along-+n (see
            // cube_all_faces_offset).
            assert!(!inward);
        }
    }

    #[test]
    fn outward_direction_sets_inward_false() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let result = classify_faces(&brep, 0.1, &[], ShellDirection::Outward, &[]).unwrap();
        for &(_, _, inward) in &result.offset_faces {
            // Material-Outward on the inward-folded cube = along −n.
            assert!(inward);
        }
    }

    #[test]
    fn per_face_thickness_override() {
        let brep = BRep::unit_cube(tol()).unwrap();
        let first_face = brep.faces.keys().next().unwrap();
        let result = classify_faces(
            &brep,
            0.1,
            &[(first_face, 0.5)],
            ShellDirection::Inward,
            &[],
        )
        .unwrap();
        // Find the overridden face.
        let overridden = result
            .offset_faces
            .iter()
            .find(|(fh, _, _)| *fh == first_face)
            .unwrap();
        assert!((overridden.1 - 0.5).abs() < 1e-12);
        // Other faces should have default thickness.
        for &(fh, dist, _) in &result.offset_faces {
            if fh != first_face {
                assert!((dist - 0.1).abs() < 1e-12);
            }
        }
    }
}
