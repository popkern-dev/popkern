//! Geometry-type-specific healing: specialized repair strategies that
//! preserve invariants for fillet, offset, and swept surfaces.
//!
//! Since all complex surfaces become NURBS after construction, surface
//! type is inferred from geometric properties rather than stored metadata:
//! - Fillet: G1-continuous with both adjacent faces (tangent blend)
//! - Offset: approximately constant distance from an adjacent face
//! - Swept: iso-parameter curves in one direction have similar shape

#[path = "type_specific_helpers.rs"]
mod helpers;

use crate::geometry::surface::SurfaceGeometry;
use crate::healing::continuity::enforcement;
use crate::healing::quality::ContinuityLevel;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle};

use helpers::{
    adjacent_faces, adjust_offset_control_points, adjust_swept_control_points,
    edges_shared_with_neighbors, find_best_offset_base, measure_cross_section_uniformity,
    measure_offset_consistency, try_detect_fillet, try_detect_offset, try_detect_swept,
};

/// Detected geometric type of a face's surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceType {
    /// A blend/fillet surface tangent to both adjacent base surfaces.
    Fillet,
    /// A surface at approximately constant distance from an adjacent face.
    Offset,
    /// A surface whose cross-sections are approximately congruent.
    Swept,
    /// An ordinary surface with no detected special type.
    General,
}

/// Result of type-specific surface detection.
#[derive(Debug, Clone)]
pub struct TypeDetectionResult {
    /// The face analyzed.
    pub face: FaceHandle,
    /// Detected surface type.
    pub surface_type: SurfaceType,
    /// Confidence of detection (0.0 to 1.0).
    pub confidence: f64,
    /// Description of what was detected.
    pub description: String,
}

/// Result of type-specific repair.
#[derive(Debug, Clone)]
pub struct TypeSpecificRepairResult {
    /// The face repaired.
    pub face: FaceHandle,
    /// The surface type that was preserved.
    pub surface_type: SurfaceType,
    /// Maximum displacement introduced.
    pub max_displacement: f64,
    /// Whether the invariant was preserved.
    pub invariant_preserved: bool,
    /// Description of what was done.
    pub description: String,
}

// ---------------------------------------------------------------------------
// Detection
// ---------------------------------------------------------------------------

/// Detect the geometric type of a face's surface.
///
/// Analyzes the face's geometry and its relationship to adjacent faces
/// to determine if it's a fillet, offset, swept, or general surface.
pub fn detect_surface_type(
    brep: &BRep,
    face: FaceHandle,
    tolerance: &Tolerance,
) -> TypeDetectionResult {
    let general = TypeDetectionResult {
        face,
        surface_type: SurfaceType::General,
        confidence: 1.0,
        description: "no special type detected".to_string(),
    };

    let Some(face_data) = brep.faces.get(face) else {
        return general;
    };
    let Some(ref sof) = face_data.surface else {
        return general;
    };
    let SurfaceGeometry::Nurbs(ref nurbs) = sof.geometry else {
        return general;
    };

    let neighbors = adjacent_faces(brep, face);

    if let Some(result) = try_detect_fillet(brep, face, &neighbors, tolerance) {
        return result;
    }

    if let Some(result) = try_detect_offset(brep, face, &sof.geometry, &neighbors, tolerance) {
        return result;
    }

    if let Some(result) = try_detect_swept(face, nurbs, tolerance) {
        return result;
    }

    general
}

/// Detect surface types for all faces in a BRep.
pub fn detect_all_surface_types(brep: &BRep, tolerance: &Tolerance) -> Vec<TypeDetectionResult> {
    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
    face_handles
        .into_iter()
        .map(|fh| detect_surface_type(brep, fh, tolerance))
        .collect()
}

// ---------------------------------------------------------------------------
// Type-specific repair
// ---------------------------------------------------------------------------

/// Repair a fillet surface by restoring tangency with adjacent faces.
///
/// For a face detected as a fillet: measures tangency deviation at
/// both boundary curves, adjusts the fillet surface's boundary-row
/// control points to restore G1 continuity with both base surfaces.
pub fn repair_fillet(
    brep: &mut BRep,
    face: FaceHandle,
    tolerance: &Tolerance,
    max_displacement: f64,
) -> TypeSpecificRepairResult {
    let detection = detect_surface_type(brep, face, tolerance);
    if detection.surface_type != SurfaceType::Fillet {
        return TypeSpecificRepairResult {
            face,
            surface_type: SurfaceType::Fillet,
            max_displacement: 0.0,
            invariant_preserved: true,
            description: "face is not a fillet; no repair needed".to_string(),
        };
    }

    let shared_edges = edges_shared_with_neighbors(brep, face);
    let mut total_disp = 0.0_f64;
    let mut all_ok = true;

    for edge in &shared_edges {
        match enforcement::enforce_g1(brep, *edge, tolerance, max_displacement) {
            Ok(result) => {
                total_disp = total_disp.max(result.max_displacement);
                if result.achieved_level < ContinuityLevel::G1 {
                    all_ok = false;
                }
            }
            Err(_) => {
                all_ok = false;
            }
        }
    }

    TypeSpecificRepairResult {
        face,
        surface_type: SurfaceType::Fillet,
        max_displacement: total_disp,
        invariant_preserved: all_ok,
        description: format!(
            "fillet repair: {} shared edges processed, max disp {total_disp:.6e}",
            shared_edges.len()
        ),
    }
}

/// Repair an offset surface by restoring constant distance from base.
///
/// For a face detected as offset: measures the distance variation
/// between this surface and the base (adjacent face). Adjusts control
/// points to reduce distance variance while preserving the offset
/// distance at the median.
pub fn repair_offset(
    brep: &mut BRep,
    face: FaceHandle,
    tolerance: &Tolerance,
    max_displacement: f64,
) -> TypeSpecificRepairResult {
    let detection = detect_surface_type(brep, face, tolerance);
    if detection.surface_type != SurfaceType::Offset {
        return TypeSpecificRepairResult {
            face,
            surface_type: SurfaceType::Offset,
            max_displacement: 0.0,
            invariant_preserved: true,
            description: "face is not an offset; no repair needed".to_string(),
        };
    }

    let neighbors = adjacent_faces(brep, face);
    let sof = match brep.faces.get(face).and_then(|f| f.surface.clone()) {
        Some(s) => s,
        None => {
            return TypeSpecificRepairResult {
                face,
                surface_type: SurfaceType::Offset,
                max_displacement: 0.0,
                invariant_preserved: false,
                description: "face missing surface geometry".to_string(),
            };
        }
    };

    let base_face = find_best_offset_base(brep, &sof.geometry, &neighbors);
    let base_geom = base_face.and_then(|bf| {
        brep.faces
            .get(bf)
            .and_then(|f| f.surface.as_ref().map(|s| s.geometry.clone()))
    });

    let Some(base_geom) = base_geom else {
        return TypeSpecificRepairResult {
            face,
            surface_type: SurfaceType::Offset,
            max_displacement: 0.0,
            invariant_preserved: false,
            description: "no base surface found for offset repair".to_string(),
        };
    };

    let SurfaceGeometry::Nurbs(ref nurbs) = sof.geometry else {
        return TypeSpecificRepairResult {
            face,
            surface_type: SurfaceType::Offset,
            max_displacement: 0.0,
            invariant_preserved: true,
            description: "non-NURBS offset surface; skipped".to_string(),
        };
    };

    let (median_dist, _, _) = measure_offset_consistency(&sof.geometry, &base_geom, 16);
    if median_dist < tolerance.linear() {
        return TypeSpecificRepairResult {
            face,
            surface_type: SurfaceType::Offset,
            max_displacement: 0.0,
            invariant_preserved: true,
            description: "offset distance near zero; skipped".to_string(),
        };
    }

    let mut new_nurbs = nurbs.clone();
    let disp =
        adjust_offset_control_points(&mut new_nurbs, &base_geom, median_dist, max_displacement);

    if disp > 0.0
        && let Some(face_mut) = brep.faces.get_mut(face)
        && let Some(ref mut s) = face_mut.surface
    {
        s.geometry = SurfaceGeometry::Nurbs(new_nurbs);
    }

    let sof2 = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry);
    let invariant_ok = if let Some(geom) = sof2 {
        let (_, max_dev, _) = measure_offset_consistency(geom, &base_geom, 16);
        max_dev < tolerance.linear() * 10.0
    } else {
        false
    };

    TypeSpecificRepairResult {
        face,
        surface_type: SurfaceType::Offset,
        max_displacement: disp,
        invariant_preserved: invariant_ok,
        description: format!(
            "offset repair: target distance {median_dist:.6e}, max disp {disp:.6e}"
        ),
    }
}

/// Repair a swept surface by restoring cross-section uniformity.
///
/// For a face detected as swept: measures the deviation of cross-sections
/// from the reference (first) cross-section. Adjusts control points to
/// reduce cross-section variation.
pub fn repair_swept(
    brep: &mut BRep,
    face: FaceHandle,
    tolerance: &Tolerance,
    max_displacement: f64,
) -> TypeSpecificRepairResult {
    let detection = detect_surface_type(brep, face, tolerance);
    if detection.surface_type != SurfaceType::Swept {
        return TypeSpecificRepairResult {
            face,
            surface_type: SurfaceType::Swept,
            max_displacement: 0.0,
            invariant_preserved: true,
            description: "face is not a swept surface; no repair needed".to_string(),
        };
    }

    let sof = match brep.faces.get(face).and_then(|f| f.surface.clone()) {
        Some(s) => s,
        None => {
            return TypeSpecificRepairResult {
                face,
                surface_type: SurfaceType::Swept,
                max_displacement: 0.0,
                invariant_preserved: false,
                description: "face missing surface geometry".to_string(),
            };
        }
    };

    let SurfaceGeometry::Nurbs(ref nurbs) = sof.geometry else {
        return TypeSpecificRepairResult {
            face,
            surface_type: SurfaceType::Swept,
            max_displacement: 0.0,
            invariant_preserved: true,
            description: "non-NURBS swept surface; skipped".to_string(),
        };
    };

    let mut new_nurbs = nurbs.clone();
    let disp = adjust_swept_control_points(&mut new_nurbs, max_displacement);

    if disp > 0.0
        && let Some(face_mut) = brep.faces.get_mut(face)
        && let Some(ref mut s) = face_mut.surface
    {
        s.geometry = SurfaceGeometry::Nurbs(new_nurbs);
    }

    let sof2 = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry);
    let invariant_ok = if let Some(SurfaceGeometry::Nurbs(n)) = sof2 {
        let (var_ratio, _) = measure_cross_section_uniformity(n, 5);
        var_ratio < 0.05
    } else {
        false
    };

    TypeSpecificRepairResult {
        face,
        surface_type: SurfaceType::Swept,
        max_displacement: disp,
        invariant_preserved: invariant_ok,
        description: format!("swept repair: max disp {disp:.6e}"),
    }
}

/// Run type-specific repair on all detected special surfaces.
pub fn repair_all_type_specific(
    brep: &mut BRep,
    tolerance: &Tolerance,
    max_displacement: f64,
) -> Vec<TypeSpecificRepairResult> {
    let detections = detect_all_surface_types(brep, tolerance);
    let mut results = Vec::new();

    for det in &detections {
        let result = match det.surface_type {
            SurfaceType::Fillet => repair_fillet(brep, det.face, tolerance, max_displacement),
            SurfaceType::Offset => repair_offset(brep, det.face, tolerance, max_displacement),
            SurfaceType::Swept => repair_swept(brep, det.face, tolerance, max_displacement),
            SurfaceType::General => continue,
        };
        results.push(result);
    }

    results
}

#[cfg(test)]
#[path = "type_specific_tests.rs"]
mod tests;
