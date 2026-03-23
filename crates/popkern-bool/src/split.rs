//! Stage 2: Split faces of both solids along intersection curves.
//!
//! Groups intersection records by face, then calls `BRep::split_face` on each
//! affected face to produce face fragments. Tracks provenance (which original
//! face each fragment came from) and which edges were created by splitting.

use std::collections::{HashMap, HashSet};

use popkern_core::{
    tolerance::Tolerance,
    topology::{BRep, EdgeHandle, FaceHandle},
};

use crate::error::{BoolDiagnostics, BoolError, Stage};
use crate::intersect::IntersectionRecord;

/// Provenance map tracking which original face each fragment came from,
/// plus which edges were created by the splitting process.
#[derive(Debug, Clone)]
pub(crate) struct FragmentMap {
    /// Maps each new fragment face handle back to the original face it was split from.
    pub origin: HashMap<FaceHandle, FaceHandle>,
    /// The set of edges created during splitting (trim curve edges).
    pub intersection_edges: HashSet<EdgeHandle>,
}

/// Result of splitting both solids along intersection curves.
#[derive(Debug)]
pub(crate) struct SplitResult {
    /// The modified solid A (with faces split).
    pub a: BRep,
    /// The modified solid B (with faces split).
    pub b: BRep,
    /// Fragment provenance for solid A.
    pub fragments_a: FragmentMap,
    /// Fragment provenance for solid B.
    pub fragments_b: FragmentMap,
}

/// Split the faces of two BRep solids along intersection curves.
///
/// For each intersection record, splits the corresponding face of A using
/// `curve_uv1` (the UV curve on surface A) and the face of B using
/// `curve_uv2` (the UV curve on surface B).
///
/// Only the first intersection curve per face is applied. Faces with no
/// intersection curves map to themselves in the fragment map.
///
/// # Errors
///
/// Returns [`BoolDiagnostics`] if any `split_face` call fails.
pub(crate) fn split_along_curves(
    mut a: BRep,
    mut b: BRep,
    records: &[IntersectionRecord],
    tolerance: &Tolerance,
) -> Result<SplitResult, BoolDiagnostics> {
    let mut fragments_a = FragmentMap {
        origin: HashMap::new(),
        intersection_edges: HashSet::new(),
    };
    let mut fragments_b = FragmentMap {
        origin: HashMap::new(),
        intersection_edges: HashSet::new(),
    };

    // Group records by face_a and face_b (take only the first curve per face).
    let mut first_curve_a: HashMap<FaceHandle, usize> = HashMap::new();
    let mut first_curve_b: HashMap<FaceHandle, usize> = HashMap::new();

    for (idx, record) in records.iter().enumerate() {
        first_curve_a.entry(record.face_a).or_insert(idx);
        first_curve_b.entry(record.face_b).or_insert(idx);
    }

    // Split faces in solid A.
    for (&face_handle, &record_idx) in &first_curve_a {
        // Skip faces that no longer exist (e.g. stale handles).
        if !a.faces.contains_key(face_handle) {
            continue;
        }
        let record = &records[record_idx];
        let curve = &record.curve;

        // For face A: curve_uv1 is the UV parameterization on surface A.
        // split_face expects (face, trim_curve_uv, trim_curve_3d, tolerance).
        match a.split_face(face_handle, &curve.curve_uv1, &curve.curve_3d, tolerance) {
            Ok(result) => {
                fragments_a.origin.insert(result.face_a, face_handle);
                fragments_a.origin.insert(result.face_b, face_handle);
                for &edge in &result.new_edges {
                    fragments_a.intersection_edges.insert(edge);
                }
            }
            Err(core_err) => {
                return Err(BoolDiagnostics {
                    error: Some(BoolError::SplitFailed {
                        face: face_handle,
                        reason: core_err.to_string(),
                    }),
                    intersection_curves: records.iter().map(|r| r.curve.clone()).collect(),
                    face_classifications: Vec::new(),
                    failed_stage: Some(Stage::Split),
                });
            }
        }
    }

    // Split faces in solid B.
    for (&face_handle, &record_idx) in &first_curve_b {
        // Skip faces that no longer exist (e.g. stale handles).
        if !b.faces.contains_key(face_handle) {
            continue;
        }
        let record = &records[record_idx];
        let curve = &record.curve;

        // For face B: curve_uv2 is the UV parameterization on surface B.
        // It becomes "uv1" relative to face B's surface in the split_face call.
        match b.split_face(face_handle, &curve.curve_uv2, &curve.curve_3d, tolerance) {
            Ok(result) => {
                fragments_b.origin.insert(result.face_a, face_handle);
                fragments_b.origin.insert(result.face_b, face_handle);
                for &edge in &result.new_edges {
                    fragments_b.intersection_edges.insert(edge);
                }
            }
            Err(core_err) => {
                return Err(BoolDiagnostics {
                    error: Some(BoolError::SplitFailed {
                        face: face_handle,
                        reason: core_err.to_string(),
                    }),
                    intersection_curves: records.iter().map(|r| r.curve.clone()).collect(),
                    face_classifications: Vec::new(),
                    failed_stage: Some(Stage::Split),
                });
            }
        }
    }

    // All unsplit faces in A map to themselves.
    for face_handle in a.faces.keys() {
        fragments_a.origin.entry(face_handle).or_insert(face_handle);
    }

    // All unsplit faces in B map to themselves.
    for face_handle in b.faces.keys() {
        fragments_b.origin.entry(face_handle).or_insert(face_handle);
    }

    Ok(SplitResult {
        a,
        b,
        fragments_a,
        fragments_b,
    })
}

#[cfg(test)]
#[path = "split_tests.rs"]
mod tests;
