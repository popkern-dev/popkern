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

use crate::error::BoolDiagnostics;
use crate::intersect::IntersectionRecord;

#[path = "split_helpers.rs"]
mod helpers;
use helpers::{find_fragment_for_curve, split_one_face};

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

/// Apply all intersection curves for one face of solid A or B.
///
/// Processes curves in order. The first curve targets `original_face`.
/// Subsequent curves locate the correct fragment via UV midpoint test
/// and split that fragment, with graceful fallback on failure.
#[allow(clippy::too_many_arguments)]
fn split_face_curves(
    brep: &mut BRep,
    original_face: FaceHandle,
    record_indices: &[usize],
    records: &[IntersectionRecord],
    tolerance: &Tolerance,
    fragments: &mut FragmentMap,
    use_uv1: bool,
) -> Result<(), BoolDiagnostics> {
    for (i, &record_idx) in record_indices.iter().enumerate() {
        let record = &records[record_idx];
        let curve_uv = if use_uv1 {
            &record.curve.curve_uv1
        } else {
            &record.curve.curve_uv2
        };

        let target_face = if i == 0 {
            original_face
        } else {
            // Multi-curve splitting: locate the fragment containing this
            // curve's UV midpoint and split that fragment.
            match find_fragment_for_curve(brep, original_face, curve_uv, fragments, tolerance) {
                Some(frag) => frag,
                None => continue,
            }
        };

        if !brep.faces.contains_key(target_face) {
            continue;
        }

        split_one_face(
            brep,
            target_face,
            &record.curve.curve_3d,
            curve_uv,
            &record.curve,
            tolerance,
            records,
            fragments,
        )?;
    }
    Ok(())
}

/// Split the faces of two BRep solids along intersection curves.
///
/// For each intersection record, splits the corresponding face of A using
/// `curve_uv1` (the UV curve on surface A) and the face of B using
/// `curve_uv2` (the UV curve on surface B).
///
/// All intersection curves per face are processed in sequence. After each
/// split, subsequent curves target whichever fragment contains the curve's
/// midpoint. Faces with no intersection curves map to themselves in the
/// fragment map.
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

    // Group records by face_a and face_b (collect all curves per face).
    let mut curves_a: HashMap<FaceHandle, Vec<usize>> = HashMap::new();
    let mut curves_b: HashMap<FaceHandle, Vec<usize>> = HashMap::new();

    for (idx, record) in records.iter().enumerate() {
        curves_a.entry(record.face_a).or_default().push(idx);
        curves_b.entry(record.face_b).or_default().push(idx);
    }

    // Split faces in solid A.
    for (&original_face, record_indices) in &curves_a {
        if !a.faces.contains_key(original_face) {
            continue;
        }
        split_face_curves(
            &mut a,
            original_face,
            record_indices,
            records,
            tolerance,
            &mut fragments_a,
            true,
        )?;
    }

    // Split faces in solid B.
    for (&original_face, record_indices) in &curves_b {
        if !b.faces.contains_key(original_face) {
            continue;
        }
        split_face_curves(
            &mut b,
            original_face,
            record_indices,
            records,
            tolerance,
            &mut fragments_b,
            false,
        )?;
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

#[cfg(test)]
#[path = "split_extra_tests.rs"]
mod extra_tests;

#[cfg(test)]
#[path = "split_multicurve_tests.rs"]
mod multicurve_tests;
