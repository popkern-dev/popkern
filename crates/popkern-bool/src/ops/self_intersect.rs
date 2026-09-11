//! Global self-intersection scan for a single body (PK-327, design §5
//! Stage 4 tier 2).
//!
//! Finds transversal intersection curves between NON-ADJACENT faces of
//! one BRep: BVH broad-phase over face AABBs, adjacency exemption
//! (faces sharing an edge through any loop), analytic-coincidence
//! exemption (same-surface pairs are the local-veto / spindle arms'
//! job — a surface cannot cross itself transversally), then
//! surface-surface intersection clipped to both trimmed domains.
//! Degenerate clip remnants (touch-at-a-vertex pairs) are dropped by a
//! 3D length gate.
//!
//! The scan is detection only: it never mutates the body. The trim
//! stage consumes the records.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::intersection::surface_surface::{
    SurfaceSurfaceIntersection, intersect_surfaces,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::clip::clip_curve_to_faces;
use crate::error::BoolError;
use crate::intersect::bounded_parameter_range;

/// One transversal self-intersection curve between two faces of the
/// scanned body, clipped to both trimmed domains.
#[derive(Debug, Clone)]
pub struct SelfIntersectionCurve {
    /// First face (scan order: `face_a < face_b` by slotmap order).
    pub face_a: FaceHandle,
    /// Second face.
    pub face_b: FaceHandle,
    /// The clipped curve: 3D plus pcurves on `face_a` (`curve_uv1`)
    /// and `face_b` (`curve_uv2`).
    pub curve: CurveOnSurface,
}

/// Result of [`find_self_intersections`].
#[derive(Debug, Clone, Default)]
pub struct SelfIntersectionScan {
    /// All transversal self-intersection curves found.
    pub curves: Vec<SelfIntersectionCurve>,
    /// Broad-phase candidate pairs that survived the exemptions and
    /// were handed to the surface-surface solver (observability).
    pub pairs_checked: usize,
    /// Candidate pairs skipped because their surfaces are analytically
    /// coincident (same-surface folding is the spindle / distance-map
    /// arms' responsibility, recorded here so those arms can see the
    /// hand-off).
    pub coincident_pairs_skipped: usize,
}

impl SelfIntersectionScan {
    /// True when the scan found no transversal self-intersections.
    #[must_use]
    pub fn is_clean(&self) -> bool {
        self.curves.is_empty()
    }
}

/// Scan `body` for transversal intersection curves between
/// non-adjacent faces.
///
/// # Errors
///
/// Returns [`BoolError::IntersectionFailed`] when the surface-surface
/// solver fails on a candidate pair (the pair is named).
pub fn find_self_intersections(
    body: &BRep,
    tolerance: &Tolerance,
) -> Result<SelfIntersectionScan, BoolError> {
    let _ = tolerance;
    let mut scan = SelfIntersectionScan::default();
    let face_keys: Vec<FaceHandle> = body.faces.keys().collect();
    // Detection authority: the validator's own pair machinery
    // (adjacency exemption, BVH pruning, trim-overlap gates, pair
    // budget). Trimming exactly what the validator would reject keeps
    // the serve gate coherent — a scan stricter than validate_all
    // false-positives on grazing corner constructions the validator
    // deliberately clears (measured on the PK-326 corner slivers).
    let pairs = body.collect_self_intersecting_pairs_among(&face_keys, true, true);
    scan.pairs_checked = pairs.len();
    // Curves shorter than this are touch-at-a-vertex clip remnants,
    // not trim regions.
    let min_curve_length = body.tolerance.linear() * 10.0;

    // Per-pair scans are pure reads of `body` (plain SlotMap storage,
    // no interior mutability), so the pair list dispatches in parallel
    // with an ordered collect (PK-277 discipline; PK-329): results are
    // folded in the caller's pair order, so curves, counters, and the
    // first-failing-pair error are positionally identical to the
    // sequential loop at any thread count.
    let per_pair = |&(face_a, face_b): &(FaceHandle, FaceHandle)| -> Result<
        (Vec<SelfIntersectionCurve>, usize),
        BoolError,
    > {
        let (Some(sa), Some(sb)) = (
            body.faces[face_a].surface.as_ref(),
            body.faces[face_b].surface.as_ref(),
        ) else {
            return Ok((Vec::new(), 0));
        };
        let (ua, va) = bounded_parameter_range(body, face_a, &sa.geometry);
        let (ub, vb) = bounded_parameter_range(body, face_b, &sb.geometry);
        let result =
            intersect_surfaces(&sa.geometry, ua, va, &sb.geometry, ub, vb, &body.tolerance)
                .map_err(|e| {
                    BoolError::IntersectionFailed(format!(
                        "self-intersection scan pair {face_a:?}×{face_b:?}: {e}"
                    ))
                })?;

        let mut curves = Vec::new();
        let mut coincident_skipped = 0usize;
        for intersection in result.intersections {
            let SurfaceSurfaceIntersection::Curve(curve) = intersection else {
                // Whole-region coincidence: the same-surface arms'
                // responsibility (spindle / distance map).
                coincident_skipped += 1;
                continue;
            };
            for clipped in clip_curve_to_faces(body, face_a, body, face_b, &curve, None) {
                if curve_length_estimate(&clipped) < min_curve_length {
                    continue;
                }
                curves.push(SelfIntersectionCurve {
                    face_a,
                    face_b,
                    curve: clipped,
                });
            }
        }
        Ok((curves, coincident_skipped))
    };

    #[cfg(feature = "parallel")]
    let outcomes: Vec<Result<(Vec<SelfIntersectionCurve>, usize), BoolError>> = {
        use rayon::prelude::*;
        pairs.par_iter().map(per_pair).collect()
    };
    #[cfg(not(feature = "parallel"))]
    let outcomes: Vec<Result<(Vec<SelfIntersectionCurve>, usize), BoolError>> =
        pairs.iter().map(per_pair).collect();

    for outcome in outcomes {
        let (curves, coincident_skipped) = outcome?;
        scan.coincident_pairs_skipped += coincident_skipped;
        scan.curves.extend(curves);
    }
    Ok(scan)
}

/// Chord-sum estimate of a curve's 3D length (16 samples — a length
/// gate, not a measurement).
fn curve_length_estimate(curve: &CurveOnSurface) -> f64 {
    let (t0, t1) = curve.curve_3d.parameter_range();
    let mut length = 0.0;
    let mut prev = curve.curve_3d.point_at(t0);
    for i in 1..=16 {
        let t = t0 + (t1 - t0) * (i as f64) / 16.0;
        let p = curve.curve_3d.point_at(t);
        length += (p - prev).length();
        prev = p;
    }
    length
}
