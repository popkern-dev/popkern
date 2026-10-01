//! Measurement helpers for `validation_sweep.rs`.
//!
//! Extracted to keep the main sweep file under the 500-line
//! structural limit. Will be reorganised further in the commercial-gate reorganisation
//! (validation sweep + commercial-grade gate) when additional
//! acceptance gates join.
//!
//! Reachable from `validation_sweep.rs` via `#[path] mod scoring`.

#![allow(dead_code)]

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

pub(crate) fn gregory_patch_of(brep: &BRep, face: FaceHandle) -> Option<GregoryPatch> {
    // `.get()` not indexing: sliver-absorb and cap-insert passes
    // occasionally delete a corner face after the pipeline has
    // recorded it in `corner_patches` (cube_extreme_40 + N ≥ 5
    // folded-sub-patch cases). Matches the pattern in
    // `kp_gate_helpers.rs::patch_of`.
    brep.faces.get(face)?.surface.as_ref().and_then(|sof| {
        if let SurfaceGeometry::Gregory(gp) = &sof.geometry {
            Some(gp.clone())
        } else {
            None
        }
    })
}

pub(crate) struct Scorecard {
    pub n_samples: usize,
    pub min_cos: f64,
    pub mean_cos: f64,
    pub max_cos: f64,
    pub max_angle_rad: f64,
    pub max_position_gap: f64,
}

pub(crate) struct RungResult {
    pub label: &'static str,
    pub outcome: Outcome,
}

pub(crate) enum Outcome {
    FilletFailed(String),
    NoGregorySubpatches,
    Success {
        n_sub_patches: usize,
        variant_counts: VariantCounts,
        brep_valid: bool,
        scorecard: Scorecard,
    },
}

#[derive(Default)]
pub(crate) struct VariantCounts {
    pub deg5: usize,
    pub biquintic: usize,
    pub biseptic: usize,
    pub two_stripe: usize,
    pub binonic: usize,
    pub cubic: usize,
    pub triangular: usize,
}

pub(crate) fn count_variants(patches: &[GregoryPatch]) -> VariantCounts {
    let mut vc = VariantCounts::default();
    for p in patches {
        match p {
            GregoryPatch::Rectangular(_) => vc.cubic += 1,
            GregoryPatch::RectangularDeg5(_) => vc.deg5 += 1,
            GregoryPatch::RectangularBiquinticFull(_) => vc.biquintic += 1,
            GregoryPatch::RectangularBiseptic(_) => vc.biseptic += 1,
            GregoryPatch::RectangularTwoStripe(_) => vc.two_stripe += 1,
            GregoryPatch::RectangularBiNonic(_)
            | GregoryPatch::RectangularBiNonicDualTwist(_)
            | GregoryPatch::RectangularBiNonicDualJet(_) => vc.binonic += 1,
            GregoryPatch::Triangular(_) => vc.triangular += 1,
        }
    }
    vc
}
