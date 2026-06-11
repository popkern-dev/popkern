//! Measurement helpers for `validation_sweep.rs`.
//!
//! Extracted to keep the main sweep file under the 500-line
//! structural limit. Will be reorganised further in Phase H.7
//! (validation sweep + commercial-grade gate) when additional
//! acceptance gates join.
//!
//! Reachable from `validation_sweep.rs` via `#[path] mod scoring`.

#![allow(dead_code)]

use popkern_core::geometry::Vec3;
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

pub(crate) fn unit_normal(gp: &GregoryPatch, u: f64, v: f64) -> Option<Vec3> {
    gp.partial_u(u, v).cross(&gp.partial_v(u, v)).normalized()
}

pub(crate) struct Scorecard {
    pub n_samples: usize,
    pub min_cos: f64,
    pub mean_cos: f64,
    pub max_cos: f64,
    pub max_angle_rad: f64,
    pub max_position_gap: f64,
}

pub(crate) fn sample_scorecard(patches: &[GregoryPatch], n_samples: usize) -> Scorecard {
    let np = patches.len();
    if np < 2 {
        return Scorecard {
            n_samples: 0,
            min_cos: 1.0,
            mean_cos: 1.0,
            max_cos: 1.0,
            max_angle_rad: 0.0,
            max_position_gap: 0.0,
        };
    }
    let ts: Vec<f64> = (1..=n_samples)
        .map(|k| (k as f64) / ((n_samples + 1) as f64))
        .collect();
    let mut cos_samples: Vec<f64> = Vec::new();
    let mut max_gap = 0.0_f64;
    for k in 0..np {
        let a = &patches[k];
        let b = &patches[(k + 1) % np];
        for &t in &ts {
            let t_b = 1.0 - t;
            let Some(na) = unit_normal(a, 1.0, t) else {
                continue;
            };
            let Some(nb) = unit_normal(b, t_b, 0.0) else {
                continue;
            };
            let cos = (na.x * nb.x + na.y * nb.y + na.z * nb.z).abs();
            cos_samples.push(cos);
            let pa = a.point_at(1.0, t);
            let pb = b.point_at(t_b, 0.0);
            let gap =
                ((pa.x - pb.x).powi(2) + (pa.y - pb.y).powi(2) + (pa.z - pb.z).powi(2)).sqrt();
            max_gap = max_gap.max(gap);
        }
    }
    let n = cos_samples.len();
    let min_cos = cos_samples.iter().copied().fold(f64::INFINITY, f64::min);
    let max_cos = cos_samples
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let mean_cos = if n > 0 {
        cos_samples.iter().sum::<f64>() / (n as f64)
    } else {
        0.0
    };
    let max_angle_rad = min_cos.clamp(-1.0, 1.0).acos();
    Scorecard {
        n_samples: n,
        min_cos,
        mean_cos,
        max_cos,
        max_angle_rad,
        max_position_gap: max_gap,
    }
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
            GregoryPatch::RectangularBiNonic(_) | GregoryPatch::RectangularBiNonicDualTwist(_) => {
                vc.binonic += 1
            }
            GregoryPatch::Triangular(_) => vc.triangular += 1,
        }
    }
    vc
}
