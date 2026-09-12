//! Mixed-material-sense boolean rescue (PK-478 stage 4).
//!
//! Orchestrates [`super::execute_boolean_pipeline`]: the primary run
//! is untouched, a failed MIXED-sense pair is retried once with its
//! inward operand reversed, and the rescue is accepted only when the
//! served body passes the interior-divergence tear check. See the
//! function doc and docs/architecture/pk478-sense-aware-merge.md.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::execute_boolean_pipeline;
use crate::cancellation::CancellationToken;
use crate::error::{BoolError, BoolOperation, PhaseMetric};
use crate::progress::ProgressReporter;

/// Execute a boolean with full pipeline knobs: telemetry, cancellation,
/// and progress reporting.
///
/// Appends a [`PhaseMetric`] after each stage that runs to completion,
/// polls `cancel` between stages, and notifies `progress` (if any)
/// after each stage with a [`ProgressEvent`]. Returns
/// [`BoolError::Cancelled`] as soon as a cancellation is observed;
/// the input model is unchanged.
pub(crate) fn execute_boolean_full(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tolerance: &Tolerance,
    metrics: &mut Vec<PhaseMetric>,
    cancel: Option<&CancellationToken>,
    progress: Option<&mut dyn ProgressReporter>,
) -> Result<BRep, BoolError> {
    // PK-478 stage 4: MIXED-sense rescue retry. A mixed operand pair
    // (one majority-inward — every blend-rebuilt body; one
    // majority-outward — sweep/extrude output) has walks that
    // traverse every shared edge in the same 3D direction; loops
    // mixing shared and operand-own edges then wrong-twin and the
    // assembly gapes or strands (measured: one-pitch crest-band
    // gaps, half-claimed strands). The primary run stays UNTOUCHED —
    // matched pairs and every historically-serving mixed pair are
    // byte-identical — and only a FAILED mixed pair is retried with
    // its inward operand reversed to outward
    // (`reverse_body_orientation`, pcurve-correct since PK-475;
    // registration is loop-covariant since the audit-first
    // graph_build, which is what dissolved the measured reversal
    // falsification). The rescue is accepted only when the served
    // body passes the interior-divergence tear check on top of the
    // pipeline's own validation: at degraded tolerance a rescued
    // assembly can serve a validator-clean TORN body (measured:
    // torus-box at 1e-2 — 9 faces, 193 open mesh edges, 22% of the
    // true volume, 4 torn edges at the calibrated 5x-linear
    // threshold); the honest primary decline is served instead.
    // `component_material_sense` is a majority statistic seeded per
    // body; `None` (degenerate sampling) disables the rescue.
    let sense_of = |body: &BRep| {
        body.faces
            .iter()
            .next()
            .and_then(|(fh, _)| body.component_material_sense(fh))
    };
    let (sa, sb) = (sense_of(&a), sense_of(&b));
    let mixed = {
        use popkern_core::validation::MaterialSense::{FoldedInward, FoldedOutward};
        matches!(
            (sa, sb),
            (Some(FoldedInward), Some(FoldedOutward)) | (Some(FoldedOutward), Some(FoldedInward))
        )
    };
    let rescue_operands = mixed.then(|| (a.clone(), b.clone()));

    // `progress` is consumed by the primary run; the rescue run is
    // internal and reports no progress.
    let primary = execute_boolean_pipeline(op, a, b, tolerance, metrics, cancel, progress);
    let err = match primary {
        Ok(result) => return Ok(result),
        Err(BoolError::Cancelled) => return Err(BoolError::Cancelled),
        Err(e) => e,
    };
    let Some((mut ra, mut rb)) = rescue_operands else {
        return Err(err);
    };
    {
        use popkern_core::validation::MaterialSense;
        if sa == Some(MaterialSense::FoldedInward) {
            ra.reverse_body_orientation();
        } else {
            rb.reverse_body_orientation();
        }
    }
    if popkern_core::trace::trace_enabled("boolean") {
        eprintln!(
            "  [exec] mixed-sense pair declined ({err}); retrying with the inward operand reversed"
        );
    }
    match execute_boolean_pipeline(op, ra, rb, tolerance, &mut Vec::new(), cancel, None) {
        Ok(rescued) => {
            // Tear bar: a divergence at 1% of the body extent is a
            // feature-scale hole, not an approximation band
            // (measured: the coarse torus-box lie reads 2.7% of its
            // diagonal on all four torn edges, while the rescued
            // thread cut's legitimate helical SSI band peaks at 0.1%
            // of its body). The 5x-linear floor keeps tiny bodies
            // from noise rejections and alone rejects the measured
            // coarse lie.
            let diag = {
                let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
                for (_, v) in &rescued.vertices {
                    let p = v.position;
                    for (i, c) in [p.x, p.y, p.z].into_iter().enumerate() {
                        lo[i] = lo[i].min(c);
                        hi[i] = hi[i].max(c);
                    }
                }
                if rescued.vertices.is_empty() {
                    0.0
                } else {
                    ((hi[0] - lo[0]).powi(2) + (hi[1] - lo[1]).powi(2) + (hi[2] - lo[2]).powi(2))
                        .sqrt()
                }
            };
            let tear_threshold = (0.01 * diag).max(5.0 * tolerance.linear().max(1e-9));
            let tears = rescued.find_tessellation_non_manifold_edges(16, tear_threshold);
            if tears.is_empty() {
                if popkern_core::trace::trace_enabled("boolean") {
                    eprintln!("  [exec] mixed-sense rescue serves (tear check clean)");
                }
                Ok(rescued)
            } else {
                if popkern_core::trace::trace_enabled("boolean") {
                    let worst = tears
                        .iter()
                        .map(|r| r.max_divergence)
                        .fold(0.0_f64, f64::max);
                    eprintln!(
                        "  [exec] mixed-sense rescue REJECTED: {} torn edge(s), worst divergence {worst:.3e}; serving the primary decline",
                        tears.len()
                    );
                }
                Err(err)
            }
        }
        Err(_) => Err(err),
    }
}
