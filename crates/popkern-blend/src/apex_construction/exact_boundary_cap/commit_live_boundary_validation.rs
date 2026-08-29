//! Live cap/support pcurve agreement validation for exact-boundary commit.

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

const LIVE_PAIR_SAMPLE_COUNT: usize = 17;

#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct LivePcurvePairResidual {
    pub residual: f64,
    pub sample_parameter: f64,
}

pub(super) fn live_pcurve_pair_residual(
    cap_pcurve: &Pcurve,
    cap_surface: &SurfaceGeometry,
    support_pcurve: &Pcurve,
    support_surface: &SurfaceGeometry,
) -> Option<LivePcurvePairResidual> {
    let (cap_start, cap_end) = cap_pcurve.parameter_range();
    let (support_start, support_end) = support_pcurve.parameter_range();
    if !cap_start.is_finite()
        || !cap_end.is_finite()
        || !support_start.is_finite()
        || !support_end.is_finite()
        || cap_start == cap_end
        || support_start == support_end
    {
        return None;
    }

    let forward = sampled_pair_residual(
        cap_pcurve,
        cap_surface,
        support_pcurve,
        support_surface,
        false,
    )?;
    let reversed = sampled_pair_residual(
        cap_pcurve,
        cap_surface,
        support_pcurve,
        support_surface,
        true,
    )?;
    Some(if reversed.residual < forward.residual {
        reversed
    } else {
        forward
    })
}

fn sampled_pair_residual(
    cap_pcurve: &Pcurve,
    cap_surface: &SurfaceGeometry,
    support_pcurve: &Pcurve,
    support_surface: &SurfaceGeometry,
    reverse_support: bool,
) -> Option<LivePcurvePairResidual> {
    let (cap_start, cap_end) = cap_pcurve.parameter_range();
    let (support_start, support_end) = support_pcurve.parameter_range();
    let mut worst = LivePcurvePairResidual {
        residual: 0.0,
        sample_parameter: 0.0,
    };
    for index in 0..LIVE_PAIR_SAMPLE_COUNT {
        let alpha = index as f64 / (LIVE_PAIR_SAMPLE_COUNT - 1) as f64;
        let cap_t = cap_start + (cap_end - cap_start) * alpha;
        let support_alpha = if reverse_support { 1.0 - alpha } else { alpha };
        let support_t = support_start + (support_end - support_start) * support_alpha;
        let cap_uv = cap_pcurve.point_at(cap_t);
        let support_uv = support_pcurve.point_at(support_t);
        let cap_point = cap_surface.point_at(cap_uv.x, cap_uv.y);
        let support_point = support_surface.point_at(support_uv.x, support_uv.y);
        let residual = cap_point.distance_to(&support_point);
        if !residual.is_finite() {
            return None;
        }
        if residual > worst.residual {
            worst = LivePcurvePairResidual {
                residual,
                sample_parameter: alpha,
            };
        }
    }
    Some(worst)
}
