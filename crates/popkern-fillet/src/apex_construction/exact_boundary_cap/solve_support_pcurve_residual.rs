//! Boundary residual checks for post-pcurve support-normal refinement.

use popkern_core::geometry::NurbsSurface;
use popkern_core::geometry::surface::SurfaceEval;

use super::super::pcurves::{PreparedExactBoundaryPcurvePair, PreparedExactBoundaryPcurves};
use super::super::ribbon::{PreparedExactBoundaryRibbon, PreparedExactBoundaryRibbons};
use super::ExactBoundarySolvedSurfaces;

pub(super) fn worst_cap_support_position_residual(
    solved: &ExactBoundarySolvedSurfaces,
    pcurves: &PreparedExactBoundaryPcurves,
    ribbons: &PreparedExactBoundaryRibbons,
) -> f64 {
    pcurves
        .pairs
        .iter()
        .filter_map(|pair| {
            let surface = solved.surfaces.iter().find(|surface| {
                surface.sector_index == pair.sector_index && surface.sector_name == pair.sector_name
            })?;
            let ribbon = ribbon_for_pair(ribbons, pair)?;
            Some(pair_cap_support_position_residual(
                &surface.surface,
                pair,
                ribbon,
            ))
        })
        .fold(0.0, f64::max)
}

fn pair_cap_support_position_residual(
    surface: &NurbsSurface,
    pair: &PreparedExactBoundaryPcurvePair,
    ribbon: &PreparedExactBoundaryRibbon,
) -> f64 {
    ribbon
        .samples
        .iter()
        .map(|sample| {
            let cap_uv = pair.cap_pcurve.point_at(
                pair.cap_pcurve.t_start
                    + (pair.cap_pcurve.t_end - pair.cap_pcurve.t_start) * sample.parameter,
            );
            surface
                .point_at(cap_uv.x, cap_uv.y)
                .distance_to(&sample.support_point)
        })
        .filter(|residual| residual.is_finite())
        .fold(0.0, f64::max)
}

fn ribbon_for_pair<'a>(
    ribbons: &'a PreparedExactBoundaryRibbons,
    pair: &PreparedExactBoundaryPcurvePair,
) -> Option<&'a PreparedExactBoundaryRibbon> {
    ribbons.ribbons.iter().find(|ribbon| {
        ribbon.sector_index == pair.sector_index
            && ribbon.side == pair.side
            && ribbon.edge == pair.edge
            && ribbon.support_half_edge == pair.support_half_edge
    })
}
