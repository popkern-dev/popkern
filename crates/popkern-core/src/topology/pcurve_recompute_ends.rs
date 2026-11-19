//! An edge end at a singular row of its face's chart, for the pcurve
//! recompute's free-azimuth snap (split out of `pcurve_recompute.rs` for
//! the 500-line ratchet).

use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

/// Whether a chart point lies at a singular row collapsing u, asked of the
/// chart service in space within a reach of 1e-4 of the edge's interior
/// `|S_u|` (its distance from the axis on a surface of revolution: a
/// projected end lands ~1e-6 of it off the row, PK-860), never under the
/// linear tolerance. `quarter_samples` are the edge's interior samples
/// `(t, u, v)`. Read by the linear tolerance alone, a ball's edge ending
/// 1e-4 off its pole kept an arbitrary azimuth and the loop tore.
pub(super) fn end_at_singular_row(
    surface: &SurfaceGeometry,
    tolerance: &Tolerance,
    quarter_samples: &[(f64, f64, f64)],
) -> impl Fn(f64, f64) -> bool {
    let du_in = quarter_samples
        .iter()
        .map(|&(_, u, v)| surface.partial_u(u, v).length())
        .fold(0.0_f64, f64::max);
    let reach = (1.0e-4 * du_in).max(tolerance.linear());
    let descriptor = crate::geometry::chart::ChartDescriptor::of(surface, tolerance);
    let surface = surface.clone();
    move |u, v| {
        crate::geometry::chart::ops::singular_point_near(
            &descriptor,
            &surface.point_at(u, v),
            reach,
        )
        .is_some_and(|(_, collapses_u)| collapses_u)
    }
}
