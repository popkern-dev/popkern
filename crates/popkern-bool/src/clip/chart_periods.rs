//! Chart periods for point-in-face classification (split per
//! file-size rule; see [`classification_chart_periods`]).

/// Chart periods for point-in-face classification: analytic
/// periodicity, PLUS the geometric closure of a NURBS chart. A lofted
/// or swept closed profile's wall is one G0-closed NURBS surface whose
/// seam pcurves alias in UV, so the even-odd contour is degenerate
/// without period-image handling — every interior point of a
/// propeller-blade wall classified OUTSIDE and its whole hub
/// intersection loop was clipped away (PK-440). Only chart-space
/// consumers may use this: evaluation is NOT periodic beyond the
/// chart, so `periodicity()` itself stays `None` for NURBS. The
/// detection lives in `popkern_core::geometry::chart_periods`, shared
/// with the marching wrap and the fragment walk (PK-443).
pub(crate) fn classification_chart_periods(
    geometry: &popkern_core::geometry::SurfaceGeometry,
) -> (Option<f64>, Option<f64>) {
    popkern_core::geometry::chart_periods::geometric_chart_periods(geometry)
}
