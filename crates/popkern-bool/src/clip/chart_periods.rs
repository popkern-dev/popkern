//! Chart periods for point-in-face classification (split per
//! file-size rule; see [`classification_chart_periods`]).

/// Chart periods for point-in-face classification: analytic
/// periodicity, PLUS the geometric closure of a NURBS chart. A lofted
/// or swept closed profile's wall is one G0-closed NURBS surface whose
/// seam pcurves alias in UV, so the even-odd contour is degenerate
/// without period-image handling — every interior point of a
/// propeller-blade wall classified OUTSIDE and its whole hub
/// intersection loop was clipped away (PK-440). Only classification
/// consumers may use this: evaluation is NOT periodic beyond the
/// chart, so `periodicity()` itself stays `None` for NURBS.
pub(crate) fn classification_chart_periods(
    geometry: &popkern_core::geometry::SurfaceGeometry,
) -> (Option<f64>, Option<f64>) {
    use popkern_core::geometry::SurfaceGeometry;
    use popkern_core::geometry::surface::SurfaceEval as _;
    if let SurfaceGeometry::Nurbs(ns) = geometry {
        let ((u0, u1), (v0, v1)) = ns.parameter_range();
        let (nu, nv) = (ns.num_u, ns.num_v);
        if nu < 2 || nv < 2 || ns.control_points.len() != nu * nv {
            return (None, None);
        }
        let cp = |i: usize, j: usize| ns.control_points[i * nv + j];
        let w = |i: usize, j: usize| ns.weights[i * nv + j];
        let diag = {
            let mut lo = cp(0, 0);
            let mut hi = cp(0, 0);
            for p in &ns.control_points {
                lo = popkern_core::geometry::Point3::new(
                    lo.x.min(p.x),
                    lo.y.min(p.y),
                    lo.z.min(p.z),
                );
                hi = popkern_core::geometry::Point3::new(
                    hi.x.max(p.x),
                    hi.y.max(p.y),
                    hi.z.max(p.z),
                );
            }
            lo.distance_to(&hi)
        };
        let eps = 1.0e-9 * (1.0 + diag);
        let u_closed = (0..nv).all(|j| {
            cp(0, j).distance_to(&cp(nu - 1, j)) < eps && (w(0, j) - w(nu - 1, j)).abs() < 1.0e-9
        });
        let v_closed = (0..nu).all(|i| {
            cp(i, 0).distance_to(&cp(i, nv - 1)) < eps && (w(i, 0) - w(i, nv - 1)).abs() < 1.0e-9
        });
        return (u_closed.then_some(u1 - u0), v_closed.then_some(v1 - v0));
    }
    geometry.periodicity()
}
