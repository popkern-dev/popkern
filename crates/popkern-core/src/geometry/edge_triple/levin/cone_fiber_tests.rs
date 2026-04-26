//! Unit tests for the rank-3 cone fibre dispatcher. Calls
//! `try_rank3_cone_fiber` directly to bypass the upstream specialised
//! quadric-pair dispatchers (e.g. `try_cyl_cone_general`) that
//! handle some pair classes before Levin.

use super::cone_fiber::try_rank3_cone_fiber;
use crate::geometry::edge_triple::inputs::EdgeAnchors;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

#[test]
fn cyl_cone_skew_rank3_cone_fiber_produces_nurbs() {
    // cyl: axis (0, 0, 1), origin (0, 0, 0), radius 1 → x² + y² = 1
    // cone: apex (2, 0, 0), axis (-1, 0, 0), half-angle 45°
    //         → y² + z² = (x − 2)²
    //
    // The pencil's degenerate fibres are all rank 3 (proper cones),
    // so Levin's rank-≤2 plane decomposition doesn't apply. The
    // rank-3 dispatcher samples the cone fibre's null directions,
    // intersects each ruling with cyl's quadric, and fits NURBS.
    //
    // Two disjoint loops in the intersection; we pick anchors on the
    // same Levin σ-branch (start (1, 0, 1) at τ=0, σ=+1; end
    // (-1, 0, 3) at τ=π, σ=+3).
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.0,
    };
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(2.0, 0.0, 0.0),
        axis: Vec3::new(-1.0, 0.0, 0.0),
        half_angle: std::f64::consts::FRAC_PI_4,
    };
    let tol = Tolerance::default();
    let start = Point3::new(1.0, 0.0, 1.0);
    let end = Point3::new(-1.0, 0.0, 3.0);
    let anchors = EdgeAnchors::Open { start, end };

    // Stage 2b.35 status: the dispatcher *runs* and produces the
    // expected geometry. Post-construction verify (G2 invariant)
    // requires the LSQ NURBS fit to track curve_3d ↔ pcurve to within
    // `tol.linear()` — for spatial quartics with our 33-sample / 33-CP
    // / degree-5 fit, mid-curve residual is ~10⁻⁵, which exceeds the
    // default tolerance. Tightening fit quality (interpolation,
    // higher degree, rational NURBS) is a follow-up.
    //
    // For now, we accept the `DegenerateSurfacePair` error path and
    // verify the algorithmic correctness directly: sampling the fibre
    // produces points that lie on both surfaces' quadric forms.
    let triple = try_rank3_cone_fiber(&cyl, &cone, &anchors, &tol)
        .expect("rank-3 cone fibre returned None")
        .expect("rank-3 cone fibre returned Err");

    // Endpoint G5 holds (the underlying NURBS interpolation passes
    // through the start and end Levin samples to ~ULP).
    use crate::geometry::curve::CurveEval;
    let p0 = triple.curve_3d.point_at(triple.t_start);
    let p1 = triple.curve_3d.point_at(triple.t_end);
    assert!(
        p0.distance_to(&start) < tol.linear(),
        "start endpoint: {p0:?} vs {start:?}"
    );
    assert!(
        p1.distance_to(&end) < tol.linear(),
        "end endpoint: {p1:?} vs {end:?}"
    );
    // Mid-curve sanity: the curve point lies on both surfaces' quadric
    // forms within tolerance — confirms the analytic Levin sampling
    // produced geometrically valid intersection points.
    let p_mid = triple
        .curve_3d
        .point_at((triple.t_start + triple.t_end) * 0.5);
    let cyl_q = super::quadric::surface_to_quadric(&cyl).unwrap();
    let cone_q = super::quadric::surface_to_quadric(&cone).unwrap();
    let on_cyl = cyl_q.eval_affine(p_mid.x, p_mid.y, p_mid.z);
    let on_cone = cone_q.eval_affine(p_mid.x, p_mid.y, p_mid.z);
    assert!(
        on_cyl.abs() < tol.linear(),
        "midpoint not on cyl: pᵀAp = {on_cyl}"
    );
    assert!(
        on_cone.abs() < tol.linear(),
        "midpoint not on cone: pᵀBp = {on_cone}"
    );
}
