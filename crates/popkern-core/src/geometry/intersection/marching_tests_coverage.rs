//! Additional coverage tests for the marching module.
//!
//! Split from `marching_tests.rs` to keep each file under 500 lines.
#[cfg(test)]
use super::*;
use crate::geometry::intersection::degeneracy::{
    DefaultDegeneracyHandler, RobustDegeneracyHandler,
};
use crate::geometry::surface::SurfaceEval;
use std::f64::consts::TAU;

// ---------------------------------------------------------------------------
// Helper: fully-degenerate NURBS surface (all control points at the same spot).
// Its normal_at is zero everywhere — useful for singularity-path testing.
// ---------------------------------------------------------------------------
fn degenerate_nurbs_surface(origin: Point3) -> SurfaceGeometry {
    use crate::geometry::nurbs_surface::NurbsSurface;
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        // All four control points at the same location ⇒ partial_u = partial_v = 0.
        control_points: vec![origin, origin, origin, origin],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0, 1.0, 1.0, 1.0],
    })
}

/// `estimate_normal_near_singularity` should return `None` when all four
/// offset positions also produce degenerate (zero) normals.
///
/// A fully-degenerate NURBS surface (all control points at the origin) has
/// zero normal everywhere, so even the offset samples are degenerate.
#[test]
fn estimate_normal_near_singularity_all_degenerate_returns_none() {
    let surf = degenerate_nurbs_surface(Point3::new(0.0, 0.0, 0.0));
    let result = super::estimate_normal_near_singularity(&surf, 0.5, 0.5, 0.1);
    assert!(
        result.is_none(),
        "all-degenerate surface should return None from estimate_normal_near_singularity"
    );
}

/// `tangent_basis` should not panic on a zero-length normal; it falls back to
/// `Vec3::new(0,0,1)` via `unwrap_or` and returns two perpendicular tangent vectors.
#[test]
fn tangent_basis_with_zero_normal_uses_fallback() {
    let zero = Vec3::new(0.0, 0.0, 0.0);
    // Should not panic — uses the unwrap_or fallback in normalized().
    let (tu, tv) = super::tangent_basis(&zero);
    // The fallback normal is (0,0,1); tu and tv should still be finite unit-ish vectors.
    assert!(tu.length() > 0.5, "tu should be non-zero after fallback");
    assert!(tv.length() > 0.5, "tv should be non-zero after fallback");
}

/// `find_tangent_march_direction` should return `None` when the normal at the
/// query point is degenerate (< 1e-12 length).
///
/// The degenerate NURBS surface always has a zero normal, so the early-exit
/// guard `if n1.length() < 1e-12 { return None; }` is exercised.
#[test]
fn find_tangent_march_direction_degenerate_normal_returns_none() {
    let degenerate = degenerate_nurbs_surface(Point3::new(0.0, 0.0, 0.0));
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let point = Point3::new(0.0, 0.0, 0.0);
    let tol = Tolerance::default().with_linear(1e-3);

    let dir =
        super::find_tangent_march_direction(&degenerate, 0.5, 0.5, &plane, 0.0, 0.0, &point, &tol);
    assert!(
        dir.is_none(),
        "degenerate surface normal should make find_tangent_march_direction return None"
    );
}

/// Marching from a point that immediately hits a singularity on surface 1 with
/// the `DefaultDegeneracyHandler` (which returns `true` ⇒ stop) should break
/// after the first step and return only the start sample.
///
/// The fully-degenerate NURBS surface has zero normal everywhere, so the
/// `n1.length() < epsilon` branch is taken on the very first iteration and
/// `handler.on_singularity` returns `true`, triggering the `break`.
#[test]
fn march_singularity_default_handler_stops_immediately() {
    let degenerate = degenerate_nurbs_surface(Point3::new(0.0, 0.0, 0.0));
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let start = MarchSample {
        u1: 0.5,
        v1: 0.5,
        u2: 0.0,
        v2: 0.0,
        point: Point3::new(0.0, 0.0, 0.0),
    };
    let tol = Tolerance::default().with_linear(1e-3);

    let (samples, is_closed) = march(
        &degenerate,
        (0.0, 1.0),
        (0.0, 1.0),
        &plane,
        (-2.0, 2.0),
        (-2.0, 2.0),
        &start,
        &tol,
        &DefaultDegeneracyHandler,
    );
    // DefaultDegeneracyHandler stops at singularity — only the start point is kept.
    assert_eq!(
        samples.len(),
        1,
        "DefaultDegeneracyHandler should stop at first singularity; got {} samples",
        samples.len()
    );
    assert!(!is_closed);
}

/// Marching from a point on a fully-degenerate NURBS surface with the
/// `RobustDegeneracyHandler` exercises the singularity retry loop: the handler
/// returns `false`, `estimate_normal_near_singularity` also returns `None`
/// (all offsets are degenerate), the halved-step retry loop exhausts all 3
/// attempts, and the `!recovered` branch breaks out of the march loop.
///
/// The march should terminate gracefully with only the start sample rather
/// than panicking.
#[test]
fn march_singularity_retry_loop_exhausted() {
    let degenerate = degenerate_nurbs_surface(Point3::new(0.0, 0.0, 0.0));
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let start = MarchSample {
        u1: 0.5,
        v1: 0.5,
        u2: 0.0,
        v2: 0.0,
        point: Point3::new(0.0, 0.0, 0.0),
    };
    let tol = Tolerance::default().with_linear(1e-3);

    let (samples, is_closed) = march(
        &degenerate,
        (0.0, 1.0),
        (0.0, 1.0),
        &plane,
        (-2.0, 2.0),
        (-2.0, 2.0),
        &start,
        &tol,
        &RobustDegeneracyHandler,
    );
    // Retry loop exhausted → only the start sample is kept.
    assert_eq!(
        samples.len(),
        1,
        "exhausted singularity retry should produce only the start sample; got {}",
        samples.len()
    );
    assert!(!is_closed);
}

/// Same as `march_singularity_retry_loop_exhausted` but with the degenerate
/// surface as surface 2 (not surface 1).  This exercises the second singularity
/// block in `march_one_direction` (lines for surface 2 normal check).
#[test]
fn march_singularity_surface2_retry_loop_exhausted() {
    let degenerate = degenerate_nurbs_surface(Point3::new(0.0, 0.0, 0.0));
    // Surface 1 is a normal plane so its normal is well-defined.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let start = MarchSample {
        u1: 0.0,
        v1: 0.0,
        u2: 0.5,
        v2: 0.5,
        point: Point3::new(0.0, 0.0, 0.0),
    };
    let tol = Tolerance::default().with_linear(1e-3);

    let (samples, is_closed) = march(
        &plane,
        (-2.0, 2.0),
        (-2.0, 2.0),
        &degenerate,
        (0.0, 1.0),
        (0.0, 1.0),
        &start,
        &tol,
        &RobustDegeneracyHandler,
    );
    assert_eq!(
        samples.len(),
        1,
        "exhausted singularity retry on surface 2 should produce only the start sample; got {}",
        samples.len()
    );
    assert!(!is_closed);
}

/// Same as `march_singularity_default_handler_stops_immediately` but with the
/// degenerate surface as surface 2.  Exercises the second singularity block
/// with the `DefaultDegeneracyHandler`.
#[test]
fn march_singularity_surface2_default_handler_stops() {
    let degenerate = degenerate_nurbs_surface(Point3::new(0.0, 0.0, 0.0));
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let start = MarchSample {
        u1: 0.0,
        v1: 0.0,
        u2: 0.5,
        v2: 0.5,
        point: Point3::new(0.0, 0.0, 0.0),
    };
    let tol = Tolerance::default().with_linear(1e-3);

    let (samples, is_closed) = march(
        &plane,
        (-2.0, 2.0),
        (-2.0, 2.0),
        &degenerate,
        (0.0, 1.0),
        (0.0, 1.0),
        &start,
        &tol,
        &DefaultDegeneracyHandler,
    );
    assert_eq!(
        samples.len(),
        1,
        "DefaultDegeneracyHandler should stop at singularity on surface 2; got {} samples",
        samples.len()
    );
    assert!(!is_closed);
}

/// Two locally-tangent planes (same normal) marching with the
/// `DefaultDegeneracyHandler` should stop at the tangency (`on_tangent_during_march`
/// returns `true`).
///
/// This exercises the `handler.on_tangent_during_march` → `break` path
/// (line 302–303 in marching.rs).
#[test]
fn march_tangent_during_march_default_handler_stops() {
    // Two identical planes: normals are parallel so N1×N2 = 0 immediately.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let start = MarchSample {
        u1: 0.0,
        v1: 0.0,
        u2: 0.0,
        v2: 0.0,
        point: Point3::new(0.0, 0.0, 0.0),
    };
    let tol = Tolerance::default().with_linear(1e-3);

    let (samples, is_closed) = march(
        &plane,
        (-2.0, 2.0),
        (-2.0, 2.0),
        &plane,
        (-2.0, 2.0),
        (-2.0, 2.0),
        &start,
        &tol,
        &DefaultDegeneracyHandler,
    );
    // DefaultDegeneracyHandler stops at tangency immediately — only the start.
    assert_eq!(
        samples.len(),
        1,
        "DefaultDegeneracyHandler should stop at tangent march; got {} samples",
        samples.len()
    );
    assert!(!is_closed);
}

/// Two identical planes with `RobustDegeneracyHandler`: the tangent march path
/// invokes `find_tangent_march_direction`.  Since the normals are parallel but
/// well-defined (non-degenerate), the perturbation finder returns a direction
/// and the march proceeds until Newton projection fails (the two coincident
/// planes have no unique intersection point), exhausting retries.
///
/// This exercises the `RobustDegeneracyHandler.on_tangent_during_march` → `false`
/// → `find_tangent_march_direction` → `Some(dir)` path.
#[test]
fn march_tangent_during_march_robust_handler_recovers() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let start = MarchSample {
        u1: 0.0,
        v1: 0.0,
        u2: 0.0,
        v2: 0.0,
        point: Point3::new(0.0, 0.0, 0.0),
    };
    let tol = Tolerance::default().with_linear(1e-3);

    // Should not panic regardless of how many samples are produced.
    let (samples, _is_closed) = march(
        &plane,
        (-2.0, 2.0),
        (-2.0, 2.0),
        &plane,
        (-2.0, 2.0),
        (-2.0, 2.0),
        &start,
        &tol,
        &RobustDegeneracyHandler,
    );
    // At minimum the start sample is always included.
    assert!(
        !samples.is_empty(),
        "march should always include at least the start sample"
    );
}

/// Marching on a torus-plane intersection exercises the periodic wrapping in
/// both u and v (the torus is periodic in both directions) and produces a closed
/// curve, verifying the closure detection path inside `march_one_direction`.
#[test]
fn march_torus_plane_closed_curve() {
    // Torus with major radius 2, minor radius 0.5, axis along Z.
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 2.0,
        minor_radius: 0.5,
    };
    // Horizontal plane at z=0 cuts through the outer equator of the torus.
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    // At torus parameters (u=0, v=0): outer equator point is (2.5, 0, 0).
    let p0 = torus.point_at(0.0, 0.0);
    let start = MarchSample {
        u1: 0.0,
        v1: 0.0,
        u2: p0.x,
        v2: p0.y,
        point: p0,
    };
    let tol = Tolerance::default().with_linear(1e-3);
    let (samples, is_closed) = march(
        &torus,
        (0.0, TAU),
        (0.0, TAU),
        &plane,
        (-4.0, 4.0),
        (-4.0, 4.0),
        &start,
        &tol,
        &DefaultDegeneracyHandler,
    );
    assert!(
        samples.len() > 3,
        "torus-plane intersection should produce multiple samples, got {}",
        samples.len()
    );
    assert!(
        is_closed,
        "torus-plane intersection should be a closed curve"
    );
    // All samples should be approximately at z=0.
    for s in &samples {
        assert!(
            s.point.z.abs() < 0.05,
            "torus-plane intersection point should be at z~0, got z={}",
            s.point.z
        );
    }
}

/// `march_one_direction` closes the curve early (before the outer `march` closure
/// check) when the new sample returns to within `tolerance.linear` of the start.
///
/// A small torus with tight tolerance exercises this inner closure detection path:
/// the inner march detects the return while `samples.len() > 3` within
/// `march_one_direction` itself and pushes the closing sample then returns early.
#[test]
fn march_one_direction_inner_closure_detection() {
    // Use a small torus so the march wraps around quickly.
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 1.0,
        minor_radius: 0.3,
    };
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let p0 = torus.point_at(0.0, 0.0);
    let start = MarchSample {
        u1: 0.0,
        v1: 0.0,
        u2: p0.x,
        v2: p0.y,
        point: p0,
    };
    // Use a slightly coarser tolerance so the march completes the full circle
    // quickly and the inner closure check fires.
    let tol = Tolerance::default().with_linear(5e-3);
    let (samples, is_closed) = march(
        &torus,
        (0.0, TAU),
        (0.0, TAU),
        &plane,
        (-3.0, 3.0),
        (-3.0, 3.0),
        &start,
        &tol,
        &DefaultDegeneracyHandler,
    );
    assert!(
        samples.len() > 3,
        "should produce multiple samples, got {}",
        samples.len()
    );
    assert!(
        is_closed,
        "small torus equator should be detected as closed"
    );
}
