//! Stage 2b.47: verification harness for the `try_levin_pair_nested`
//! dispatch path.
//!
//! These tests pin two properties:
//!
//! 1. **Reachability** — the nested path is wired through dispatch.rs
//!    and produces a `Some(Ok(_))` triple on a representative quadric
//!    pair. Without this, the wiring landed in 2b.46 could regress
//!    silently (e.g., if a future refactor moved the helper but
//!    forgot to keep the chain wired).
//! 2. **Correctness** — when the nested path produces a triple, the
//!    triple satisfies G5 (anchor accuracy) and G2 (positional
//!    invariant on both surfaces) at `tol.linear()`. The symbolic
//!    plane normals must round-trip through the f64 boundary
//!    conversion without losing precision below tolerance.
//!
//! Tests target `try_levin_pair_nested` directly (rather than going
//! through `compute_edge_geometry`) to avoid the upstream dispatch
//! short-circuiting via `try_tangent_pair` or `try_levin_pair`. End-
//! to-end coverage through the public entry point lands once we have
//! a CAD-realistic input that *only* the nested path can solve.

use crate::geometry::edge_triple::inputs::EdgeAnchors;
use crate::geometry::edge_triple::levin::dispatch::{try_levin_pair, try_levin_pair_nested};
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

fn z_axis_cyl(radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius,
    }
}

fn x_axis_cyl(radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(1.0, 0.0, 0.0),
        radius,
    }
}

/// Steinmetz: two unit cylinders, axes z and x, meeting at the origin.
/// Pencil `λA + B` at `λ = −1` gives `D = diag(1, 0, −1, 0)`, rank-2
/// mixed-sign with eigenvalues `±1`. Plane normals: `e_x ± e_z`,
/// i.e., the canonical Steinmetz planes `x + z = 0` and `x − z = 0`.
/// The nested path's `Q[√d, √μ_+]` arithmetic with `α = 1` reduces
/// to ordinary `Q` arithmetic (structurally non-trivial — the `√α`
/// component carries the rational — but numerically equivalent).
#[test]
fn nested_path_produces_triple_on_perpendicular_cylinders() {
    let cyl_z = z_axis_cyl(1.0);
    let cyl_x = x_axis_cyl(1.0);
    // Anchor on the plane x + z = 0, on both cylinders. With unit
    // cylinders sharing the origin and perpendicular axes, the
    // intersection has 4 ellipse pieces; pick endpoints of one piece
    // on the upper-half (+y) plane.
    // For y² = 1 − x², z = −x: parametrise x = cos θ, y = sin θ,
    // z = −cos θ. At θ = 0: (1, 0, −1). At θ = π/2: (0, 1, 0).
    let start = Point3::new(1.0, 0.0, -1.0);
    let end = Point3::new(0.0, 1.0, 0.0);
    let anchors = EdgeAnchors::Open { start, end };
    let tol = Tolerance::default();

    // Sanity: anchors satisfy x + z = 0 and lie on both cylinders.
    assert!((start.x + start.z).abs() < 1e-12);
    assert!((end.x + end.z).abs() < 1e-12);
    assert!((start.x.powi(2) + start.y.powi(2) - 1.0).abs() < 1e-12);
    assert!((start.y.powi(2) + start.z.powi(2) - 1.0).abs() < 1e-12);

    let result = try_levin_pair_nested(&cyl_z, &cyl_x, &anchors, &tol)
        .expect("nested path returns Some on perpendicular cylinders");
    let triple = result.expect("nested path returns Ok on perpendicular cylinders");

    // G5: endpoints exact at tol.linear().
    use crate::geometry::curve::CurveEval;
    let p_start = triple.curve_3d.point_at(triple.t_start);
    let p_end = triple.curve_3d.point_at(triple.t_end);
    assert!(
        p_start.distance_to(&start) <= tol.linear(),
        "start: {:?} vs {:?}",
        p_start,
        start
    );
    assert!(
        p_end.distance_to(&end) <= tol.linear(),
        "end: {:?} vs {:?}",
        p_end,
        end
    );

    // G2: positional invariant on both surfaces at every sample.
    let q = triple
        .verify(&cyl_z, &cyl_x, &tol, 32)
        .expect("verify on nested-path triple");
    assert!(
        q.max_positional_residual <= tol.linear(),
        "G2 residual {} > tol {}",
        q.max_positional_residual,
        tol.linear()
    );
}

/// Anchors not on the rank-2 fibre's plane → `try_levin_pair_nested`
/// returns `None` rather than producing a wrong triple. The pencil
/// for two perpendicular cylinders has only two real planes; an
/// anchor pair off both planes can't be fitted.
#[test]
fn nested_path_returns_none_when_anchors_off_plane() {
    let cyl_z = z_axis_cyl(1.0);
    let cyl_x = x_axis_cyl(1.0);
    // Anchors that are NOT on either Steinmetz plane (x = ±z). Pick
    // one on the +z arm of cyl_z that's not at z = ±x.
    let start = Point3::new(1.0, 0.0, 0.5); // x + z = 1.5 ≠ 0; x − z = 0.5 ≠ 0
    let end = Point3::new(0.0, 1.0, 0.5);
    let anchors = EdgeAnchors::Open { start, end };
    let tol = Tolerance::default();

    let result = try_levin_pair_nested(&cyl_z, &cyl_x, &anchors, &tol);
    assert!(
        result.is_none(),
        "nested path must return None when anchors miss every fibre plane"
    );
}

/// Non-quadric input (torus) → `try_levin_pair_nested` returns
/// `None` because `surface_to_quadric_rat` rejects it. The dispatch
/// then falls through to the next chain step.
#[test]
fn nested_path_returns_none_on_non_quadric_input() {
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 5.0,
        minor_radius: 1.0,
    };
    let cyl = z_axis_cyl(1.0);
    let anchors = EdgeAnchors::Open {
        start: Point3::new(1.0, 0.0, 0.0),
        end: Point3::new(0.0, 1.0, 0.0),
    };
    let tol = Tolerance::default();

    assert!(try_levin_pair_nested(&torus, &cyl, &anchors, &tol).is_none());
    assert!(try_levin_pair_nested(&cyl, &torus, &anchors, &tol).is_none());
}

/// **Equivalence**: when both `try_levin_pair` (f64) and
/// `try_levin_pair_nested` (symbolic) succeed on the same input,
/// they must agree on the recovered intersection.
///
/// This test exercises the symbolic-numeric consistency of the
/// nested path: same anchors, same surfaces, same tolerance ⇒ both
/// paths produce a `EdgeTriple` whose curve passes through the same
/// 3-D points (within tolerance) and verifies G2 on both surfaces.
///
/// We do NOT assert the curves are bitwise-identical: f64 Jacobi and
/// symbolic Q[√d, √α] follow different code paths, so the underlying
/// `CurveGeometry` representations may have different parameter
/// orientations or basis conventions. The geometric content is what's
/// invariant.
#[test]
fn nested_and_f64_levin_paths_agree_on_perpendicular_cylinders() {
    let cyl_z = z_axis_cyl(1.0);
    let cyl_x = x_axis_cyl(1.0);
    let start = Point3::new(1.0, 0.0, -1.0);
    let end = Point3::new(0.0, 1.0, 0.0);
    let anchors = EdgeAnchors::Open { start, end };
    let tol = Tolerance::default();

    let f64_triple = try_levin_pair(&cyl_z, &cyl_x, &anchors, &tol)
        .expect("f64 Levin Some")
        .expect("f64 Levin Ok");
    let nested_triple = try_levin_pair_nested(&cyl_z, &cyl_x, &anchors, &tol)
        .expect("nested Levin Some")
        .expect("nested Levin Ok");

    use crate::geometry::curve::CurveEval;
    // Endpoints reach the same anchors — both paths are anchored.
    let f64_start = f64_triple.curve_3d.point_at(f64_triple.t_start);
    let nested_start = nested_triple.curve_3d.point_at(nested_triple.t_start);
    assert!(f64_start.distance_to(&start) <= tol.linear());
    assert!(nested_start.distance_to(&start) <= tol.linear());

    let f64_end = f64_triple.curve_3d.point_at(f64_triple.t_end);
    let nested_end = nested_triple.curve_3d.point_at(nested_triple.t_end);
    assert!(f64_end.distance_to(&end) <= tol.linear());
    assert!(nested_end.distance_to(&end) <= tol.linear());

    // Midpoint sample — both curves should pass through approximately
    // the same point. They're parameterised potentially differently,
    // so we test both midpoints lie on both surfaces.
    let f64_mid = f64_triple
        .curve_3d
        .point_at((f64_triple.t_start + f64_triple.t_end) * 0.5);
    let nested_mid = nested_triple
        .curve_3d
        .point_at((nested_triple.t_start + nested_triple.t_end) * 0.5);

    let on_cyl_z = |p: &Point3| (p.x * p.x + p.y * p.y - 1.0).abs();
    let on_cyl_x = |p: &Point3| (p.y * p.y + p.z * p.z - 1.0).abs();
    assert!(on_cyl_z(&f64_mid) <= tol.linear() && on_cyl_x(&f64_mid) <= tol.linear());
    assert!(on_cyl_z(&nested_mid) <= tol.linear() && on_cyl_x(&nested_mid) <= tol.linear());

    // G2 on both surfaces for both triples.
    let f64_q = f64_triple
        .verify(&cyl_z, &cyl_x, &tol, 32)
        .expect("f64 verify");
    let nested_q = nested_triple
        .verify(&cyl_z, &cyl_x, &tol, 32)
        .expect("nested verify");
    assert!(f64_q.max_positional_residual <= tol.linear());
    assert!(nested_q.max_positional_residual <= tol.linear());
}

/// Cost-bound regression: `try_levin_pair_nested` returns quickly
/// (microseconds) on inputs whose rational quadrics have huge
/// numerators/denominators after f64-roundtrip. Without the bounded
/// `find_small_rational_fiber`, the older `find_rational_degenerate_fiber`
/// would trial-divide through `divisors_bigint` to ~2³² and hang.
/// We assert "returns within wall-clock 200 ms" as a smoke test.
#[test]
fn nested_path_does_not_hang_on_irrational_axis() {
    // Cylinder with axis (1/√2, 0, 1/√2) — irrational axis components,
    // f64-roundtrip gives ~50-bit numerators in the rational quadric.
    let inv_root2 = 1.0 / 2.0_f64.sqrt();
    let cyl_skew = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(inv_root2, 0.0, inv_root2),
        radius: 1.0,
    };
    let cyl_z = z_axis_cyl(1.0);
    let anchors = EdgeAnchors::Open {
        start: Point3::new(1.0, 0.0, 0.0),
        end: Point3::new(0.0, 1.0, 0.0),
    };
    let tol = Tolerance::default();

    let t0 = std::time::Instant::now();
    let _ = try_levin_pair_nested(&cyl_z, &cyl_skew, &anchors, &tol);
    let elapsed = t0.elapsed();
    assert!(
        elapsed.as_millis() < 200,
        "try_levin_pair_nested took {} ms — bounded fiber search regressed",
        elapsed.as_millis()
    );
}
