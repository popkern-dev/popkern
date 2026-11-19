//! Tests for `interval_curvature` (split out of the module file).
use super::*;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};

#[test]
fn sphere_bounds_contain_exact_curvature() {
    let r = 2.0;
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: r,
        axis: crate::geometry::Vec3::Z,
    };
    let exact_k = 1.0 / (r * r);

    // Sample away from poles.
    let bounds = guaranteed_curvature_bounds(
        &sphere,
        (0.5, 2.5),
        (
            std::f64::consts::FRAC_PI_4,
            3.0 * std::f64::consts::FRAC_PI_4,
        ),
        0.1,
        3,
    );

    assert!(
        bounds.gaussian.contains(exact_k),
        "Gaussian bounds [{}, {}] should contain exact K={exact_k}",
        bounds.gaussian.lo,
        bounds.gaussian.hi
    );
}

#[test]
fn plane_bounds_contain_zero() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };

    let bounds = guaranteed_curvature_bounds(&plane, (0.0, 1.0), (0.0, 1.0), 0.01, 3);

    assert!(
        bounds.gaussian.contains(0.0),
        "Plane Gaussian bounds [{}, {}] should contain K=0",
        bounds.gaussian.lo,
        bounds.gaussian.hi
    );
    assert!(
        bounds.mean.contains(0.0),
        "Plane mean bounds [{}, {}] should contain H=0",
        bounds.mean.lo,
        bounds.mean.hi
    );
}

#[test]
fn subdivision_tightens_bounds() {
    let r = 3.0;
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: r,
        axis: crate::geometry::Vec3::Z,
    };

    let coarse = guaranteed_curvature_bounds(
        &sphere,
        (0.5, 2.5),
        (
            std::f64::consts::FRAC_PI_4,
            3.0 * std::f64::consts::FRAC_PI_4,
        ),
        10.0, // Very loose target -> no subdivision
        0,
    );

    let fine = guaranteed_curvature_bounds(
        &sphere,
        (0.5, 2.5),
        (
            std::f64::consts::FRAC_PI_4,
            3.0 * std::f64::consts::FRAC_PI_4,
        ),
        0.01,
        4,
    );

    // Fine bounds should be tighter (or equal) than coarse.
    assert!(
        fine.gaussian_width <= coarse.gaussian_width + 1e-10,
        "Fine width {} should be <= coarse width {}",
        fine.gaussian_width,
        coarse.gaussian_width
    );
}

#[test]
fn cylinder_bounds_contain_zero_gaussian() {
    let r = 5.0;
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        radius: r,
    };

    let bounds =
        guaranteed_curvature_bounds(&cyl, (0.0, std::f64::consts::TAU), (-1.0, 1.0), 0.1, 3);

    // Cylinder has K = 0 everywhere.
    assert!(
        bounds.gaussian.contains(0.0),
        "Cylinder Gaussian bounds [{}, {}] should contain K=0",
        bounds.gaussian.lo,
        bounds.gaussian.hi
    );
}

#[test]
fn guaranteed_bounds_results_have_fields() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };

    let bounds = guaranteed_curvature_bounds(&plane, (0.0, 1.0), (0.0, 1.0), 0.01, 2);

    assert!(bounds.gaussian_width >= 0.0);
    assert!(bounds.gaussian.lo <= bounds.gaussian.hi);
    assert!(bounds.mean.lo <= bounds.mean.hi);
}

/// PK-510 counterexample: the retired sampled-Lipschitz estimate
/// (five curvature samples per cell, Lipschitz constant from their
/// differences, 2× safety) against the rigorous enclosure, on a torus
/// strip over the minor angle `v ∈ [−π/2, 3π/2]`.
/// `K = cos v / (r (R + r cos v))` is 0 at the root cell's corners and
/// centre (v = −π/2, π/2, 3π/2). The estimate therefore collapses to
/// [0, 0], while the true curvature spans `−1/(r(R−r))` to
/// `1/(r(R+r))`. The rigorous bound contains both extremes.
#[test]
fn the_rigorous_bound_contains_extremes_the_sampled_estimate_missed() {
    use crate::healing::curvature::surface_curvature::gaussian_curvature;
    let (major, minor) = (2.0, 0.5);
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: major,
        minor_radius: minor,
    };
    let (u_range, v_range) = (
        (0.0, 1.0),
        (-std::f64::consts::FRAC_PI_2, 1.5 * std::f64::consts::PI),
    );
    let (true_hi, true_lo) = (
        1.0 / (minor * (major + minor)),
        -1.0 / (minor * (major - minor)),
    );
    // The retired estimate on the root cell.
    let (um, vm) = (0.5 * (u_range.0 + u_range.1), 0.5 * (v_range.0 + v_range.1));
    let points = [
        (u_range.0, v_range.0),
        (u_range.1, v_range.0),
        (u_range.0, v_range.1),
        (u_range.1, v_range.1),
        (um, vm),
    ];
    let samples: Vec<(f64, f64, f64)> = points
        .iter()
        .filter_map(|&(u, v)| gaussian_curvature(&torus, u, v).map(|k| (u, v, k)))
        .collect();
    let mut lipschitz = 0.0_f64;
    for a in &samples {
        for b in &samples {
            let d = ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt();
            if d > 1e-15 {
                lipschitz = lipschitz.max((a.2 - b.2).abs() / d);
            }
        }
    }
    let diagonal = ((u_range.1 - u_range.0).powi(2) + (v_range.1 - v_range.0).powi(2)).sqrt();
    let slack = 2.0 * lipschitz * diagonal;
    let sampled_hi = samples
        .iter()
        .map(|s| s.2)
        .fold(f64::NEG_INFINITY, f64::max)
        + slack;
    let sampled_lo = samples.iter().map(|s| s.2).fold(f64::INFINITY, f64::min) - slack;
    assert!(
        true_hi > sampled_hi && true_lo < sampled_lo,
        "negative control: the sampled estimate [{sampled_lo:e}, {sampled_hi:e}] misses [{true_lo}, {true_hi}]"
    );
    let bounds = guaranteed_curvature_bounds(&torus, u_range, v_range, 0.5, 6);
    assert_eq!(bounds.degenerate_cells, 0);
    assert!(
        bounds.gaussian.lo <= true_lo && true_hi <= bounds.gaussian.hi,
        "rigorous [{}, {}] must contain [{true_lo}, {true_hi}]",
        bounds.gaussian.lo,
        bounds.gaussian.hi
    );
}

/// A chart singularity (the sphere pole) cannot be enclosed; the
/// bound is then the whole line, never a guess.
#[test]
fn a_chart_singularity_reads_unbounded() {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
        axis: crate::geometry::Vec3::Z,
    };
    let pole = crate::geometry::surface::SurfaceEval::parameter_range(&sphere).1;
    let bounds = guaranteed_curvature_bounds(&sphere, (0.0, 1.0), (pole.1 - 0.3, pole.1), 0.1, 3);
    assert!(bounds.degenerate_cells > 0);
    assert!(bounds.gaussian.lo.is_infinite() && bounds.gaussian.hi.is_infinite());
}
