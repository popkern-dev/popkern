//! Gates for the variational cap. Fixtures are plain geometry (a
//! spherical n-gon whose neighbour is the sphere itself), so the
//! expected cap is known: the sphere continued across the opening.

use super::chart::CharacteristicMap;
use super::linalg::lstsq;
use super::space::{SectorSpace, apply};
use super::spoke::spoke_rows;
use super::*;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::{Point3, Vec3};
use faer::Mat;

/// Rational quadratic great-circle arc between two unit vectors.
fn great_arc(a: Vec3, b: Vec3) -> NurbsCurve {
    let half = 0.5 * a.dot(&b).clamp(-1.0, 1.0).acos();
    let w = half.cos();
    let mid = (a + b) * (0.5 / (w * w));
    NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(a.x, a.y, a.z),
            Point3::new(mid.x, mid.y, mid.z),
            Point3::new(b.x, b.y, b.z),
        ],
        weights: vec![1.0, w, 1.0],
    }
}

/// A spherical n-gon at polar angle `theta` about +z on the unit
/// sphere, sides counter-clockwise about the outward normal.
fn spherical_ngon(n: usize, theta: f64) -> Vec<VariationalCapSide> {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let corner = |k: usize| {
        let a = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
        Vec3::new(theta.sin() * a.cos(), theta.sin() * a.sin(), theta.cos())
    };
    // Orient: the sphere's natural normal against the outward radius.
    let p = corner(0);
    let (u, v) = sphere
        .closest_point(&Point3::new(p.x, p.y, p.z))
        .expect("projection");
    let reversed = sphere.normal_at(u, v).dot(&p) < 0.0;
    (0..n)
        .map(|k| VariationalCapSide {
            curve: great_arc(corner(k), corner(k + 1)),
            neighbour: sphere.clone(),
            neighbour_reversed: reversed,
        })
        .collect()
}

#[test]
fn pulled_back_polynomials_satisfy_the_spoke_rows() {
    for n in 3..=8 {
        for knots_inside in [false, true] {
            let mut space = SectorSpace::bezier(n, 7);
            if knots_inside {
                let mut k = space.knots.clone();
                k.insert(8, 0.3);
                k.insert(9, 0.65);
                space = SectorSpace::new(n, 7, k);
            }
            let phi = CharacteristicMap::new(n);
            let x = pulled_back_cubic(&space, &phi);
            let rows = spoke_rows(&space, phi.c);
            let worst = rows
                .iter()
                .flat_map(|r| x.iter().map(move |c| apply(r, c).abs()))
                .fold(0.0, f64::max);
            assert!(worst < 1e-10, "n {n} knots {knots_inside}: {worst:.3e}");
        }
    }
}

/// Coefficients of `F(r_n^s φ_n(u, v))` per sector, `F` a fixed cubic
/// map R² → R³, fitted exactly by least squares (the pulled-back cubic
/// is bicubic, inside the bi-7 space).
fn pulled_back_cubic(space: &SectorSpace, phi: &CharacteristicMap) -> [Vec<f64>; 3] {
    let n = space.n;
    let m = space.m;
    let c = phi.c;
    let s = (2.0 * std::f64::consts::PI / n as f64).sin();
    let tau = (std::f64::consts::PI / n as f64).tan();
    let f = |x: f64, y: f64| {
        [
            x + 0.3 * y * y,
            y - 0.2 * x * y + 0.1 * x * x * x,
            0.7 * x * x - 0.4 * x * y + 0.9 * y * y + 0.2 * y * y * y,
        ]
    };
    let grid = 2 * m;
    let mut out: [Vec<f64>; 3] = std::array::from_fn(|_| vec![0.0; space.total()]);
    for sector in 0..n {
        let rot = 2.0 * std::f64::consts::PI * sector as f64 / n as f64;
        let mut a = Mat::<f64>::zeros((grid + 1) * (grid + 1), m * m);
        let mut b = Mat::<f64>::zeros((grid + 1) * (grid + 1), 3);
        for i in 0..=grid {
            for j in 0..=grid {
                let (u, v) = (i as f64 / grid as f64, j as f64 / grid as f64);
                let row = i * (grid + 1) + j;
                let fx = u + c * v * (1.0 - u);
                let fy = v * (s * (1.0 - u) + tau * u);
                let (x, y) = (
                    rot.cos() * fx - rot.sin() * fy,
                    rot.sin() * fx + rot.cos() * fy,
                );
                let val = f(x, y);
                for (k, vk) in val.iter().enumerate() {
                    b.write(row, k, *vk);
                }
                for &(idx, w) in &space.functional(sector, u, v, 0, 0) {
                    a.write(row, idx - sector * m * m, w);
                }
            }
        }
        let sol = lstsq(a, b, 1e-13).expect("full rank");
        for (k, coord) in out.iter_mut().enumerate() {
            for l in 0..m * m {
                coord[sector * m * m + l] = sol.read(l, k);
            }
        }
    }
    out
}

#[test]
fn spherical_ngons_meet_the_targets_with_exact_spokes() {
    let targets = VariationalCapTargets::g2(1e-6, 1e-4, 5e-3);
    for n in [3, 4, 5, 6, 8] {
        let sides = spherical_ngon(n, 0.6);
        let cap = fit_variational_cap(&sides, &targets, &VariationalCapOptions::default(), 1e-7)
            .unwrap_or_else(|e| panic!("n {n}: {e}"));
        let r = &cap.report;
        eprintln!(
            "n {n}: levels {} unknowns {} gn {} worst rim {:?} spoke {:?} reg {:.3e}",
            r.refinement_levels,
            r.unknowns,
            r.gauss_newton_iterations,
            r.worst_rim(&targets),
            r.worst_spoke,
            r.min_regularity
        );
        assert!(r.meets(&targets), "n {n}");
        assert!(
            r.worst_spoke.0 < 1e-10 && r.worst_spoke.1 < 1e-9,
            "n {n}: {:?}",
            r.worst_spoke
        );
        assert!(r.worst_spoke.2 < 1e-6, "n {n}: {:?}", r.worst_spoke);
        // The cap stays near the sphere it continues. The bending
        // energy has the spherical cap as a critical point; the residual
        // sag (2e-3 at n = 4 to 8e-3 at n = 6, 2026-09-28) comes from
        // the parameterisation regulariser. The coordinate-wise
        // parametric thin plate sagged 6e-3 to 1.1e-1.
        for sector in &cap.sectors {
            for i in 0..=8 {
                for j in 0..=8 {
                    let p = sector.point_at(i as f64 / 8.0, j as f64 / 8.0);
                    let radius = Vec3::new(p.x, p.y, p.z).length();
                    assert!((radius - 1.0).abs() < 1e-2, "n {n}: radius {radius}");
                }
            }
        }
    }
}

#[test]
fn an_unreachable_target_declines_with_the_worst_sample() {
    let sides = spherical_ngon(5, 0.6);
    let targets = VariationalCapTargets::g2(1e-14, 1e-14, 1e-14);
    let options = VariationalCapOptions {
        max_refinement_levels: 1,
        ..Default::default()
    };
    match fit_variational_cap(&sides, &targets, &options, 1e-7) {
        Err(VariationalCapError::BudgetExhausted { worst, levels, .. }) => {
            assert!(levels <= 1);
            assert!(worst.parameter.is_finite());
        }
        other => panic!("expected BudgetExhausted, got {other:?}"),
    }
}

#[test]
fn an_open_boundary_and_a_degenerate_side_decline() {
    let mut sides = spherical_ngon(4, 0.6);
    let targets = VariationalCapTargets::g1(1e-6, 1e-4);
    sides[1].curve.control_points[0].x += 1e-3;
    assert!(matches!(
        fit_variational_cap(&sides, &targets, &VariationalCapOptions::default(), 1e-7),
        Err(VariationalCapError::BoundaryNotClosed { side: 0, .. })
    ));
    let sides = spherical_ngon(2, 0.6);
    assert!(matches!(
        fit_variational_cap(&sides, &targets, &VariationalCapOptions::default(), 1e-7),
        Err(VariationalCapError::UnsupportedSideCount { got: 2 })
    ));
}

/// Worst radial deviation of a cap from the unit sphere.
fn sphere_sag(cap: &VariationalCap) -> f64 {
    let mut sag = 0.0f64;
    for sector in &cap.sectors {
        for i in 0..=8 {
            for j in 0..=8 {
                let p = sector.point_at(i as f64 / 8.0, j as f64 / 8.0);
                sag = sag.max((Vec3::new(p.x, p.y, p.z).length() - 1.0).abs());
            }
        }
    }
    sag
}

/// E3 (2026-09-28): the coordinate-wise parametric thin plate is not
/// geometric under the non-affine characteristic map, and a spherical
/// n-gon's interior sagged 6e-3 (n = 4) to 1.1e-1 (n = 8) of the
/// radius. The bending energy (frozen reference normal) has the
/// spherical cap as a critical point: measured 2e-3 to 8e-3.
#[test]
fn the_bending_energy_follows_the_sphere_the_parametric_one_sags() {
    let targets = VariationalCapTargets::g2(1e-6, 1e-4, 5e-3);
    for n in [5, 8] {
        let sides = spherical_ngon(n, 0.6);
        let bending =
            fit_variational_cap(&sides, &targets, &VariationalCapOptions::default(), 1e-7)
                .expect("bending fit");
        let parametric = fit_variational_cap(
            &sides,
            &targets,
            &VariationalCapOptions {
                energy: CapEnergy::CharacteristicThinPlate,
                ..Default::default()
            },
            1e-7,
        )
        .expect("parametric fit");
        let (b, p) = (sphere_sag(&bending), sphere_sag(&parametric));
        assert!(b < 1e-2, "n {n}: bending sag {b:.3e}");
        assert!(b < 0.25 * p, "n {n}: bending {b:.3e} vs parametric {p:.3e}");
    }
}
