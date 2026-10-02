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
pub(super) fn spherical_ngon(n: usize, theta: f64) -> Vec<VariationalCapSide> {
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
        // The cap stays near the sphere it continues (the energy is zero
        // on it); the octagon's strongly non-affine chart leaves the
        // constrained degree-7 space 1.7e-2 of sag (measured).
        for sector in &cap.sectors {
            for i in 0..=8 {
                for j in 0..=8 {
                    let p = sector.point_at(i as f64 / 8.0, j as f64 / 8.0);
                    let radius = Vec3::new(p.x, p.y, p.z).length();
                    assert!(
                        (radius - 1.0).abs() < if n == 8 { 2.5e-2 } else { 5e-3 },
                        "n {n}: radius {radius}"
                    );
                }
            }
        }
    }
}

#[test]
fn an_unreachable_target_returns_the_best_level_flagged() {
    let sides = spherical_ngon(5, 0.6);
    let targets = VariationalCapTargets::g2(1e-14, 1e-14, 1e-14);
    let options = VariationalCapOptions {
        max_refinement_levels: 1,
        ..Default::default()
    };
    let cap = fit_variational_cap(&sides, &targets, &options, 1e-7).expect("best effort");
    assert!(!cap.report.targets_met);
    assert!(!cap.report.meets(&targets));
    let worst = cap.report.worst_rim(&targets).expect("rim");
    assert!(worst.parameter.is_finite());
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

/// E3 (2026-09-28/29): the coordinate-wise parametric thin plate is not
/// geometric under the non-affine characteristic map, and a spherical
/// n-gon's interior sagged 6e-3 (n = 4) to 1.1e-1 (n = 8) of the
/// radius; a frozen-normal bending energy was not stationary at the
/// sphere either (5e-3..3e-2). The mean-curvature-variation energy is
/// zero on the sphere: measured 1.8e-4 (n = 4) and 6.9e-4 (n = 5).
#[test]
fn mean_curvature_variation_follows_the_sphere() {
    let targets = VariationalCapTargets::g2(1e-6, 1e-4, 5e-3);
    for (n, bar) in [(4, 5e-4), (5, 2e-3)] {
        let sides = spherical_ngon(n, 0.6);
        let fair = fit_variational_cap(&sides, &targets, &VariationalCapOptions::default(), 1e-7)
            .expect("fit");
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
        let (f, p) = (sphere_sag(&fair), sphere_sag(&parametric));
        assert!(f < bar, "n {n}: sag {f:.3e}");
        assert!(f < 0.25 * p, "n {n}: {f:.3e} vs parametric {p:.3e}");
    }
}

#[test]
fn warm_start_transfer_reproduces_the_previous_level_exactly() {
    let coarse = SectorSpace::bezier(5, 7);
    let total = coarse.total();
    let coords: [Vec<f64>; 3] = [
        (0..total).map(|k| (k as f64 * 0.31).sin()).collect(),
        (0..total).map(|k| (k as f64 * 0.17).cos()).collect(),
        (0..total).map(|k| (k as f64 * 0.07).sin() * 2.0).collect(),
    ];
    let mut knots = coarse.knots.clone();
    knots.extend([0.5, 0.75, 0.875]);
    knots.sort_by(f64::total_cmp);
    let fine = SectorSpace::new(5, 7, knots);
    let moved = super::fit::transfer(&coarse, &coords, &fine).expect("transfer");
    for s in 0..5 {
        let a = coarse.sector_surface(&coords, s);
        let b = fine.sector_surface(&moved, s);
        for &(u, v) in &[(0.1, 0.9), (0.6, 0.3), (0.99, 0.8), (1.0, 1.0)] {
            assert!(a.point_at(u, v).distance_to(&b.point_at(u, v)) < 1e-12);
        }
    }
}

#[test]
fn a_thin_sector_serves_only_inside_the_curvature_limit() {
    // A spherical pentagon with one short side (0.012 rad of azimuth
    // against ~1.6): the short side's sectors are thin (regularity
    // 1.1e-2, under the 2e-2 served floor), and the fit bends beside it
    // to curvature 18 on the unit sphere (measured).
    let sides = spherical_polygon(&[0.0, 0.012, 1.6, 3.2, 4.7], 0.6);
    let targets = VariationalCapTargets::g1(1e-6, 1e-4);
    let with_limit = |limit| VariationalCapOptions {
        thin_sector_curvature: limit,
        ..Default::default()
    };
    // Without a limit the served floor alone decides.
    assert!(matches!(
        fit_variational_cap(&sides, &targets, &with_limit(None), 1e-7),
        Err(VariationalCapError::Folded { .. })
    ));
    // Under a limit the cap stays inside, the thin cap serves.
    let cap = fit_variational_cap(&sides, &targets, &with_limit(Some(100.0)), 1e-7)
        .expect("a thin cap inside the limit serves");
    let regularity = cap.report.min_regularity;
    assert!((5e-3..2e-2).contains(&regularity), "{regularity}");
    // Past a limit under its curvature it declines typed.
    match fit_variational_cap(&sides, &targets, &with_limit(Some(4.0)), 1e-7) {
        Err(VariationalCapError::ThinSectorBends {
            curvature, limit, ..
        }) => assert!(curvature > limit && limit == 4.0, "{curvature}"),
        other => panic!("expected ThinSectorBends: {other:?}"),
    }
    // Under the hard floor no limit admits the cap.
    let sides = spherical_polygon(&[0.0, 0.003, 1.6, 3.2, 4.7], 0.6);
    assert!(matches!(
        fit_variational_cap(&sides, &targets, &with_limit(Some(1e9)), 1e-7),
        Err(VariationalCapError::Folded { .. })
    ));
}

#[test]
fn the_fit_takes_up_to_twenty_four_sides() {
    // The bridged setback border of an eight-strip apex has 16 sides.
    let targets = VariationalCapTargets::g1(1e-6, 1e-4);
    let options = VariationalCapOptions {
        max_refinement_levels: 0,
        ..Default::default()
    };
    let cap =
        fit_variational_cap(&spherical_ngon(16, 0.6), &targets, &options, 1e-7).expect("16 sides");
    assert_eq!(cap.sectors.len(), 16);
    assert!(cap.report.min_regularity >= 2e-2);
    assert!(matches!(
        fit_variational_cap(
            &spherical_ngon(MAX_CAP_SIDES + 1, 0.6),
            &targets,
            &options,
            1e-7
        ),
        Err(VariationalCapError::UnsupportedSideCount { got: 25 })
    ));
}

/// A spherical polygon at polar angle `theta` about +z with corners at
/// the azimuths `angles` (counter-clockwise).
pub(super) fn spherical_polygon(angles: &[f64], theta: f64) -> Vec<VariationalCapSide> {
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, 0.0),
        radius: 1.0,
    };
    let n = angles.len();
    let corner = |k: usize| {
        let a = angles[k % n];
        Vec3::new(theta.sin() * a.cos(), theta.sin() * a.sin(), theta.cos())
    };
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
