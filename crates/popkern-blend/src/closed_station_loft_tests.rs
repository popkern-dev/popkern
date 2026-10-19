// This file is only compiled under #[cfg(test)] (see closed_station_loft.rs).
use super::*;
use popkern_core::geometry::SurfaceEval;

/// Exact rational quarter-arc rows around a tilted, non-circular
/// closed spine: `n` stations of a tube of radius 0.5 about an ellipse.
fn station_rows(n: usize) -> (Vec<Vec<Point3>>, Vec<Vec<f64>>) {
    use popkern_core::geometry::PortableMath;
    let w_mid = std::f64::consts::FRAC_1_SQRT_2;
    let mut rows = Vec::new();
    let mut weights = Vec::new();
    for i in 0..n {
        let a = std::f64::consts::TAU * i as f64 / n as f64;
        let (s, c) = (a.psin(), a.pcos());
        let centre = Point3::new(5.0 * c, 3.0 * s, 0.4 * c);
        let out = popkern_core::geometry::Vec3::new(c, s, 0.0);
        let up = popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0);
        rows.push(vec![
            centre + out * 0.5,
            centre + (out + up) * 0.5,
            centre + up * 0.5,
        ]);
        weights.push(vec![1.0, w_mid, 1.0]);
    }
    (rows, weights)
}

#[test]
fn the_loft_is_periodic_and_interpolates_every_station() {
    let (rows, weights) = station_rows(24);
    let knots_v = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let loft = closed_station_loft(&rows, &weights, 2, &knots_v).expect("loft");
    assert!((loft.u_period().expect("periodic") - 1.0).abs() < 1e-15);
    assert_eq!(loft.parameter_range().0, (0.0, 1.0));
    // Station i sits on knot i + 3 (chord-length parameters).
    for (i, row) in rows.iter().enumerate() {
        let u = loft.knots_u[i + 3];
        for (v, p) in [(0.0, row[0]), (1.0, row[2])] {
            let gap = loft.point_at(u, v).distance_to(&p);
            assert!(gap < 1e-12, "station {i} v {v}: {gap:e}");
        }
    }
    // C2 across the seam, and the rails close.
    let (p0, du0, _) = loft.point_and_first_partials(0.0, 0.5);
    let (p1, du1, _) = loft.point_and_first_partials(1.0, 0.5);
    assert!(p0.distance_to(&p1) < 1e-12);
    assert!((du0 - du1).length() < 1e-9 * du0.length());
}

#[test]
fn too_few_or_ragged_stations_decline() {
    let (rows, weights) = station_rows(2);
    let knots_v = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    assert!(closed_station_loft(&rows, &weights, 2, &knots_v).is_none());
    let (mut rows, weights) = station_rows(8);
    rows[3].pop();
    assert!(closed_station_loft(&rows, &weights, 2, &knots_v).is_none());
}
