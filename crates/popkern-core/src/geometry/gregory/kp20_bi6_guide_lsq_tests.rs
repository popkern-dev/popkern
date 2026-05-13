use super::*;

#[test]
fn solves_planar_valence_five_tensor_border() {
    let tensor = planar_tensor_border(5);
    let sectors = solve_kp20_guide_lsq(&tensor, Point3::new(0.0, 0.0, 0.0)).unwrap();
    assert_eq!(sectors.len(), 5);
    for sector in sectors {
        for cp in sector.cps {
            assert!(cp.x.is_finite());
            assert!(cp.y.is_finite());
            assert!(cp.z.is_finite());
        }
    }
}

#[test]
fn rejects_kp18_excluded_regular_valence_four() {
    let tensor = planar_tensor_border(4);
    let err = solve_kp20_guide_lsq(&tensor, Point3::new(0.0, 0.0, 0.0)).unwrap_err();
    assert_eq!(err, Kp20GuideLsqError::UnsupportedValence(4));
}

#[test]
fn reduced_constraints_leave_kp18_lsq_free_count() {
    for n in [5, 6, 7, 8, 9, 10] {
        let vars = n * TRI6_CP_COUNT;
        let rows = constraint_rows(n).unwrap();
        let independent = independent_rows(vars, &rows);
        assert_eq!(vars - independent.len() - 1, 12 * n + 9);
    }
}

fn planar_tensor_border(n: usize) -> TensorBorderBi3 {
    let mut sectors = Vec::new();
    for s in 0..n {
        let a = core::f64::consts::TAU * s as f64 / n as f64;
        let b = core::f64::consts::TAU * (s + 1) as f64 / n as f64;
        let p0 = Point3::new(a.cos(), a.sin(), 0.0);
        let p3 = Point3::new(b.cos(), b.sin(), 0.0);
        let mut rows = [[Point3::new(0.0, 0.0, 0.0); 4]; 3];
        for i in 0..4 {
            let t = i as f64 / 3.0;
            let x = (1.0 - t) * p0.x + t * p3.x;
            let y = (1.0 - t) * p0.y + t * p3.y;
            rows[0][i] = Point3::new(x, y, 0.0);
            rows[1][i] = Point3::new(0.8 * x, 0.8 * y, 0.0);
            rows[2][i] = Point3::new(0.6 * x, 0.6 * y, 0.0);
        }
        sectors.push(rows);
    }
    TensorBorderBi3 { sectors }
}
