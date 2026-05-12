use super::*;

const EPS: f64 = 1e-9;

#[test]
fn tri6_basis_partitions_unity() {
    for (u, v) in [(0.0, 0.0), (0.2, 0.3), (0.7, 0.1), (0.0, 1.0)] {
        let sum: f64 = bernstein_tri6_all(u, v).iter().sum();
        assert!((sum - 1.0).abs() < EPS, "sum={sum} at ({u}, {v})");
    }
}

#[test]
fn tri6_first_partials_sum_to_zero() {
    let u = 0.23;
    let v = 0.31;
    let du: f64 = bernstein_tri6_partial_u(u, v).iter().sum();
    let dv: f64 = bernstein_tri6_partial_v(u, v).iter().sum();
    assert!(du.abs() < EPS, "du sum {du}");
    assert!(dv.abs() < EPS, "dv sum {dv}");
}

#[test]
fn tri6_partials_match_finite_differences() {
    let u = 0.21;
    let v = 0.27;
    let h = 1e-6;
    let b = bernstein_tri6_all(u, v);
    let bu = bernstein_tri6_all(u + h, v);
    let bv = bernstein_tri6_all(u, v + h);
    let du = bernstein_tri6_partial_u(u, v);
    let dv = bernstein_tri6_partial_v(u, v);
    for i in 0..TRI6_CP_COUNT {
        assert!(((bu[i] - b[i]) / h - du[i]).abs() < 2e-4);
        assert!(((bv[i] - b[i]) / h - dv[i]).abs() < 2e-4);
    }
}

#[test]
fn tri6_third_partials_sum_to_zero() {
    let u = 0.18;
    let v = 0.22;
    for values in [
        bernstein_tri6_third_partial_uuu(u, v),
        bernstein_tri6_third_partial_uuv(u, v),
        bernstein_tri6_third_partial_uvv(u, v),
        bernstein_tri6_third_partial_vvv(u, v),
    ] {
        let sum: f64 = values.iter().sum();
        assert!(sum.abs() < 1e-8, "third partial sum {sum}");
    }
}
