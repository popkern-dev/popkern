use super::tests_support::*;
use super::*;

// ---- Corner CP identities ----

#[test]
fn bi5_star_corner_equals_tri5_star() {
    // At bi-5 CP (0, 0) the value is the patch value at (0, 0),
    // which in simplex coords is the star vertex (u_s = 0, v_s = 0,
    // w_s = 1). So bi5[0][0] = cps[IDX_500].
    let mut cps = [Point3::new(0.0, 0.0, 0.0); TRI5_CP_COUNT];
    for i in 0..TRI5_CP_COUNT {
        cps[i] = Point3::new(i as f64 + 1.0, i as f64 * 2.0, -(i as f64));
    }
    let bi5 = tri5_to_bi5_degree_raise(&cps);
    assert!(approx_p(bi5[0][0], cps[IDX_500], 1e-12));
}

#[test]
fn bi5_mk_corner_equals_tri5_mk() {
    // At bi-5 CP (5, 0): u = 1, v = 0 → simplex coords (1, 0, 0) = m_k.
    let mut cps = [Point3::new(0.0, 0.0, 0.0); TRI5_CP_COUNT];
    for i in 0..TRI5_CP_COUNT {
        cps[i] = Point3::new(i as f64, -(i as f64), i as f64 / 2.0);
    }
    let bi5 = tri5_to_bi5_degree_raise(&cps);
    // IDX_050 = (i = 5, j = 0) is the m_k corner in this codebase's
    // convention (paper label "050" = barycentric (0, 5, 0)).
    assert!(approx_p(bi5[5][0], cps[IDX_050], 1e-12));
}

#[test]
fn bi5_mk1_corner_equals_tri5_mk1() {
    // At bi-5 CP (0, 5): u = 0, v = 1 → simplex (0, 1, 0) = m_{k+1}.
    let mut cps = [Point3::new(0.0, 0.0, 0.0); TRI5_CP_COUNT];
    for i in 0..TRI5_CP_COUNT {
        cps[i] = Point3::new(-(i as f64), i as f64, i as f64 + 1.0);
    }
    let bi5 = tri5_to_bi5_degree_raise(&cps);
    // IDX_005 = (i = 0, j = 5) is the m_{k+1} corner.
    assert!(approx_p(bi5[0][5], cps[IDX_005], 1e-12));
}
