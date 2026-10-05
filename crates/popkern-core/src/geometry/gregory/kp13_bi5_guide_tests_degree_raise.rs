use super::tests_support::*;
use super::*;
use popkern_math::PortableMath;

// ---- tri5_to_bi5_degree_raise: polynomial equivalence ----

#[test]
fn degree_raise_preserves_polynomial_inside_simplex() {
    // Build a deterministic tri-5 CP vector and verify that the
    // degree-raised bi-5 patch evaluates to the same value as the
    // tri-5 patch, at interior points of the simplex.
    let mut cps = [Point3::new(0.0, 0.0, 0.0); TRI5_CP_COUNT];
    for i in 0..TRI5_CP_COUNT {
        let t = i as f64 * 0.37;
        cps[i] = Point3::new(t.psin(), t.pcos(), 0.5 * t.psin() * t.pcos());
    }
    let bi5 = tri5_to_bi5_degree_raise(&cps);
    let sector = Bi5GuideSector { cps: bi5 };

    for &(u, v) in &[(0.1, 0.1), (0.25, 0.5), (0.3, 0.3), (0.5, 0.25), (0.7, 0.2)] {
        let tri_val = {
            let b = bernstein_tri5_all(u, v);
            let mut out = Point3::new(0.0, 0.0, 0.0);
            for i in 0..TRI5_CP_COUNT {
                out.x += b[i] * cps[i].x;
                out.y += b[i] * cps[i].y;
                out.z += b[i] * cps[i].z;
            }
            out
        };
        let bi_val = sector.eval(u, v);
        assert!(
            approx_p(tri_val, bi_val, 1e-10),
            "(u, v) = ({u}, {v}): tri={tri_val:?}, bi={bi_val:?}"
        );
    }
}

#[test]
fn degree_raise_preserves_polynomial_outside_simplex() {
    // Key observation: g△ is polynomial, so the identity
    //   tri5_polynomial(u, v) = bi5_polynomial(u, v)
    // holds for ALL (u, v) ∈ R², including (u, v) ∈ [0, 1]² where
    // u + v > 1 (the v_k region).
    let mut cps = [Point3::new(0.0, 0.0, 0.0); TRI5_CP_COUNT];
    for i in 0..TRI5_CP_COUNT {
        let t = i as f64 * 0.123;
        cps[i] = Point3::new(t.pcos(), -t.psin(), t.psin() + 0.2);
    }
    let bi5 = tri5_to_bi5_degree_raise(&cps);
    let sector = Bi5GuideSector { cps: bi5 };

    for &(u, v) in &[(0.8, 0.7), (1.0, 1.0), (0.5, 0.9), (0.95, 0.6)] {
        let tri_val = eval_tri5_extended(&cps, u, v);
        let bi_val = sector.eval(u, v);
        assert!(
            approx_p(tri_val, bi_val, 1e-9),
            "(u, v) = ({u}, {v}) in v_k region: tri={tri_val:?}, bi={bi_val:?}"
        );
    }
}

#[test]
fn degree_raise_linearity() {
    // The degree-raise is linear in the tri-5 CPs.
    let mut cps_a = [Point3::new(0.0, 0.0, 0.0); TRI5_CP_COUNT];
    let mut cps_b = [Point3::new(0.0, 0.0, 0.0); TRI5_CP_COUNT];
    for i in 0..TRI5_CP_COUNT {
        cps_a[i] = Point3::new(i as f64, 0.0, 0.0);
        cps_b[i] = Point3::new(0.0, i as f64, (-1.0_f64).powi(i as i32) * 0.5);
    }
    let bi_a = tri5_to_bi5_degree_raise(&cps_a);
    let bi_b = tri5_to_bi5_degree_raise(&cps_b);
    // Sum:
    let mut cps_sum = cps_a;
    for i in 0..TRI5_CP_COUNT {
        cps_sum[i].x += cps_b[i].x;
        cps_sum[i].y += cps_b[i].y;
        cps_sum[i].z += cps_b[i].z;
    }
    let bi_sum = tri5_to_bi5_degree_raise(&cps_sum);
    for p in 0..6 {
        for q in 0..6 {
            let expected_x = bi_a[p][q].x + bi_b[p][q].x;
            let expected_y = bi_a[p][q].y + bi_b[p][q].y;
            let expected_z = bi_a[p][q].z + bi_b[p][q].z;
            assert!(
                (bi_sum[p][q].x - expected_x).abs() < 1e-9
                    && (bi_sum[p][q].y - expected_y).abs() < 1e-9
                    && (bi_sum[p][q].z - expected_z).abs() < 1e-9,
                "linearity fails at ({p}, {q})"
            );
        }
    }
}
