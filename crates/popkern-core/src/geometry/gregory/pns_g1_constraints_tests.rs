//! Tests for [`pns_g1_constraints`]: eq (2) bi-4 and eq (2³) bi-3.

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::Point3;

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn approx_pt(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

#[test]
fn bi4_zero_inputs_give_zero_outputs() {
    let zero = p(0.0, 0.0, 0.0);
    let row1 = [zero; 5];
    for n in [5_usize, 6, 7, 8] {
        let r = bi4_g1_row_from_row1(zero, row1, row1, n).expect("determined valence must solve");
        assert!(approx_pt(r.cur10, zero, 1e-12), "n={n}");
        assert!(approx_pt(r.cur20, zero, 1e-12), "n={n}");
        assert!(approx_pt(r.cur30, zero, 1e-12), "n={n}");
        assert!(approx_pt(r.cur40, zero, 1e-12), "n={n}");
        assert!(approx_pt(r.nxt21, zero, 1e-12), "n={n}");
    }
}

#[test]
fn bi4_n6_matches_hand_calculation() {
    // For n=6: c = cos(π/3) = 0.5.
    // Simple inputs:
    //   cur00 = (1, 2, 3)
    //   cur_row1 = [(1,0,0), (0,1,0), (0,0,1), (2,0,0), (0,2,0)]
    //   nxt_row1 = [(0,0,1), (1,0,0), (0,1,0), (0,2,0), (2,0,0)]
    let cur00 = p(1.0, 2.0, 3.0);
    let cur_row1 = [
        p(1.0, 0.0, 0.0),
        p(0.0, 1.0, 0.0),
        p(0.0, 0.0, 1.0),
        p(2.0, 0.0, 0.0),
        p(0.0, 2.0, 0.0),
    ];
    let nxt_row1 = [
        p(0.0, 0.0, 1.0),
        p(1.0, 0.0, 0.0),
        p(0.0, 1.0, 0.0),
        p(0.0, 2.0, 0.0),
        p(2.0, 0.0, 0.0),
    ];
    let r =
        bi4_g1_row_from_row1(cur00, cur_row1, nxt_row1, 6).expect("determined valence must solve");

    // cur10 = (2(c-1)·cur00 + cur01 + nxt01) / (2c)
    //       = (-1·(1,2,3) + (1,0,0) + (0,0,1)) / 1
    //       = (-1+1+0, -2+0+0, -3+0+1) = (0, -2, -2)
    assert!(approx_pt(r.cur10, p(0.0, -2.0, -2.0), 1e-12));

    // cur40 = (cur41 + nxt41)/2 = ((0,2,0)+(2,0,0))/2 = (1, 1, 0)
    assert!(approx_pt(r.cur40, p(1.0, 1.0, 0.0), 1e-12));

    // cur30 = [4(cur31+nxt31) - c(cur41+nxt41)] / (2(4-c))
    //       = [4((2,0,0)+(0,2,0)) - 0.5((0,2,0)+(2,0,0))] / 7
    //       = [(8,8,0) - (1,1,0)] / 7
    //       = (7, 7, 0) / 7 = (1, 1, 0)
    assert!(approx_pt(r.cur30, p(1.0, 1.0, 0.0), 1e-12));

    // cur20 = [4c(cur11+nxt11) + (3c-4)(cur01+nxt01+2(c-1)cur00)] / (6c²)
    // 3c-4 = -2.5, 2(c-1) = -1, 6c² = 1.5
    // inner = (1,0,0)+(0,0,1) + (-1)(1,2,3) = (0,-2,-2)
    // = [4·0.5·((0,1,0)+(1,0,0)) + (-2.5)(0,-2,-2)] / 1.5
    // = [2·(1,1,0) + (0, 5, 5)] / 1.5
    // = [(2,1,0) + (0,5,5)] / 1.5   (wait, 2·(1,1,0) = (2,2,0))
    // = [(2,2,0) + (0,5,5)] / 1.5
    // = (2, 7, 5) / 1.5 = (4/3, 14/3, 10/3)
    assert!(approx_pt(
        r.cur20,
        p(4.0 / 3.0, 14.0 / 3.0, 10.0 / 3.0),
        1e-12
    ));

    // nxt21 = -cur21 + (2-c)·cur20 + c·cur30
    //       = -(0,0,1) + 1.5·(4/3, 14/3, 10/3) + 0.5·(1, 1, 0)
    //       = (0,0,-1) + (2, 7, 5) + (0.5, 0.5, 0)
    //       = (2.5, 7.5, 4.0)
    assert!(approx_pt(r.nxt21, p(2.5, 7.5, 4.0), 1e-12));
}

/// Construct a bi-4 row-0 `[cur_{00,10,20,30,40}]` via eq (2) from
/// a given row-1 + central, then evaluate the G¹ identity
/// `nxt_v(u,0) + cur_v(u,0) = 2c(1-u) cur_u(u,0)` at several u
/// values. This is the strong test: it verifies the transcription
/// matches the paper's underlying polynomial identity.
#[test]
fn bi4_g1_polynomial_identity_holds() {
    // Set up two adjacent bi-4 patches with their j=1 row + shared corners.
    // `cur` has full 5×5 grid; we initialise j=1 arbitrarily, and
    // j=0 via eq (2).
    let cur00 = p(0.0, 0.0, 0.0);
    let cur_row1 = [
        p(0.3, 0.1, 0.0),
        p(0.5, 0.3, 0.1),
        p(0.7, 0.4, 0.2),
        p(0.9, 0.2, 0.0),
        p(1.0, 0.0, 0.2),
    ];
    let nxt_row1 = [
        p(0.3, 0.1, 0.0), // = cur_row1[0] — they share the "nxt_{01}" = sector k+1's cur01
        p(0.4, 0.5, 0.1),
        p(0.5, 0.6, 0.3),
        p(0.7, 0.4, 0.1),
        p(1.0, 0.0, 0.2), // = cur_row1[4] — shared corner
    ];
    // Note: for these test inputs we don't enforce the shared-boundary
    // constraint cur_row1[0] = nxt_row1[0] and cur_row1[4] = nxt_row1[4]
    // (which would be the case in a real cap — the corners of adjacent
    // sectors coincide). The test inputs satisfy this.

    for n in [5_usize, 6, 7, 8] {
        let r = bi4_g1_row_from_row1(cur00, cur_row1, nxt_row1, n)
            .expect("determined valence must solve");
        let c = (2.0 * std::f64::consts::PI / (n as f64)).cos();

        // Assemble cur j=0 row:
        //   [cur00, cur10, cur20, cur30, cur40]
        let cur_row0 = [cur00, r.cur10, r.cur20, r.cur30, r.cur40];

        // Assemble nxt j=1 row (with nxt21 replaced):
        let mut nxt_row1_final = nxt_row1;
        nxt_row1_final[2] = r.nxt21;

        // G¹ eq (1): nxt_v(u, 0) + cur_v(u, 0) = 2c(1 − u) · cur_u(u, 0)
        //
        // For a bi-4 patch p(u, v) = Σ_{i,j} p_{ij} B_i^4(u) B_j^4(v),
        //   p_u(u, 0) = 4 · Σ_{i=0}^3 (p_{i+1,0} − p_{i,0}) B_i^3(u)
        //   p_v(u, 0) = 4 · Σ_{i=0}^4 (p_{i,1} − p_{i,0}) B_i^4(u)
        //
        // Note: `nxt_v(u, 0)` here means the v-derivative of the
        // next patch at its shared-boundary (which is also v=0 in
        // its own local frame). Because adjacent sectors share the
        // boundary row `cur_row0 == nxt_row0` (same curve), the
        // relevant first differences are `nxt_{i,1} − nxt_{i,0}`
        // with `nxt_{i,0} = cur_{i,0}` (shared).
        let nxt_row0 = cur_row0; // shared boundary

        for &u in &[0.0, 0.2, 0.5, 0.7, 1.0] {
            let cur_u_at_u0 = deriv_u_bi4(&cur_row0, u);
            let cur_v_at_u0 = deriv_v_bi4(&cur_row0, &cur_row1, u);
            let nxt_v_at_u0 = deriv_v_bi4(&nxt_row0, &nxt_row1_final, u);

            let lhs_x = nxt_v_at_u0.x + cur_v_at_u0.x;
            let lhs_y = nxt_v_at_u0.y + cur_v_at_u0.y;
            let lhs_z = nxt_v_at_u0.z + cur_v_at_u0.z;
            let fac = 2.0 * c * (1.0 - u);
            let rhs_x = fac * cur_u_at_u0.x;
            let rhs_y = fac * cur_u_at_u0.y;
            let rhs_z = fac * cur_u_at_u0.z;

            let eps = 1e-10;
            assert!(
                (lhs_x - rhs_x).abs() < eps,
                "n={n}, u={u}: lhs.x={lhs_x} rhs.x={rhs_x}"
            );
            assert!(
                (lhs_y - rhs_y).abs() < eps,
                "n={n}, u={u}: lhs.y={lhs_y} rhs.y={rhs_y}"
            );
            assert!(
                (lhs_z - rhs_z).abs() < eps,
                "n={n}, u={u}: lhs.z={lhs_z} rhs.z={rhs_z}"
            );
        }
    }
}

/// Degree-4 boundary u-derivative at v=0:
/// `p_u(u, 0) = 4 · Σ_{i=0}^3 (p_{i+1,0} − p_{i,0}) B_i^3(u)`
fn deriv_u_bi4(row0: &[Point3; 5], u: f64) -> Point3 {
    let b3 = [
        (1.0 - u).powi(3),
        3.0 * u * (1.0 - u).powi(2),
        3.0 * u * u * (1.0 - u),
        u.powi(3),
    ];
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..4 {
        let dx = row0[i + 1].x - row0[i].x;
        let dy = row0[i + 1].y - row0[i].y;
        let dz = row0[i + 1].z - row0[i].z;
        out.x += 4.0 * b3[i] * dx;
        out.y += 4.0 * b3[i] * dy;
        out.z += 4.0 * b3[i] * dz;
    }
    out
}

/// Degree-4 boundary v-derivative at v=0:
/// `p_v(u, 0) = 4 · Σ_{i=0}^4 (p_{i,1} − p_{i,0}) B_i^4(u)`
fn deriv_v_bi4(row0: &[Point3; 5], row1: &[Point3; 5], u: f64) -> Point3 {
    let b4 = {
        let omu = 1.0 - u;
        [
            omu.powi(4),
            4.0 * u * omu.powi(3),
            6.0 * u * u * omu * omu,
            4.0 * u.powi(3) * omu,
            u.powi(4),
        ]
    };
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..5 {
        let dx = row1[i].x - row0[i].x;
        let dy = row1[i].y - row0[i].y;
        let dz = row1[i].z - row0[i].z;
        out.x += 4.0 * b4[i] * dx;
        out.y += 4.0 * b4[i] * dy;
        out.z += 4.0 * b4[i] * dz;
    }
    out
}

#[test]
fn bi3_row20_n3_matches_hand_calc() {
    // For n=3: c = -1/2.
    // cur_{20} = [3(cur21+nxt21) - c(cur31+nxt31)] / (2(3-c))
    //         = [3((1,0,0)+(0,1,0)) - (-0.5)((0,0,1)+(1,1,0))] / 7
    //         = [(3,3,0) + 0.5·(1,1,1)] / 7
    //         = (3.5, 3.5, 0.5) / 7
    //         = (0.5, 0.5, 1.0/14)
    let cur21 = p(1.0, 0.0, 0.0);
    let cur31 = p(0.0, 0.0, 1.0);
    let nxt21 = p(0.0, 1.0, 0.0);
    let nxt31 = p(1.0, 1.0, 0.0);
    let got = bi3_g1_row20(cur21, cur31, nxt21, nxt31, 3);
    assert!(approx_pt(got, p(0.5, 0.5, 1.0 / 14.0), 1e-12));
}

#[test]
fn bi3_row20_n5_zero_inputs() {
    let zero = p(0.0, 0.0, 0.0);
    let got = bi3_g1_row20(zero, zero, zero, zero, 5);
    assert!(approx_pt(got, zero, 1e-12));
}
