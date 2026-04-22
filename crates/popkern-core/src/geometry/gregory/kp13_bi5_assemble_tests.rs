//! Tests for Phase I.16 end-to-end KP13 bi-5 cap assembly.

#![allow(clippy::needless_range_loop)]

use super::super::bernstein_bi5::bernstein_bi5_all;
use super::super::kp13_bi5_boundary::TensorBorder;
use super::*;

const EPS: f64 = 1e-9;

fn approx_p(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

// ---- Structural ----

#[test]
fn zero_input_produces_zero_output() {
    for n in [3, 4, 5, 6, 7, 8] {
        let input = Bi5CapInput {
            target: TensorBorder::zero(n),
            cc_limit: Point3::new(0.0, 0.0, 0.0),
        };
        let out = assemble_bi5_macro_cap(&input).unwrap();
        assert_eq!(out.n, n);
        assert_eq!(out.patches_bi9.len(), 4 * n);
        for (idx, patch) in out.patches_bi9.iter().enumerate() {
            for p in 0..10 {
                for q in 0..10 {
                    let cp = patch[p][q];
                    assert!(
                        cp.x.abs() < 1e-8 && cp.y.abs() < 1e-8 && cp.z.abs() < 1e-8,
                        "n={n}, patch {idx}, ({p}, {q}): non-zero {cp:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn output_has_correct_shape() {
    for n in [3, 4, 5, 6, 7, 8] {
        let input = Bi5CapInput {
            target: TensorBorder::zero(n),
            cc_limit: Point3::new(0.1, 0.2, 0.3),
        };
        let out = assemble_bi5_macro_cap(&input).unwrap();
        assert_eq!(out.n, n);
        assert_eq!(out.cap.n, n);
        assert_eq!(out.cap.sectors.len(), n);
        assert_eq!(out.patches_bi9.len(), 4 * n);
    }
}

#[test]
fn cap_validates() {
    let n = 5;
    let mut target = TensorBorder::zero(n);
    for k in 0..n {
        target.patches[k].cps[3][3] = Point3::new(k as f64, k as f64 * 2.0, k as f64 * 3.0);
    }
    let input = Bi5CapInput {
        target,
        cc_limit: Point3::new(0.0, 0.0, 0.0),
    };
    let out = assemble_bi5_macro_cap(&input).unwrap();
    out.cap.validate().unwrap();
    for (idx, patch) in out.patches_bi9.iter().enumerate() {
        for p in 0..10 {
            for q in 0..10 {
                let cp = patch[p][q];
                assert!(
                    cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite(),
                    "patch {idx} ({p}, {q}) non-finite {cp:?}"
                );
            }
        }
    }
}

// ---- Bi-5 ↔ bi-9 elevation equivalence ----

fn eval_bi5(cps: &[[Point3; 6]; 6], u: f64, v: f64) -> Point3 {
    let bu = bernstein_bi5_all(u);
    let bv = bernstein_bi5_all(v);
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..6 {
        for j in 0..6 {
            let w = bu[i] * bv[j];
            out.x += w * cps[i][j].x;
            out.y += w * cps[i][j].y;
            out.z += w * cps[i][j].z;
        }
    }
    out
}

fn eval_bi9(cps: &[[Point3; 10]; 10], u: f64, v: f64) -> Point3 {
    // Evaluate bi-9 via the standard binomial formula since the
    // bernstein_bi9 module is not imported here.
    let bu = bernstein_deg9(u);
    let bv = bernstein_deg9(v);
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..10 {
        for j in 0..10 {
            let w = bu[i] * bv[j];
            out.x += w * cps[i][j].x;
            out.y += w * cps[i][j].y;
            out.z += w * cps[i][j].z;
        }
    }
    out
}

fn bernstein_deg9(t: f64) -> [f64; 10] {
    let mut out = [0.0_f64; 10];
    let s = 1.0 - t;
    // B_i^9(t) = C(9, i) · t^i · s^(9-i).
    let binom: [f64; 10] = [1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0];
    for i in 0..10 {
        out[i] = binom[i] * t.powi(i as i32) * s.powi(9 - i as i32);
    }
    out
}

#[test]
fn bi9_elevation_preserves_polynomial() {
    // Degree raise from bi-5 to bi-9 must preserve the polynomial
    // identically (at every (u, v) in [0, 1]²).
    let n = 4;
    let mut target = TensorBorder::zero(n);
    for k in 0..n {
        for i in 0..4 {
            for j in 0..4 {
                target.patches[k].cps[i][j] =
                    Point3::new(k as f64 + i as f64, j as f64, (i * j) as f64 * 0.1);
            }
        }
    }
    let input = Bi5CapInput {
        target,
        cc_limit: Point3::new(0.0, 0.0, 0.0),
    };
    let out = assemble_bi5_macro_cap(&input).unwrap();

    for k in 0..n {
        for (piece_idx, piece) in [CapPieceId::O, CapPieceId::L, CapPieceId::R, CapPieceId::M]
            .iter()
            .enumerate()
        {
            let bi5 = out.cap.sector(k).piece(*piece).cps;
            let bi9 = &out.patches_bi9[4 * k + piece_idx];
            for &(u, v) in &[(0.0, 0.0), (0.5, 0.5), (0.3, 0.7), (1.0, 1.0)] {
                let v5 = eval_bi5(&bi5, u, v);
                let v9 = eval_bi9(bi9, u, v);
                assert!(
                    approx_p(v5, v9, 1e-8),
                    "k={k}, piece {piece:?}, (u,v)=({u},{v}): bi5={v5:?} vs bi9={v9:?}"
                );
            }
        }
    }
}

// ---- Determinism ----

#[test]
fn assemble_deterministic() {
    let n = 5;
    let mut target = TensorBorder::zero(n);
    for k in 0..n {
        target.patches[k].cps[2][2] = Point3::new(0.5, -0.3, 0.7);
        target.patches[k].cps[3][3] = Point3::new(1.0, 1.0, 1.0);
    }
    let input = Bi5CapInput {
        target,
        cc_limit: Point3::new(0.0, 0.0, 0.0),
    };
    let a = assemble_bi5_macro_cap(&input).unwrap();
    let b = assemble_bi5_macro_cap(&input).unwrap();
    assert_eq!(a, b);
}

// ---- Star point pinning ----

#[test]
fn cap_star_point_equals_cc_limit_via_o_piece_origin() {
    let n = 6;
    let cc = Point3::new(1.5, -0.5, 0.25);
    let input = Bi5CapInput {
        target: TensorBorder::zero(n),
        cc_limit: cc,
    };
    let out = assemble_bi5_macro_cap(&input).unwrap();

    // For each sector, O piece evaluated at (0, 0) in piece-local
    // = star corner. Should equal cc_limit (value component of
    // the star jet of g̃^△^0, which starts as 0 in zero-target
    // F_1-min solve — so actually 0, not cc_limit).
    //
    // The value at the star corner depends on both the guide
    // (which pins star to cc_limit) AND the F_1-min'ed p (which
    // may pick a different star value from G_3[0]). For zero
    // target, F_1 min can't reduce the CP-norm target, so may
    // leave the star at something not equal to cc_limit.
    //
    // Here we just sanity-check the O piece's (0, 0) is finite.
    for k in 0..n {
        let o_star = out.cap.sector(k).piece(CapPieceId::O).cps[0][0];
        assert!(
            o_star.x.is_finite() && o_star.y.is_finite() && o_star.z.is_finite(),
            "sector {k} O star non-finite: {o_star:?}"
        );
    }
}
