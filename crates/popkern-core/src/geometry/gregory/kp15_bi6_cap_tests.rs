//! Tests for the bi-6 polar cap scaffolding.
//!
//! Strategy: verify the orchestrator wires together I.1–I.9
//! pieces correctly by running on synthetic symmetric inputs and
//! checking the outputs satisfy basic invariants:
//! - Star point matches independent Catmull-Clark.
//! - Row-0 CPs (shared-edge) match the sector inputs.
//! - Circulant-solver outputs match the independently-solved
//!   ring values.
//! - bi-9 grids are valid-shape with finite CPs.

// Paper-notation test module: loop indices, helper arities, and tuple shapes
// mirror KP15/KP20/PR08 subscripts and equation forms; keep the index form.
#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::Point3;
use crate::geometry::gregory::kp15_bi6_rules::rule_5_reduce_edge;

fn pt(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

fn approx_eq(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

/// Build a deterministic pentagon test input. Geometry is rough
/// — we're testing the wiring, not commercial-grade output.
fn pentagon_input() -> Bi6CapInput {
    let n = 5;
    let params = Bi6Reparams::for_valence(n);
    let c7 = pt(0.0, 0.0, 1.0);
    let mut c5s = Vec::with_capacity(n);
    let mut c6s = Vec::with_capacity(n);
    let mut sectors = Vec::with_capacity(n);
    for k in 0..n {
        let theta = (k as f64) * 2.0 * std::f64::consts::PI / (n as f64);
        let dx = theta.cos();
        let dy = theta.sin();
        c5s.push(pt(dx * 0.3, dy * 0.3, 0.9));
        c6s.push(pt(dx * 0.7, dy * 0.7, 0.6));
        sectors.push(CirculantSectorRhs {
            p_00: c7,
            p_10: pt(dx * 0.1, dy * 0.1, 0.95),
            p_20: pt(dx * 0.25, dy * 0.25, 0.85),
            p_40: pt(dx * 0.6, dy * 0.6, 0.4),
            p_50: pt(dx * 0.8, dy * 0.8, 0.25),
            p_60: pt(dx, dy, 0.1),
            p_01: pt(-dy * 0.1, dx * 0.1, 0.95),
            p_11: pt((dx - dy) * 0.1, (dy + dx) * 0.1, 0.85),
            p_02: pt(-dy * 0.2, dx * 0.2, 0.88),
            p_42: pt(-dy * 0.5, dx * 0.5, 0.35),
            p_52: pt(-dy * 0.6, dx * 0.6, 0.2),
            p_62: pt(-dy * 0.7, dx * 0.7, 0.05),
        });
    }
    Bi6CapInput {
        n,
        params,
        c7,
        c5s,
        c6s,
        sectors,
        p0_22: pt(0.0, 0.0, 0.95),
        p0_21_pinned: None,
        interior_filler: None,
    }
}

#[test]
fn pentagon_assembles_without_error() {
    let input = pentagon_input();
    let out = assemble_bi6_polar_cap(&input).expect("pentagon should assemble");
    assert_eq!(out.patches_bi6.len(), 5);
    assert_eq!(out.patches_bi9.len(), 5);
}

#[test]
fn star_point_matches_catmull_clark() {
    let input = pentagon_input();
    let out = assemble_bi6_polar_cap(&input).unwrap();
    let expected = catmull_clark_star_point(input.n, input.c7, &input.c5s, &input.c6s);
    assert!(approx_eq(out.star_point, expected, 1e-12));
    // Every patch's (0, 0) CP equals the star point.
    for patch in &out.patches_bi6 {
        assert!(approx_eq(patch[0][0], expected, 1e-12));
    }
}

#[test]
fn row0_cps_match_sector_inputs() {
    // The row-0 BB-CPs (i, 0) for i ∈ {1, 2, 4, 5, 6} come
    // straight from the caller-supplied sector data. The
    // (3, 0) CP comes from rule (5) edge reduction.
    let input = pentagon_input();
    let out = assemble_bi6_polar_cap(&input).unwrap();
    for (k, s) in input.sectors.iter().enumerate() {
        let g = &out.patches_bi6[k];
        assert!(approx_eq(g[1][0], s.p_10, 1e-12));
        assert!(approx_eq(g[2][0], s.p_20, 1e-12));
        assert!(approx_eq(g[4][0], s.p_40, 1e-12));
        assert!(approx_eq(g[5][0], s.p_50, 1e-12));
        assert!(approx_eq(g[6][0], s.p_60, 1e-12));
        // (3, 0) is the rule-5 degree reduction of the six
        // edge CPs.
        let expected = rule_5_reduce_edge(s.p_00, s.p_10, s.p_20, s.p_40, s.p_50, s.p_60);
        assert!(approx_eq(g[3][0], expected, 1e-12));
    }
}

#[test]
fn row2_tensor_border_cps_match_sector_inputs() {
    let input = pentagon_input();
    let out = assemble_bi6_polar_cap(&input).unwrap();
    for (k, s) in input.sectors.iter().enumerate() {
        let g = &out.patches_bi6[k];
        assert!(approx_eq(g[4][2], s.p_42, 1e-12));
        assert!(approx_eq(g[5][2], s.p_52, 1e-12));
        assert!(approx_eq(g[6][2], s.p_62, 1e-12));
    }
}

#[test]
fn star_row1_row2_via_s1_inputs() {
    let input = pentagon_input();
    let out = assemble_bi6_polar_cap(&input).unwrap();
    for (k, s) in input.sectors.iter().enumerate() {
        let g = &out.patches_bi6[k];
        // (0, 1), (1, 1), (0, 2) come from the sector's own
        // `p_01, p_11, p_02` — the (S1) rule output stored
        // on the previous sector becomes the (0, 1) of the
        // next sector, but our scaffold stores the per-sector
        // CP directly. Here we just verify we see them.
        assert!(approx_eq(g[0][1], s.p_01, 1e-12));
        assert!(approx_eq(g[1][1], s.p_11, 1e-12));
        assert!(approx_eq(g[0][2], s.p_02, 1e-12));
    }
}

#[test]
fn p22_at_sector_0_equals_pinned_free_param() {
    let input = pentagon_input();
    let out = assemble_bi6_polar_cap(&input).unwrap();
    assert!(approx_eq(out.patches_bi6[0][2][2], input.p0_22, 1e-12));
}

#[test]
fn bi9_elevation_corner_cps_match_bi6() {
    // Degree elevation preserves the 4 corner BB-CPs exactly.
    let input = pentagon_input();
    let out = assemble_bi6_polar_cap(&input).unwrap();
    for k in 0..input.n {
        let bi6 = &out.patches_bi6[k];
        let bi9 = &out.patches_bi9[k];
        assert!(approx_eq(bi9[0][0], bi6[0][0], 1e-10));
        assert!(approx_eq(bi9[0][9], bi6[0][6], 1e-10));
        assert!(approx_eq(bi9[9][0], bi6[6][0], 1e-10));
        assert!(approx_eq(bi9[9][9], bi6[6][6], 1e-10));
    }
}

#[test]
fn all_output_cps_finite() {
    let input = pentagon_input();
    let out = assemble_bi6_polar_cap(&input).unwrap();
    for patch in &out.patches_bi6 {
        for row in patch.iter() {
            for cp in row.iter() {
                assert!(cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite());
            }
        }
    }
    for patch in &out.patches_bi9 {
        for row in patch.iter() {
            for cp in row.iter() {
                assert!(cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite());
            }
        }
    }
}

#[test]
fn hexagon_requires_p0_21_pin() {
    let mut input = pentagon_input();
    input.n = 6;
    input.params = Bi6Reparams::for_valence(6);
    // Rebuild sectors, c5s, c6s for n = 6.
    let n = 6;
    input.c5s.clear();
    input.c6s.clear();
    input.sectors.clear();
    for k in 0..n {
        let theta = (k as f64) * 2.0 * std::f64::consts::PI / (n as f64);
        let dx = theta.cos();
        let dy = theta.sin();
        input.c5s.push(pt(dx * 0.3, dy * 0.3, 0.9));
        input.c6s.push(pt(dx * 0.7, dy * 0.7, 0.6));
        input.sectors.push(CirculantSectorRhs {
            p_00: input.c7,
            p_10: pt(dx * 0.1, dy * 0.1, 0.95),
            p_20: pt(dx * 0.25, dy * 0.25, 0.85),
            p_40: pt(dx * 0.6, dy * 0.6, 0.4),
            p_50: pt(dx * 0.8, dy * 0.8, 0.25),
            p_60: pt(dx, dy, 0.1),
            p_01: pt(-dy * 0.1, dx * 0.1, 0.95),
            p_11: pt((dx - dy) * 0.1, (dy + dx) * 0.1, 0.85),
            p_02: pt(-dy * 0.2, dx * 0.2, 0.88),
            p_42: pt(-dy * 0.5, dx * 0.5, 0.35),
            p_52: pt(-dy * 0.6, dx * 0.6, 0.2),
            p_62: pt(-dy * 0.7, dx * 0.7, 0.05),
        });
    }
    // Without pin at n=6: should fail.
    assert!(assemble_bi6_polar_cap(&input).is_err());
    // With pin: succeeds.
    input.p0_21_pinned = Some(pt(0.0, 0.0, 0.9));
    let out = assemble_bi6_polar_cap(&input).expect("n=6 with pin");
    assert_eq!(out.patches_bi6.len(), 6);
}

#[test]
fn rejects_bad_inputs() {
    let mut input = pentagon_input();
    input.n = 2;
    assert!(assemble_bi6_polar_cap(&input).is_err());

    let mut input = pentagon_input();
    input.sectors.pop();
    assert!(assemble_bi6_polar_cap(&input).is_err());

    let mut input = pentagon_input();
    input.c5s.pop();
    assert!(assemble_bi6_polar_cap(&input).is_err());
}

#[test]
fn interior_filler_is_applied() {
    // Supply an interior filler at (3, 3) = specific value;
    // verify the output uses it.
    let mut input = pentagon_input();
    let sentinel = pt(99.0, 42.0, -7.0);
    let mut filler: Vec<[[Option<Point3>; 7]; 7]> = vec![[[None; 7]; 7]; input.n];
    filler[0][3][3] = Some(sentinel);
    input.interior_filler = Some(filler);
    let out = assemble_bi6_polar_cap(&input).unwrap();
    assert!(approx_eq(out.patches_bi6[0][3][3], sentinel, 1e-12));
}

#[test]
fn via_table_pentagon_assembles_without_error() {
    let input = pentagon_input();
    let out = super::assemble_bi6_polar_cap_via_table(&input, 5, 0.0).unwrap();
    assert_eq!(out.patches_bi6.len(), input.n);
    assert_eq!(out.patches_bi9.len(), input.n);
    for grid in &out.patches_bi9 {
        for row in grid {
            for cp in row {
                assert!(cp.x.is_finite() && cp.y.is_finite() && cp.z.is_finite());
            }
        }
    }
}

#[test]
fn via_table_satisfies_g1_g2_machine_zero_and_lowers_f5_vs_existing() {
    // Production-side wiring check: the new `assemble_bi6_polar_cap_via_table`
    // public entry point must produce caps that satisfy KP15 eq (2)
    // and eq (3) to machine zero at every cross-sector edge AND
    // achieve strictly lower F_5 than the existing
    // `assemble_bi6_polar_cap` path on the same input. These are
    // the two empirical claims of the Phase E ADR.
    use crate::geometry::gregory::kp15_bi6_f5::f_kappa_quadratic_form;
    use crate::geometry::gregory::kp15_bi6_local_constraint::{
        eq2_g1_residual_bb_degree7, eq3_g2_residual_phi_multiplied,
    };

    let input = pentagon_input();
    let n = input.n;
    let q = f_kappa_quadratic_form(5);
    let f5 = |patches: &[[[Point3; 7]; 7]]| -> f64 {
        let mut acc = 0.0;
        for grid in patches {
            for r in 0..49 {
                let (i_r, j_r) = (r / 7, r % 7);
                for c in 0..49 {
                    let (i_c, j_c) = (c / 7, c % 7);
                    let p_r = grid[i_r][j_r];
                    let p_c = grid[i_c][j_c];
                    let qq = q[r][c];
                    if qq == 0.0 {
                        continue;
                    }
                    acc += qq * (p_r.x * p_c.x + p_r.y * p_c.y + p_r.z * p_c.z);
                }
            }
        }
        acc
    };

    let out_legacy = assemble_bi6_polar_cap(&input).unwrap();
    let out_table = super::assemble_bi6_polar_cap_via_table(&input, 5, 0.0).unwrap();

    let f5_legacy = f5(&out_legacy.patches_bi6);
    let f5_table = f5(&out_table.patches_bi6);

    let mut worst_eq2 = 0.0_f64;
    let mut worst_eq3 = 0.0_f64;
    for k in 0..n {
        let km1 = (k + n - 1) % n;
        let bb2 = eq2_g1_residual_bb_degree7(
            &out_table.patches_bi6[km1],
            &out_table.patches_bi6[k],
            &input.params,
        );
        for v in bb2 {
            worst_eq2 = worst_eq2.max(v.length());
        }
        let bb3 = eq3_g2_residual_phi_multiplied(
            &out_table.patches_bi6[km1],
            &out_table.patches_bi6[k],
            &input.params,
        );
        for v in &bb3 {
            worst_eq3 = worst_eq3.max(v.length());
        }
    }

    eprintln!(
        "Via-table pentagon: F_5 legacy = {f5_legacy:.3e}, F_5 table = {f5_table:.3e}, \
         worst eq(2) = {worst_eq2:.3e}, worst eq(3) = {worst_eq3:.3e}"
    );

    assert!(
        worst_eq2 < 1e-8,
        "table path: worst eq (2) residual = {worst_eq2:.3e}"
    );
    assert!(
        worst_eq3 < 1e-8,
        "table path: worst eq (3) residual = {worst_eq3:.3e}"
    );
    assert!(
        f5_table < f5_legacy,
        "table F_5 = {f5_table:.3e} ≥ legacy F_5 = {f5_legacy:.3e}"
    );
}
