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
