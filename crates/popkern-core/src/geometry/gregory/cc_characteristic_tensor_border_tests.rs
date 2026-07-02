use super::*;
use crate::geometry::gregory::cc_characteristic_ring::cc_characteristic_ring;

const TOL: f64 = 1.0e-12;

fn close(a: Point3, b: [f64; 2], label: &str) {
    assert!(
        (a.x - b[0]).abs() < TOL && (a.y - b[1]).abs() < TOL && a.z.abs() < TOL,
        "{label}: got ({}, {}, {}), expected ({}, {})",
        a.x,
        a.y,
        a.z,
        b[0],
        b[1],
    );
}

/// Uniform cubic B-spline basis on one span.
fn bspline_basis(t: f64) -> [f64; 4] {
    let s = 1.0 - t;
    [
        s * s * s / 6.0,
        (3.0 * t * t * t - 6.0 * t * t + 4.0) / 6.0,
        (-3.0 * t * t * t + 3.0 * t * t + 3.0 * t + 1.0) / 6.0,
        t * t * t / 6.0,
    ]
}

fn bb_eval_cubic(row: &[Point3; 4], t: f64) -> [f64; 2] {
    let s = 1.0 - t;
    let w = [s * s * s, 3.0 * s * s * t, 3.0 * s * t * t, t * t * t];
    let mut acc = [0.0_f64; 2];
    for (wi, cp) in w.iter().zip(row) {
        acc[0] += wi * cp.x;
        acc[1] += wi * cp.y;
    }
    acc
}

#[test]
fn spoke_endpoints_land_on_unit_ngon_vertices_counterclockwise() {
    for n in [3, 4, 5, 6, 7, 8, 10] {
        let border = cc_characteristic_tensor_border(n).unwrap();
        assert_eq!(border.sectors.len(), n);
        let alpha = core::f64::consts::TAU / n as f64;
        for (k, sector) in border.sectors.iter().enumerate() {
            let a0 = k as f64 * alpha;
            let a1 = (k + 1) as f64 * alpha;
            close(
                sector.near_c.rows[0][0],
                [a0.cos(), a0.sin()],
                &format!("n={n} sector {k} m^k"),
            );
            close(
                sector.near_m.rows[0][3],
                [a1.cos(), a1.sin()],
                &format!("n={n} sector {k} m^(k+1)"),
            );
            // The two pieces share the sector corner c.
            let c_end = sector.near_c.rows[0][3];
            let c_start = sector.near_m.rows[0][0];
            close(c_start, [c_end.x, c_end.y], &format!("n={n} sector {k} c"));
        }
    }
}

#[test]
fn sector_corner_matches_self_similar_ring_value() {
    // c = chi(1/2, 1/2) = lambda * chi(1, 1): the corner must equal
    // the normalized lambda * f(1,1,0), an independent path through
    // the Omega_2 grid (self-similarity of the characteristic map).
    for n in [3, 5, 6, 7, 8] {
        let border = cc_characteristic_tensor_border(n).unwrap();
        let ring = cc_characteristic_ring(n).unwrap();

        // Recompute the normalization exactly as the extractor does.
        let mask = [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0];
        let corner_of = |grid: &[[[f64; 2]; 4]; 4]| -> [f64; 2] {
            let mut acc = [0.0_f64; 2];
            for (r, wr) in mask.iter().enumerate() {
                for (c, wc) in mask.iter().enumerate() {
                    acc[0] += wr * wc * grid[r][c][0];
                    acc[1] += wr * wc * grid[r][c][1];
                }
            }
            acc
        };
        let m_u = corner_of(&ring.grids[0]);
        let m_v = corner_of(&ring.grids[2]);
        let norm_sq = m_u[0] * m_u[0] + m_u[1] * m_u[1];
        let zeta = [m_u[0] / norm_sq, -m_u[1] / norm_sq];
        let mul = |a: [f64; 2], b: [f64; 2]| -> [f64; 2] {
            [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
        };
        let conj = mul(zeta, m_v)[1] < 0.0;
        let mut expected = mul(zeta, ring.f110);
        expected = [ring.lambda * expected[0], ring.lambda * expected[1]];
        if conj {
            expected[1] = -expected[1];
        }
        assert!(
            (border.corner[0] - expected[0]).abs() < TOL
                && (border.corner[1] - expected[1]).abs() < TOL,
            "n={n}: corner {:?} != lambda*f110 normalized {:?}",
            border.corner,
            expected,
        );
    }
}

#[test]
fn border_pieces_reproduce_ring_jets_across_the_boundary() {
    // Decisive scale/orientation gate: at border samples, the bi-3
    // block's value, first, and second depth-derivatives must equal
    // the ring surface's B-spline jets across the edge (with the
    // cap-side depth pointing OPPOSITE the ring-side parameter, so
    // odd orders flip sign).
    for n in [3, 5, 6, 7, 8] {
        let border = cc_characteristic_tensor_border(n).unwrap();
        let ring = cc_characteristic_ring(n).unwrap();

        // Recompute normalization (as in the extractor).
        let mask = [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0];
        let corner_of = |grid: &[[[f64; 2]; 4]; 4]| -> [f64; 2] {
            let mut acc = [0.0_f64; 2];
            for (r, wr) in mask.iter().enumerate() {
                for (c, wc) in mask.iter().enumerate() {
                    acc[0] += wr * wc * grid[r][c][0];
                    acc[1] += wr * wc * grid[r][c][1];
                }
            }
            acc
        };
        let m_u = corner_of(&ring.grids[0]);
        let m_v = corner_of(&ring.grids[2]);
        let norm_sq = m_u[0] * m_u[0] + m_u[1] * m_u[1];
        let zeta = [m_u[0] / norm_sq, -m_u[1] / norm_sq];
        let mul = |a: [f64; 2], b: [f64; 2]| -> [f64; 2] {
            [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
        };
        let conj = mul(zeta, m_v)[1] < 0.0;
        let normalize = |v: [f64; 2]| -> [f64; 2] {
            let s = mul(zeta, v);
            if conj { [s[0], -s[1]] } else { s }
        };

        // Ring-side jets on Omega_1's s = 0 edge (E_u), at along
        // parameter t: value/d/ds/d2/ds2 masks over the grid columns.
        let value_mask = [1.0 / 6.0, 4.0 / 6.0, 1.0 / 6.0, 0.0];
        let dds_mask = [-0.5, 0.0, 0.5, 0.0];
        let d2ds2_mask = [1.0, -2.0, 1.0, 0.0];
        let piece = &border.sectors[0].near_c;
        for t in [0.2, 0.5, 0.8] {
            let basis = bspline_basis(t);
            let mut jet = [[0.0_f64; 2]; 3];
            for (r, br) in basis.iter().enumerate() {
                for c in 0..4 {
                    let v = normalize(ring.grids[0][r][c]);
                    jet[0][0] += br * value_mask[c] * v[0];
                    jet[0][1] += br * value_mask[c] * v[1];
                    jet[1][0] += br * dds_mask[c] * v[0];
                    jet[1][1] += br * dds_mask[c] * v[1];
                    jet[2][0] += br * d2ds2_mask[c] * v[0];
                    jet[2][1] += br * d2ds2_mask[c] * v[1];
                }
            }

            let r0 = bb_eval_cubic(&piece.rows[0], t);
            let r1 = bb_eval_cubic(&piece.rows[1], t);
            let r2 = bb_eval_cubic(&piece.rows[2], t);
            // Cap-side bi-3 depth derivatives at the boundary.
            let d1 = [3.0 * (r1[0] - r0[0]), 3.0 * (r1[1] - r0[1])];
            let d2 = [
                6.0 * (r2[0] - 2.0 * r1[0] + r0[0]),
                6.0 * (r2[1] - 2.0 * r1[1] + r0[1]),
            ];
            for axis in 0..2 {
                assert!(
                    (r0[axis] - jet[0][axis]).abs() < TOL,
                    "n={n} t={t}: border value != ring value",
                );
                assert!(
                    (d1[axis] + jet[1][axis]).abs() < TOL,
                    "n={n} t={t}: first depth derivative must be the ring's, reflected",
                );
                assert!(
                    (d2[axis] - jet[2][axis]).abs() < TOL,
                    "n={n} t={t}: second depth derivative must match the ring's",
                );
            }
        }
    }
}

#[test]
fn valence_four_border_is_the_exact_identity_lattice() {
    // n = 4: chi is the identity map, so in normalized coordinates the
    // near_c piece of sector 0 is exactly rows[j][a] = (1 - j/3, a/3).
    let border = cc_characteristic_tensor_border(4).unwrap();
    let piece = &border.sectors[0].near_c;
    for (j, row) in piece.rows.iter().enumerate() {
        for (a, cp) in row.iter().enumerate() {
            close(
                *cp,
                [1.0 - j as f64 / 3.0, a as f64 / 3.0],
                &format!("n=4 near_c[{j}][{a}]"),
            );
        }
    }
    let piece = &border.sectors[0].near_m;
    for (j, row) in piece.rows.iter().enumerate() {
        for (a, cp) in row.iter().enumerate() {
            close(
                *cp,
                [1.0 - a as f64 / 3.0, 1.0 - j as f64 / 3.0],
                &format!("n=4 near_m[{j}][{a}]"),
            );
        }
    }
}

#[test]
fn sectors_are_rotations_and_pieces_are_diagonal_mirrors() {
    for n in [3, 5, 6, 8] {
        let border = cc_characteristic_tensor_border(n).unwrap();
        let alpha = core::f64::consts::TAU / n as f64;
        let s0 = &border.sectors[0];
        for (k, sector) in border.sectors.iter().enumerate() {
            let rot = [(k as f64 * alpha).cos(), (k as f64 * alpha).sin()];
            for j in 0..3 {
                for a in 0..4 {
                    let base = s0.near_c.rows[j][a];
                    let expect = [
                        base.x * rot[0] - base.y * rot[1],
                        base.x * rot[1] + base.y * rot[0],
                    ];
                    close(
                        sector.near_c.rows[j][a],
                        expect,
                        &format!("n={n} sector {k} rotation"),
                    );
                }
            }
        }
        // Diagonal symmetry: reflecting across the sector bisector
        // (z -> w * conj(z), w = e^{2*pi*i/n}) maps the near_c piece
        // reversed onto near_m.
        let w = [alpha.cos(), alpha.sin()];
        for j in 0..3 {
            for a in 0..4 {
                let base = s0.near_c.rows[j][3 - a];
                let reflected = [w[0] * base.x + w[1] * base.y, w[1] * base.x - w[0] * base.y];
                close(
                    s0.near_m.rows[j][a],
                    reflected,
                    &format!("n={n} diagonal mirror [{j}][{a}]"),
                );
            }
        }
    }
}

#[test]
fn rejects_invalid_valence() {
    assert!(cc_characteristic_tensor_border(2).is_err());
}
