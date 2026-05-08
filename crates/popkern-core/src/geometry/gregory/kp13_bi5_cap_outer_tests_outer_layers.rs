#![cfg(test)]
use super::tests_support::*;
use super::*;

#[test]
fn full_layer_transfer_zero_target_leaves_cap_zero() {
    for n in [3, 4, 5, 6, 7, 8] {
        let mut cap = Bi5Cap::zero(n);
        let target = TensorBorder::zero(n);
        populate_outer_layers_full_identity(&mut cap, &target);
        cap.validate().unwrap();
        for k in 0..n {
            for id in CapPieceId::all() {
                let piece = cap.sector(k).piece(id);
                for p in 0..6 {
                    for q in 0..6 {
                        let cp = piece.cps[p][q];
                        assert!(
                            cp.x.abs() < 1e-10 && cp.y.abs() < 1e-10 && cp.z.abs() < 1e-10,
                            "n={n}, k={k}, id={id:?}, ({p}, {q}): non-zero {cp:?}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn full_layer_transfer_leaves_o_piece_zero() {
    let n = 5;
    let mut cap = Bi5Cap::zero(n);
    let mut target = TensorBorder::zero(n);
    for k in 0..n {
        for i in 0..4 {
            for j in 0..4 {
                target.patches[k].cps[i][j] =
                    Point3::new(k as f64 + i as f64 * 0.1, j as f64 * 0.2, 0.0);
            }
        }
    }
    populate_outer_layers_full_identity(&mut cap, &target);

    for k in 0..n {
        let o = cap.sector(k).piece(CapPieceId::O);
        for p in 0..6 {
            for q in 0..6 {
                let cp = o.cps[p][q];
                assert!(
                    cp.x.abs() < 1e-10 && cp.y.abs() < 1e-10 && cp.z.abs() < 1e-10,
                    "n={n}, k={k}, O piece ({p}, {q}): non-zero {cp:?}"
                );
            }
        }
    }
}

#[test]
fn full_layer_transfer_sets_exactly_outer_two_layers() {
    // For each non-o piece, layers adjacent to the outer boundary
    // (i ∈ {4, 5} for r/m, j ∈ {4, 5} for l/m) should be set; other
    // CPs remain zero.
    let n = 6;
    let mut cap = Bi5Cap::zero(n);
    let mut target = TensorBorder::zero(n);
    for k in 0..n {
        for i in 0..4 {
            for j in 0..4 {
                target.patches[k].cps[i][j] =
                    Point3::new(1.0 + i as f64 + j as f64, 2.0 * (i + j) as f64, k as f64);
            }
        }
    }
    populate_outer_layers_full_identity(&mut cap, &target);

    for k in 0..n {
        // R piece: set at i ∈ {4, 5} (any j). Other i stay zero.
        let r = cap.sector(k).piece(CapPieceId::R);
        for i in 0..4 {
            for j in 0..6 {
                let cp = r.cps[i][j];
                assert!(
                    cp.x.abs() < 1e-10 && cp.y.abs() < 1e-10 && cp.z.abs() < 1e-10,
                    "k={k}, R ({i}, {j}): should be zero, got {cp:?}"
                );
            }
        }
        // L piece: set at j ∈ {4, 5}. Other j stay zero.
        let l = cap.sector(k).piece(CapPieceId::L);
        for i in 0..6 {
            for j in 0..4 {
                let cp = l.cps[i][j];
                assert!(
                    cp.x.abs() < 1e-10 && cp.y.abs() < 1e-10 && cp.z.abs() < 1e-10,
                    "k={k}, L ({i}, {j}): should be zero, got {cp:?}"
                );
            }
        }
        // M piece: set at i ∈ {4, 5} OR j ∈ {4, 5}. Others stay zero.
        let m = cap.sector(k).piece(CapPieceId::M);
        for i in 0..4 {
            for j in 0..4 {
                let cp = m.cps[i][j];
                assert!(
                    cp.x.abs() < 1e-10 && cp.y.abs() < 1e-10 && cp.z.abs() < 1e-10,
                    "k={k}, M ({i}, {j}): should be zero, got {cp:?}"
                );
            }
        }
    }
}

#[test]
fn full_layer_transfer_matches_bi5_evaluation_at_outer_points() {
    // Under β = identity, each non-o piece's outer boundary should
    // evaluate to the same 3D point as the tensor-border b^k at
    // the sector-local coordinate. Verify at outer-edge corners.
    use super::super::super::kp13_bi5_boundary::Bi3TensorBorder;

    let n = 5;
    let mut cap = Bi5Cap::zero(n);
    let mut target = TensorBorder::zero(n);
    // Set up a non-trivial bi-3 target.
    for k in 0..n {
        for i in 0..4 {
            for j in 0..4 {
                target.patches[k].cps[i][j] = Point3::new(
                    (i as f64) + (j as f64) * 0.5 + (k as f64) * 0.1,
                    (i as f64) * 0.3 - (j as f64) * 0.2,
                    (k as f64) * (i + j) as f64 * 0.05,
                );
            }
        }
    }
    populate_outer_layers_full_identity(&mut cap, &target);

    for k in 0..n {
        let sector = cap.sector(k);
        // r piece at piece-local (1, 0) = sector (1, 0) = m_k.
        let r = sector.piece(CapPieceId::R);
        let r10 = r.eval(1.0, 0.0);
        let b_mk = eval_bi3(&target.patches[k], 1.0, 0.0);
        assert!(
            approx_p(r10, b_mk, 1e-10),
            "k={k}: r(1,0)={r10:?} vs b(1,0)={b_mk:?}"
        );
        // m piece at piece-local (1, 1) = sector (1, 1) = v_k.
        let m = sector.piece(CapPieceId::M);
        let m11 = m.eval(1.0, 1.0);
        let b_vk = eval_bi3(&target.patches[k], 1.0, 1.0);
        assert!(
            approx_p(m11, b_vk, 1e-10),
            "k={k}: m(1,1)={m11:?} vs b(1,1)={b_vk:?}"
        );
        // l piece at piece-local (0, 1) = sector (0, 1) = m_{k+1}.
        let l = sector.piece(CapPieceId::L);
        let l01 = l.eval(0.0, 1.0);
        let b_mk1 = eval_bi3(&target.patches[k], 0.0, 1.0);
        assert!(
            approx_p(l01, b_mk1, 1e-10),
            "k={k}: l(0,1)={l01:?} vs b(0,1)={b_mk1:?}"
        );
    }

    fn eval_bi3(p: &Bi3TensorBorder, u: f64, v: f64) -> Point3 {
        let bu = [
            (1.0 - u).powi(3),
            3.0 * u * (1.0 - u).powi(2),
            3.0 * u * u * (1.0 - u),
            u * u * u,
        ];
        let bv = [
            (1.0 - v).powi(3),
            3.0 * v * (1.0 - v).powi(2),
            3.0 * v * v * (1.0 - v),
            v * v * v,
        ];
        let mut out = Point3::new(0.0, 0.0, 0.0);
        for i in 0..4 {
            for j in 0..4 {
                let w = bu[i] * bv[j];
                out.x += w * p.cps[i][j].x;
                out.y += w * p.cps[i][j].y;
                out.z += w * p.cps[i][j].z;
            }
        }
        out
    }
}

#[test]
fn full_layer_transfer_matches_bi3_along_u1_edge_of_r_piece() {
    // The r piece's outer boundary u=1 corresponds to the tensor
    // border's edge u=1 (sector-local). Sample along this edge.
    use super::super::super::kp13_bi5_boundary::Bi3TensorBorder;

    let n = 3;
    let mut cap = Bi5Cap::zero(n);
    let mut target = TensorBorder::zero(n);
    for k in 0..n {
        for i in 0..4 {
            for j in 0..4 {
                target.patches[k].cps[i][j] =
                    Point3::new((i * j) as f64, (i + j) as f64, k as f64 * 0.5);
            }
        }
    }
    populate_outer_layers_full_identity(&mut cap, &target);

    for k in 0..n {
        let r = cap.sector(k).piece(CapPieceId::R);
        let b = &target.patches[k];
        // piece-local (1, v) corresponds to sector (1, v/2) = b(1, v/2).
        for &v in &[0.0_f64, 0.25, 0.5, 0.75, 1.0] {
            let r_val = r.eval(1.0, v);
            let b_val = eval_bi3_local(b, 1.0, v / 2.0);
            assert!(
                approx_p(r_val, b_val, 1e-9),
                "k={k}, v={v}: r(1, v)={r_val:?} vs b(1, v/2)={b_val:?}"
            );
        }
    }

    fn eval_bi3_local(p: &Bi3TensorBorder, u: f64, v: f64) -> Point3 {
        let bu = [
            (1.0 - u).powi(3),
            3.0 * u * (1.0 - u).powi(2),
            3.0 * u * u * (1.0 - u),
            u * u * u,
        ];
        let bv = [
            (1.0 - v).powi(3),
            3.0 * v * (1.0 - v).powi(2),
            3.0 * v * v * (1.0 - v),
            v * v * v,
        ];
        let mut out = Point3::new(0.0, 0.0, 0.0);
        for i in 0..4 {
            for j in 0..4 {
                let w = bu[i] * bv[j];
                out.x += w * p.cps[i][j].x;
                out.y += w * p.cps[i][j].y;
                out.z += w * p.cps[i][j].z;
            }
        }
        out
    }
}

#[test]
fn full_layer_transfer_deterministic() {
    let n = 5;
    let mut target = TensorBorder::zero(n);
    for k in 0..n {
        target.patches[k].cps[2][3] = Point3::new(1.0, 2.0, 3.0);
        target.patches[k].cps[3][2] = Point3::new(-1.0, 0.5, 2.5);
    }
    let mut cap_a = Bi5Cap::zero(n);
    let mut cap_b = Bi5Cap::zero(n);
    populate_outer_layers_full_identity(&mut cap_a, &target);
    populate_outer_layers_full_identity(&mut cap_b, &target);
    assert_eq!(cap_a, cap_b);
}
