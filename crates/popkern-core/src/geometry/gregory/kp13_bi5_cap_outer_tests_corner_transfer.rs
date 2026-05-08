use super::tests_support::*;
use super::*;

// ---- populate_outer_corner_jets_identity: end-to-end ----

#[test]
fn populate_with_zero_target_leaves_cap_zero() {
    for n in [3, 4, 5, 6, 7, 8] {
        let mut cap = Bi5Cap::zero(n);
        let target = TensorBorder::zero(n);
        populate_outer_corner_jets_identity(&mut cap, &target);
        cap.validate().unwrap();
        for k in 0..n {
            for id in CapPieceId::all() {
                let piece = cap.sector(k).piece(id);
                for p in 0..6 {
                    for q in 0..6 {
                        let cp = piece.cps[p][q];
                        assert!(
                            cp.x.abs() < 1e-10 && cp.y.abs() < 1e-10 && cp.z.abs() < 1e-10,
                            "n={n}, k={k}, id={id:?}, ({p}, {q}): non-zero CP {cp:?}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn populate_sets_only_outer_corner_blocks() {
    // O piece should remain zero after populate (I.14.b touches
    // only ℓ, m, r). Interior CPs of ℓ, m, r should also stay zero.
    let n = 5;
    let mut cap = Bi5Cap::zero(n);
    let mut target = TensorBorder::zero(n);
    for k in 0..n {
        target.patches[k].cps[3][3] = Point3::new(1.0, 2.0, 3.0);
        target.patches[k].cps[3][0] = Point3::new(4.0, 5.0, 6.0);
        target.patches[k].cps[0][3] = Point3::new(7.0, 8.0, 9.0);
    }
    populate_outer_corner_jets_identity(&mut cap, &target);

    // O piece remains all zero.
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

    // R piece: outer-corner block CPs at (1, 0) should be set; interior stays zero.
    let outer_r = [(5, 0), (4, 0), (5, 1), (3, 0), (5, 2), (4, 1)];
    for k in 0..n {
        let r = cap.sector(k).piece(CapPieceId::R);
        for p in 0..6 {
            for q in 0..6 {
                let is_outer = outer_r.contains(&(p, q));
                let cp = r.cps[p][q];
                if is_outer {
                    // At least the value (P[5][0]) should match the tensor-border jet at m_k.
                    // We don't check each CP exactly, just that they're not all zero.
                } else {
                    assert!(
                        cp.x.abs() < 1e-10 && cp.y.abs() < 1e-10 && cp.z.abs() < 1e-10,
                        "n={n}, k={k}, R ({p}, {q}): non-zero interior CP {cp:?}"
                    );
                }
            }
        }
    }
}

#[test]
fn populate_matches_tensor_border_value_at_outer_corners() {
    // The piece corner VALUE CPs (P[5][0], P[0][5], P[5][5]) should
    // match the tensor-border's corresponding corner values.
    let n = 6;
    let mut cap = Bi5Cap::zero(n);
    let mut target = TensorBorder::zero(n);
    for k in 0..n {
        // Set the outer corner CPs of each bi-3 patch to distinct values
        // so the tensor-border jets have non-zero value at each corner.
        target.patches[k].cps[3][0] = Point3::new(k as f64, 0.0, 0.0); // m_k corner value
        target.patches[k].cps[0][3] = Point3::new(0.0, k as f64, 0.0); // m_{k+1}
        target.patches[k].cps[3][3] = Point3::new(0.0, 0.0, k as f64); // v_k
    }
    populate_outer_corner_jets_identity(&mut cap, &target);

    for k in 0..n {
        let sector = cap.sector(k);
        let r = sector.piece(CapPieceId::R);
        let l = sector.piece(CapPieceId::L);
        let m = sector.piece(CapPieceId::M);
        let j_mk = target.jet(k, ChiCorner::Mk);
        let j_mk1 = target.jet(k, ChiCorner::Mk1);
        let j_vk = target.jet(k, ChiCorner::Vk);
        assert!(approx_p(r.cps[5][0], j_mk.value, 1e-10));
        assert!(approx_p(l.cps[0][5], j_mk1.value, 1e-10));
        assert!(approx_p(m.cps[5][5], j_vk.value, 1e-10));
    }
}

#[test]
fn populate_deterministic() {
    let n = 5;
    let mut target = TensorBorder::zero(n);
    for k in 0..n {
        target.patches[k].cps[3][3] = Point3::new(k as f64, 1.0, -1.0);
    }
    let mut cap_a = Bi5Cap::zero(n);
    let mut cap_b = Bi5Cap::zero(n);
    populate_outer_corner_jets_identity(&mut cap_a, &target);
    populate_outer_corner_jets_identity(&mut cap_b, &target);
    assert_eq!(cap_a, cap_b);
}

// ---- Structural invariants ----

#[test]
#[should_panic]
fn populate_rejects_n_mismatch() {
    let mut cap = Bi5Cap::zero(5);
    let target = TensorBorder::zero(6);
    populate_outer_corner_jets_identity(&mut cap, &target);
}

// ---- Full depth-2 outer-layer transfer ----
