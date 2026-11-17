//! The shared second-form evaluation against the per-pair one (PK-1223).

use super::*;

fn series(seed: f64) -> Series {
    let mut c = [Iv::point(0.0); SEAM_TERMS];
    for (k, slot) in c.iter_mut().enumerate() {
        let mid = (seed * (k as f64 + 1.3)).sin() * 2.0 + if k == 0 { 3.0 } else { 0.0 };
        *slot = Iv::new(mid - 1e-9 * (k as f64 + 1.0), mid + 2e-9 * (k as f64 + 1.0));
    }
    Taylor::new(c)
}

fn vector(seed: f64) -> SeriesVec {
    TaylorVec3 {
        x: series(seed),
        y: series(seed * 1.7 + 0.4),
        z: series(seed * 2.3 - 0.9),
    }
}

fn bits(s: &Series) -> Vec<(u64, u64)> {
    s.c.iter()
        .map(|i| (i.lo.to_bits(), i.hi.to_bits()))
        .collect()
}

#[test]
fn shared_second_forms_are_the_per_pair_forms_bit_for_bit() {
    for seed in [0.31, 1.7, 4.05] {
        let side = SideState {
            p: vector(seed),
            su: vector(seed + 0.11),
            sv: vector(seed + 2.9),
            suu: vector(seed + 0.5),
            suv: vector(seed + 1.4),
            svv: vector(seed + 3.3),
            tangent: vector(seed + 0.8),
            umbilic: None,
        };
        let (unit, e1, e2) = (vector(seed + 5.0), vector(seed + 6.1), vector(seed + 7.2));
        let shared = second_forms(&side, unit, e1, e2);
        let each = [
            second_form(&side, unit, e1, e1),
            second_form(&side, unit, e1, e2),
            second_form(&side, unit, e2, e2),
        ];
        for (a, b) in shared.iter().zip(&each) {
            assert_eq!(bits(a), bits(b));
        }
    }
}
