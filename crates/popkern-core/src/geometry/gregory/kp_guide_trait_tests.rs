//! Tests for the [`super::KpGuide`] trait: the bi-quintic impl
//! derived from a degree-elevated bicubic must produce
//! identical trait outputs at every sample point. If this holds,
//! `derive_row1_frame` / `derive_row2_frame` called against the
//! bi-quintic produce the same bidegree-9 patch CPs as when
//! called against the bicubic — the "no regression" invariant
//! of step 2.

use super::*;
use crate::geometry::{Point3, Vec3};

fn bumpy_bicubic() -> KpGuidePatch {
    let mut g = KpGuidePatch::zeros();
    for i in 0..4 {
        for j in 0..4 {
            let s = i as f64 / 3.0;
            let t = j as f64 / 3.0;
            g.grid[i][j] = Point3::new(
                s,
                t,
                0.4 * (s * (1.0 - s)) * (t * (1.0 - t)) + 0.15 * (s - 0.5) * (t - 0.5),
            );
        }
    }
    g
}

fn vec3_max(a: Vec3, b: Vec3) -> f64 {
    ((a.x - b.x).abs())
        .max((a.y - b.y).abs())
        .max((a.z - b.z).abs())
}

fn cps10_max(a: [Vec3; 10], b: [Vec3; 10]) -> f64 {
    (0..10)
        .map(|i| vec3_max(a[i], b[i]))
        .fold(0.0_f64, f64::max)
}

#[test]
fn trait_outputs_match_between_bicubic_and_elevated_biquintic() {
    let bi = bumpy_bicubic();
    let bq = KpGuidePatchBiQuintic::from_bicubic_elevated(&bi);

    let uvs = [(0.0, 0.0), (0.3, 0.4), (0.5, 0.5), (0.9, 0.8), (1.0, 1.0)];
    for &(u, v) in &uvs {
        let ba = KpGuide::partial_u(&bi, u, v);
        let bqa = KpGuide::partial_u(&bq, u, v);
        assert!(vec3_max(ba, bqa) < 1e-12, "partial_u at ({u},{v})");
        let ba = KpGuide::partial_v(&bi, u, v);
        let bqa = KpGuide::partial_v(&bq, u, v);
        assert!(vec3_max(ba, bqa) < 1e-12);
        let ba = KpGuide::second_partial_uu(&bi, u, v);
        let bqa = KpGuide::second_partial_uu(&bq, u, v);
        assert!(vec3_max(ba, bqa) < 1e-11);
        let ba = KpGuide::second_partial_uv(&bi, u, v);
        let bqa = KpGuide::second_partial_uv(&bq, u, v);
        assert!(vec3_max(ba, bqa) < 1e-11);
        let ba = KpGuide::second_partial_vv(&bi, u, v);
        let bqa = KpGuide::second_partial_vv(&bq, u, v);
        assert!(vec3_max(ba, bqa) < 1e-11);
    }

    // Cross-derivative CPs at degree 9 must match. Both paths
    // elevate to the same 10-CP representation of the same curve.
    assert!(
        cps10_max(
            KpGuide::cross_deriv_cps_u0_degree9(&bi),
            KpGuide::cross_deriv_cps_u0_degree9(&bq)
        ) < 1e-10
    );
    assert!(
        cps10_max(
            KpGuide::cross_deriv_cps_u1_degree9(&bi),
            KpGuide::cross_deriv_cps_u1_degree9(&bq)
        ) < 1e-10
    );
    assert!(
        cps10_max(
            KpGuide::cross_deriv_cps_v0_degree9(&bi),
            KpGuide::cross_deriv_cps_v0_degree9(&bq)
        ) < 1e-10
    );
    assert!(
        cps10_max(
            KpGuide::cross_deriv_cps_v1_degree9(&bi),
            KpGuide::cross_deriv_cps_v1_degree9(&bq)
        ) < 1e-10
    );
    assert!(
        cps10_max(
            KpGuide::cross_second_deriv_cps_uu_u0_degree9(&bi),
            KpGuide::cross_second_deriv_cps_uu_u0_degree9(&bq)
        ) < 1e-10
    );
    assert!(
        cps10_max(
            KpGuide::cross_second_deriv_cps_uu_u1_degree9(&bi),
            KpGuide::cross_second_deriv_cps_uu_u1_degree9(&bq)
        ) < 1e-10
    );
    assert!(
        cps10_max(
            KpGuide::cross_second_deriv_cps_vv_v0_degree9(&bi),
            KpGuide::cross_second_deriv_cps_vv_v0_degree9(&bq)
        ) < 1e-10
    );
    assert!(
        cps10_max(
            KpGuide::cross_second_deriv_cps_vv_v1_degree9(&bi),
            KpGuide::cross_second_deriv_cps_vv_v1_degree9(&bq)
        ) < 1e-10
    );
}
