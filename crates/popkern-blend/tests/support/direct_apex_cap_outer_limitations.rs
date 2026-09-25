//! Recorded-limitation arm split out of the direct-apex Gregory
//! cap-outer set (see `direct_apex_validation_matrix.rs`).

use crate::direct_apex_validate_helpers::{build_under_direct_apex, native_cap_faces};
use crate::helpers::measure_cap_outer_boundary;

/// pyramid4_mixed's cap-outer arm of the Gregory set, separated so the
/// cube/pyramid3 coverage stays gated while this class's G1 limitation
/// is recorded.
#[test]
fn direct_apex_5_cap_outer_pyramid4_mixed_g1_limitation() {
    let result = build_under_direct_apex(true, &[0.04, 0.06, 0.05, 0.07]);
    let cap_faces = native_cap_faces(&result).into_iter().collect::<Vec<_>>();
    let m = measure_cap_outer_boundary(&result.brep, &cap_faces, 100);
    assert!(
        m.worst_g0_gap <= 5.0e-3 && m.worst_g1_normal <= 5.0e-2,
        "pyramid4_mixed cap_outer g0={:.3e} g1={:.3e} vs gates 5e-3/5e-2",
        m.worst_g0_gap,
        m.worst_g1_normal,
    );
}
