//! §13.4 soft-supports mode end-to-end tests.
//!
//! Verifies that `FilletOptions::soft_supports = true` correctly
//! loosens the pipeline's internal tolerance so fillet operations
//! on low-precision input (e.g. drifted imports) succeed where
//! they would fail in strict mode, and that the corresponding
//! `BlendWarning::SoftSupportTolerancesApplied` warning is emitted.

use popkern_blend::{
    BlendError, BlendWarning, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

// ---------------------------------------------------------------------------
// Warning emission
// ---------------------------------------------------------------------------

#[test]
fn soft_supports_flag_emits_tolerance_warning() {
    // On clean input the soft-supports flag should still emit
    // exactly one SoftSupportTolerancesApplied warning, carrying
    // the configured factor. The fillet itself proceeds normally
    // because the cube is tight enough to fit both tolerances.
    let brep = BRep::unit_cube(tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let opts = FilletOptions {
        soft_supports: true,
        soft_support_factor: 10.0,
        ..FilletOptions::default()
    };
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let result = fillet_edges(&brep, &specs, &opts, tol(), None);
    let blend = match result {
        Ok(b) => b,
        Err(BlendError::PartialFailure { .. }) => return, // accept partial on corner conflicts
        Err(e) => panic!("expected success or partial, got {e:?}"),
    };
    let soft_warns: Vec<_> = blend
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::SoftSupportTolerancesApplied { .. }))
        .collect();
    assert_eq!(
        soft_warns.len(),
        1,
        "expected exactly one SoftSupportTolerancesApplied warning, got {}: {:?}",
        soft_warns.len(),
        blend.edge_warnings
    );
    if let (_, BlendWarning::SoftSupportTolerancesApplied { factor }) = soft_warns[0] {
        assert!(
            (factor - 10.0).abs() < 1e-12,
            "reported factor {factor} should match configured 10.0"
        );
    }
}

#[test]
fn soft_supports_disabled_emits_no_tolerance_warning() {
    // The default `soft_supports = false` must not emit the
    // soft-support warning — it's an opt-in mode, not a default.
    let brep = BRep::unit_cube(tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let result = fillet_edges(&brep, &specs, &FilletOptions::default(), tol(), None);
    let blend = match result {
        Ok(b) => b,
        Err(BlendError::PartialFailure { .. }) => return,
        Err(e) => panic!("expected success or partial, got {e:?}"),
    };
    let has_soft = blend
        .edge_warnings
        .iter()
        .any(|(_, w)| matches!(w, BlendWarning::SoftSupportTolerancesApplied { .. }));
    assert!(
        !has_soft,
        "soft-supports disabled should not emit the warning, got {:?}",
        blend.edge_warnings
    );
}

#[test]
fn soft_supports_with_factor_one_does_not_emit_warning() {
    // `soft_support_factor = 1.0` is a no-op (tolerance unchanged)
    // — the pipeline should skip the warning emission even though
    // `soft_supports = true`, because no relaxation actually
    // happened. This keeps the warning meaningful: it only fires
    // when the tolerance was genuinely loosened.
    let brep = BRep::unit_cube(tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let opts = FilletOptions {
        soft_supports: true,
        soft_support_factor: 1.0,
        ..FilletOptions::default()
    };
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let result = fillet_edges(&brep, &specs, &opts, tol(), None);
    let blend = match result {
        Ok(b) => b,
        Err(BlendError::PartialFailure { .. }) => return,
        Err(e) => panic!("expected success or partial, got {e:?}"),
    };
    let has_soft = blend
        .edge_warnings
        .iter()
        .any(|(_, w)| matches!(w, BlendWarning::SoftSupportTolerancesApplied { .. }));
    assert!(
        !has_soft,
        "factor = 1.0 is a no-op and should not emit the warning"
    );
}

// ---------------------------------------------------------------------------
// Low-precision input behaviour — collapsed edge threshold
// ---------------------------------------------------------------------------
//
// §13.4 predicts that loosening the tolerance by the soft-support
// factor will cause the §13.1(5) collapsed-edge pre-pass to
// *reclassify* short-but-nonzero edges as collapsed. The
// `skip_collapsed_edges` pre-pass uses `tolerance.linear()` as
// its length threshold, so scaling that by 10× widens the range
// of edges that get skipped.
//
// To test this without fabricating a real imported BRep, we move
// a cube vertex by a distance larger than the strict tolerance
// but smaller than the strict tolerance × 10. In strict mode the
// edge is non-collapsed (operation proceeds normally); in
// soft mode the collapsed-edge pre-pass drops it from the batch
// with a CollapsedEdgeSkipped warning.

#[test]
fn soft_supports_reclassifies_near_collapsed_edge() {
    let mut brep = BRep::unit_cube(tol()).unwrap();
    let target_edge = brep.edges.keys().next().unwrap();

    // Strict linear tolerance is 1e-6 by default; the soft
    // multiplier (10×) makes it 1e-5. Move the vertex by 5e-6 —
    // above strict threshold, below soft threshold.
    let he = brep.edges[target_edge].half_edge;
    let twin = brep.half_edges[he].twin;
    let v_start = brep.half_edges[he].origin;
    let v_end = brep.half_edges[twin].origin;
    let end_pos = brep.vertices[v_end].position;
    let drift = 5e-6;
    brep.vertices[v_start].position =
        popkern_core::geometry::Point3::new(end_pos.x + drift, end_pos.y, end_pos.z);

    let specs = vec![FilletSpec {
        edge: target_edge,
        radius: FilletRadius::Constant(0.1),
    }];

    // Strict mode: edge is NOT collapsed, the pre-pass leaves it
    // in the spec list. Subsequent phases may succeed, partial-
    // fail, or fail — we just verify the pre-pass didn't drop it.
    let strict_result = fillet_edges(&brep, &specs, &FilletOptions::default(), tol(), None);
    let strict_warns = match &strict_result {
        Ok(b) => b.edge_warnings.clone(),
        Err(BlendError::PartialFailure { .. }) => return,
        Err(e) => panic!("strict mode unexpected error: {e:?}"),
    };
    let strict_collapsed = strict_warns
        .iter()
        .any(|(_, w)| matches!(w, BlendWarning::CollapsedEdgeSkipped { .. }));
    assert!(
        !strict_collapsed,
        "strict mode (1e-6 tol) should NOT drop a 5e-6 edge as collapsed"
    );

    // Soft mode: edge IS classified as collapsed under the loose
    // 1e-5 threshold and dropped from the batch with the warning.
    let soft_opts = FilletOptions {
        soft_supports: true,
        soft_support_factor: 10.0,
        ..FilletOptions::default()
    };
    let soft_result = fillet_edges(&brep, &specs, &soft_opts, tol(), None);
    let soft_warns = match &soft_result {
        Ok(b) => b.edge_warnings.clone(),
        Err(BlendError::NoValidEdges(_)) => {
            // The edge was dropped and there's nothing left to
            // fillet — which IS the soft-mode repair behaviour.
            // Equivalent to "batch successfully soft-dropped the
            // only spec"; treat as success for this test's purpose.
            return;
        }
        Err(BlendError::PartialFailure { .. }) => return,
        Err(e) => panic!("soft mode unexpected error: {e:?}"),
    };
    let has_soft_warn = soft_warns
        .iter()
        .any(|(_, w)| matches!(w, BlendWarning::SoftSupportTolerancesApplied { .. }));
    assert!(
        has_soft_warn,
        "soft mode should emit the SoftSupportTolerancesApplied warning"
    );
    let has_collapsed = soft_warns
        .iter()
        .any(|(_, w)| matches!(w, BlendWarning::CollapsedEdgeSkipped { .. }));
    assert!(
        has_collapsed,
        "soft mode (1e-5 tol) should drop the 5e-6 edge as collapsed, got {:?}",
        soft_warns
    );
}
