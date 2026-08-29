use super::*;

fn unit_plane_patch() -> Kp20Bi6Patch {
    let mut patch = Kp20Bi6Patch::zero();
    for row in 0..7 {
        for column in 0..7 {
            patch.cps[row][column] = Point3::new(row as f64 / 6.0, column as f64 / 6.0, 0.0);
        }
    }
    patch
}

fn projection_at(side: SupportCellPatchSide, ts: &[f64]) -> Vec<Point3> {
    let patch = unit_plane_patch();
    ts.iter()
        .map(|&t| super::boundary_point(&patch, side, t))
        .collect()
}

#[test]
fn none_caller_returns_uniform_no_projection_status() {
    let patch = unit_plane_patch();
    let outcome = cap_vertex_splits_along_axis(&patch, None, 4);
    assert_eq!(
        outcome.diagnostic.status,
        SupportCellSplitAxisStatus::UniformNoProjection
    );
    assert_eq!(outcome.params, vec![1.0 / 4.0, 2.0 / 4.0, 3.0 / 4.0]);
}

#[test]
fn uniform_vertex_set_routes_through_equivalence_shortcut() {
    let patch = unit_plane_patch();
    let vertices = projection_at(
        SupportCellPatchSide::VMin,
        &[1.0 / 4.0, 2.0 / 4.0, 3.0 / 4.0],
    );
    let projection = CapVertexAxisProjection {
        side: SupportCellPatchSide::VMin,
        vertices: &vertices,
    };
    let outcome = cap_vertex_splits_along_axis(&patch, Some(&projection), 4);
    assert_eq!(
        outcome.diagnostic.status,
        SupportCellSplitAxisStatus::UniformEquivalent
    );
    assert_eq!(outcome.params, vec![1.0 / 4.0, 2.0 / 4.0, 3.0 / 4.0]);
    assert_eq!(outcome.diagnostic.projected_count, 3);
    assert_eq!(outcome.diagnostic.uniform_fill_count, 0);
}

#[test]
fn non_uniform_vertex_set_produces_vertex_aware_sequence() {
    let patch = unit_plane_patch();
    let vertices = projection_at(SupportCellPatchSide::VMin, &[0.1, 0.4, 0.85]);
    let projection = CapVertexAxisProjection {
        side: SupportCellPatchSide::VMin,
        vertices: &vertices,
    };
    let outcome = cap_vertex_splits_along_axis(&patch, Some(&projection), 4);
    assert_eq!(
        outcome.diagnostic.status,
        SupportCellSplitAxisStatus::VertexAware
    );
    assert_eq!(outcome.params.len(), 3);
    assert!((outcome.params[0] - 0.1).abs() < 1.0e-9);
    assert!((outcome.params[1] - 0.4).abs() < 1.0e-9);
    assert!((outcome.params[2] - 0.85).abs() < 1.0e-9);
}

#[test]
fn undercount_fills_largest_gaps_uniformly() {
    let patch = unit_plane_patch();
    // Two cap-vertices but cells = 4 requires 3 interior splits.
    let vertices = projection_at(SupportCellPatchSide::VMin, &[0.2, 0.9]);
    let projection = CapVertexAxisProjection {
        side: SupportCellPatchSide::VMin,
        vertices: &vertices,
    };
    let outcome = cap_vertex_splits_along_axis(&patch, Some(&projection), 4);
    assert_eq!(
        outcome.diagnostic.status,
        SupportCellSplitAxisStatus::VertexAware
    );
    assert_eq!(outcome.diagnostic.projected_count, 2);
    assert_eq!(outcome.diagnostic.uniform_fill_count, 1);
    assert_eq!(outcome.params.len(), 3);
    // Cap-vertices (0.2, 0.9); largest pre-fill gap is (0.2, 0.9);
    // midpoint = 0.55 is inserted. Newton-projected parameters
    // carry tiny FP drift, so compare within a generous tolerance.
    let approx = |target: f64| outcome.params.iter().any(|t| (*t - target).abs() < 1.0e-9);
    assert!(approx(0.2), "missing 0.2 in {:?}", outcome.params);
    assert!(approx(0.55), "missing 0.55 in {:?}", outcome.params);
    assert!(approx(0.9), "missing 0.9 in {:?}", outcome.params);
}

#[test]
fn overcount_falls_back_to_uniform_with_status() {
    let patch = unit_plane_patch();
    // Five cap-vertices on V-axis, but cells = 2 only takes 1 interior split.
    let vertices = projection_at(SupportCellPatchSide::VMin, &[0.1, 0.3, 0.5, 0.7, 0.9]);
    let projection = CapVertexAxisProjection {
        side: SupportCellPatchSide::VMin,
        vertices: &vertices,
    };
    let outcome = cap_vertex_splits_along_axis(&patch, Some(&projection), 2);
    assert_eq!(
        outcome.diagnostic.status,
        SupportCellSplitAxisStatus::UniformFallback
    );
    assert_eq!(
        outcome.diagnostic.fallback_reason,
        Some(SupportCellSplitFallbackReason::Overcount)
    );
    assert_eq!(outcome.params, vec![0.5]);
}

#[test]
fn min_arc_guard_rejects_thin_subcells() {
    let patch = unit_plane_patch();
    // Two vertices clustered to leave a sub-cell well under the
    // 5%-of-bbox-diagonal minimum (bbox_diag ≈ √2 ≈ 1.414;
    // 5% ≈ 0.0707).
    let vertices = projection_at(SupportCellPatchSide::VMin, &[0.4, 0.41]);
    let projection = CapVertexAxisProjection {
        side: SupportCellPatchSide::VMin,
        vertices: &vertices,
    };
    let outcome = cap_vertex_splits_along_axis(&patch, Some(&projection), 3);
    assert_eq!(
        outcome.diagnostic.status,
        SupportCellSplitAxisStatus::UniformFallback
    );
    assert_eq!(
        outcome.diagnostic.fallback_reason,
        Some(SupportCellSplitFallbackReason::MinArc)
    );
}
