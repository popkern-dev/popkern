//! Mixed-radius direct-apex dispatch policy.

use popkern_core::topology::VertexHandle;

use super::EdgeGeometry;

/// Per-apex promotion set, computed BEFORE Phase A so anchor policy
/// and dispatch share one decision (PK-127).
///
/// The promotion used to be decided only at dispatch time, after
/// Phase A had already placed LEGACY (trim-crossing, no planner
/// override) anchors around the apex — the repeated-apex partition
/// plan then met anchor positions it cannot consume ("cannot
/// duplicate stale residual strip seam" / "no existing BRep vertex
/// within 1e-6") and every default-options promoted apex fell back
/// to the unmeasured legacy cap. Phase A consults this set so edges
/// incident to a promoted apex get the direct-apex anchor
/// convention (perp-walk positions + planner-prepass overrides),
/// exactly as if `use_direct_apex_construction` were set.
pub(super) fn promoted_mixed_radius_apexes_for_mode(
    vertex_mode: &super::VertexMode,
    edge_geometries: &[EdgeGeometry],
    use_rational_gregory: bool,
    use_kp_binonic_gregory: bool,
) -> std::collections::BTreeSet<VertexHandle> {
    // Vertex blends dispatch only in Smooth mode with 2+ fillets —
    // the same gate the dispatch loop applies.
    if *vertex_mode != super::VertexMode::Smooth || edge_geometries.len() < 2 {
        return std::collections::BTreeSet::new();
    }
    let mut vertex_edges: std::collections::BTreeMap<VertexHandle, Vec<usize>> =
        std::collections::BTreeMap::new();
    for (index, geometry) in edge_geometries.iter().enumerate() {
        vertex_edges
            .entry(geometry.analysis.vertex_start)
            .or_default()
            .push(index);
        vertex_edges
            .entry(geometry.analysis.vertex_end)
            .or_default()
            .push(index);
    }
    vertex_edges
        .into_iter()
        .filter(|(_, edge_indices)| {
            edge_indices.len() >= 2
                && should_promote_mixed_radius_gregory_apex(
                    edge_indices,
                    edge_geometries,
                    use_rational_gregory,
                    use_kp_binonic_gregory,
                )
        })
        .map(|(vertex, _)| vertex)
        .collect()
}

pub(super) fn should_promote_mixed_radius_gregory_apex(
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    use_rational_gregory: bool,
    use_kp_binonic_gregory: bool,
) -> bool {
    // N=3 (trihedral) and N=4 mixed-radius apexes promote. N=4 joined
    // when PK-099's orientation-honest convexity exposed that
    // same-side-convex quad apexes never reach the rational three-quad
    // dispatch on the legacy route (the overlap-split corner loop is
    // not 2n-dispatch-ready), so legacy shipped an unmeasured fallback
    // cap (pyramid4_mixed: cap/support g0 = 1.34). The direct-apex
    // partition construction closes the same apex at 1.5e-14, and its
    // snapshot rollback restores legacy behavior on any failure.
    if !use_rational_gregory || !matches!(edge_indices.len(), 3 | 4) {
        return false;
    }
    let radii = edge_indices
        .iter()
        .filter_map(|&index| {
            edge_geometries
                .get(index)
                .map(|geometry| geometry.effective_radius)
        })
        .collect::<Vec<_>>();
    let Some((&first, rest)) = radii.split_first() else {
        return false;
    };
    if rest.iter().all(|radius| (*radius - first).abs() <= 1.0e-9) {
        return false;
    }
    if use_kp_binonic_gregory {
        let min_radius = radii.iter().copied().fold(f64::INFINITY, f64::min);
        let max_radius = radii.iter().copied().fold(0.0_f64, f64::max);
        if min_radius <= 0.0 {
            return false;
        }
        // Exact-boundary direct-apex remains the preferred mixed-radius
        // construction when the radius band is modest. Explicit K-P/binonic
        // requests own high-asymmetry trihedra so the legacy K-P ladder stays
        // measurable instead of being masked by automatic direct-apex routing.
        return max_radius / min_radius <= 2.0;
    }
    true
}
