//! Direct-apex dispatch promotion policy.
//!
//! Decides, per apex, whether the direct-apex construction runs by
//! default (as if `use_direct_apex_construction` were set). Two
//! promoted classes:
//!
//! * **Mixed-radius apexes at N ∈ {3, 4}** (PK-099/PK-127): the
//!   legacy route ships an unmeasured fallback cap on same-side
//!   convex quad apexes (pyramid4_mixed: cap/support g0 = 1.34); the
//!   direct-apex partition closes the same apex at 1.5e-14.
//! * **Uniform-radius cylindrical N = 4 apexes** (PK-229): the ONLY
//!   class whose legacy N-quad poke commit reached serving output
//!   under default options, at sub-commercial quality (cap_outer
//!   g0 = 2.18e-1 — the rolling-ball sphere that exists at regular
//!   pyramid apexes is never found because `try_spherical_blend`'s
//!   circumcenter solve is singular on co-circular arc samples, and
//!   the n < 5 poke acceptance ships the Gregory fallback
//!   unmeasured). The direct-apex exact-boundary route serves the
//!   same apex at g0 = 3.7e-8 under pure defaults (measured on the
//!   PK-228 pin fixture). Restricted to all-`Cylinder` strips so
//!   conic-profile (NURBS-strip) corners keep their current serve.
//!
//! Both classes rely on the direct-apex entry-snapshot rollback: any
//! construction failure restores the corner and the legacy dispatch
//! serves exactly as before (`DirectApexFallbackToLegacy` recorded).

use popkern_core::topology::VertexHandle;

use super::EdgeGeometry;

/// Surface-family flags that alter which apex classes promote.
#[derive(Clone, Copy)]
pub(super) struct PromotionFlags {
    pub use_rational_gregory: bool,
    pub use_kp_binonic_gregory: bool,
    /// True when any explicit legacy surface-family flag is set
    /// (biquintic-full, biseptic, two-stripe, K-P binonic, K-P
    /// polar — all default-off). Explicit family requests keep the
    /// legacy promotion ladder measurable on its reaching fixture
    /// (uniform pyramid apexes), so the uniform-arm promotion stands
    /// down — the same precedent as the binonic high-asymmetry
    /// carve-out on the mixed arm.
    pub explicit_legacy_surface_family: bool,
}

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
pub(super) fn promoted_apexes_for_mode(
    vertex_mode: &super::VertexMode,
    edge_geometries: &[EdgeGeometry],
    flags: PromotionFlags,
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
                && should_promote_direct_apex(edge_indices, edge_geometries, flags)
        })
        .map(|(vertex, _)| vertex)
        .collect()
}

pub(super) fn should_promote_direct_apex(
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    flags: PromotionFlags,
) -> bool {
    // N=3 (trihedral) and N=4 apexes only. N=4 joined when PK-099's
    // orientation-honest convexity exposed that same-side-convex quad
    // apexes never reach the rational three-quad dispatch on the
    // legacy route (the overlap-split corner loop is not
    // 2n-dispatch-ready), so legacy shipped an unmeasured fallback
    // cap. The direct-apex partition construction closes the same
    // apex at 1.5e-14, and its snapshot rollback restores legacy
    // behavior on any failure.
    if !flags.use_rational_gregory || !matches!(edge_indices.len(), 3 | 4) {
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
    if radii.len() != edge_indices.len() {
        return false;
    }
    let Some((&first, rest)) = radii.split_first() else {
        return false;
    };
    let uniform_effective = rest.iter().all(|radius| (*radius - first).abs() <= 1.0e-9);
    if uniform_effective {
        // PK-229 uniform arm: constant-radius cylindrical N = 4
        // apexes (pyramid apex slant blends). N = 3 uniform corners
        // stay legacy: cube-class corners serve the exact
        // rolling-ball sphere there, and the oblique
        // (non-orthogonal) trihedral apex class needs its own
        // geometric discriminator before it can promote — see the
        // PK-229 ledger follow-ups. Variable-radius edges (start ≠
        // end) and non-cylindrical (conic-profile NURBS) strips are
        // not the measured class and keep their current serve.
        return edge_indices.len() == 4
            && !flags.explicit_legacy_surface_family
            && edge_indices.iter().all(|&index| {
                edge_geometries.get(index).is_some_and(|g| {
                    (g.radius_at_start - g.radius_at_end).abs() <= 1.0e-9
                        && matches!(
                            g.surface,
                            popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. }
                        )
                })
            });
    }
    if flags.use_kp_binonic_gregory {
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
