//! Direct-apex dispatch promotion policy.
//!
//! Decides, per apex, whether the direct-apex construction runs by
//! default (as if `use_direct_apex_construction` were set). Three
//! promoted classes:
//!
//! * **Mixed-radius apexes at N ∈ {3, 4}** (PK-099/PK-127): the
//!   vertex-blend route ships an unmeasured fallback cap on same-side
//!   convex quad apexes (pyramid4_mixed: cap/support g0 = 1.34); the
//!   direct-apex partition closes the same apex at 1.5e-14.
//! * **Uniform-radius cylindrical N = 4 apexes** (PK-229): the ONLY
//!   class whose N-quad poke commit reached serving output
//!   under default options, at sub-commercial quality (cap_outer
//!   g0 = 2.18e-1 — the rolling-ball sphere that exists at regular
//!   pyramid apexes is never found because `try_spherical_blend`'s
//!   circumcenter solve is singular on co-circular arc samples, and
//!   the n < 5 poke acceptance ships the Gregory fallback
//!   unmeasured). The direct-apex exact-boundary route serves the
//!   same apex at g0 = 3.7e-8 under pure defaults (measured on the
//!   PK-228 pin fixture). Restricted to all-`Cylinder` strips so
//!   conic-profile (NURBS-strip) corners keep their current serve.
//! * **Uniform-radius cylindrical OBLIQUE N = 3 apexes** (PK-232):
//!   before the kef `outgoing` repair the class shipped NO smooth
//!   cap at all (anchor collection walked from stale handles and
//!   the blend fell to Mitred); post-repair the N-quad poke serves
//!   g0 = 7.2e-2 with recorded declines, while direct-apex serves
//!   the same apex at g0 = 4.9e-8 under pure defaults. Gated on the
//!   axis-obliqueness discriminator: mutually-orthogonal
//!   (cube-class) corners stay vertex-blend for their exact rolling-ball
//!   sphere serve.
//!
//! All classes rely on the direct-apex entry-snapshot rollback: any
//! construction failure restores the corner and the vertex-blend dispatch
//! serves exactly as before (`DirectApexFallbackToVertexBlend` recorded).
//!
//! Every class is additionally gated on the vertex classifier:
//! only `SameSideConvex` apexes promote (PK-231 — the machinery's
//! measured domain is all-convex; an OppositeSide saddle promoted
//! through the mixed arm planned the Open rung and its vacuous
//! no-op "success" deleted the pre-ring saddle cap entirely).

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, VertexHandle};

use super::EdgeGeometry;

/// Surface-family flags that alter which apex classes promote.
#[derive(Clone, Copy)]
pub(super) struct PromotionFlags {
    pub use_rational_gregory: bool,
    pub use_kp_binonic_gregory: bool,
    /// True when any explicit sub-patch family flag is set
    /// (biquintic-full, biseptic, two-stripe, K-P binonic, K-P
    /// polar — all default-off). Explicit family requests keep the
    /// flag-gated promotion ladder measurable on its reaching fixture
    /// (uniform pyramid apexes), so the uniform-arm promotion stands
    /// down — the same precedent as the binonic high-asymmetry
    /// carve-out on the mixed arm.
    pub explicit_subpatch_family_requested: bool,
}

/// Per-apex promotion set, computed BEFORE anchor splitting so anchor policy
/// and dispatch share one decision (PK-127).
///
/// The promotion used to be decided only at dispatch time, after
/// the anchor-split stage had already placed trim-crossing (no planner
/// override) anchors around the apex — the repeated-apex partition
/// plan then met anchor positions it cannot consume ("cannot
/// duplicate stale residual strip seam" / "no existing BRep vertex
/// within 1e-6") and every default-options promoted apex fell back
/// to the unmeasured vertex-blend cap. the anchor-split stage consults this set so edges
/// incident to a promoted apex get the direct-apex anchor
/// convention (perp-walk positions + planner-prepass overrides),
/// exactly as if `use_direct_apex_construction` were set.
pub(super) fn promoted_apexes_for_mode(
    brep: &BRep,
    vertex_mode: &super::VertexMode,
    edge_geometries: &[EdgeGeometry],
    flags: PromotionFlags,
    tolerance: &Tolerance,
) -> std::collections::BTreeSet<VertexHandle> {
    // Vertex blends dispatch only in Smooth mode with 2+ fillets —
    // the same gate the dispatch loop applies.
    if *vertex_mode != super::VertexMode::Smooth || edge_geometries.len() < 2 {
        return std::collections::BTreeSet::new();
    }
    // Closed edges count ONCE at their closure vertex (PK-296) —
    // mirrors the dispatch loop's map so promotion can never target
    // a closure vertex the dispatch loop will not visit.
    let vertex_edges = crate::analysis::vertex_incidences(
        edge_geometries
            .iter()
            .map(|eg| (eg.analysis.vertex_start, eg.analysis.vertex_end)),
    );
    vertex_edges
        .into_iter()
        .filter(|(vertex, edge_indices)| {
            edge_indices.len() >= 2
                && should_promote_direct_apex(edge_indices, edge_geometries, flags)
                && apex_is_same_side_convex(brep, *vertex, edge_indices, edge_geometries, tolerance)
        })
        .map(|(vertex, _)| vertex)
        .collect()
}

/// Every promoted class was measured on all-convex apexes (pyramid
/// slant blends, cube corners) — the direct-apex planner models the
/// rolling ball on ONE side of every support, so a concave strip is
/// invisible to its footprint search (PK-231: the OppositeSide
/// notched-cube saddle planned "only 2 fillets" out of 3, landed the
/// Open rung, and the vacuous no-op success swallowed the vertex-blend
/// saddle machinery entirely). Only `SameSideConvex` promotes;
/// saddle/concave/tangent/mixed apexes keep their vertex-blend dispatch,
/// which owns the class-specific construction strategies.
fn apex_is_same_side_convex(
    brep: &BRep,
    vertex: VertexHandle,
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> bool {
    let edges: Vec<_> = edge_indices
        .iter()
        .filter_map(|&index| edge_geometries.get(index).map(|g| g.analysis.edge))
        .collect();
    let report = crate::vertex_classifier::classify_vertex_blend(brep, vertex, &edges, tolerance);
    matches!(
        report.classification,
        crate::vertex_classifier::VertexBlendClassification::SameSideConvex
    )
}

/// True when every pair of the corner's cylinder-strip axes is
/// orthogonal — the cube-corner (prismatic trihedral) signature. The
/// oblique pyramid apex measures pairwise |cos| = 8.5e-2 on the
/// flattest ledger fixture (n_gonal_pyramid(3, 1.0, 0.8)), four
/// orders above the 1e-4 band, so the discriminator has no
/// borderline class among known fixtures. Non-cylinder strips never
/// reach this test (the uniform arm requires all-Cylinder first).
pub(super) fn cylinder_axes_mutually_orthogonal(
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
) -> bool {
    let axes: Vec<popkern_core::geometry::Vec3> = edge_indices
        .iter()
        .filter_map(|&index| {
            let g = edge_geometries.get(index)?;
            match g.surface {
                popkern_core::geometry::surface::SurfaceGeometry::Cylinder { axis, .. } => {
                    axis.normalized()
                }
                _ => None,
            }
        })
        .collect();
    if axes.len() != edge_indices.len() {
        // Missing/degenerate axis: treat as the protected class so
        // the corner keeps its current serve.
        return true;
    }
    axes_mutually_orthogonal(&axes)
}

/// Pairwise |cos| ≤ 1e-4 for every axis pair. The band is angular
/// and scale-free (axes are unit), sized to absorb constructor
/// round-off on genuinely orthogonal corners while sitting four
/// orders below the flattest known oblique fixture (8.5e-2).
fn axes_mutually_orthogonal(axes: &[popkern_core::geometry::Vec3]) -> bool {
    axes.iter()
        .enumerate()
        .all(|(i, a)| axes[i + 1..].iter().all(|b| a.dot(b).abs() <= 1.0e-4))
}

pub(super) fn should_promote_direct_apex(
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    flags: PromotionFlags,
) -> bool {
    // N=3 (trihedral) and N=4 apexes only. N=4 joined when PK-099's
    // orientation-honest convexity exposed that same-side-convex quad
    // apexes never reach the rational three-quad dispatch on the
    // vertex-blend route (the overlap-split corner loop is not
    // 2n-dispatch-ready), so vertex-blend shipped an unmeasured fallback
    // cap. The direct-apex partition construction closes the same
    // apex at 1.5e-14, and its snapshot rollback restores the vertex-blend serve
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
        // Uniform arm: constant-radius cylindrical apexes.
        // Variable-radius edges (start ≠ end) and non-cylindrical
        // (conic-profile NURBS) strips are not the measured class
        // and keep their current serve; explicit sub-patch family
        // requests keep the flag-gated ladder measurable on its
        // reaching fixture.
        let cylindrical_constant = edge_indices.iter().all(|&index| {
            edge_geometries.get(index).is_some_and(|g| {
                (g.radius_at_start - g.radius_at_end).abs() <= 1.0e-9
                    && matches!(
                        g.surface,
                        popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. }
                    )
            })
        });
        if flags.explicit_subpatch_family_requested || !cylindrical_constant {
            return false;
        }
        // N = 4 (PK-229): pyramid apex slant blends, measured
        // g0 2.18e-1 N-quad poke -> 3.7e-8 direct-apex.
        if edge_indices.len() == 4 {
            return true;
        }
        // N = 3 (PK-232): only the OBLIQUE trihedral class promotes
        // (measured g0 7.2e-2 N-quad poke -> 4.9e-8 direct-apex on
        // n_gonal_pyramid(3, 1.0, 0.8) r=0.08 under pure defaults).
        // Mutually-orthogonal (cube-class) corners stay vertex-blend:
        // their carve produces the arc+bridge boundary ring that
        // serves the EXACT rolling-ball sphere via surface-fix
        // promotion, which no fitted cap can beat.
        return !cylinder_axes_mutually_orthogonal(edge_indices, edge_geometries);
    }
    if flags.use_kp_binonic_gregory {
        let min_radius = radii.iter().copied().fold(f64::INFINITY, f64::min);
        let max_radius = radii.iter().copied().fold(0.0_f64, f64::max);
        if min_radius <= 0.0 {
            return false;
        }
        // Exact-boundary direct-apex remains the preferred mixed-radius
        // construction when the radius band is modest. Explicit K-P/binonic
        // requests own high-asymmetry trihedra so the flag-gated K-P ladder stays
        // measurable instead of being masked by automatic direct-apex routing.
        return max_radius / min_radius <= 2.0;
    }
    true
}

#[cfg(test)]
mod tests {
    use popkern_core::geometry::Vec3;

    use super::axes_mutually_orthogonal;

    #[test]
    fn cube_corner_axes_are_orthogonal() {
        let axes = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        assert!(axes_mutually_orthogonal(&axes));
    }

    #[test]
    fn oblique_pyramid_axes_are_not_orthogonal() {
        // Slant-edge directions of n_gonal_pyramid(3, 1.0, 0.8):
        // pairwise cos = (cos 120° + h²) / (1 + h²) ≈ 8.54e-2 — the
        // flattest ledger fixture, four orders above the band.
        let h = 0.8_f64;
        let axes: Vec<Vec3> = (0..3)
            .map(|i| {
                let theta = 2.0 * std::f64::consts::PI * f64::from(i) / 3.0;
                Vec3::new(theta.cos(), theta.sin(), -h)
                    .normalized()
                    .expect("slant direction is non-degenerate")
            })
            .collect();
        assert!(!axes_mutually_orthogonal(&axes));
        let worst = (0..3)
            .flat_map(|i| ((i + 1)..3).map(move |j| (i, j)))
            .map(|(i, j)| axes[i].dot(&axes[j]).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            (worst - 8.54e-2).abs() < 1e-3,
            "expected the documented 8.5e-2 signature, got {worst:.3e}"
        );
    }

    #[test]
    fn round_off_orthogonal_axes_stay_in_band() {
        // A numerically-constructed cube corner: 1e-7-level cross
        // terms must still classify as the protected orthogonal
        // class.
        let axes = [
            Vec3::new(1.0, 1.0e-7, -3.0e-8),
            Vec3::new(-2.0e-8, 1.0, 5.0e-8),
            Vec3::new(4.0e-8, -1.0e-7, 1.0),
        ];
        assert!(axes_mutually_orthogonal(&axes));
    }
}
