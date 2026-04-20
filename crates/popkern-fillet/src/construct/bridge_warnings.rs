//! §5.3 item #6 & #7 — per-vertex UX warnings emitted by
//! `apply_blend_geometry` after `insert_vertex_blend` returns.
//!
//! * Item #6: [`emit_bridge_fallback_warnings`] dedups
//!   [`BlendWarning::SurfaceBridgeFallback`] across sub-patches.
//! * Item #7: [`emit_tangent_corner_warning`] fires a single
//!   [`BlendWarning::TangentCornerCap`] when the tangent-plane
//!   corner method was used.

use popkern_core::topology::{EdgeHandle, VertexHandle};

use crate::analysis::EdgeAnalysis;
use crate::construct::EdgeGeometry;
use crate::error::BlendWarning;
use crate::types::{CornerMethod, CornerPatchInfo};

/// Push one deduped
/// [`BlendWarning::SurfaceBridgeFallback`] per unique surface
/// kind found across the vertex's corner patches. Attaches the
/// warning to the first incident edge (the per-edge warning
/// vector shape mandated by `ApplyBlendResult::warnings`).
pub(crate) fn emit_bridge_fallback_warnings(
    infos: &[CornerPatchInfo],
    vertex: VertexHandle,
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    warnings: &mut Vec<(EdgeHandle, BlendWarning)>,
) {
    let _ = std::marker::PhantomData::<EdgeAnalysis>;
    let mut seen: std::collections::BTreeSet<&'static str> = std::collections::BTreeSet::new();
    for info in infos {
        for kind in &info.bridge_fallbacks {
            if !seen.insert(kind) {
                continue;
            }
            let edge_idx = edge_indices[0];
            if edge_idx < edge_geometries.len() {
                warnings.push((
                    edge_geometries[edge_idx].analysis.edge,
                    BlendWarning::SurfaceBridgeFallback {
                        vertex,
                        surface_kind: kind,
                    },
                ));
            }
        }
    }
}

/// Push a single
/// [`BlendWarning::TangentCornerCap`] when any of the corner's
/// sub-patches used [`CornerMethod::TangentPlane`]. One warning
/// per vertex regardless of sub-patch count.
pub(crate) fn emit_tangent_corner_warning(
    infos: &[CornerPatchInfo],
    vertex: VertexHandle,
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    warnings: &mut Vec<(EdgeHandle, BlendWarning)>,
) {
    let used_tangent_plane = infos
        .iter()
        .any(|info| info.method == CornerMethod::TangentPlane);
    if !used_tangent_plane {
        return;
    }
    let Some(first) = infos.first() else {
        return;
    };
    let edge_idx = edge_indices[0];
    if edge_idx < edge_geometries.len() {
        warnings.push((
            edge_geometries[edge_idx].analysis.edge,
            BlendWarning::TangentCornerCap {
                vertex,
                n_sides: first.n_sides,
            },
        ));
    }
}

/// §5.3 item #10 — emit one [`BlendWarning::CornerRemediationOverlap`]
/// when any sub-patch reports `overlap_remediated = true`. One
/// warning per vertex regardless of sub-patch count (a single
/// overlap remediation collapses the corner to one polar patch
/// so this only fires once anyway, but the dedup is defensive).
pub(crate) fn emit_overlap_remediation_warning(
    infos: &[CornerPatchInfo],
    vertex: VertexHandle,
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    warnings: &mut Vec<(EdgeHandle, BlendWarning)>,
) {
    let remediated = infos.iter().any(|info| info.overlap_remediated);
    if !remediated {
        return;
    }
    let Some(first) = infos.first() else {
        return;
    };
    let edge_idx = edge_indices[0];
    if edge_idx < edge_geometries.len() {
        warnings.push((
            edge_geometries[edge_idx].analysis.edge,
            BlendWarning::CornerRemediationOverlap {
                vertex,
                n_sides: first.n_sides,
            },
        ));
    }
}

/// Phase H.8 step-D — emit one
/// [`BlendWarning::CornerContinuityRelaxed`] when any sub-patch's
/// `relaxed_continuity_coplanarity` is `Some`. One warning per
/// vertex; the worst coplanarity across sub-patches is reported.
pub(crate) fn emit_corner_continuity_relaxed_warning(
    infos: &[CornerPatchInfo],
    vertex: VertexHandle,
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    warnings: &mut Vec<(EdgeHandle, BlendWarning)>,
) {
    let worst = infos
        .iter()
        .filter_map(|info| info.relaxed_continuity_coplanarity)
        .fold(f64::NEG_INFINITY, f64::max);
    if !worst.is_finite() {
        return;
    }
    let Some(first) = infos.first() else {
        return;
    };
    let edge_idx = edge_indices[0];
    if edge_idx < edge_geometries.len() {
        warnings.push((
            edge_geometries[edge_idx].analysis.edge,
            BlendWarning::CornerContinuityRelaxed {
                vertex,
                n_sides: first.n_sides,
                worst_coplanarity: worst,
            },
        ));
    }
}
