//! §5.3 item #6 — emit deduped surface-bridge-fallback warnings.
//!
//! Extracted from `construct/mod.rs` to keep that file under the
//! 500-line structural cap. A corner with multiple sub-patches
//! (e.g. three-quad Gregory dispatch at a trihedral vertex) can
//! have the same surface-fallback name replicated across every
//! sub-patch's `CornerPatchInfo::bridge_fallbacks`; we emit one
//! warning per UNIQUE surface kind so UI consumers see each
//! distinct issue exactly once.

use popkern_core::topology::{EdgeHandle, VertexHandle};

use crate::analysis::EdgeAnalysis;
use crate::construct::EdgeGeometry;
use crate::error::BlendWarning;
use crate::types::CornerPatchInfo;

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
    let _ = std::marker::PhantomData::<EdgeAnalysis>; // doc link
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
