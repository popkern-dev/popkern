use std::collections::{BTreeSet, HashMap};

use crate::topology::{
    BRep, BrepExportQualityScope, BrepSewingEdgePair, EdgeHandle, FaceHandle, HalfEdgeHandle,
    LoopHandle, VertexHandle,
};

/// Pre-repair topology scope affected by PnS boundary welding candidates.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PolyhedralBoundaryWeldScope {
    pub candidate_pairs: usize,
    pub components: Vec<usize>,
    pub edges: Vec<EdgeHandle>,
    pub vertices: Vec<VertexHandle>,
    pub half_edges: Vec<HalfEdgeHandle>,
    pub loops: Vec<LoopHandle>,
    pub faces: Vec<FaceHandle>,
}

impl PolyhedralBoundaryWeldScope {
    pub fn export_quality_scope(&self) -> BrepExportQualityScope {
        BrepExportQualityScope {
            edges: self.edges.clone(),
            vertices: self.vertices.clone(),
            half_edges: self.half_edges.clone(),
            loops: self.loops.clone(),
            faces: self.faces.clone(),
        }
    }
}

pub(super) fn build_boundary_weld_scope(
    brep: &BRep,
    candidates: &[BrepSewingEdgePair],
    edge_component: &HashMap<EdgeHandle, usize>,
) -> PolyhedralBoundaryWeldScope {
    let mut scope = PolyhedralBoundaryWeldScope {
        candidate_pairs: candidates.len(),
        ..Default::default()
    };
    let mut components = BTreeSet::new();
    for candidate in candidates {
        for edge in [candidate.first, candidate.second] {
            if let Some(component) = edge_component.get(&edge) {
                components.insert(*component);
            }
            push_edge_scope(brep, edge, &mut scope);
        }
    }
    scope.components = components.into_iter().collect();
    sort_scope_handles(&mut scope);
    scope
}

fn push_edge_scope(brep: &BRep, edge: EdgeHandle, scope: &mut PolyhedralBoundaryWeldScope) {
    push_unique(&mut scope.edges, edge);
    let Some(edge_data) = brep.edges.get(edge) else {
        return;
    };
    push_half_edge_scope(brep, edge_data.half_edge, scope);
    if let Some(half_edge) = brep.half_edges.get(edge_data.half_edge) {
        push_half_edge_scope(brep, half_edge.twin, scope);
    }
}

fn push_half_edge_scope(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    scope: &mut PolyhedralBoundaryWeldScope,
) {
    let Some(half_edge_data) = brep.half_edges.get(half_edge) else {
        return;
    };
    push_unique(&mut scope.half_edges, half_edge);
    push_unique(&mut scope.vertices, half_edge_data.origin);
    if let Some(next) = brep.half_edges.get(half_edge_data.next) {
        push_unique(&mut scope.vertices, next.origin);
    }
    push_unique(&mut scope.loops, half_edge_data.loop_);
    if let Some(loop_data) = brep.loops.get(half_edge_data.loop_) {
        push_unique(&mut scope.faces, loop_data.face);
    }
}

fn push_unique<T: Copy + Eq>(items: &mut Vec<T>, item: T) {
    if !items.contains(&item) {
        items.push(item);
    }
}

fn sort_scope_handles(scope: &mut PolyhedralBoundaryWeldScope) {
    scope.edges.sort_by_key(|handle| format!("{handle:?}"));
    scope.vertices.sort_by_key(|handle| format!("{handle:?}"));
    scope.half_edges.sort_by_key(|handle| format!("{handle:?}"));
    scope.loops.sort_by_key(|handle| format!("{handle:?}"));
    scope.faces.sort_by_key(|handle| format!("{handle:?}"));
}
