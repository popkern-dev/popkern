//! The result edges and vertices a complex split plan names (PK-987).
//!
//! The assembler builds one result edge per graph edge and one vertex
//! per graph node, and pairs the fragments' half-edges first come,
//! first served: a third fragment on an edge mints a duplicate with a
//! free twin. Under `NonManifoldPolicy::SplitToManifold` the complex
//! pairs them instead (`complex::split_plan`): each selected wedge of
//! a non-manifold edge is its own result edge, each fan of a split
//! vertex its own vertex, and every fragment half-edge the plan names
//! wires to its planned edge. Everything else keeps the default.

use std::collections::{BTreeMap, HashMap};

use popkern_core::constructors::builder::BRepBuilder;
use popkern_core::topology::{HalfEdgeHandle, VertexHandle};

use crate::complex::split_plan::SplitPlan;
use crate::merge::graph::IntersectionGraph;

/// The planned edges, built.
pub(super) struct SplitBuild<'a> {
    plan: Option<&'a SplitPlan>,
    /// Half-edge pair of each planned edge: (start → end, end → start).
    edges: Vec<(HalfEdgeHandle, HalfEdgeHandle)>,
}

impl<'a> SplitBuild<'a> {
    /// Create the plan's vertex copies and edges (copy 0 of a node is
    /// its default vertex).
    pub(super) fn new(
        plan: Option<&'a SplitPlan>,
        builder: &mut BRepBuilder,
        graph: &IntersectionGraph,
        node_to_vertex: &HashMap<usize, VertexHandle>,
        carry: &mut super::assemble_far_carry::Carry<'_>,
        carried: &mut crate::serve_region::Carried,
    ) -> Self {
        let mut edges = Vec::new();
        if let Some(plan) = plan {
            let mut copies: BTreeMap<(usize, usize), VertexHandle> = BTreeMap::new();
            for (&node, &n) in &plan.vertex_copies {
                for k in 1..n {
                    copies.insert((node, k), builder.add_vertex(graph.nodes[node].position));
                }
            }
            let vertex = |(node, k): (usize, usize)| match k {
                0 => node_to_vertex.get(&node).copied(),
                _ => copies.get(&(node, k)).copied(),
            };
            for e in &plan.edges {
                let (Some(vs), Some(ve)) = (vertex(e.start), vertex(e.end)) else {
                    // The node is no kept fragment's: no wiring names it.
                    edges.push((HalfEdgeHandle::default(), HalfEdgeHandle::default()));
                    continue;
                };
                let (edge_h, fwd, twin) = builder.add_edge(vs, ve);
                carry.edge(&graph.half_edges[e.curve_half_edge].source, edge_h, carried);
                super::assemble_edge_geometry::attach_edge_curve(
                    builder,
                    graph,
                    e.curve_half_edge,
                    edge_h,
                    e.start.0,
                    e.end.0,
                );
                edges.push((fwd, twin));
            }
        }
        Self { plan, edges }
    }

    /// The planned half-edge wiring `fragment`'s loop half-edge
    /// `listed` (as the fragment lists it), traversed against its
    /// listed direction when `reversed`; `None` when the plan leaves it
    /// to the default pairing.
    pub(super) fn half_edge(
        &self,
        fragment: usize,
        listed: usize,
        reversed: bool,
    ) -> Option<HalfEdgeHandle> {
        let &(edge, forward) = self.plan?.uses.get(&(fragment, listed))?;
        let (fwd, twin) = *self.edges.get(edge)?;
        Some(if forward != reversed { fwd } else { twin })
    }
}
