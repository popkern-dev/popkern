//! Synthetic bridge targets for the gap closure: exact boundary-edge
//! × surface junctions that never became graph nodes (split from
//! `graph_gap_closure.rs` under the 500-line structural budget).

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::BRep;

use super::graph::{IntersectionGraph, Operand};
/// A bridge target: an existing graph node, or an exact boundary
/// junction not yet represented by one (minted on certification).
#[derive(Clone, Copy, Debug)]
pub(super) enum Target {
    Node(usize),
    Junction(Point3),
}

/// Collect exact boundary-edge × other-surface junctions near the tip
/// as bridge targets. A march that stops short of its own boundary
/// crossing (the grazing exhaust) leaves its TRUE terminus with
/// nothing to be split at — the junction must be minted and the
/// boundary edge split there afterward (the caller re-runs the
/// node-driven edge split pass).
#[allow(clippy::too_many_arguments)]
pub(super) fn push_synthetic_junction_targets(
    cands: &mut Vec<(f64, Target)>,
    graph: &IntersectionGraph,
    original_a: &BRep,
    original_b: &BRep,
    face_a: popkern_core::topology::FaceHandle,
    face_b: popkern_core::topology::FaceHandle,
    p: &Point3,
    outward: &Vec3,
    cap: f64,
    lin: f64,
) {
    let pair_faces = [
        (Operand::A, original_a, face_a),
        (Operand::B, original_b, face_b),
    ];
    for (op, brep, face) in pair_faces {
        let other_surf = match op {
            Operand::A => original_b
                .faces
                .get(face_b)
                .and_then(|f| f.surface.as_ref())
                .map(|sf| &sf.geometry),
            Operand::B => original_a
                .faces
                .get(face_a)
                .and_then(|f| f.surface.as_ref())
                .map(|sf| &sf.geometry),
        };
        let Some(other_surf) = other_surf else {
            continue;
        };
        let (other_brep, other_face) = match op {
            Operand::A => (original_b, face_b),
            Operand::B => (original_a, face_a),
        };
        let ow = crate::intersect::bounded_parameter_range(other_brep, other_face, other_surf);
        for eh in crate::clip_endpoint_refine::boundary_edges(brep, face) {
            let Some(edge) = brep.edges.get(eh) else {
                continue;
            };
            let Some(c) = edge.curve.as_ref() else {
                continue;
            };
            let range = (c.t_start.min(c.t_end), c.t_start.max(c.t_end));
            let geometry = crate::clip_endpoint_refine::clone_curve_geometry(&c.geometry);
            let Ok(result) = popkern_core::geometry::intersection::intersect_surface_curve(
                other_surf,
                ow.0,
                ow.1,
                &geometry,
                range,
                &original_a.tolerance,
            ) else {
                continue;
            };
            for isect in &result.intersections {
                let popkern_core::geometry::intersection::SurfaceCurveIntersection::Point {
                    point,
                    ..
                } = isect
                else {
                    continue;
                };
                let d = point.distance_to(p);
                if d < lin * 10.0 || d > cap {
                    continue;
                }
                // A junction already represented by a node
                // (within snap) is the node's business.
                if graph.find_node(point, lin * 20.0).is_some() {
                    continue;
                }
                let Some(dir) = (*point - *p).normalized() else {
                    continue;
                };
                if outward.dot(&dir) < -0.5 {
                    continue;
                }
                cands.push((d, Target::Junction(*point)));
            }
        }
    }
    cands.sort_by(|a, b| a.0.total_cmp(&b.0));
    cands.truncate(8);
}

/// Dead-end tips, judged PER FACE — the same criterion as the walk's
/// graze filter: a node whose incidence ON ONE FACE'S chart is a
/// single shared record pair (no boundary continuation there)
/// dead-ends that face even when the node carries records of other
/// pairs (measured: the nut's runout-chain ends are two-record
/// junction nodes globally, yet each face sees exactly one). Maps
/// tip node -> its single record half-edge.
pub(super) fn per_face_dead_end_tips(
    graph: &IntersectionGraph,
) -> std::collections::BTreeMap<usize, usize> {
    let is_shared = |h: usize, graph: &IntersectionGraph| {
        matches!(
            graph.half_edges[h].source,
            super::graph::EdgeSource::Shared { .. }
        )
    };

    let mut per_face_tips: std::collections::BTreeMap<usize, usize> =
        std::collections::BTreeMap::new();
    let faces: Vec<_> = graph.face_edges.keys().cloned().collect();
    for f in faces {
        let hes = &graph.face_edges[&f];
        let mut node_inc: std::collections::BTreeMap<usize, Vec<usize>> =
            std::collections::BTreeMap::new();
        for &h in hes {
            let a = graph.half_edges[h].origin;
            let b = graph.half_edges[graph.half_edges[h].twin].origin;
            node_inc.entry(a).or_default().push(h);
            node_inc.entry(b).or_default().push(h);
        }
        for (n, inc) in node_inc {
            let records: std::collections::BTreeSet<usize> = inc
                .iter()
                .filter(|&&h| is_shared(h, graph))
                .map(|&h| h.min(graph.half_edges[h].twin))
                .collect();
            let has_boundary = inc.iter().any(|&h| !is_shared(h, graph));
            if records.len() == 1 && !has_boundary {
                let rh = *inc
                    .iter()
                    .find(|&&h| is_shared(h, graph))
                    .unwrap_or(&inc[0]);
                per_face_tips.entry(n).or_insert(rh);
            }
        }
    }
    per_face_tips
}
