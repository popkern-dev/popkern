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
