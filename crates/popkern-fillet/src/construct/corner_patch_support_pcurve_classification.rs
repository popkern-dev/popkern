//! External cap/support pcurve classification diagnostics.

use std::collections::BTreeSet;

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

#[derive(Debug, Clone, Copy, Default)]
pub(super) struct ExternalCapSupportPcurveClassification {
    pub(super) unrecognized_cap_boundary_pcurve_edges: usize,
    pub(super) unrecognized_cap_boundary_hint_edges: usize,
    pub(super) recognized_cap_boundary_hint_edges: usize,
}

pub(super) fn classify_external_cap_support_pcurves(
    brep: &BRep,
    cap_faces: &BTreeSet<FaceHandle>,
) -> ExternalCapSupportPcurveClassification {
    let mut out = ExternalCapSupportPcurveClassification::default();
    for &cap_face in cap_faces {
        let Some(half_edges) = face_outer_half_edges(brep, cap_face) else {
            continue;
        };
        for half_edge in half_edges {
            if !is_external_cap_support_half_edge(brep, cap_faces, half_edge) {
                continue;
            }
            let Some(pcurve) = brep.half_edges[half_edge].pcurve.as_ref() else {
                continue;
            };
            let hinted = brep.surface_boundary_hint(half_edge).is_some();
            if canonical_boundary_span(pcurve) {
                out.recognized_cap_boundary_hint_edges += hinted as usize;
            } else {
                out.unrecognized_cap_boundary_pcurve_edges += 1;
                out.unrecognized_cap_boundary_hint_edges += hinted as usize;
            }
        }
    }
    out
}

fn face_outer_half_edges(brep: &BRep, face: FaceHandle) -> Option<Vec<HalfEdgeHandle>> {
    let start = brep.loops.get(brep.faces.get(face)?.outer_loop)?.half_edge;
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        out.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            return Some(out);
        }
    }
    None
}

fn is_external_cap_support_half_edge(
    brep: &BRep,
    cap_faces: &BTreeSet<FaceHandle>,
    half_edge: HalfEdgeHandle,
) -> bool {
    let Some(twin_face) = brep
        .half_edges
        .get(half_edge)
        .and_then(|half_edge| brep.half_edges.get(half_edge.twin))
        .and_then(|twin| brep.loops.get(twin.loop_))
        .map(|loop_| loop_.face)
    else {
        return false;
    };
    !cap_faces.contains(&twin_face)
}

fn canonical_boundary_span(pcurve: &Pcurve) -> bool {
    let (t0, t1) = pcurve.parameter_range();
    let start = pcurve.point_at(t0);
    let end = pcurve.point_at(t1);
    (close(start.x, 0.0) && close(end.x, 0.0) && usable_span(start.y, end.y))
        || (close(start.x, 1.0) && close(end.x, 1.0) && usable_span(start.y, end.y))
        || (close(start.y, 0.0) && close(end.y, 0.0) && usable_span(start.x, end.x))
        || (close(start.y, 1.0) && close(end.y, 1.0) && usable_span(start.x, end.x))
}

fn usable_span(a: f64, b: f64) -> bool {
    (a - b).abs() > 1.0e-9 && (0.0..=1.0).contains(&a) && (0.0..=1.0).contains(&b)
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1.0e-9
}
