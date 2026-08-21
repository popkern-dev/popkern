//! HONEST-DECLINE strand gate for graph-merge assembly (PK-304/305).
//!
//! A stranded half is a region no kept fragment claimed — the
//! self-twin masking downstream keeps the body validate()-clean while
//! the SHAPE silently loses the unclaimed region (measured: the
//! rotated torus × torus union served with 8 masked strands and 5.5%
//! of the union volume missing). A boolean that cannot claim both
//! sides of every edge has mis-partitioned an operand; declining
//! transactionally is the contract. POPKERN_ALLOW_STRANDS=1 restores
//! the old masking for diagnosis.

use std::collections::HashMap;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

use crate::error::BoolDiagnostics;
use crate::error::{BoolError, Stage};
use crate::merge::graph::IntersectionGraph;

/// Decline when a lethal strand class is present. Benign strand
/// classes, calibrated against the tuned corpus:
/// 1. Open-sheet operands (or synthetic test fragments with empty
///    source breps) legitimately leave boundary twins unwired.
/// 2. Fully-unclaimed edges drop whole with two discarded regions (a
///    subtraction's record between two discarded sides).
/// 3. Zero-length record residue and original/seam half-wiring (the
///    cone-apex seam piece class) mask harmlessly.
/// 4. A strand covered by a wired PARALLEL DUPLICATE in the same
///    direction (duplicate record copies each claimed once — the
///    sphere-cube subtraction's mirrored rim records).
///
/// Lethal = half-claimed, non-degenerate, SHARED-record, uncovered.
#[allow(clippy::too_many_arguments)]
pub(super) fn decline_on_lethal_strands(
    brep: &BRep,
    graph: &IntersectionGraph,
    he_map: &HashMap<usize, HalfEdgeHandle>,
    dangling_hes: &[HalfEdgeHandle],
    brep_a: &BRep,
    brep_b: &BRep,
    tolerance: &Tolerance,
) -> Result<(), BoolDiagnostics> {
    let default_loop = LoopHandle::default();
    let operands_closed = !brep_a.faces.is_empty()
        && !brep_b.faces.is_empty()
        && !brep_a.has_boundary_edges()
        && !brep_b.has_boundary_edges();
    if !operands_closed || std::env::var_os("POPKERN_ALLOW_STRANDS").is_some_and(|v| v == "1") {
        return Ok(());
    }
    let rev_map: HashMap<HalfEdgeHandle, usize> = he_map.iter().map(|(g, h)| (*h, *g)).collect();
    let lethal: Vec<HalfEdgeHandle> = dangling_hes
        .iter()
        .copied()
        .filter(|&heh| {
            let Some(he) = brep.half_edges.get(heh) else {
                return false;
            };
            let twin_wired = brep
                .half_edges
                .get(he.twin)
                .is_some_and(|t| t.loop_ != default_loop);
            if !twin_wired {
                return false;
            }
            let is_shared = rev_map.get(&heh).is_some_and(|&g| {
                matches!(
                    graph.half_edges[g].source,
                    crate::merge::graph::EdgeSource::Shared { .. }
                )
            });
            if !is_shared {
                return false;
            }
            let a = brep.vertices.get(he.origin).map(|v| v.position);
            let b = brep
                .half_edges
                .get(he.twin)
                .and_then(|t| brep.vertices.get(t.origin))
                .map(|v| v.position);
            let (Some(a), Some(b)) = (a, b) else {
                return false;
            };
            let gate = tolerance.linear() * 10.0;
            if a.distance_to(&b) <= gate {
                return false;
            }
            let covered = brep.half_edges.iter().any(|(hh, cand)| {
                if hh == heh || cand.loop_ == default_loop {
                    return false;
                }
                let ca = brep.vertices.get(cand.origin).map(|v| v.position);
                let cb = brep
                    .half_edges
                    .get(cand.twin)
                    .and_then(|t| brep.vertices.get(t.origin))
                    .map(|v| v.position);
                match (ca, cb) {
                    (Some(ca), Some(cb)) => ca.distance_to(&a) < gate && cb.distance_to(&b) < gate,
                    _ => false,
                }
            });
            !covered
        })
        .collect();
    if lethal.is_empty() {
        return Ok(());
    }
    if std::env::var_os("PK_BOOL_TRACE").is_some() {
        for &heh in &lethal {
            let Some(he) = brep.half_edges.get(heh) else {
                continue;
            };
            let a = brep.vertices.get(he.origin).map(|v| v.position);
            let b = brep
                .half_edges
                .get(he.twin)
                .and_then(|t| brep.vertices.get(t.origin))
                .map(|v| v.position);
            if let (Some(a), Some(b)) = (a, b) {
                eprintln!(
                    "  [lethal-strand] {heh:?} ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4}) ghe {:?}",
                    a.x,
                    a.y,
                    a.z,
                    b.x,
                    b.y,
                    b.z,
                    rev_map.get(&heh)
                );
            }
        }
    }
    Err(BoolDiagnostics::from_error(
        BoolError::MergeFailed {
            stage: format!("{:?}", Stage::Merge),
            reason: format!(
                "assembly stranded {} half-claimed intersection-record side(s) \
                 ({} dangling total); declining rather than serving a body \
                 with missing regions",
                lethal.len(),
                dangling_hes.len()
            ),
        },
        Stage::Merge,
    ))
}

/// Mask benign dangling halves: remove each unwired half, self-twin
/// its surviving twin, strip the edge's curve (the watertight check
/// would read a self-pointing twin), retire fully-dead edges, and
/// drop orphaned vertices.
pub(super) fn mask_dangling_halves(
    brep: &mut BRep,
    graph: &IntersectionGraph,
    he_map: &HashMap<usize, HalfEdgeHandle>,
    dangling_hes: &[HalfEdgeHandle],
) {
    // Collect twin info before mutating.
    let twin_fixups: Vec<(HalfEdgeHandle, HalfEdgeHandle)> = dangling_hes
        .iter()
        .filter_map(|&heh| brep.half_edges.get(heh).map(|he| (heh, he.twin)))
        .collect();

    if std::env::var_os("PK_BOOL_TRACE").is_some() {
        let rev: std::collections::HashMap<HalfEdgeHandle, usize> =
            he_map.iter().map(|(g, h)| (*h, *g)).collect();
        for (heh, _) in &twin_fixups {
            if let Some(he) = brep.half_edges.get(*heh) {
                let a = brep.vertices.get(he.origin).map(|v| v.position);
                let b = brep
                    .half_edges
                    .get(he.twin)
                    .and_then(|t| brep.vertices.get(t.origin))
                    .map(|v| v.position);
                let src = rev
                    .get(heh)
                    .map(|&g| format!("ghe {g} {:?}", graph.half_edges[g].source));
                if let (Some(a), Some(b)) = (a, b) {
                    eprintln!(
                        "  [assemble-strand] {heh:?} ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}) {src:?}",
                        a.x, a.y, a.z, b.x, b.y, b.z
                    );
                }
            }
        }
    }
    for (heh, twin_h) in &twin_fixups {
        // Point the surviving twin's twin to itself so it doesn't
        // reference a removed handle. A self-twinned survivor whose
        // edge later degenerates to zero length is retired by the
        // spur/self-twin postprocess pass (PK-293) — the walk-level
        // half-wiring that strands these halves is the recorded
        // follow-up.
        if let Some(twin) = brep.half_edges.get_mut(*twin_h) {
            twin.twin = *twin_h;
        }
        // Remove curve geometry from the edge — the watertight check
        // evaluates curves at endpoints and compares to vertex positions,
        // but with a self-pointing twin the check would fail.
        if let Some(surviving) = brep.half_edges.get(*twin_h) {
            let edge_h = surviving.edge;
            if let Some(edge) = brep.edges.get_mut(edge_h) {
                edge.curve = None;
            }
        }
        brep.half_edges.remove(*heh);
    }

    // Remove edges where BOTH half-edges were removed.  If either
    // half-edge survives, update `edge.half_edge` to point to the
    // surviving one so validation finds a valid handle.
    let mut edges_to_remove: Vec<EdgeHandle> = Vec::new();
    let mut edges_to_fixup: Vec<(EdgeHandle, HalfEdgeHandle)> = Vec::new();

    for (eh, edge) in brep.edges.iter() {
        let canonical_alive = brep.half_edges.contains_key(edge.half_edge);
        if canonical_alive {
            continue; // Edge is fine.
        }
        // Canonical was removed — check if any other half-edge references
        // this edge.
        let survivor = brep
            .half_edges
            .iter()
            .find(|(_, he)| he.edge == eh)
            .map(|(heh, _)| heh);
        match survivor {
            Some(surviving_he) => edges_to_fixup.push((eh, surviving_he)),
            None => edges_to_remove.push(eh),
        }
    }

    for (eh, surviving_he) in &edges_to_fixup {
        if let Some(edge) = brep.edges.get_mut(*eh) {
            edge.half_edge = *surviving_he;
        }
    }

    for eh in &edges_to_remove {
        brep.edges.remove(*eh);
    }

    // Fix up vertex outgoing handles that point to removed half-edges.
    let outgoing_map: HashMap<VertexHandle, HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .map(|(heh, he)| (he.origin, heh))
        .collect();
    for (vh, vertex) in &mut brep.vertices {
        if !brep.half_edges.contains_key(vertex.outgoing) {
            if let Some(&heh) = outgoing_map.get(&vh) {
                vertex.outgoing = heh;
            } else {
                vertex.outgoing = HalfEdgeHandle::default();
            }
        }
    }

    // Remove orphaned vertices.
    let orphan_verts: Vec<VertexHandle> = brep
        .vertices
        .iter()
        .filter(|(_, v)| v.outgoing == HalfEdgeHandle::default())
        .map(|(vh, _)| vh)
        .collect();

    for vh in &orphan_verts {
        brep.vertices.remove(*vh);
    }
}
