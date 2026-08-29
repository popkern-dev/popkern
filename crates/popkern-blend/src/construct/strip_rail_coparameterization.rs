//! Rail co-parameterization for swept (NURBS) blend strips (PK-139).
//!
//! Variable-radius sweeps assemble rail pcurves whose parameterization
//! diverges from the rail edge curves on curved rails (chord-length
//! edge params vs the skin's uniform stations — measured 1.755e-2
//! same-fraction deviation at a by-points bulge apex). The images
//! coincide but every consumer that samples curve and pcurve at the
//! same fraction (the mesher's boundary conformity above all) reads a
//! crack. `fix_stale_pcurves` cannot serve this: its staleness gate
//! scales with the body's accumulated imprecision, which such bodies
//! exceed.
//!
//! This scoped pass (same pattern as `settle_terminal_arc_faces`)
//! re-fits the offending half-edge pcurves through the projection
//! recompute — co-parameterized with the edge curve by construction —
//! under an improve-or-revert contract per half-edge.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Re-fit strip-face pcurves whose same-fraction deviation exceeds
/// the conformity gate. Returns the number of half-edges improved.
pub(crate) fn coparameterize_strip_rails(
    result: &mut BRep,
    faces: &[FaceHandle],
    tolerance: &Tolerance,
) -> usize {
    let gate = tolerance.linear().max(1.0e-9) * 10.0;
    let trace = popkern_core::trace::trace_enabled("terminal-arc");
    let mut improved = 0_usize;
    // Measurement phase: per-face offender detection is read-only and
    // face-local (a refit rewrites only that face's own half-edge
    // pcurves), so it runs in parallel under the `parallel` feature
    // (PK-302), returning each offending half-edge with its measured
    // deviation — the exact values the serial walk would compute,
    // since no mutation precedes them. The refit loop below stays
    // serial in the unchanged face order and consumes the recorded
    // deviations instead of re-measuring.
    let offenders: Vec<Vec<(popkern_core::topology::HalfEdgeHandle, f64)>> = {
        let result_ref: &BRep = result;
        let per_face = |&face: &FaceHandle| -> Vec<(popkern_core::topology::HalfEdgeHandle, f64)> {
            face_loop_half_edges(result_ref, face)
                .map(|hes| {
                    hes.iter()
                        .filter_map(|&he| {
                            let report = result_ref.measure_half_edge_pcurve_deviation(he)?;
                            (report.max_same_fraction_deviation > gate)
                                .then_some((he, report.max_same_fraction_deviation))
                        })
                        .collect()
                })
                .unwrap_or_default()
        };
        #[cfg(feature = "parallel")]
        {
            use rayon::prelude::*;
            faces.par_iter().map(per_face).collect()
        }
        #[cfg(not(feature = "parallel"))]
        {
            faces.iter().map(per_face).collect()
        }
    };
    for (&face, face_offenders) in faces.iter().zip(offenders) {
        if face_offenders.is_empty() {
            continue;
        }
        // Improve-or-revert per half-edge: clear the offenders, refit
        // the whole face through the projection recompute, and keep
        // each replacement only when it measures strictly better.
        let mut cleared = Vec::new();
        for (he, deviation) in face_offenders {
            let Some(original) = result.half_edges.get(he).and_then(|h| h.pcurve.clone()) else {
                continue;
            };
            result.half_edges[he].pcurve = None;
            cleared.push((he, original, deviation));
        }
        if cleared.is_empty() {
            continue;
        }
        let _ = result.recompute_face_pcurves(face);
        for (he, original, original_deviation) in cleared {
            let refit_deviation = result
                .half_edges
                .get(he)
                .and_then(|h| h.pcurve.as_ref())
                .and_then(|_| result.measure_half_edge_pcurve_deviation(he))
                .map(|r| r.max_same_fraction_deviation);
            match refit_deviation {
                Some(d) if d < original_deviation => {
                    improved += 1;
                    if trace {
                        eprintln!(
                            "[rail-copar] {he:?} on {face:?}: {original_deviation:.3e} -> {d:.3e}"
                        );
                    }
                }
                _ => {
                    if trace {
                        eprintln!(
                            "[rail-copar] {he:?} on {face:?}: refit {refit_deviation:?} not \
                             better than {original_deviation:.3e} — reverted"
                        );
                    }
                    if let Some(h) = result.half_edges.get_mut(he) {
                        h.pcurve = Some(original);
                    }
                }
            }
        }
    }
    improved
}

/// The outer-loop half-edges of a NURBS-surfaced face, or `None` when
/// the face is absent, unsurfaced, non-NURBS, or loop-broken. Shared
/// by the parallel pre-filter and the serial refit loop so both walk
/// the same sequence.
fn face_loop_half_edges(
    brep: &BRep,
    face: FaceHandle,
) -> Option<Vec<popkern_core::topology::HalfEdgeHandle>> {
    let is_nurbs = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Nurbs(_)));
    if !is_nurbs {
        return None;
    }
    let loop_handle = brep.faces.get(face)?.outer_loop;
    let start = brep.loops.get(loop_handle)?.half_edge;
    let mut hes = Vec::new();
    let mut cursor = start;
    for _ in 0..brep.half_edges.len() {
        hes.push(cursor);
        // A broken loop keeps the partial walk (matching the historic
        // serial measurement) rather than skipping the face.
        let Some(data) = brep.half_edges.get(cursor) else {
            break;
        };
        cursor = data.next;
        if cursor == start {
            break;
        }
    }
    Some(hes)
}
