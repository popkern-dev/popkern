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
    let trace = std::env::var_os("PK_TERMINAL_ARC_TRACE").is_some();
    let mut improved = 0_usize;
    for &face in faces {
        let is_nurbs = result
            .faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Nurbs(_)));
        if !is_nurbs {
            continue;
        }
        let Some(loop_handle) = result.faces.get(face).map(|f| f.outer_loop) else {
            continue;
        };
        let Some(start) = result.loops.get(loop_handle).map(|l| l.half_edge) else {
            continue;
        };
        let mut hes = Vec::new();
        let mut cursor = start;
        for _ in 0..result.half_edges.len() {
            hes.push(cursor);
            let Some(data) = result.half_edges.get(cursor) else {
                break;
            };
            cursor = data.next;
            if cursor == start {
                break;
            }
        }
        let mut needs_refit = false;
        for &he in &hes {
            let Some(report) = result.measure_half_edge_pcurve_deviation(he) else {
                continue;
            };
            if report.max_same_fraction_deviation <= gate {
                continue;
            }
            needs_refit = true;
        }
        if !needs_refit {
            continue;
        }
        // Improve-or-revert per half-edge: clear the offenders, refit
        // the whole face through the projection recompute, and keep
        // each replacement only when it measures strictly better.
        let mut cleared = Vec::new();
        for &he in &hes {
            let Some(report) = result.measure_half_edge_pcurve_deviation(he) else {
                continue;
            };
            if report.max_same_fraction_deviation <= gate {
                continue;
            }
            let Some(original) = result.half_edges.get(he).and_then(|h| h.pcurve.clone()) else {
                continue;
            };
            result.half_edges[he].pcurve = None;
            cleared.push((he, original, report.max_same_fraction_deviation));
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
