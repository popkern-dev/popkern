//! §13.1(1) cusp-subdivide repair.
//!
//! Spec text: *"Cusp in spine/surface (vanishing tangent): subdivide,
//! reparameterize."* Interior cusps on a fillet edge curve break the
//! rolling-ball contact solver's Newton iteration (the local frame
//! goes singular as the tangent vanishes), so the spec's repair
//! strategy is to split the edge at the cusp and fillet each
//! sub-segment independently. The downstream vertex-blend pass then
//! joins the two sub-fillets at the new vertex.
//!
//! This module implements the subdivision as a **pre-pass** that
//! runs at the top of [`crate::fillet_edges_internal`], before
//! Phase 1 (Analyze). It clones the input BRep, scans every
//! requested edge for interior cusps, calls
//! [`popkern_core::topology::BRep::split_edge`] at each cusp
//! parameter, and replaces each subdivided spec with two sub-specs
//! pointing at the new edges.
//!
//! ## Detector
//!
//! [`find_interior_cusp`] samples the edge curve's un-normalized
//! first derivative at N interior parameter points. The `derivative_at`
//! return preserves magnitude (unlike `tangent_at`, which is
//! normalized and always unit-length or exactly zero). A sample
//! whose magnitude falls below `tolerance.linear()` indicates a
//! cusp; the detector returns the normalized parameter `t ∈ (0, 1)`
//! suitable for [`BRep::split_edge`]. If multiple samples qualify,
//! returns the one with the smallest magnitude (closest to the
//! cusp point).
//!
//! ## Iterative application
//!
//! A single edge may contain multiple interior cusps. The pre-pass
//! applies one split per iteration and re-scans the growing spec
//! list until no more cusps are found or a safety cap is reached.
//! Each split halves the remaining cusps in the worst case, so the
//! iteration count is O(log) in practice.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

use crate::types::FilletSpec;

/// Maximum pre-pass iterations before giving up. Each iteration
/// applies at most one cusp split per spec; a single edge rarely
/// has more than 2–3 interior cusps in practice, so this cap is
/// a generous safety net rather than a normal-case limit.
const MAX_SUBDIVIDE_ITERATIONS: usize = 32;

/// Detect an interior cusp on the given edge and return the
/// normalized parameter `t ∈ (0, 1)` where it was found, or `None`
/// when no interior cusp is detected.
///
/// Samples the edge curve's un-normalized first derivative at 20
/// interior points (excluding both endpoints, which are handled
/// separately by
/// [`crate::degeneracy::check_edge_tangent_cusp`]) and returns the
/// parameter whose derivative magnitude is both smallest **and**
/// below `tolerance.linear()`. This detector is intentionally
/// conservative — it does not flag wide flat regions where the
/// derivative merely happens to be small at one sample, only
/// genuine cusps where the derivative vanishes.
pub fn find_interior_cusp(brep: &BRep, edge: EdgeHandle, tolerance: &Tolerance) -> Option<f64> {
    let edge_data = brep.edges.get(edge)?;
    let curve_on_edge = edge_data.curve.as_ref()?;
    let curve = &curve_on_edge.geometry;
    let t_start = curve_on_edge.t_start;
    let t_end = curve_on_edge.t_end;
    let span = t_end - t_start;
    if span.abs() < 1e-15 {
        return None;
    }

    // Sample the interior — i in 1..N (skip endpoints 0 and N).
    // N = 20 gives a step of 5% of the parameter range, which
    // catches any cusp wider than that and comfortably most
    // narrower cusps via the continuous approach to zero.
    const N: usize = 20;
    let mut best: Option<(f64, f64)> = None; // (t_normalized, magnitude)
    let tol = tolerance.linear();
    for i in 1..N {
        let frac = (i as f64) / (N as f64);
        let t_native = t_start + span * frac;
        let d = curve.derivative_at(t_native, 1);
        let mag = d.length();
        if mag < tol {
            let better = best.map(|(_, bm)| mag < bm).unwrap_or(true);
            if better {
                best = Some((frac, mag));
            }
        }
    }
    best.map(|(t, _)| t)
}

/// §13.1(1) pre-pass: subdivide fillet edges at every interior cusp.
///
/// Returns the (cloned, possibly mutated) BRep and a new spec list
/// in which each original spec either survives unchanged (no cusp)
/// or is replaced by two sub-specs pointing at the split edges.
///
/// Iterates until no more cusps are found or
/// [`MAX_SUBDIVIDE_ITERATIONS`] is reached. On a healthy input
/// (no cusps), the first iteration bails out with the original
/// spec list and a trivially-cloned BRep, so the pre-pass is
/// essentially free in the common case — the only work is the
/// `N=20`-sample scan per spec.
///
/// Errors from `split_edge` are swallowed: the detector short-
/// circuits on any split that fails (returning the original spec
/// list) and lets the downstream pipeline report whatever failure
/// naturally surfaces. This keeps the pre-pass's contract simple:
/// *"if a clean subdivision is possible, apply it; otherwise don't
/// interfere."*
pub fn subdivide_cusp_edges(
    brep: &BRep,
    specs: &[FilletSpec],
    tolerance: &Tolerance,
) -> (BRep, Vec<FilletSpec>) {
    let mut current_brep = brep.clone();
    let mut current_specs: Vec<FilletSpec> = specs.to_vec();

    for _ in 0..MAX_SUBDIVIDE_ITERATIONS {
        let mut next_specs: Vec<FilletSpec> = Vec::with_capacity(current_specs.len());
        let mut did_split = false;
        for spec in &current_specs {
            match find_interior_cusp(&current_brep, spec.edge, tolerance) {
                Some(t_normalized) => {
                    match current_brep.split_edge(spec.edge, t_normalized) {
                        Ok(result) => {
                            next_specs.push(FilletSpec {
                                edge: result.edge_a,
                                radius: spec.radius.clone(),
                            });
                            next_specs.push(FilletSpec {
                                edge: result.edge_b,
                                radius: spec.radius.clone(),
                            });
                            did_split = true;
                        }
                        Err(_) => {
                            // split_edge failed (handle stale,
                            // t out of range, etc.). Keep the
                            // original spec — the downstream
                            // detectors will report whatever the
                            // real problem is.
                            next_specs.push(spec.clone());
                        }
                    }
                }
                None => next_specs.push(spec.clone()),
            }
        }
        current_specs = next_specs;
        if !did_split {
            break;
        }
    }

    (current_brep, current_specs)
}

#[cfg(test)]
#[path = "cusp_subdivide_tests.rs"]
mod tests;
