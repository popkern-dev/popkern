//! Cap-vertex-aware admissibility gate for Phase 5 split_edge call sites.
//!
//! Stage B1 of `docs/decisions/2026-06-21-kp20-cap-vertex-aware-
//! support-cell-split.md`. Wraps `brep.split_edge(edge, t)` call sites
//! in `construct::kp20_polar_*` with a check that the resulting vertex
//! 3D position is reachable by the adjacent face's existing boundary
//! geometry within the per-row displacement cap (`MAX_BOUNDARY_ROW_
//! DISPLACEMENT = 0.1`).
//!
//! The four sites surveyed in
//! `docs/reports/2026-06-21-kp20-split-edge-sites.md` are wired in
//! Stage B2 (one site per substage). Stage B1 ships the helper +
//! unit tests; the helper is not yet called by any production code.
//!
//! Algorithm shape:
//!   1. Evaluate the edge curve at `candidate_t` to get the 3D
//!      position `p` where the new vertex would land.
//!   2. For each face incident on the edge, find the boundary
//!      curve's parameter that lands closest to `p` using
//!      [`super::support_cell_boundary_projection::project_to_boundary`].
//!   3. Compute the 3D displacement from `p` to the closest on-curve
//!      point. The maximum across both faces is the displacement the
//!      boundary-solve would have to absorb.
//!   4. If the worst displacement is within `cap`, the candidate is
//!      admissible — return [`SplitParameterDecision::UseCandidate`].
//!   5. Else search along the edge for a `t` whose evaluation lands
//!      within cap of both faces' boundaries. Bisection between
//!      [0, candidate_t] and [candidate_t, 1] up to a small
//!      iteration cap. Return [`SplitParameterDecision::UseProjected`]
//!      on success.
//!   6. If no such `t` exists, return
//!      [`SplitParameterDecision::SkipUnreachable`] with the best
//!      observed displacement; the caller falls through to the
//!      completion ladder's next candidate.
//!
//! Stage B1 ships a structural implementation that fills in step 4
//! (admissibility check on the candidate) and steps 5–6 fall back to
//! `SkipUnreachable` without searching — Stage B2 expands the search
//! once empirical evidence on hexagon6_variable_route tells us what
//! search shape is needed.

use popkern_core::topology::{BRep, EdgeHandle};

/// Outcome of the admissibility check at a `brep.split_edge` call site.
#[derive(Debug, Clone, Copy, PartialEq)]
#[allow(dead_code)] // first caller lands in Stage B2
pub(crate) enum SplitParameterDecision {
    /// Caller's candidate parameter is admissible — proceed with
    /// `brep.split_edge(edge, candidate_t)`.
    UseCandidate(f64),
    /// A different parameter is closer to a reachable on-curve point;
    /// caller should split at this projected value instead.
    UseProjected(f64),
    /// No parameter exists at which the resulting vertex would be
    /// reachable by both faces' boundaries within `cap`. Caller should
    /// skip the split and let the completion ladder advance.
    SkipUnreachable {
        candidate_t: f64,
        best_displacement: f64,
    },
}

/// `MAX_BOUNDARY_ROW_DISPLACEMENT` from the boundary fit
/// (`kp20_polar_cap_support_boundary_fit.rs:7`). Mirrored here to avoid
/// a cross-module borrow at the call site.
#[allow(dead_code)]
pub(crate) const SPLIT_SITE_DISPLACEMENT_CAP: f64 = 0.1;

/// Check whether splitting `edge` at `candidate_t` would land a vertex
/// reachable by the existing surfaces of both incident faces.
///
/// Stage B1 wires the displacement-driven decision branches; the
/// projection/search branches return `SkipUnreachable` until Stage B2
/// supplies empirical evidence for the search shape on
/// hexagon6_variable_route.
#[allow(dead_code)]
pub(crate) fn cap_vertex_aware_split_parameter(
    brep: &BRep,
    edge: EdgeHandle,
    candidate_t: f64,
    cap: f64,
) -> SplitParameterDecision {
    let displacement = displacement_at(brep, edge, candidate_t).unwrap_or(0.0);
    decide_from_displacement(candidate_t, displacement, cap)
}

/// Pure decision logic separated from the BRep lookup so the three
/// outcome branches are unit-testable without synthetic BRep
/// fixtures.
#[allow(dead_code)]
fn decide_from_displacement(
    candidate_t: f64,
    displacement: f64,
    cap: f64,
) -> SplitParameterDecision {
    if displacement <= cap {
        SplitParameterDecision::UseCandidate(candidate_t)
    } else {
        // Stage B1 always falls through to skip — Stage B2 expands
        // here once empirical evidence tells us what to search.
        SplitParameterDecision::SkipUnreachable {
            candidate_t,
            best_displacement: displacement,
        }
    }
}

/// Worst-case 3D distance from the edge-curve evaluation at `t` to
/// the closest on-curve point of either incident face's existing
/// boundary surface. Returns `None` if the displacement cannot be
/// estimated from the available BRep data (e.g. an edge curve is
/// missing or a face has no surface).
///
/// Stage B1 surface lookup is intentionally absent — the candidate
/// site's face geometry will be threaded in by Stage B2. We exercise
/// the edge-curve evaluation path so the BRep-shape assumptions are
/// pinned by tests, and return `Some(0.0)` so the decision logic
/// routes through `UseCandidate`.
#[allow(dead_code)]
fn displacement_at(brep: &BRep, edge: EdgeHandle, _t: f64) -> Option<f64> {
    // Validate the BRep shape so Stage B2 can rely on these
    // invariants when it threads the surface lookup in: the edge
    // must be live and carry a curve representative.
    let _ = brep.edges.get(edge)?.curve.as_ref()?;
    Some(0.0)
}

#[cfg(test)]
#[path = "kp20_polar_support_cell_cap_vertex_split_site_tests.rs"]
mod tests;
