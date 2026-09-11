//! PK-478 E1 diagnostic: closed-edge side-witness census.
//!
//! For every CLOSED, non-seam edge of a body, measures the three
//! candidate registration witnesses side by side so their agreement
//! with the body's stored-loop convention can be asserted per
//! population (the E1 exploration of
//! `docs/architecture/pk478-sense-aware-merge.md`):
//!
//! - **W1** — the shipping witness's basis: material-left of the
//!   RAY-PARITY outward normal (geometry-only, reversal-invariant).
//! - **W2** — material-left of the FLAG-FOLDED normal (covariant
//!   with a body reversal on honest faces).
//! - **PC** — the half-edge pcurve's claimed traversal
//!   (`pcurve_claims_forward`), the stored-loop truth wherever
//!   pcurves are honestly authored (fresh constructor/sweep/blend
//!   output; NOT trustworthy on legacy or re-registered bodies).
//!
//! This is a promoted diagnostic (structured rows a test can
//! assert), not a trace. It changes no behavior.

use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use super::graph_build_closed_side::{ClosedSideContext, census_probe};

/// One closed edge's witness measurements, per half-edge side.
#[derive(Debug, Clone)]
pub struct ClosedSideCensusRow {
    /// The closed edge.
    pub edge: EdgeHandle,
    /// Face of the edge's canonical half-edge.
    pub canonical_face: FaceHandle,
    /// Face of the twin half-edge.
    pub twin_face: FaceHandle,
    /// W1 vote for the canonical half-edge's face (material-left of
    /// ray-parity outward; positive = that face would traverse the
    /// stored curve forward). `None` = unmeasurable.
    pub w1_canonical: Option<f64>,
    /// W1 vote for the twin's face.
    pub w1_twin: Option<f64>,
    /// W2 vote for the canonical half-edge's face (material-left of
    /// the flag-folded normal).
    pub w2_canonical: Option<f64>,
    /// W2 vote for the twin's face.
    pub w2_twin: Option<f64>,
    /// Pcurve claim for the canonical half-edge: `Some(true)` = its
    /// pcurve traces the stored curve forward.
    pub pc_canonical: Option<bool>,
    /// Pcurve claim for the twin half-edge.
    pub pc_twin: Option<bool>,
    /// Winding-honesty verdicts for the two faces (canonical, twin).
    pub honest: (Option<bool>, Option<bool>),
}

/// Measure every closed, non-seam edge of `brep`.
#[doc(hidden)]
pub fn closed_side_census(brep: &BRep) -> Vec<ClosedSideCensusRow> {
    let honesty: std::collections::HashMap<FaceHandle, Option<bool>> = brep
        .audit_winding_honesty()
        .into_iter()
        .map(|v| (v.face, v.honest))
        .collect();
    let mut ctx = ClosedSideContext::new(brep);
    let mut rows = Vec::new();
    for (eh, edge) in &brep.edges {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        if he.origin != twin.origin {
            continue; // open edge
        }
        let face_of = |lh| brep.loops.get(lh).map(|l| l.face);
        let (Some(f_can), Some(f_twin)) = (face_of(he.loop_), face_of(twin.loop_)) else {
            continue;
        };
        if f_can == f_twin {
            continue; // seam edge: exempt in the shipping witness too
        }
        let probe = census_probe(&mut ctx, edge.half_edge);
        rows.push(ClosedSideCensusRow {
            edge: eh,
            canonical_face: f_can,
            twin_face: f_twin,
            w1_canonical: probe.w1_own,
            w1_twin: probe.w1_twin,
            w2_canonical: probe.w2_own,
            w2_twin: probe.w2_twin,
            pc_canonical: probe.pc_own,
            pc_twin: probe.pc_twin,
            honest: (
                honesty.get(&f_can).copied().flatten(),
                honesty.get(&f_twin).copied().flatten(),
            ),
        });
    }
    rows
}
