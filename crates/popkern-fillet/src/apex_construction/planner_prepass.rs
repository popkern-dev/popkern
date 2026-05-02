//! Phase 6c-A2 §28.H Session 1 — planner pre-pass.
//!
//! Runs the A2 planner BEFORE Phase A's perp-walk anchor placement,
//! using only edge geometries (radii, axes, host face normals from
//! `EdgeAnalysis` + sphere-fit termination math).  The planner
//! itself reads from BRep state and `FilletSurfaceResult`, both of
//! which are available pre-Phase-A:
//!
//! * **BRep state**: the pre-Phase-A snapshot (cube faces with
//!   original edges; no perp-anchors yet).
//! * **FilletSurfaceResult**: built directly from `EdgeGeometry`
//!   (the surface and trim curves are computed in
//!   `fillet_surface::build_fillet_surface` BEFORE Phase A).
//!
//! The output is a `HashMap<VertexHandle, Vec<Point3>>` mapping
//! each apex vertex (cube/pyramid corner with ≥3 incident filleted
//! edges) to the planner's emitted `apex_vertex` target positions
//! (2N entries per apex, alternating FilletEnd / FaceEdge boundary
//! endpoints).
//!
//! Subsequent §28.H sessions (2-3) consume this map to position
//! Phase A's anchors at planner targets, so post-Phase-A snap finds
//! a zero-distance match instead of needing §28.C relocation.
//!
//! See ADR `docs/decisions/2026-04-30-phase-6c-a2-direct-apex-
//! construction.md` §36.3.

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, VertexHandle};

use super::plan::plan_apex_resolution_lenient;
use crate::construct::EdgeGeometry;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{FilletRadius, FilletSpec};

/// Run the A2 planner against the pre-Phase-A BRep state.
///
/// For each vertex `V` with ≥3 incident filleted edges, build the
/// per-apex `(FilletSpec[], FilletSurfaceResult[])` lists from
/// `edge_geometries`, then call
/// [`plan_apex_resolution_lenient`].  Successful plans (whose
/// rung is Direct or Setback — not Open) contribute their
/// `target_topology.apex_vertices` to the returned map.
///
/// Planner errors are silently swallowed: a vertex that can't be
/// planned (e.g. degenerate sphere-cylinder intersection on
/// extreme radii) just doesn't appear in the output.  Phase A
/// then runs unmodified for that vertex's incident edges, falling
/// back to perp-walk-math anchor positions.
///
/// `brep` should be the pre-Phase-A snapshot — calling this
/// function after Phase A produces the same output but wastes
/// work (the planner reads only original edge endpoints, not
/// perp-anchors).
pub(crate) fn plan_apex_targets_pre_phase_a(
    brep: &BRep,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> HashMap<VertexHandle, Vec<Point3>> {
    let trace = std::env::var_os("POPKERN_A2_PREPASS_TRACE").is_some();
    let mut out: HashMap<VertexHandle, Vec<Point3>> = HashMap::new();

    // Build per-vertex incidence: for each V, list the
    // `edge_geometries` indices whose edge has V as an endpoint.
    let mut v_to_fillet_indices: HashMap<VertexHandle, Vec<usize>> = HashMap::new();
    for (i, eg) in edge_geometries.iter().enumerate() {
        let edge = eg.analysis.edge;
        let Some(edge_data) = brep.edges.get(edge) else {
            continue;
        };
        let primary = edge_data.half_edge;
        let Some(he) = brep.half_edges.get(primary) else {
            continue;
        };
        let Some(twin_he) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let v_a = he.origin;
        let v_b = twin_he.origin;
        v_to_fillet_indices.entry(v_a).or_default().push(i);
        v_to_fillet_indices.entry(v_b).or_default().push(i);
    }

    // For each apex candidate (V with ≥3 incident filleted edges),
    // run the planner.
    for (vh, fillet_indices) in &v_to_fillet_indices {
        if fillet_indices.len() < 3 {
            continue;
        }

        let specs: Vec<FilletSpec> = fillet_indices
            .iter()
            .map(|&i| build_fillet_spec(&edge_geometries[i]))
            .collect();
        let surfaces: Vec<FilletSurfaceResult> = fillet_indices
            .iter()
            .map(|&i| build_fillet_surface(&edge_geometries[i]))
            .collect();

        match plan_apex_resolution_lenient(brep, *vh, &specs, &surfaces, tolerance, None) {
            Ok(plan) => {
                if plan.cap_kind.is_open() {
                    if trace {
                        eprintln!("[a2-prepass] V={vh:?}: planner returned Open rung; skipping");
                    }
                    continue;
                }
                let targets = plan.target_topology.apex_vertices;
                if trace {
                    eprintln!(
                        "[a2-prepass] V={vh:?}@{:?}: {} planner targets",
                        brep.vertices.get(*vh).map(|v| v.position),
                        targets.len()
                    );
                    for (k, t) in targets.iter().enumerate() {
                        eprintln!("[a2-prepass]   target[{k}]={t:?}");
                    }
                }
                out.insert(*vh, targets);
            }
            Err(e) => {
                if trace {
                    eprintln!(
                        "[a2-prepass] V={vh:?}: planner failed ({e:?}); falling back to \
                         perp-walk-math at this vertex"
                    );
                }
            }
        }
    }
    out
}

/// Construct a [`FilletSpec`] from one [`EdgeGeometry`].  Mirrors
/// the spec construction in `crate::construct::a2_dispatch::
/// try_a2_apex_construction` so the planner sees the same data
/// pre- and post-Phase-A.
fn build_fillet_spec(eg: &EdgeGeometry) -> FilletSpec {
    let radius = if (eg.radius_at_start - eg.radius_at_end).abs() < 1e-12 {
        FilletRadius::Constant(eg.radius_at_start)
    } else {
        FilletRadius::Variable {
            start: eg.radius_at_start,
            end: eg.radius_at_end,
        }
    };
    FilletSpec {
        edge: eg.analysis.edge,
        radius,
    }
}

/// Construct a [`FilletSurfaceResult`] from one [`EdgeGeometry`].
/// Mirrors the surface construction in `crate::construct::a2_
/// dispatch::try_a2_apex_construction`.
fn build_fillet_surface(eg: &EdgeGeometry) -> FilletSurfaceResult {
    let downgrade_reason = eg.continuity_downgrade.as_ref().map(|(_, s)| s.clone());
    FilletSurfaceResult {
        surface: eg.surface.clone(),
        trim_curve_1: eg.trim_3d_1.clone(),
        trim_curve_2: eg.trim_3d_2.clone(),
        continuity_downgrade: downgrade_reason,
    }
}
