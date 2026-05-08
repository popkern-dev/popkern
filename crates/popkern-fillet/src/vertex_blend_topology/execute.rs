//! Phase 6c-C.2b — apex-region deletion *executor*.
//!
//! Consumes a `DeletionPlan` (from `plan.rs`) and rewrites the
//! BRep to match the target topology.
//!
//! # Status
//!
//! * **Implemented:** the safe `mef` step that inserts the new
//!   cap-tangent arc edges across the cap face's outer loop.
//!   Each `mef` splits the cap face into a (smaller real cap)
//!   and a (lens sub-face containing the obsolete generators
//!   and the rolling-ball-footprint vertex between them).
//!
//! * **Deferred to a follow-up:** the lens-collapse step.
//!   Naive `kef` of one of the lens's generator edges would
//!   merge the lens into one of the two adjacent cylinder
//!   faces.  But the *other* generator of the same lens
//!   borders a *different* cylinder face, so a second `kef`
//!   would merge incompatible cylinder surfaces.  The right
//!   sequence requires either:
//!   1. A non-Euler direct slot-map mutation (delete the lens
//!      face + 2 generators + footprint vertex; splice the
//!      cylinder face boundaries to skip them).
//!   2. A multi-face combined op that handles the cross-
//!      cylinder junction at each footprint atomically.
//!
//!   Both paths require careful re-routing of the cylinder
//!   faces' outer-loop half-edge `.next` / `.prev` chains and
//!   are best implemented after a proof-of-concept with
//!   explicit before/after BRep snapshots.
//!
//! Until the lens-collapse step lands, this executor returns a
//! `PartialExecution` report describing how far it got.  The
//! cap surface stays a Sphere primitive (correct centre after
//! the centroid fix in `surface_fix.rs`); the topology now
//! includes 3 extra edges (the cap-tangent arcs) on the cap
//! face's outer loop, alongside the existing 6 generators.
//!
//! Net empirical effect: G0 unchanged from the centroid-fix
//! baseline (still 2.93e-2 on cube_uniform) because the lens
//! sub-faces still exist and the metric reads the cap face's
//! union of all sub-faces.  Once the lens-collapse step lands,
//! G0 should drop to near-zero.

use popkern_core::geometry::CurveEval;
use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use super::execute_geometry::{
    build_cap_tangent_half_circle, find_loop_he_pair, find_shared_cylinder, fit_pcurve_along_arc,
    he_after_origin_in_loop,
};
use super::plan::DeletionPlan;

/// Outcome of running an executor against a plan.
///
/// `PartialExecution` is the only success variant today; the
/// `Complete` variant is reserved for when the lens-collapse
/// step is implemented.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub(crate) struct ExecutionReport {
    /// New cap-tangent arc edges inserted via `mef`.
    pub(crate) inserted_arcs: Vec<EdgeHandle>,
    /// Lens sub-faces created by the `mef` calls — each
    /// contains a rolling-ball-footprint vertex and 2 obsolete
    /// generator edges.  These need a follow-up collapse pass.
    pub(crate) lens_faces: Vec<FaceHandle>,
    /// Vertex handles of the rolling-ball footprints that
    /// remain (one per lens).  Their incident-edge count is
    /// still 2 (the 2 obsolete generators).
    pub(crate) remaining_footprints: Vec<VertexHandle>,
}

/// Errors from `execute_deletion_plan`.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub(crate) enum ExecuteError {
    /// Plan didn't pass `is_actionable`.  Caller should bail
    /// without mutating.
    PlanNotActionable,
    /// Couldn't find a half-edge originating from a kept
    /// vertex in the cap face's outer loop.  Means the cap
    /// face has been mutated since the plan was built.
    HalfEdgeLookupFailed { vertex: VertexHandle },
    /// `mef` failed during arc insertion.  Reports the BRep's
    /// error message verbatim.
    MefFailed(String),
    /// Couldn't determine which fillet cylinder hosts the
    /// cap-tangent arc for a given `(V_a, V_b)` pair.  Means
    /// the BRep doesn't have the expected structure (e.g. the
    /// kept vertices aren't both on a single cylinder).
    CylinderLookupFailed {
        v_a: VertexHandle,
        v_b: VertexHandle,
    },
    /// Actionable plan was missing a field required by execution.
    MissingPlanField(&'static str),
}

/// Run the executor.  See module docs for what's implemented.
#[allow(dead_code)]
pub(crate) fn execute_deletion_plan(
    brep: &mut BRep,
    plan: &DeletionPlan,
) -> Result<ExecutionReport, ExecuteError> {
    if !plan.is_actionable() {
        return Err(ExecuteError::PlanNotActionable);
    }
    let cap_face = plan
        .cap_face
        .ok_or(ExecuteError::MissingPlanField("cap_face"))?;
    let (cap_centre, cap_radius) = plan
        .cap_sphere
        .ok_or(ExecuteError::MissingPlanField("cap_sphere"))?;

    let mut inserted_arcs: Vec<EdgeHandle> = Vec::new();
    let mut lens_faces: Vec<FaceHandle> = Vec::new();
    let mut remaining_footprints: Vec<VertexHandle> = Vec::new();

    for &(v_a, v_b, _) in &plan.insert_cap_tangent_arcs {
        let cyl_surf = find_shared_cylinder(brep, v_a, v_b)
            .ok_or(ExecuteError::CylinderLookupFailed { v_a, v_b })?;

        let arc_curve =
            build_cap_tangent_half_circle(brep, cap_centre, cap_radius, &cyl_surf, v_a, v_b);

        // Walk the cap face's outer loop. The two half-edges
        // we want to mef between are: (a) the half-edge
        // originating at v_a in the cap face's loop and (b) the
        // half-edge originating at v_b in the cap face's loop.
        // These are typically separated by exactly ONE vertex
        // in the loop (the footprint between them); after the
        // earlier mef calls in this same execution the
        // separation may shrink to zero (loop has been split
        // already).
        let (he_a, he_b) = find_loop_he_pair(brep, cap_face, v_a, v_b)
            .ok_or(ExecuteError::HalfEdgeLookupFailed { vertex: v_a })?;

        // Capture the footprint vertex between v_a and v_b
        // before mef — it's the loop's vertex sandwiched between
        // v_a and v_b.
        let footprint = he_after_origin_in_loop(brep, he_a).filter(|&v| v != v_b);

        let mef = brep
            .mef(he_a, he_b)
            .map_err(|e| ExecuteError::MefFailed(format!("{e:?}")))?;
        let (t_start, t_end) = arc_curve.parameter_range();
        // Capture cap surface for pcurve fitting before mutating
        // the edge.
        let cap_surf_clone = brep
            .faces
            .get(cap_face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| s.geometry.clone());
        if let Some(edge_data) = brep.edges.get_mut(mef.edge) {
            edge_data.curve = Some(CurveOnEdge {
                geometry: popkern_core::geometry::CurveGeometry::Nurbs(arc_curve.clone()),
                t_start,
                t_end,
            });
        }
        // Phase 6c-C.2b/2d: fit pcurves on both half-edges of the
        // new arc.  Without these, downstream consumers (mesh,
        // bool, the G0/G1 metric) read `pcurve = None` and either
        // skip the edge or fall back to a stale value.  Both pcurves
        // are NURBS fits of `project_to_uv(arc_3d_sample)` against
        // the respective surface.
        // mef.half_edge is on the OLD face (cap), mef.half_edge_twin
        // is on the NEW face (lens).  Per `MefResult` docs.
        let cap_he = mef.half_edge;
        let lens_he = mef.half_edge_twin;
        if let Some(cap_surf) = cap_surf_clone {
            // cap_he runs v_a → v_b (the arc's natural direction);
            // fit forward.  lens_he is the twin running v_b → v_a;
            // fit reversed so its pcurve(t_lo) lands on v_b.
            let cap_pcurve = fit_pcurve_along_arc(&cap_surf, &arc_curve, false);
            let cyl_pcurve = fit_pcurve_along_arc(&cyl_surf, &arc_curve, true);
            if let Some(pc) = cap_pcurve {
                brep.half_edges[cap_he].pcurve = Some(pc.into());
            }
            if let Some(pc) = cyl_pcurve {
                brep.half_edges[lens_he].pcurve = Some(pc.into());
            }
        }
        inserted_arcs.push(mef.edge);
        lens_faces.push(mef.face);
        if let Some(fp) = footprint {
            remaining_footprints.push(fp);
        }
    }

    Ok(ExecutionReport {
        inserted_arcs,
        lens_faces,
        remaining_footprints,
    })
}
