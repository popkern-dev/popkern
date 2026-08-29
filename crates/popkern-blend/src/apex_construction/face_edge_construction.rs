//! Face-edge construction for constructive apex building.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use super::constructive_geometry::attach_curve_to_edge_aligned;
use super::constructive_topology::find_outgoing_half_edge;
use super::error::ApexBuildError;
use super::plan::ApexPlan;

pub(super) fn construct_sphere_cap_face_from_v_incident_edges(
    brep: &mut BRep,
    plan: &ApexPlan,
) -> Result<FaceHandle, ApexBuildError> {
    crate::construct::blend_insert::kef_cross_face_v_edges(brep, plan.apex_vertex).map_err(|e| {
        ApexBuildError::PlanMismatch {
            reason: format!("construct_sphere_cap_face_from_v_incident_edges: {e}"),
        }
    })
}

/// Create the FaceEdge bridges between adjacent sector corners, one
/// `mef` per BRIDGED corner. Welded corners (shared-corner contract,
/// PK-219) have no bridge — their slot holds `None`, and the two
/// sectors share the corner vertex directly.
pub(super) fn create_gregory_face_edges(
    brep: &mut BRep,
    plan: &ApexPlan,
    apex_vertex_handles: &[VertexHandle],
) -> Result<Vec<Option<EdgeHandle>>, ApexBuildError> {
    let n = plan.cap_kind.n_fillets();
    if apex_vertex_handles.len() != 2 * n {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "create_gregory_face_edges: expected 2N={} apex_vertex_handles, got {}",
                2 * n,
                apex_vertex_handles.len()
            ),
        });
    }

    let trace = crate::direct_apex_trace::enabled();
    let mut face_edges: Vec<Option<EdgeHandle>> = Vec::with_capacity(n);

    for k in 0..n {
        let next_k = (k + 1) % n;
        let t_b_k = apex_vertex_handles[2 * k + 1];
        let t_a_next = apex_vertex_handles[2 * next_k];
        if plan.cap_kind.corner_is_welded(k) {
            if t_b_k != t_a_next {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "create_gregory_face_edges: corner {k} is declared welded but slots \
                         hold distinct vertices {t_b_k:?}/{t_a_next:?}"
                    ),
                });
            }
            if trace {
                eprintln!(
                    "[direct-apex-trace §14.7/B] corner {k} welded at {t_b_k:?} — no FaceEdge"
                );
            }
            face_edges.push(None);
            continue;
        }

        let he_origin =
            find_outgoing_half_edge(brep, t_b_k).ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "create_gregory_face_edges: T_b[{k}]={t_b_k:?} has no outgoing half-edge"
                ),
            })?;
        let he_dest = find_outgoing_half_edge(brep, t_a_next).ok_or_else(|| {
            ApexBuildError::PlanMismatch {
                reason: format!(
                    "create_gregory_face_edges: T_a[{next_k}]={t_a_next:?} has no outgoing half-edge"
                ),
            }
        })?;

        let loop_origin = brep.half_edges[he_origin].loop_;
        let loop_dest = brep.half_edges[he_dest].loop_;
        if loop_origin != loop_dest {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "create_gregory_face_edges: T_b[{k}] and T_a[{next_k}] are not co-located \
                     (loops {loop_origin:?} vs {loop_dest:?}); apex vertex host-plane \
                     dispatch was wrong for this fillet pair"
                ),
            });
        }

        let mef_result =
            brep.mef(he_dest, he_origin)
                .map_err(|e| ApexBuildError::EulerOpFailed {
                    op: "mef",
                    step: 100 + k,
                    cause: format!("FaceEdge[{k}]: {e}"),
                })?;

        let curve =
            crate::apex_construction::constructive_geometry::planned_face_edge_curve(plan, k)?;
        attach_curve_to_edge_aligned(brep, mef_result.edge, curve)?;

        if trace {
            let p = curve.degree;
            let n_cp = curve.control_points.len() - 1;
            let c0 = curve.point_at(curve.knots[p]);
            let c1 = curve.point_at(curve.knots[n_cp + 1]);
            eprintln!(
                "[direct-apex-trace §14.7/B] mef FaceEdge[{k}] T_b[{k}]={t_b_k:?} -> T_a[{next_k}]={t_a_next:?} \
                 -> edge={:?}, new face={:?}, curve endpoints=({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                mef_result.edge, mef_result.face, c0.x, c0.y, c0.z, c1.x, c1.y, c1.z
            );
        }
        face_edges.push(Some(mef_result.edge));
    }

    Ok(face_edges)
}
