//! Boundary and sector helpers for DirectApex Gregory cap assembly.

use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::error::ApexBuildError;
use super::gregory_cap::GregoryCapSector;
use super::plan::ApexPlan;
use super::target_topology::{TargetCapFace, TargetEdgeRole};
use crate::blend::boundary::{BlendBoundary, BoundarySegment};

pub(super) fn validate_physical_sector_handles(
    sectors: &[GregoryCapSector],
    n: usize,
    apex_edge_handles: &[Option<EdgeHandle>],
) -> Result<(), ApexBuildError> {
    if sectors.len() != n {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "Gregory assembly expected {n} physical sectors, got {}",
                sectors.len()
            ),
        });
    }
    for sector in sectors {
        let fillet_slot = 2 * sector.sector_index;
        let face_slot = fillet_slot + 1;
        let expected_fillet = apex_edge_handles.get(fillet_slot).copied().flatten();
        let expected_face = apex_edge_handles.get(face_slot).copied().flatten();
        if expected_fillet != Some(sector.fillet_end_edge) || expected_face != sector.face_edge {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "Gregory sector {} does not match physical apex-edge handles: \
                     sector FilletEnd={:?} FaceEdge={:?}, apex slots FilletEnd={:?} FaceEdge={:?}",
                    sector.sector_index,
                    sector.fillet_end_edge,
                    sector.face_edge,
                    expected_fillet,
                    expected_face
                ),
            });
        }
    }
    Ok(())
}

pub(super) fn validate_planned_gregory_loop_roles(
    plan: &ApexPlan,
    outer_indices: &[usize],
    n: usize,
) -> Result<(), ApexBuildError> {
    for (k, &edge_idx) in outer_indices.iter().enumerate() {
        let target_edge = plan
            .target_topology
            .apex_edges
            .get(edge_idx)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "outer_loop_edge_indices[{k}]={edge_idx} out of range for \
                     apex_edges (len {})",
                    plan.target_topology.apex_edges.len(),
                ),
            })?;
        match target_edge.role {
            TargetEdgeRole::CapVsCyl { fillet_idx } if k % 2 == 0 && fillet_idx < n => {}
            TargetEdgeRole::CapVsCyl { fillet_idx } => {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "Gregory cap loop has CapVsCyl fillet {fillet_idx} at invalid \
                         outer-loop position {k} for N={n}"
                    ),
                });
            }
            TargetEdgeRole::CapVsPlane { .. } if k % 2 == 1 => {}
            TargetEdgeRole::CapVsPlane { .. } => {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "Gregory cap loop expects FaceEdge at odd index; found CapVsPlane \
                         at outer-loop position {k}"
                    ),
                });
            }
            TargetEdgeRole::InnerSpoke => {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "Gregory cap outer loop should not contain InnerSpoke at index {k}"
                    ),
                });
            }
        }
    }
    Ok(())
}

pub(super) fn oriented_cap_edge_curve(
    brep: &BRep,
    cap_face: FaceHandle,
    edge: EdgeHandle,
    role: &str,
) -> Result<NurbsCurve, ApexBuildError> {
    let edge_data = brep
        .edges
        .get(edge)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!("Gregory {role} edge {edge:?} is not present in the BRep"),
        })?;
    let curve_on_edge = edge_data
        .curve
        .as_ref()
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!("Gregory {role} edge {edge:?} has no physical curve"),
        })?;
    let cap_half_edge = half_edge_for_face_on_edge(brep, edge, cap_face).ok_or_else(|| {
        ApexBuildError::PlanMismatch {
            reason: format!("Gregory {role} edge {edge:?} is not bounded by cap face {cap_face:?}"),
        }
    })?;
    let curve = curve_on_edge
        .geometry
        .to_nurbs_in_range((curve_on_edge.t_start, curve_on_edge.t_end));
    if cap_half_edge == edge_data.half_edge {
        Ok(curve)
    } else {
        Ok(reverse_nurbs_curve(&curve))
    }
}

pub(super) fn half_edge_for_face_on_edge(
    brep: &BRep,
    edge: EdgeHandle,
    face: FaceHandle,
) -> Option<HalfEdgeHandle> {
    let edge_data = brep.edges.get(edge)?;
    let primary = edge_data.half_edge;
    let primary_face = brep
        .half_edges
        .get(primary)
        .and_then(|he| brep.loops.get(he.loop_))
        .map(|l| l.face);
    if primary_face == Some(face) {
        return Some(primary);
    }
    let twin = brep.half_edges.get(primary)?.twin;
    let twin_face = brep
        .half_edges
        .get(twin)
        .and_then(|he| brep.loops.get(he.loop_))
        .map(|l| l.face);
    (twin_face == Some(face)).then_some(twin)
}

fn reverse_nurbs_curve(curve: &NurbsCurve) -> NurbsCurve {
    let first = *curve.knots.first().unwrap_or(&0.0);
    let last = *curve.knots.last().unwrap_or(&1.0);
    let mut knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|k| first + last - *k)
        .collect();
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut control_points = curve.control_points.clone();
    control_points.reverse();
    let mut weights = curve.weights.clone();
    weights.reverse();
    NurbsCurve {
        degree: curve.degree,
        knots,
        control_points,
        weights,
    }
}

pub(super) fn missing_adjacent_face_error(
    brep: &BRep,
    cap_face: FaceHandle,
    edge: EdgeHandle,
    role: &str,
    fillet_index: usize,
) -> ApexBuildError {
    let edge_alive = brep.edges.contains_key(edge);
    let adj = if edge_alive {
        let primary = brep.edges[edge].half_edge;
        let pa = brep
            .half_edges
            .get(primary)
            .and_then(|he| brep.loops.get(he.loop_))
            .map(|l| l.face);
        let twin = brep.half_edges.get(primary).map(|he| he.twin);
        let pb = twin
            .and_then(|t| brep.half_edges.get(t))
            .and_then(|he| brep.loops.get(he.loop_))
            .map(|l| l.face);
        format!("primary_face={pa:?}, twin_face={pb:?}, cap_face={cap_face:?}")
    } else {
        "edge no longer exists in BRep".to_string()
    };
    ApexBuildError::PlanMismatch {
        reason: format!(
            "Gregory {role}[{fillet_index}] BRep edge {edge:?} has no twin face on \
             the non-cap side ({adj})"
        ),
    }
}

pub(super) fn inject_physical_hex_boundary_cps(
    boundary: &mut BlendBoundary,
) -> Result<(), ApexBuildError> {
    let n = boundary.segments.len() / 2;
    let mut hex_boundary_cps_deg9 = Vec::with_capacity(n);
    for k in 0..n {
        let row0_curve = match &boundary.segments[2 * k] {
            BoundarySegment::FilletEnd { curve, .. } => curve,
            _ => {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "Gregory physical boundary injection expected FilletEnd at segment {}",
                        2 * k
                    ),
                });
            }
        };
        let col9_curve = match &boundary.segments[2 * k + 1] {
            BoundarySegment::FaceEdge { curve, .. } => curve,
            _ => {
                return Err(ApexBuildError::PlanMismatch {
                    reason: format!(
                        "Gregory physical boundary injection expected FaceEdge at segment {}",
                        2 * k + 1
                    ),
                });
            }
        };
        let row0 = sample_and_fit_deg9(row0_curve).ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "failed to fit degree-9 physical FilletEnd boundary for Gregory sector {k}"
            ),
        })?;
        let col9 = sample_and_fit_deg9(col9_curve).ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "failed to fit degree-9 physical FaceEdge boundary for Gregory sector {k}"
            ),
        })?;
        hex_boundary_cps_deg9.push((row0, col9));
    }
    boundary.hex_boundary_cps_deg9 = hex_boundary_cps_deg9;
    Ok(())
}

fn sample_and_fit_deg9(curve: &NurbsCurve) -> Option<[Point3; 10]> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::fitting::fit_curve;

    let (t0, t1) = curve.parameter_range();
    if (t1 - t0).abs() < 1e-15 {
        return None;
    }
    let n_samples = 40usize;
    let mut samples = Vec::with_capacity(n_samples);
    let mut params = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let s = i as f64 / (n_samples - 1) as f64;
        let t = t0 + s * (t1 - t0);
        samples.push(curve.point_at(t));
        params.push(s);
    }
    let fitted = fit_curve(&samples, &params, 9, 10).ok()?;
    if fitted.control_points.len() != 10 {
        return None;
    }
    let mut out = [Point3::new(0.0, 0.0, 0.0); 10];
    for (i, cp) in fitted.control_points.iter().enumerate() {
        out[i] = *cp;
    }
    Some(out)
}

pub(super) fn gregory_sectors_from_plan_and_apex_edge_handles(
    plan: &ApexPlan,
    apex_edge_handles: &[Option<EdgeHandle>],
) -> Result<Vec<GregoryCapSector>, ApexBuildError> {
    let outer_indices = match &plan.target_topology.cap_face {
        TargetCapFace::Gregory {
            outer_loop_edge_indices,
        } => outer_loop_edge_indices,
        other => {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "Gregory sector derivation needs Gregory target topology, got {other:?}"
                ),
            });
        }
    };
    let _ = outer_indices;
    let n = plan.cap_kind.n_fillets();
    if apex_edge_handles.len() < 2 * n {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "Gregory sector derivation expected at least {} apex-edge slots, got {}",
                2 * n,
                apex_edge_handles.len()
            ),
        });
    }

    let mut sectors = Vec::with_capacity(n);
    for sector_index in 0..n {
        let fillet_slot = 2 * sector_index;
        let face_slot = fillet_slot + 1;
        let fillet_plan_index = plan
            .target_topology
            .fillet_end_edge_index(sector_index)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!("Gregory sector {sector_index} has no planned CapVsCyl apex_edge"),
            })?;
        let fillet_index = gregory_loop_fillet_index(plan, fillet_plan_index, n)?;
        let fillet_end_edge = apex_edge_handles
            .get(fillet_slot)
            .copied()
            .flatten()
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "Gregory sector {sector_index} missing physical FilletEnd edge \
                         at apex-edge slot {fillet_slot}"
                ),
            })?;
        let face_edge = apex_edge_handles.get(face_slot).copied().flatten();
        if face_edge.is_none() && !plan.cap_kind.corner_is_welded(sector_index) {
            return Err(ApexBuildError::PlanMismatch {
                reason: format!(
                    "Gregory sector {sector_index} missing physical FaceEdge edge \
                         at apex-edge slot {face_slot} (corner not declared welded)"
                ),
            });
        }
        sectors.push(GregoryCapSector {
            sector_index,
            fillet_index,
            fillet_end_edge,
            face_edge,
        });
    }
    Ok(sectors)
}

fn gregory_loop_fillet_index(
    plan: &ApexPlan,
    apex_edge_index: usize,
    n: usize,
) -> Result<usize, ApexBuildError> {
    let edge = plan
        .target_topology
        .apex_edges
        .get(apex_edge_index)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "Gregory sector references apex_edge[{apex_edge_index}], but apex_edges len is {}",
                plan.target_topology.apex_edges.len()
            ),
        })?;
    match edge.role {
        TargetEdgeRole::CapVsCyl { fillet_idx } if fillet_idx < n => Ok(fillet_idx),
        TargetEdgeRole::CapVsCyl { fillet_idx } => Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "Gregory sector references out-of-range fillet index {fillet_idx} for N={n}"
            ),
        }),
        ref other => Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "Gregory sector expected CapVsCyl apex edge at index {apex_edge_index}, got {other:?}"
            ),
        }),
    }
}

pub(super) fn find_face_across_edge(
    brep: &BRep,
    edge: EdgeHandle,
    cap_face: FaceHandle,
) -> Option<FaceHandle> {
    let edge_data = brep.edges.get(edge)?;
    let primary = edge_data.half_edge;
    let primary_he = brep.half_edges.get(primary)?;
    let twin_he = brep.half_edges.get(primary_he.twin)?;
    let face_a = brep.loops.get(primary_he.loop_)?.face;
    let face_b = brep.loops.get(twin_he.loop_)?.face;
    if face_a == cap_face {
        Some(face_b)
    } else if face_b == cap_face {
        Some(face_a)
    } else {
        None
    }
}
