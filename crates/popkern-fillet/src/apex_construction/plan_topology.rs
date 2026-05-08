//! Phase 6c-A2.1.d/e — cap kind decision and target topology assembly.
//!
//! `decide_cap_kind` (A2.1.d): given the per-fillet terminations,
//! decide whether the apex needs a Sphere or Gregory cap, recording
//! reasoning in [`PlanTelemetry`].
//!
//! `build_target_topology` (A2.1.e): from the chosen [`CapKind`] +
//! terminations, build the declarative [`TargetTopology`] describing
//! the desired final apex region.

use popkern_core::geometry::{NurbsCurve, Point3};

use crate::center_curve::fit_nurbs_through_points;
use crate::types::FilletSpec;

use super::cap_kind::CapKind;
use super::error::ApexPlanError;
use super::plan_terminations::PerFilletTermination;
use super::target_topology::{TargetCapFace, TargetEdge, TargetEdgeRole, TargetTopology};
use super::telemetry::PlanTelemetry;

/// Decide whether the apex needs a Sphere cap (uniform/symmetric)
/// or an N-sided Gregory cap (variable radii / distinct
/// terminations).  Per ADR §4.1: the decision is tolerance-aware.
///
/// Sphere cap criteria (BOTH must hold):
/// 1. All fillet radii are within `tolerance.linear()` of each other.
/// 2. All adjacent FaceEdge bridges collapse to within
///    `tolerance.linear() * 10` (i.e. T_b[k] ≈ T_a[k+1]).
///
/// Otherwise → Gregory cap with distinct per-fillet terminations.
///
/// Records reasoning + diagnostic numbers in `telemetry`.
pub(crate) fn decide_cap_kind(
    terminations: &[PerFilletTermination],
    fillet_specs: &[FilletSpec],
    tolerance: &popkern_core::tolerance::Tolerance,
    telemetry: &mut PlanTelemetry,
) -> Result<CapKind, ApexPlanError> {
    let n_fil = terminations.len();
    let tol = tolerance.linear();

    // Criterion 1: radii uniformity.
    let radii: Vec<f64> = terminations.iter().map(|t| t.cap_radius).collect();
    let r_min = radii.iter().copied().fold(f64::INFINITY, f64::min);
    let r_max = radii.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let r_spread = r_max - r_min;
    let r_avg = radii.iter().sum::<f64>() / (n_fil as f64);
    let radii_uniform = r_spread <= tol;
    telemetry.set_value("radii_min", r_min);
    telemetry.set_value("radii_max", r_max);
    telemetry.set_value("radii_spread", r_spread);

    // Criterion 2: bridge-collapse check.
    //
    // Adjacent FaceEdge bridge: between T_b[k] (FilletEnd[k] end) and
    // T_a[(k+1) % n_fil] (FilletEnd[k+1] start).  Bridge collapses
    // when these two points coincide (within tol * 10).
    //
    // Gates an extra factor of 10 because Sphere caps tolerate small
    // numerical noise that wouldn't disqualify the symmetric pattern;
    // commercial kernels use a similar relaxed "convergence" tol.
    let bridge_tol = tol * 10.0;
    let mut max_bridge_gap = 0.0_f64;
    for k in 0..n_fil {
        let next_k = (k + 1) % n_fil;
        let gap = terminations[k]
            .endpoint_b
            .distance_to(&terminations[next_k].endpoint_a);
        if gap > max_bridge_gap {
            max_bridge_gap = gap;
        }
    }
    let bridges_collapse = max_bridge_gap < bridge_tol;
    telemetry.set_value("max_bridge_collapse", max_bridge_gap);
    telemetry.set_value("bridge_collapse_threshold", bridge_tol);

    // Sphere if BOTH criteria hold; else Gregory.
    if radii_uniform && bridges_collapse {
        // Compute average cap centre (per-fillet centres should
        // coincide for the symmetric case).
        let mut centre_sum = Point3::new(0.0, 0.0, 0.0);
        for t in terminations {
            centre_sum.x += t.cap_centre.x;
            centre_sum.y += t.cap_centre.y;
            centre_sum.z += t.cap_centre.z;
        }
        let n_f = n_fil as f64;
        let center = Point3::new(centre_sum.x / n_f, centre_sum.y / n_f, centre_sum.z / n_f);
        // Tangent arc endpoints: in the collapsed-bridge case, every
        // T_b[k] equals the adjacent T_a[k+1], so we collect either
        // and the result has n_fil unique points.  We use T_a[k] for
        // each k.
        let tangent_arc_endpoints: Vec<Point3> =
            terminations.iter().map(|t| t.endpoint_a).collect();
        telemetry.cap_kind_reason = format!(
            "uniform radii within {tol:.3e} (spread {r_spread:.3e}) + bridges collapse to \
             {max_bridge_gap:.3e} < {bridge_tol:.3e} — Sphere cap chosen"
        );
        Ok(CapKind::Sphere {
            center,
            radius: r_avg,
            tangent_arc_endpoints,
        })
    } else {
        // Diagnostic: which criterion failed?
        let reason = match (radii_uniform, bridges_collapse) {
            (false, false) => format!(
                "variable radii (spread {r_spread:.3e}) and non-collapsed bridges \
                 (max gap {max_bridge_gap:.3e}) — Gregory cap chosen"
            ),
            (false, true) => {
                format!("variable radii (spread {r_spread:.3e}) — Gregory cap chosen")
            }
            (true, false) => format!(
                "uniform radii but bridges don't collapse (max gap {max_bridge_gap:.3e}) — Gregory cap chosen"
            ),
            (true, true) => unreachable!(),
        };
        telemetry.cap_kind_reason = reason;
        let _ = fillet_specs; // currently unused; reserved for variable-radius variants
        Ok(gregory_cap_kind(terminations))
    }
}

pub(crate) fn gregory_cap_kind(terminations: &[PerFilletTermination]) -> CapKind {
    let per_fillet_terminations: Vec<[Point3; 2]> = terminations
        .iter()
        .map(|t| [t.endpoint_a, t.endpoint_b])
        .collect();
    let per_fillet_caps: Vec<(Point3, f64)> = terminations
        .iter()
        .map(|t| (t.cap_centre, t.cap_radius))
        .collect();
    CapKind::Gregory {
        per_fillet_terminations,
        per_fillet_caps,
    }
}

/// Build the declarative [`TargetTopology`] from the chosen
/// [`CapKind`] and the per-fillet terminations.
///
/// For **Sphere** caps:
/// * `apex_vertices`: N points, one per fillet's tangent arc endpoint
///   (T_a[k] for each k; bridges collapse so these equal T_b[k-1]).
/// * `apex_edges`: N edges, each a FilletEnd arc on a cap-vs-Cyl
///   boundary.  No FaceEdge bridges (collapsed).
/// * `cap_face`: `Sphere` with center + radius + outer-loop indices.
///
/// For **Gregory** caps:
/// * `apex_vertices`: 2N points (T_a[k] and T_b[k] per fillet).
/// * `apex_edges`: 2N edges alternating CapVsCyl (FilletEnd[k]) and
///   CapVsPlane (FaceEdge[k]).
/// * `cap_face`: `Gregory` with 2N outer-loop edge indices.
///
/// `face_boundary_updates` is left empty here — the builder
/// (A2.2) computes face splicing from the apex topology +
/// existing BRep state.  See ADR §3.3.
pub(crate) fn build_target_topology(
    cap_kind: &CapKind,
    terminations: &[PerFilletTermination],
) -> TargetTopology {
    let n_fil = terminations.len();
    match cap_kind {
        CapKind::Sphere {
            center,
            radius,
            tangent_arc_endpoints,
        } => {
            let apex_vertices: Vec<Point3> = tangent_arc_endpoints.clone();
            let mut apex_edges: Vec<TargetEdge> = Vec::with_capacity(n_fil);
            let mut outer_loop_edge_indices: Vec<usize> = Vec::with_capacity(n_fil);
            for (k, term) in terminations.iter().enumerate() {
                let next_k = (k + 1) % n_fil;
                apex_edges.push(TargetEdge {
                    origin_vertex_idx: k,
                    dest_vertex_idx: next_k,
                    curve: term.curve_3d.clone(),
                    role: TargetEdgeRole::CapVsCyl {
                        fillet_idx: term.fillet_idx,
                    },
                });
                outer_loop_edge_indices.push(k);
            }
            let cap_face = TargetCapFace::Sphere {
                center: *center,
                radius: *radius,
                outer_loop_edge_indices,
            };
            TargetTopology {
                apex_vertices,
                apex_edges,
                cap_face,
                face_boundary_updates: Vec::new(),
            }
        }
        CapKind::Open { .. } => unreachable!(
            "build_target_topology is never called with CapKind::Open — \
             plan_open_rung constructs the (empty) topology directly",
        ),
        CapKind::Gregory { .. } => {
            // 2N apex vertices: vertex layout chosen to make the
            // cap loop walk visit them in alternating order:
            //   index 2k   = T_a[k]    (FilletEnd[k] start)
            //   index 2k+1 = T_b[k]    (FilletEnd[k] end = FaceEdge[k] start)
            //
            // Cap-loop walk: T_a[0] → FilletEnd[0] → T_b[0] →
            //                FaceEdge[0] → T_a[1] → FilletEnd[1] → ...
            let mut apex_vertices: Vec<Point3> = Vec::with_capacity(2 * n_fil);
            for t in terminations {
                apex_vertices.push(t.endpoint_a);
                apex_vertices.push(t.endpoint_b);
            }
            let mut apex_edges: Vec<TargetEdge> = Vec::with_capacity(2 * n_fil);
            let mut outer_loop_edge_indices: Vec<usize> = Vec::with_capacity(2 * n_fil);
            for (k, term) in terminations.iter().enumerate() {
                let next_k = (k + 1) % n_fil;
                let v_ta_k = 2 * k;
                let v_tb_k = 2 * k + 1;
                let v_ta_next = 2 * next_k;
                // FilletEnd[k] — CapVsCyl edge from T_a[k] to T_b[k].
                apex_edges.push(TargetEdge {
                    origin_vertex_idx: v_ta_k,
                    dest_vertex_idx: v_tb_k,
                    curve: term.curve_3d.clone(),
                    role: TargetEdgeRole::CapVsCyl {
                        fillet_idx: term.fillet_idx,
                    },
                });
                outer_loop_edge_indices.push(2 * k);
                // FaceEdge[k] — CapVsPlane straight-line bridge from
                // T_b[k] to T_a[(k+1) % N].
                // Bridge between adjacent FilletEnd terminations.
                // For Gregory caps, the two endpoints are non-degenerate
                // (the dispatch check in decide_cap_kind would have
                // picked Sphere if they collapsed), so the degenerate
                // path in fit_nurbs_through_points is unreachable.
                let bridge =
                    fit_nurbs_through_points(&[term.endpoint_b, terminations[next_k].endpoint_a])
                        .unwrap_or_else(|_| NurbsCurve {
                            degree: 1,
                            knots: vec![0.0, 0.0, 1.0, 1.0],
                            control_points: vec![term.endpoint_b, terminations[next_k].endpoint_a],
                            weights: vec![1.0, 1.0],
                        });
                apex_edges.push(TargetEdge {
                    origin_vertex_idx: v_tb_k,
                    dest_vertex_idx: v_ta_next,
                    curve: bridge,
                    role: TargetEdgeRole::CapVsPlane {
                        // Shared-face index = which face this bridge
                        // sits on.  Builder maps to actual FaceHandle
                        // via `boundary.segments[2k+1].face`.
                        shared_face_idx: k,
                    },
                });
                outer_loop_edge_indices.push(2 * k + 1);
            }
            let cap_face = TargetCapFace::Gregory {
                outer_loop_edge_indices,
            };
            TargetTopology {
                apex_vertices,
                apex_edges,
                cap_face,
                face_boundary_updates: Vec::new(),
            }
        }
    }
}
