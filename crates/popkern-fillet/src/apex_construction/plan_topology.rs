//! Phase 6c-DirectApex.1.d/e — cap kind decision and target topology assembly.
//!
//! `decide_cap_kind` (DirectApex.1.d): given the per-fillet terminations,
//! decide whether the apex needs a Sphere or Gregory cap, recording
//! reasoning in [`PlanTelemetry`].
//!
//! `build_target_topology` (DirectApex.1.e): from the chosen [`CapKind`] +
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
/// Sphere cap criteria (ALL must hold):
/// 1. All fillet radii are within `tolerance.linear()` of each other.
/// 2. All adjacent FaceEdge bridges collapse to within
///    `tolerance.linear() * 10` (i.e. T_b[k] ≈ T_a[k+1]).
/// 3. The trihedral support-normal frame is mutually orthogonal. A
///    non-orthogonal pyramid apex may be topologically trihedral and
///    radius-uniform, but a single rolling-ball sphere cannot satisfy
///    the cap-to-fillet support-normal postcondition along every outer
///    seam.
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

    let single_sphere_topology_supported = n_fil == 3;
    telemetry.set_value(
        "single_sphere_topology_supported",
        if single_sphere_topology_supported {
            1.0
        } else {
            0.0
        },
    );
    let max_face_normal_dot = max_unique_face_normal_dot(terminations);
    let sphere_normal_frame_supported = max_face_normal_dot <= 1.0e-6;
    telemetry.set_value("max_face_normal_dot", max_face_normal_dot);
    telemetry.set_value(
        "single_sphere_normal_frame_supported",
        if sphere_normal_frame_supported {
            1.0
        } else {
            0.0
        },
    );

    // Sphere if all criteria hold; else Gregory. A collapsed bridge
    // loop is not sufficient for N >= 4: a single sphere can be
    // topologically closed while still failing the cap-to-cylinder
    // boundary geometry on pyramid-like apices.
    if radii_uniform
        && bridges_collapse
        && single_sphere_topology_supported
        && sphere_normal_frame_supported
    {
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
        let reason = match (
            radii_uniform,
            bridges_collapse,
            single_sphere_topology_supported,
            sphere_normal_frame_supported,
        ) {
            (false, false, _, _) => format!(
                "variable radii (spread {r_spread:.3e}) and non-collapsed bridges \
                 (max gap {max_bridge_gap:.3e}) — Gregory cap chosen"
            ),
            (false, true, _, _) => {
                format!("variable radii (spread {r_spread:.3e}) — Gregory cap chosen")
            }
            (true, false, _, _) => format!(
                "uniform radii but bridges don't collapse (max gap {max_bridge_gap:.3e}) — Gregory cap chosen"
            ),
            (true, true, false, _) => format!(
                "uniform radii and collapsed bridges, but single Sphere caps are accepted only for trihedral corners; N={n_fil} — Gregory cap chosen"
            ),
            (true, true, true, false) => format!(
                "uniform trihedral radii and collapsed bridges, but support normals are non-orthogonal \
                 (max |dot| {max_face_normal_dot:.3e}) — Gregory cap chosen"
            ),
            (true, true, true, true) => unreachable!(),
        };
        telemetry.cap_kind_reason = reason;
        let _ = fillet_specs; // currently unused; reserved for variable-radius variants
        let cap_kind = gregory_cap_kind(terminations, tolerance);
        let welded = cap_kind.welded_corner_count();
        telemetry.set_value("welded_corner_count", welded as f64);
        if welded > 0 {
            telemetry.cap_kind_reason = format!(
                "{} — {welded}/{n_fil} corners welded (shared-corner class: adjacent \
                 terminations within {bridge_tol:.3e})",
                telemetry.cap_kind_reason
            );
        }
        Ok(cap_kind)
    }
}

/// Snap a clamped NURBS curve's first/last control points to the
/// given endpoints so the curve agrees bitwise with the (possibly
/// welded) corner vertices. The displacement is bounded by the weld
/// gate (≤ 10·tol / 2 per endpoint) — below every downstream
/// endpoint-match gate; the arc interior is unchanged.
fn snap_clamped_curve_endpoints(curve: &mut NurbsCurve, start: &Point3, end: &Point3) {
    if let Some(first) = curve.control_points.first_mut() {
        *first = *start;
    }
    if let Some(last) = curve.control_points.last_mut() {
        *last = *end;
    }
}

fn max_unique_face_normal_dot(terminations: &[PerFilletTermination]) -> f64 {
    let mut normals = Vec::new();
    for termination in terminations {
        push_unique_normal(&mut normals, termination.face_a_normal);
        push_unique_normal(&mut normals, termination.face_b_normal);
    }
    let mut max_dot = 0.0_f64;
    for i in 0..normals.len() {
        for j in (i + 1)..normals.len() {
            let Some(a) = normals[i].normalized() else {
                continue;
            };
            let Some(b) = normals[j].normalized() else {
                continue;
            };
            max_dot = max_dot.max(a.dot(&b).abs());
        }
    }
    max_dot
}

fn push_unique_normal(
    normals: &mut Vec<popkern_core::geometry::Vec3>,
    normal: popkern_core::geometry::Vec3,
) {
    let Some(unit) = normal.normalized() else {
        return;
    };
    if normals.iter().any(|existing| {
        existing
            .normalized()
            .is_some_and(|other| unit.dot(&other).abs() > 1.0 - 1.0e-8)
    }) {
        return;
    }
    normals.push(unit);
}

/// Build the Gregory [`CapKind`], classifying each corner between
/// adjacent terminations as Bridged or Welded (shared-corner class,
/// PK-219).
///
/// A corner welds when the gap between fillet `k`'s T_b and fillet
/// `k+1`'s T_a is within `tolerance.linear().max(1e-9) * 10` — the
/// same gate as [`decide_cap_kind`]'s bridge-collapse criterion, so
/// the welded contract takes over exactly the set the Sphere
/// dispatch measures as collapsed. The weld point is the midpoint
/// of the two terminations (exact for planar supports: both
/// endpoints lie on the shared support plane, so their convex
/// combination does too), and both adjacent termination entries are
/// SNAPPED to it bitwise so the cap boundary closes exactly at the
/// corner. The snap displacement is ≤ half the gate — below every
/// downstream endpoint gate.
pub(crate) fn gregory_cap_kind(
    terminations: &[PerFilletTermination],
    tolerance: &popkern_core::tolerance::Tolerance,
) -> CapKind {
    let n_fil = terminations.len();
    let weld_gate = tolerance.linear().max(1.0e-9) * 10.0;
    let mut per_fillet_terminations: Vec<[Point3; 2]> = terminations
        .iter()
        .map(|t| [t.endpoint_a, t.endpoint_b])
        .collect();
    let mut corner_kinds: Vec<super::cap_kind::ApexCornerKind> = Vec::with_capacity(n_fil.max(1));
    for k in 0..n_fil {
        let next_k = (k + 1) % n_fil;
        let t_b = per_fillet_terminations[k][1];
        let t_a_next = per_fillet_terminations[next_k][0];
        let gap = t_b.distance_to(&t_a_next);
        if gap <= weld_gate && n_fil >= 2 {
            let point = Point3::new(
                0.5 * (t_b.x + t_a_next.x),
                0.5 * (t_b.y + t_a_next.y),
                0.5 * (t_b.z + t_a_next.z),
            );
            per_fillet_terminations[k][1] = point;
            per_fillet_terminations[next_k][0] = point;
            corner_kinds.push(super::cap_kind::ApexCornerKind::Welded { point });
        } else {
            corner_kinds.push(super::cap_kind::ApexCornerKind::Bridged);
        }
    }
    let per_fillet_caps: Vec<(Point3, f64)> = terminations
        .iter()
        .map(|t| (t.cap_centre, t.cap_radius))
        .collect();
    CapKind::Gregory {
        per_fillet_terminations,
        per_fillet_caps,
        corner_kinds,
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
/// (DirectApex.2) computes face splicing from the apex topology +
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
        CapKind::Gregory {
            per_fillet_terminations,
            corner_kinds,
            ..
        } => {
            // 2N apex vertices: vertex layout chosen to make the
            // cap loop walk visit them in alternating order:
            //   index 2k   = T_a[k]    (FilletEnd[k] start)
            //   index 2k+1 = T_b[k]    (FilletEnd[k] end = FaceEdge[k] start)
            //
            // Cap-loop walk: T_a[0] → FilletEnd[0] → T_b[0] →
            //                FaceEdge[0] → T_a[1] → FilletEnd[1] → ...
            //
            // Corner points come from the CapKind's terminations,
            // which [`gregory_cap_kind`] has already SNAPPED at
            // welded corners: for a welded corner k, slots 2k+1 and
            // 2(k+1) hold the same (bitwise-equal) weld point and NO
            // FaceEdge is emitted — the shared-corner contract
            // (PK-219). The historical claim that Gregory bridge
            // endpoints are always non-degenerate was false for
            // uniform pyramids (Sphere dispatch also requires
            // trihedral + orthogonal normals), which is exactly how
            // the PK-216 zero-length-edge corruption entered.
            let mut apex_vertices: Vec<Point3> = Vec::with_capacity(2 * n_fil);
            for t in per_fillet_terminations {
                apex_vertices.push(t[0]);
                apex_vertices.push(t[1]);
            }
            let mut apex_edges: Vec<TargetEdge> = Vec::with_capacity(2 * n_fil);
            let mut outer_loop_edge_indices: Vec<usize> = Vec::with_capacity(2 * n_fil);
            for (k, term) in terminations.iter().enumerate() {
                let next_k = (k + 1) % n_fil;
                let v_ta_k = 2 * k;
                let v_tb_k = 2 * k + 1;
                let v_ta_next = 2 * next_k;
                // FilletEnd[k] — CapVsCyl edge from T_a[k] to T_b[k],
                // with endpoints snapped to the (possibly welded)
                // corner points so curve and vertex agree bitwise.
                let mut curve = term.curve_3d.clone();
                snap_clamped_curve_endpoints(
                    &mut curve,
                    &apex_vertices[v_ta_k],
                    &apex_vertices[v_tb_k],
                );
                apex_edges.push(TargetEdge {
                    origin_vertex_idx: v_ta_k,
                    dest_vertex_idx: v_tb_k,
                    curve,
                    role: TargetEdgeRole::CapVsCyl {
                        fillet_idx: term.fillet_idx,
                    },
                });
                outer_loop_edge_indices.push(apex_edges.len() - 1);
                if corner_kinds.get(k).is_some_and(|kind| kind.is_welded()) {
                    // Welded corner: T_b[k] ≡ T_a[k+1]; no FaceEdge
                    // exists between them.
                    continue;
                }
                // FaceEdge[k] — CapVsPlane straight-line bridge from
                // T_b[k] to T_a[(k+1) % N] (bridged corners only).
                let bridge_start = apex_vertices[v_tb_k];
                let bridge_end = apex_vertices[v_ta_next];
                let bridge =
                    fit_nurbs_through_points(&[bridge_start, bridge_end]).unwrap_or_else(|_| {
                        NurbsCurve {
                            degree: 1,
                            knots: vec![0.0, 0.0, 1.0, 1.0],
                            control_points: vec![bridge_start, bridge_end],
                            weights: vec![1.0, 1.0],
                        }
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
                outer_loop_edge_indices.push(apex_edges.len() - 1);
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
