//! Phase 6c-A2.3.b — A2 apex-construction dispatch helper.
//!
//! Bridges the legacy `apply_blend_geometry` apex-vertex-blend
//! loop to the A2 planner + builder
//! ([`crate::apex_construction::plan_apex_resolution`] +
//! [`crate::apex_construction::construct_apex`]).
//!
//! Called from `apply_blend_geometry` only when
//! `FilletOptions::use_a2_apex_construction == true`.  Returns
//! `Ok(A2DispatchOutcome { built_faces, corner_patches })` on success (caller
//! skips the legacy `insert_vertex_blend` path, merges `built_faces` into the
//! post-loop protected-face set, and threads patch reports to diagnostics).
//! Returns
//! `Err(reason)` when A2 can't finish; caller emits
//! [`BlendWarning::A2FellBackToLegacy`](crate::error::BlendWarning::A2FellBackToLegacy)
//! and runs legacy as a backup.

use std::collections::HashMap;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use super::EdgeGeometry;

use crate::apex_construction::cap_kind::CapKind;
use crate::apex_construction::construct::construct_apex_with_surfaces;
use crate::apex_construction::constructive_builder::construct_apex_constructively;
use crate::apex_construction::plan::plan_apex_resolution_lenient;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{CornerPatchInfo, FilletRadius, FilletSpec};

/// What `try_a2_apex_construction` produced on success.
///
/// `built_faces` lists the cap face(s) `construct_apex` authored
/// (`ApexBuildResult::created_faces`).  The caller threads this
/// list into the protected-face set passed to `sliver_absorb` /
/// `cap_insert` / `late_wire_jobs` so those post-loop stages
/// don't re-author or absorb the A2-built cap.
///
/// `corner_patches` carries any cap patch reports produced by
/// the apex builder so continuity diagnostics and validation see
/// native A2 cap faces.
#[derive(Debug, Clone, Default)]
pub(crate) struct A2DispatchOutcome {
    pub built_faces: Vec<FaceHandle>,
    pub corner_patches: Vec<CornerPatchInfo>,
}

/// Try the A2 apex-construction path for one vertex.
///
/// On success, the BRep has been mutated (apex region built by
/// `construct_apex`) and `Ok(A2DispatchOutcome { built_faces })`
/// is returned.  On any A2 error (plan or build), the BRep is
/// rolled back by `construct_apex`'s own snapshot logic and
/// `Err(reason)` is returned describing why A2 couldn't finish.
pub(super) fn try_a2_apex_construction(
    brep: &mut BRep,
    vertex: VertexHandle,
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    _descendants: &HashMap<EdgeHandle, Vec<EdgeHandle>>,
    tolerance: &Tolerance,
) -> Result<A2DispatchOutcome, String> {
    // Edge handles in `edge_geometries` are the ORIGINAL pre-
    // pipeline handles.  Phase C kef's the fillet edges so
    // those handles are no longer live in the BRep — but the
    // planner's lenient entry point (`plan_apex_resolution_
    // lenient`) skips the edge-existence and edge-incidence
    // checks that would otherwise fire.  The planner reads
    // geometry from `fillet_surfaces`, not from `spec.edge`,
    // so the stale handle is harmless metadata.
    let specs: Vec<FilletSpec> = edge_indices
        .iter()
        .map(|&i| {
            let eg = &edge_geometries[i];
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
        })
        .collect();
    let surfaces: Vec<FilletSurfaceResult> = edge_indices
        .iter()
        .map(|&i| {
            let eg = &edge_geometries[i];
            // EdgeGeometry's continuity_downgrade is
            // Option<(ContinuityLevel, String)>; the surface
            // result's is Option<String> (just the reason).
            // Map by extracting the second element.
            let downgrade_reason = eg.continuity_downgrade.as_ref().map(|(_, s)| s.clone());
            FilletSurfaceResult {
                surface: eg.surface.clone(),
                trim_curve_1: eg.trim_3d_1.clone(),
                trim_curve_2: eg.trim_3d_2.clone(),
                continuity_downgrade: downgrade_reason,
            }
        })
        .collect();

    if crate::a2_trace::enabled() {
        eprintln!(
            "[a2-trace] try_a2_apex_construction: vertex={vertex:?}; \
             pre-call BRep has {} faces, {} edges, {} vertices",
            brep.faces.len(),
            brep.edges.len(),
            brep.vertices.len(),
        );
        for (fh, face) in &brep.faces {
            use popkern_core::geometry::surface::SurfaceGeometry;
            let kind = face
                .surface
                .as_ref()
                .map(|s| match &s.geometry {
                    SurfaceGeometry::Plane { origin, normal, .. } => {
                        format!("Plane(origin={origin:?}, normal={normal:?})")
                    }
                    SurfaceGeometry::Cylinder {
                        origin,
                        axis,
                        radius,
                        ..
                    } => format!("Cylinder(origin={origin:?}, axis={axis:?}, r={radius:.3e})"),
                    SurfaceGeometry::Sphere { center, radius } => {
                        format!("Sphere(center={center:?}, r={radius:.3e})")
                    }
                    SurfaceGeometry::Torus { .. } => "Torus".to_string(),
                    SurfaceGeometry::Cone { .. } => "Cone".to_string(),
                    SurfaceGeometry::Nurbs(_) => "Nurbs".to_string(),
                    SurfaceGeometry::Gregory(_) => "Gregory".to_string(),
                })
                .unwrap_or_else(|| "None".to_string());
            eprintln!("[a2-trace]   {fh:?}: {kind}");
        }
        // Dump every cylinder face's outer-loop walk so we can
        // see which edges bound which cylinder, and whether two
        // cylinder face handles share a seam edge.
        for (fh, face) in &brep.faces {
            let is_cyl = face.surface.as_ref().is_some_and(|s| {
                matches!(
                    s.geometry,
                    popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. }
                )
            });
            if !is_cyl {
                continue;
            }
            let outer_loop = face.outer_loop;
            let Some(loop_data) = brep.loops.get(outer_loop) else {
                continue;
            };
            let start = loop_data.half_edge;
            let mut cursor = start;
            let mut chain: Vec<String> = Vec::new();
            for _ in 0..(brep.half_edges.len() + 1) {
                let Some(he) = brep.half_edges.get(cursor) else {
                    break;
                };
                let origin_pos = brep.vertices[he.origin].position;
                chain.push(format!(
                    "v={:?}@({:.3},{:.3},{:.3})",
                    he.origin, origin_pos.x, origin_pos.y, origin_pos.z
                ));
                cursor = he.next;
                if cursor == start {
                    break;
                }
            }
            eprintln!(
                "[a2-trace]   {fh:?} outer loop ({} half-edges): {}",
                chain.len(),
                chain.join(" -> ")
            );
        }
    }
    let plan = plan_apex_resolution_lenient(brep, vertex, &specs, &surfaces, tolerance, None)
        .map_err(|plan_err| format!("plan_apex_resolution: {plan_err}"))?;

    // §A2.4 option 2 diagnostic: detect Phase-C cylinder-face
    // incoherence (loops with off-cyl vertices, per ADR §25).
    // Today this is purely informational under
    // `POPKERN_A2_TRACE=1` — it fires on both Sphere fixtures
    // (where A2 still succeeds because the planner finds shared
    // faces despite off-cyl loop vertices) and Gregory fixtures
    // (where A2 will fail downstream).  Future slices will use
    // this detection to drive the actual `kef`/`split_edge`
    // surgery that rebuilds coherent cylinder loops.
    if crate::a2_trace::enabled() {
        let incoherent =
            crate::apex_construction::cylinder_face_coherence::find_off_cyl_loop_vertices(brep);
        if !incoherent.is_empty() {
            eprintln!(
                "[a2-trace] cylinder-face coherence: {} face(s) with off-cyl vertices — {}",
                incoherent.len(),
                crate::apex_construction::cylinder_face_coherence::describe_incoherence(
                    &incoherent
                ),
            );
            for entry in &incoherent {
                for (vh, pos, dist) in &entry.off_cyl_vertices {
                    eprintln!(
                        "[a2-trace]   {:?} (r={:.4}) v={vh:?} pos=({:.4},{:.4},{:.4}) radial_dist={:.4}",
                        entry.face, entry.radius, pos.x, pos.y, pos.z, dist,
                    );
                }
            }
        }
    }

    // Stage 2/3/4/5 dispatch: Sphere AND Gregory both try the constructive
    // builder first; on PlanMismatch (e.g. pyramid Sphere with sphere-fit
    // T-positions that don't coincide with existing trim-trim corners,
    // any Gregory edge case the constructive path doesn't yet handle),
    // fall through to legacy A2 (§28.C-relocate path) so we don't regress
    // fixtures that were succeeding under legacy. Open is unchanged
    // (legacy short-circuits it as no-op).
    let build_result = match plan.cap_kind {
        CapKind::Sphere { .. } | CapKind::Gregory { .. } => {
            try_constructive_with_legacy_fallthrough(brep, &plan, &surfaces, tolerance)?
        }
        CapKind::Open { .. } => {
            construct_apex_with_surfaces(brep, &plan, &surfaces, tolerance, None)
                .map_err(|build_err| format!("construct_apex: {build_err}"))?
        }
    };
    Ok(A2DispatchOutcome {
        built_faces: build_result.created_faces,
        corner_patches: build_result.corner_patches,
    })
}

/// Try the constructive builder for one apex; on PlanMismatch, fall through
/// to the legacy `construct_apex_with_surfaces`. Either result is returned
/// as `Ok(ApexBuildResult)`. If both fail, the legacy error is propagated
/// upward via `Err(String)`.
///
/// The constructive builder's snapshot-rollback ensures that on its own
/// failure, the BRep is byte-equivalent to the pre-call state — so the
/// subsequent legacy call sees the same input as if the constructive path
/// were never attempted.
fn try_constructive_with_legacy_fallthrough(
    brep: &mut BRep,
    plan: &crate::apex_construction::plan::ApexPlan,
    surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
) -> Result<crate::apex_construction::construct::ApexBuildResult, String> {
    let snapshot = brep.clone();
    let execution_options =
        crate::apex_construction::constructive_builder_options::ConstructiveApexExecutionOptions::from_environment_adapter();
    match construct_apex_constructively(brep, plan, surfaces, tolerance, None, execution_options) {
        Ok(result) => {
            if matches!(plan.cap_kind, CapKind::Gregory { .. }) {
                let constructive_err =
                    if !native_gregory_has_expected_cylinder_seams(brep, &result, plan.n_fillets())
                    {
                        Some(format!(
                            "native Gregory cap has fewer than {} cap-to-cylinder seams; \
                             treating as incomplete A2 topology",
                            plan.n_fillets()
                        ))
                    } else {
                        native_gregory_path_a_error(brep)
                            .map(|reason| {
                                format!(
                                    "native Gregory cap failed structural acceptance: {reason}; \
                                 treating as incomplete A2 topology"
                                )
                            })
                            .or_else(|| {
                                native_gregory_cap_outer_g0_error(brep, &result, 5.0e-3).map(
                                |reason| {
                                    format!(
                                        "native Gregory cap failed geometric acceptance: {reason}; \
                                         treating as incomplete A2 topology"
                                    )
                                },
                            )
                            })
                    };
                if let Some(constructive_err) = constructive_err {
                    if crate::a2_trace::enabled() {
                        eprintln!(
                            "[a2-trace] constructive fell through to legacy A2: {constructive_err}"
                        );
                    }
                    *brep = snapshot;
                    return construct_apex_with_surfaces(brep, plan, surfaces, tolerance, None)
                        .map_err(|build_err| {
                            format!(
                                "construct_apex (after constructive fail: {constructive_err}): \
                                 {build_err}"
                            )
                        });
                }
            }
            Ok(result)
        }
        Err(constructive_err) => {
            if crate::a2_trace::enabled() {
                eprintln!("[a2-trace] constructive fell through to legacy A2: {constructive_err}");
            }
            construct_apex_with_surfaces(brep, plan, surfaces, tolerance, None).map_err(
                |build_err| {
                    format!(
                        "construct_apex (after constructive fail: {constructive_err}): \
                         {build_err}"
                    )
                },
            )
        }
    }
}

fn native_gregory_path_a_error(brep: &BRep) -> Option<String> {
    use popkern_core::error::ValidationError;

    brep.validate_all().into_iter().find_map(|err| match err {
        ValidationError::BrokenCycle { .. } | ValidationError::ShellGap { .. } => {
            Some(format!("{err:?}"))
        }
        ValidationError::OpenLoop { loop_, .. } => {
            if crate::a2_trace::enabled() {
                trace_loop_for_acceptance_failure(brep, loop_);
            }
            Some(format!("{err:?}"))
        }
        _ => None,
    })
}

fn trace_loop_for_acceptance_failure(brep: &BRep, loop_: popkern_core::topology::LoopHandle) {
    let Some(loop_data) = brep.loops.get(loop_) else {
        return;
    };
    let face = loop_data.face;
    let start = loop_data.half_edge;
    let mut cursor = start;
    eprintln!("[a2-trace] acceptance failure loop {loop_:?} face={face:?}");
    for index in 0..=brep.half_edges.len() {
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        let origin = brep.vertices[he.origin].position;
        let edge_curve = brep.edges.get(he.edge).and_then(|edge| edge.curve.as_ref());
        let curve_summary = edge_curve
            .map(|curve| {
                use popkern_core::geometry::curve::CurveEval;
                let p0 = curve.geometry.point_at(curve.t_start);
                let p1 = curve.geometry.point_at(curve.t_end);
                format!(
                    "curve=({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                    p0.x, p0.y, p0.z, p1.x, p1.y, p1.z
                )
            })
            .unwrap_or_else(|| "curve=<none>".to_string());
        eprintln!(
            "[a2-trace]   [{index}] he={cursor:?} edge={:?} origin={:?}@({:.3},{:.3},{:.3}) next={:?} {curve_summary}",
            he.edge, he.origin, origin.x, origin.y, origin.z, he.next
        );
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
}

fn native_gregory_cap_outer_g0_error(
    brep: &BRep,
    result: &crate::apex_construction::construct::ApexBuildResult,
    max_gap: f64,
) -> Option<String> {
    use popkern_core::geometry::surface::SurfaceEval;

    let cap_faces = result
        .corner_patches
        .iter()
        .map(|p| p.face)
        .collect::<std::collections::BTreeSet<_>>();
    if cap_faces.is_empty() {
        return Some("no native cap patch faces reported".to_string());
    }

    let mut worst_gap = 0.0f64;
    let mut worst_edge = None;
    let mut measured_edges = 0usize;
    for (edge_handle, edge) in brep.edges.iter() {
        let he_a = edge.half_edge;
        let Some(ha) = brep.half_edges.get(he_a) else {
            continue;
        };
        let Some(hb) = brep.half_edges.get(ha.twin) else {
            continue;
        };
        let face_a = brep.loops[ha.loop_].face;
        let face_b = brep.loops[hb.loop_].face;
        if face_a == face_b {
            continue;
        }
        let a_is_cap = cap_faces.contains(&face_a);
        let b_is_cap = cap_faces.contains(&face_b);
        if a_is_cap == b_is_cap {
            continue;
        }
        let (Some(pa), Some(pb)) = (ha.pcurve.as_ref(), hb.pcurve.as_ref()) else {
            continue;
        };
        let Some(sa) = brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let Some(sb) = brep.faces.get(face_b).and_then(|f| f.surface.as_ref()) else {
            continue;
        };

        let (pa_lo, pa_hi) = pa.parameter_range();
        let (pb_lo, pb_hi) = pb.parameter_range();
        let uv_a = pa.point_at(0.5 * (pa_lo + pa_hi));
        let uv_b = pb.point_at(0.5 * (pb_lo + pb_hi));
        let p_a = sa.geometry.point_at(uv_a.x, uv_a.y);
        let p_b = sb.geometry.point_at(uv_b.x, uv_b.y);
        let gap = p_a.distance_to(&p_b);
        measured_edges += 1;
        if crate::a2_trace::cap_g0_enabled() {
            trace_native_gregory_cap_seam_g0(
                brep,
                edge_handle,
                he_a,
                ha.twin,
                face_a,
                face_b,
                a_is_cap,
                gap,
            );
        }
        if gap > worst_gap {
            worst_gap = gap;
            worst_edge = Some(edge_handle);
        }
        if gap > max_gap {
            if crate::a2_trace::cap_g0_enabled() {
                continue;
            }
            return Some(format!(
                "cap_outer_g0 {gap:.3e} on {edge_handle:?} exceeds {max_gap:.3e}"
            ));
        }
    }

    if measured_edges == 0 {
        return Some("no cap-vs-noncap seams measured".to_string());
    }
    if crate::a2_trace::enabled() {
        eprintln!(
            "[a2-trace] native Gregory cap_outer_g0 postcondition: worst={worst_gap:.3e} over {measured_edges} edge(s)"
        );
    }
    if worst_gap > max_gap {
        let Some(worst_edge) = worst_edge else {
            return Some(format!(
                "cap_outer_g0 {worst_gap:.3e} exceeds {max_gap:.3e}, but no worst edge was recorded"
            ));
        };
        return Some(format!(
            "cap_outer_g0 {worst_gap:.3e} on {:?} exceeds {max_gap:.3e}",
            worst_edge
        ));
    }
    None
}

#[allow(clippy::too_many_arguments)]
fn trace_native_gregory_cap_seam_g0(
    brep: &BRep,
    edge_handle: popkern_core::topology::EdgeHandle,
    he_a: popkern_core::topology::HalfEdgeHandle,
    he_b: popkern_core::topology::HalfEdgeHandle,
    face_a: popkern_core::topology::FaceHandle,
    face_b: popkern_core::topology::FaceHandle,
    a_is_cap: bool,
    midpoint_gap: f64,
) {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

    let surface_kind = |face| {
        brep.faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| match &s.geometry {
                SurfaceGeometry::Plane { .. } => "Plane",
                SurfaceGeometry::Cylinder { .. } => "Cylinder",
                SurfaceGeometry::Cone { .. } => "Cone",
                SurfaceGeometry::Sphere { .. } => "Sphere",
                SurfaceGeometry::Torus { .. } => "Torus",
                SurfaceGeometry::Nurbs(_) => "Nurbs",
                SurfaceGeometry::Gregory(_) => "Gregory",
            })
            .unwrap_or("<none>")
    };
    let he_a_data = &brep.half_edges[he_a];
    let he_b_data = &brep.half_edges[he_b];
    eprintln!(
        "[a2-trace cap-g0] seam edge={edge_handle:?} he=({he_a:?},{he_b:?}) \
         faces=({face_a:?}:{}, {face_b:?}:{}) cap_side={} midpoint_gap={midpoint_gap:.6e}",
        surface_kind(face_a),
        surface_kind(face_b),
        if a_is_cap { "A" } else { "B" }
    );

    let Some(edge_curve) = brep
        .edges
        .get(edge_handle)
        .and_then(|edge| edge.curve.as_ref())
    else {
        eprintln!("[a2-trace cap-g0]   edge has no 3D curve");
        return;
    };
    let Some(surface_a) = brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) else {
        return;
    };
    let Some(surface_b) = brep.faces.get(face_b).and_then(|f| f.surface.as_ref()) else {
        return;
    };
    let (Some(pcurve_a), Some(pcurve_b)) = (he_a_data.pcurve.as_ref(), he_b_data.pcurve.as_ref())
    else {
        eprintln!("[a2-trace cap-g0]   missing pcurve on one side");
        return;
    };

    let edge_point = |s: f64| {
        let t = edge_curve.t_start + s * (edge_curve.t_end - edge_curve.t_start);
        edge_curve.geometry.point_at(t)
    };
    let side_point = |pcurve: &popkern_core::geometry::pcurve::Pcurve,
                      surface: &popkern_core::geometry::surface::SurfaceOnFace,
                      s: f64| {
        let (lo, hi) = pcurve.parameter_range();
        let uv = pcurve.point_at(lo + s * (hi - lo));
        surface.geometry.point_at(uv.x, uv.y)
    };
    let forward = |pcurve: &popkern_core::geometry::pcurve::Pcurve,
                   surface: &popkern_core::geometry::surface::SurfaceOnFace| {
        let p0 = side_point(pcurve, surface, 0.0);
        let p1 = side_point(pcurve, surface, 1.0);
        p0.distance_to(&edge_point(0.0)) + p1.distance_to(&edge_point(1.0))
            <= p0.distance_to(&edge_point(1.0)) + p1.distance_to(&edge_point(0.0))
    };
    let a_forward = forward(pcurve_a, surface_a);
    let b_forward = forward(pcurve_b, surface_b);

    for sample in [0.0, 0.5, 1.0] {
        let s_a = if a_forward { sample } else { 1.0 - sample };
        let s_b = if b_forward { sample } else { 1.0 - sample };
        let p_edge = edge_point(sample);
        let p_a = side_point(pcurve_a, surface_a, s_a);
        let p_b = side_point(pcurve_b, surface_b, s_b);
        eprintln!(
            "[a2-trace cap-g0]   s={sample:.1} cap_noncap={:.6e} \
             sideA_edge={:.6e} sideB_edge={:.6e}",
            p_a.distance_to(&p_b),
            p_a.distance_to(&p_edge),
            p_b.distance_to(&p_edge)
        );
    }
}

fn native_gregory_has_expected_cylinder_seams(
    brep: &BRep,
    result: &crate::apex_construction::construct::ApexBuildResult,
    n_fillets: usize,
) -> bool {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let cap_faces = result
        .corner_patches
        .iter()
        .map(|p| p.face)
        .collect::<std::collections::BTreeSet<_>>();
    if cap_faces.is_empty() {
        return false;
    }

    let mut cylinder_seams = 0usize;
    let mut cylinder_faces = std::collections::BTreeSet::new();
    for (_, edge) in brep.edges.iter() {
        let he_a = edge.half_edge;
        let Some(ha) = brep.half_edges.get(he_a) else {
            continue;
        };
        let Some(hb) = brep.half_edges.get(ha.twin) else {
            continue;
        };
        let face_a = brep.loops[ha.loop_].face;
        let face_b = brep.loops[hb.loop_].face;
        if face_a == face_b {
            continue;
        }
        let a_is_cap = cap_faces.contains(&face_a);
        let b_is_cap = cap_faces.contains(&face_b);
        if a_is_cap == b_is_cap {
            continue;
        }
        let neighbor = if a_is_cap { face_b } else { face_a };
        if matches!(
            brep.faces
                .get(neighbor)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry),
            Some(SurfaceGeometry::Cylinder { .. })
        ) {
            cylinder_seams += 1;
            cylinder_faces.insert(neighbor);
        }
    }

    if crate::a2_trace::enabled() {
        eprintln!(
            "[a2-trace] native Gregory seam postcondition: {cylinder_seams} seams across {} distinct cylinder face(s), expected {n_fillets}",
            cylinder_faces.len()
        );
    }
    cylinder_faces.len() >= n_fillets
}
