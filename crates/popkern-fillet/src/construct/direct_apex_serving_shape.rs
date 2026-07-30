//! Serving-shape defect signature + absorb/consolidate fixpoint
//! for the constructive direct-apex dispatch (PK-205/PK-206).

use popkern_core::topology::BRep;

use crate::apex_construction::constructive_builder::construct_apex_constructively;

#[path = "direct_apex_serving_shape_redo.rs"]
mod serving_shape_redo;

use super::vfree_fallback;

/// One absorb → consolidate round-trip, iterated to a bounded
/// fixpoint: the global weld is what makes ghost fragments
/// unglued/absorbable, and an absorption exposes new dead sides
/// for the weld (PK-205/PK-206).
pub(super) fn run_absorb_consolidate_fixpoint(
    brep: &mut BRep,
    result: &mut crate::apex_construction::construct::ApexBuildResult,
) {
    for round in 0..3 {
        let absorbed = crate::apex_construction::absorb_coplanar_residual_scaffold_late(brep);
        if absorbed > 0 && crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace] late carrier absorption round {round} merged {absorbed} \
                 scaffold fragment(s)"
            );
        }
        let (consolidated, diagnostic) = vfree_fallback::consolidate_served_shell(brep);
        if let Some(diagnostic) = diagnostic {
            result.diagnostics.push(diagnostic);
        }
        if absorbed == 0 && !consolidated {
            break;
        }
    }
}

/// Serving-shape defect signature: glued surfaceless residual
/// membranes or dead boundary sides in the body a blend is about
/// to serve. Both break measure/tessellation/STEP (a surfaceless
/// face has no geometric definition; a dead side is an open-sheet
/// rim), so a constructive build that ends here has an incomplete
/// support-rewrite handover.
#[derive(Clone, Copy, Debug)]
pub(super) struct ServedShapeDefect {
    pub(super) surfaceless_residuals: usize,
    pub(super) dead_sides: usize,
}

impl ServedShapeDefect {
    pub(super) fn describe(&self) -> String {
        format!(
            "{} glued surfaceless residual face(s), {} dead boundary side(s)",
            self.surfaceless_residuals, self.dead_sides
        )
    }

    /// Strictly better: no component worse, at least one better.
    pub(super) fn strictly_better_than(&self, other: &Self) -> bool {
        self.surfaceless_residuals <= other.surfaceless_residuals
            && self.dead_sides <= other.dead_sides
            && (self.surfaceless_residuals < other.surfaceless_residuals
                || self.dead_sides < other.dead_sides)
    }
}

pub(super) fn served_shape_defect(brep: &BRep) -> Option<ServedShapeDefect> {
    let surfaceless_residuals = brep
        .faces
        .iter()
        .filter(|(_, face)| {
            face.role == popkern_core::topology::TopologyRole::ConstructionResidual
                && face.surface.is_none()
        })
        .count();
    let dead_sides = brep.boundary_edges.len();
    if surfaceless_residuals == 0 && dead_sides == 0 {
        return None;
    }
    Some(ServedShapeDefect {
        surfaceless_residuals,
        dead_sides,
    })
}

pub(super) fn apply_serving_shape_gate(
    brep: &mut BRep,
    snapshot: &BRep,
    plan: &crate::apex_construction::plan::ApexPlan,
    surfaces: &[crate::fillet_surface::FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
    result: &mut crate::apex_construction::construct::ApexBuildResult,
) {
    let Some(in_place_defect) = served_shape_defect(brep) else {
        return;
    };
    // Rung 1 (PK-214): complete the in-place handover on a clone —
    // retire the glued mirror-ghost complex, run the consolidation
    // with the winding reconciliation forced past its flap gate,
    // then close the partnerless chord rings by boundary
    // completion. Accepted only on a DEFECT-FREE census plus a clean
    // validate; anything partial falls through to the aggregate
    // redo unchanged, so the rung can only improve the served body.
    {
        let mut completed = brep.clone();
        let mut completed_result = result.clone();
        let retired =
            crate::apex_construction::v_free_inplace_handover::retire_mirror_ghost_complex(
                &mut completed,
            );
        if retired > 0 {
            let consolidation = crate::apex_construction::v_free_shell_consolidation::
                consolidate_vestigial_support_shell_forcing_winding(&mut completed);
            run_absorb_consolidate_fixpoint(&mut completed, &mut completed_result);
            let rebuilt = crate::apex_construction::v_free_shell_consolidation::
                complete_partnerless_chord_rings(&mut completed);
            run_absorb_consolidate_fixpoint(&mut completed, &mut completed_result);
            let records_remap = {
                let mut sectors = std::mem::take(&mut completed_result.exact_boundary_sectors);
                let ok = crate::apex_construction::v_free_inplace_handover::
                    remap_exact_boundary_sector_records(
                        &completed,
                        &mut sectors,
                        &mut completed_result.exact_boundary_report.sectors,
                    );
                completed_result.exact_boundary_sectors = sectors;
                ok
            };
            if records_remap
                && served_shape_defect(&completed).is_none()
                && completed.validate().is_ok()
            {
                if crate::direct_apex_trace::enabled() {
                    eprintln!(
                        "[direct-apex-trace] in-place serving-shape gate: handover completed \
                         in place (ghosts {retired}, reversed {}, rebuilt rings {rebuilt}); \
                         was {in_place_defect:?}",
                        consolidation.reversed_faces
                    );
                }
                completed_result.diagnostics.push(
                    crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
                        crate::types::ExactBoundaryCapDiagnostic {
                            stage: "constructive_serving_shape".to_string(),
                            code: "constructive_serving_shape.handover_completed".to_string(),
                            detail_code: Some(format!(
                                "was {}; ghosts {retired}, reversed {}, rebuilt rings {rebuilt}",
                                in_place_defect.describe(),
                                consolidation.reversed_faces
                            )),
                            sector_name: "direct_apex".to_string(),
                            edge: None,
                            support_half_edge: None,
                            support_face: None,
                            residual: None,
                            outer_side: None,
                            sample_parameter: None,
                            tolerance: None,
                            shared_fact_ids: Vec::new(),
                        },
                    ),
                );
                *brep = completed;
                *result = completed_result;
                return;
            }
            if crate::direct_apex_trace::enabled() {
                eprintln!(
                    "[direct-apex-trace] in-place handover completion incomplete \
                     (defect {:?}); falling through to aggregate redo",
                    served_shape_defect(&completed)
                );
                for edge in &completed.boundary_edges {
                    let owner = completed
                        .half_edges
                        .iter()
                        .find(|(_, he)| he.edge == *edge)
                        .map(|(h, he)| {
                            let face = completed.loops.get(he.loop_).map(|l| l.face);
                            let a = completed
                                .vertices
                                .get(he.origin)
                                .map(|v| v.position)
                                .unwrap_or(popkern_core::geometry::Point3::new(9.9, 9.9, 9.9));
                            let b = completed
                                .half_edges
                                .get(he.next)
                                .and_then(|n| completed.vertices.get(n.origin))
                                .map(|v| v.position)
                                .unwrap_or(popkern_core::geometry::Point3::new(9.9, 9.9, 9.9));
                            format!(
                                "{h:?} face={face:?} span=({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
                                a.x, a.y, a.z, b.x, b.y, b.z
                            )
                        })
                        .unwrap_or_else(|| "unowned".to_string());
                    eprintln!(
                        "[direct-apex-trace]   completion residue: dead edge {edge:?} {owner}"
                    );
                }
            }
        }
    }
    {
        if crate::direct_apex_trace::enabled() {
            eprintln!(
                "[direct-apex-trace] in-place serving-shape gate: {in_place_defect:?}; \
                 attempting aggregate redo"
            );
        }
        // Keep-better semantics: redo the build through the
        // aggregate bridge path and serve whichever body has
        // the strictly better serving signature. The redo
        // never worsens the served shape, and the outcome is
        // recorded either way — as an aggregate_fallback
        // diagnostic when the aggregate body wins (measured:
        // the rational-Gregory family goes 12 membranes/18
        // dead sides → 0/0), or as a residual_defect_retained
        // diagnostic when the incumbent in-place body is kept
        // (measured: pyramid4_uniform keeps its 1-membrane
        // shape because the aggregate redo is worse there —
        // 9 membranes/20 dead sides; completing that handover
        // is PK-214).
        let in_place_body = brep.clone();
        let in_place_result = result.clone();
        *brep = snapshot.clone();
        let aggregate_outcome = {
            let _aggregate_only = crate::apex_construction::
                repeated_apex_side_support_bridge_operator::InPlaceRewriteDisabledGuard::activate();
            let aggregate_options =
                crate::apex_construction::constructive_builder_options::
                    ConstructiveApexExecutionOptions::from_environment_adapter();
            construct_apex_constructively(brep, plan, surfaces, tolerance, None, aggregate_options)
        };
        let serving_diagnostic = |code: &str, detail: String| {
            crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(
                crate::types::ExactBoundaryCapDiagnostic {
                    stage: "constructive_serving_shape".to_string(),
                    code: format!("constructive_serving_shape.{code}"),
                    detail_code: Some(detail),
                    sector_name: "direct_apex".to_string(),
                    edge: None,
                    support_half_edge: None,
                    support_face: None,
                    residual: None,
                    outer_side: None,
                    sample_parameter: None,
                    tolerance: None,
                    shared_fact_ids: Vec::new(),
                },
            )
        };
        let keep_in_place = |brep: &mut BRep,
                             result: &mut crate::apex_construction::construct::ApexBuildResult,
                             aggregate_note: String| {
            *brep = in_place_body.clone();
            *result = in_place_result.clone();
            result.diagnostics.push(serving_diagnostic(
                "residual_defect_retained",
                format!(
                    "in-place: {}; aggregate redo not better: {aggregate_note}",
                    in_place_defect.describe()
                ),
            ));
        };
        match aggregate_outcome {
            Ok(mut aggregate_result) => {
                run_absorb_consolidate_fixpoint(brep, &mut aggregate_result);
                let aggregate_defect = served_shape_defect(brep);
                let aggregate_wins = match aggregate_defect.as_ref() {
                    None => true,
                    Some(aggregate) => aggregate.strictly_better_than(&in_place_defect),
                };
                if aggregate_wins {
                    aggregate_result.diagnostics.push(serving_diagnostic(
                        "aggregate_fallback",
                        format!(
                            "in-place: {}; aggregate: {}",
                            in_place_defect.describe(),
                            aggregate_defect
                                .as_ref()
                                .map(ServedShapeDefect::describe)
                                .unwrap_or_else(|| "clean".to_string())
                        ),
                    ));
                    *result = aggregate_result;
                } else {
                    keep_in_place(
                        brep,
                        result,
                        aggregate_defect
                            .as_ref()
                            .map(ServedShapeDefect::describe)
                            .unwrap_or_else(|| "clean".to_string()),
                    );
                }
            }
            Err(err) => {
                keep_in_place(brep, result, format!("aggregate build failed: {err}"));
            }
        }
    }
    serving_shape_redo::apply_apex_retaining_redo_rung(
        brep, snapshot, plan, surfaces, tolerance, result,
    );
}

/// Pre-call dispatch trace: face-kind census and cylinder loop
/// walks, printed under `POPKERN_DIRECT_APEX_TRACE`.
pub(super) fn trace_pre_call_state(brep: &BRep, vertex: popkern_core::topology::VertexHandle) {
    if !crate::direct_apex_trace::enabled() {
        return;
    }
    {
        eprintln!(
            "[direct-apex-trace] try_direct_apex_construction: vertex={vertex:?}; \
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
                    SurfaceGeometry::PolyhedralSpline(_) => "PolyhedralSpline".to_string(),
                })
                .unwrap_or_else(|| "None".to_string());
            eprintln!("[direct-apex-trace]   {fh:?}: {kind}");
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
                "[direct-apex-trace]   {fh:?} outer loop ({} half-edges): {}",
                chain.len(),
                chain.join(" -> ")
            );
        }
    }
}
