//! the sphere-fit edge-wiring workirectApex apex planner — orchestration.
//!
//! `plan_apex_resolution` is the **pure** entry point that produces
//! an [`ApexPlan`] from a BRep + apex vertex + fillet specs.
//! It does NOT mutate the BRep — callers can preview, inspect
//! telemetry, or cancel.
//!
//! Algorithm (per ADR §6.2):
//! 1. Validate inputs ([`plan_validate`](super::plan_validate)).
//! 2. Compute boundary segments
//!    ([`plan_terminations::compute_apex_boundary`](super::plan_terminations::compute_apex_boundary)).
//! 3. Try the [`ConstructionRung`] ladder Direct → Setback → Open
//!    ([`plan_with_rung_ladder`]):
//!    * Direct: full per-fillet sphere-fit + Sphere/Gregory cap-kind
//!      decision + target topology assembly.
//!    * Setback ([`plan_rung_fallback::plan_setback_rung`](super::plan_rung_fallback::plan_setback_rung)):
//!      lenient per-fillet sphere-fit, ≥3 successful fillets required.
//!    * Open ([`plan_rung_fallback::plan_open_rung`](super::plan_rung_fallback::plan_open_rung)):
//!      empty apex topology, fillets terminate at rolling-ball
//!      footprints.
//! 4. Telemetry records the rung used + any fallback reason.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, VertexHandle};
use serde::{Deserialize, Serialize};

use crate::cancellation::CancellationToken;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::FilletSpec;

use super::cap_kind::CapKind;
use super::error::ApexPlanError;
use super::plan_rung_fallback::{plan_open_rung, plan_setback_rung};
use super::plan_terminations::{compute_apex_boundary, compute_per_fillet_terminations};
use super::plan_topology::{build_target_topology, decide_cap_kind, gregory_cap_kind};
use super::plan_validate::{validate_inputs, validate_inputs_lenient};
use super::target_topology::TargetTopology;
use super::telemetry::{ConstructionRung, PlanTelemetry};

/// A plan for resolving an apex vertex's fillet topology.
///
/// Produced by [`plan_apex_resolution`]; consumed by
/// `construct_apex_topology` (when DirectApex.2 lands).  Introspectable for
/// UI preview.
///
/// See ADR `docs/decisions/2026-04-30-phase-6c-a2-direct-apex-construction.md` §3.1.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ApexPlan {
    /// The apex vertex this plan resolves.
    pub apex_vertex: VertexHandle,
    /// Cap kind decision (Sphere or Gregory) with computed geometry.
    pub cap_kind: CapKind,
    /// Final topology spec — vertices, edges, faces with their geometry.
    pub target_topology: TargetTopology,
    /// Decision telemetry — explains why this plan was chosen, which
    /// rung was used, and exposes diagnostic numerics.  UI-friendly.
    pub telemetry: PlanTelemetry,
}

impl ApexPlan {
    /// Number of fillets meeting at the apex.
    pub fn n_fillets(&self) -> usize {
        self.cap_kind.n_fillets()
    }

    /// One-line, end-user-readable summary of the plan.  Suitable
    /// for status bars, tooltips, or log lines.  Includes the rung,
    /// cap-kind dispatch, and fillet count.
    ///
    /// Examples:
    /// * `"Direct: Sphere cap, 3 fillets"`
    /// * `"Setback: Gregory cap, 4 fillets (1 omitted)"`
    /// * `"Open: no cap; 4 fillets terminate at rolling-ball footprints"`
    pub fn plan_summary(&self) -> String {
        let rung = match &self.telemetry.rung {
            ConstructionRung::Direct => "Direct".to_string(),
            ConstructionRung::Setback { setback_distance } => {
                if *setback_distance > 0.0 {
                    format!("Setback ({setback_distance:.3e} m)")
                } else {
                    "Setback".to_string()
                }
            }
            ConstructionRung::Open => "Open".to_string(),
        };
        let n = self.n_fillets();
        match &self.cap_kind {
            CapKind::Sphere { .. } => format!("{rung}: Sphere cap, {n} fillets"),
            CapKind::Gregory { .. } => format!("{rung}: Gregory cap, {n} fillets"),
            CapKind::Open {
                n_fillets_attempted,
            } => format!(
                "{rung}: no cap; {n_fillets_attempted} fillets terminate at rolling-ball footprints",
            ),
        }
    }

    /// Axis-aligned bounding box of the apex region.
    ///
    /// Returns `(min_corner, max_corner)` enclosing every
    /// `apex_vertex` plus every cap-sphere centre.  Useful for
    /// camera framing in UI previews.  Returns `None` for the Open
    /// rung, which has no apex geometry to box.
    pub fn apex_aabb(
        &self,
    ) -> Option<(
        popkern_core::geometry::Point3,
        popkern_core::geometry::Point3,
    )> {
        use popkern_core::geometry::Point3;
        let mut points: Vec<Point3> = self.target_topology.apex_vertices.clone();
        match &self.cap_kind {
            CapKind::Sphere { center, .. } => points.push(*center),
            CapKind::Gregory {
                per_fillet_caps, ..
            } => {
                points.extend(per_fillet_caps.iter().map(|(c, _)| *c));
            }
            CapKind::Open { .. } => {}
        }
        if points.is_empty() {
            return None;
        }
        let mut min = points[0];
        let mut max = points[0];
        for p in &points[1..] {
            min.x = min.x.min(p.x);
            min.y = min.y.min(p.y);
            min.z = min.z.min(p.z);
            max.x = max.x.max(p.x);
            max.y = max.y.max(p.y);
            max.z = max.z.max(p.z);
        }
        Some((min, max))
    }

    /// Sample points on (or around) the cap surface for preview
    /// rendering.  Approximate — purpose is "show the user what the
    /// cap looks like before committing to construct_apex".
    ///
    /// For Sphere caps: returns the tangent-arc endpoints plus the
    /// cap centre.  Drawing line segments between adjacent endpoints
    /// approximates the cap's outer loop; the centre helps a UI
    /// orient a sphere-cap glyph.
    ///
    /// For Gregory caps: returns all 2N apex vertices (T_a + T_b
    /// for each fillet).  Drawing the polyline through these
    /// approximates the Gregory cap's outer loop.
    ///
    /// For Open: returns an empty vector.
    pub fn cap_surface_preview(&self) -> Vec<popkern_core::geometry::Point3> {
        match &self.cap_kind {
            CapKind::Sphere {
                center,
                tangent_arc_endpoints,
                ..
            } => {
                let mut out = tangent_arc_endpoints.clone();
                out.push(*center);
                out
            }
            CapKind::Gregory {
                per_fillet_terminations,
                ..
            } => per_fillet_terminations
                .iter()
                .flat_map(|[ta, tb]| [*ta, *tb])
                .collect(),
            CapKind::Open { .. } => Vec::new(),
        }
    }
}

/// Compute an [`ApexPlan`] for resolving an apex vertex's fillet
/// topology.  Pure — does not mutate the BRep.
///
/// # Inputs
///
/// * `brep` — the BRep being filleted (read-only).
/// * `apex` — the vertex where multiple fillets meet.
/// * `fillet_specs` — the fillet specs (one per incident fillet edge),
///   in the same order as the cap-loop walk would visit them.
/// * `fillet_surfaces` — fillet surface data (cylinder/torus surfaces
///   + trim curves), one per fillet, computed earlier in the
///     pipeline by `geometry_compute`.
/// * `tolerance` — geometric tolerance for cap-kind dispatch and
///   coincidence checks.
/// * `cancel` — optional cooperative cancellation token.  Polled at
///   each rung boundary; when cancelled the planner returns
///   [`ApexPlanError::Cancelled`] without further work.  Pass `None`
///   for non-cancellable use.
///
/// # Output
///
/// `Ok(ApexPlan)` on success; `Err(ApexPlanError)` describing the
/// failure mode for user-actionable feedback.
///
/// # Errors
///
/// See [`ApexPlanError`] for the full failure taxonomy.
pub fn plan_apex_resolution(
    brep: &BRep,
    apex: VertexHandle,
    fillet_specs: &[FilletSpec],
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
    cancel: Option<&CancellationToken>,
) -> Result<ApexPlan, ApexPlanError> {
    // Step 1: validate inputs (DirectApex.1.b).  Validation errors don't
    // benefit from rung fallback — invalid inputs are user errors.
    validate_inputs(brep, apex, fillet_specs)?;
    check_cancelled(cancel)?;

    // Step 2: compute boundary segments (DirectApex.1.c).  Same — boundary
    // failure is structural, not a sphere-fit failure.
    let boundary = compute_apex_boundary(brep, apex, fillet_surfaces, tolerance)?;
    check_cancelled(cancel)?;

    // Steps 3-5 with rung fallback: try Direct → Setback → Open.
    plan_with_rung_ladder(
        brep,
        apex,
        &boundary,
        fillet_specs,
        fillet_surfaces,
        tolerance,
        cancel,
    )
}

/// Mid-pipeline-friendly variant of [`plan_apex_resolution`].
///
/// Identical to the public entry point except it skips the
/// edge-existence and edge-incidence input checks, which fail
/// when the `apply_blend_geometry` pipeline calls the
/// planner mid-flight after the consumed-merge stage has kef'd the original
/// fillet edges (their handles are no longer in
/// `brep.edges`).  Radius / apex-vertex / fillet-count checks
/// still apply.
///
/// `pub(crate)` — only the fillet pipeline's DirectApex dispatch
/// (`crate::construct::direct_apex_dispatch`) calls this.  External
/// callers should use [`plan_apex_resolution`] which validates
/// the spec strictly.
pub(crate) fn plan_apex_resolution_lenient(
    brep: &BRep,
    apex: VertexHandle,
    fillet_specs: &[FilletSpec],
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
    cancel: Option<&CancellationToken>,
) -> Result<ApexPlan, ApexPlanError> {
    validate_inputs_lenient(brep, apex, fillet_specs)?;
    check_cancelled(cancel)?;

    let boundary = compute_apex_boundary(brep, apex, fillet_surfaces, tolerance)?;
    check_cancelled(cancel)?;

    plan_with_rung_ladder(
        brep,
        apex,
        &boundary,
        fillet_specs,
        fillet_surfaces,
        tolerance,
        cancel,
    )
}

/// Polls the cancellation token (if any) and returns
/// [`ApexPlanError::Cancelled`] if it has been signalled.
fn check_cancelled(cancel: Option<&CancellationToken>) -> Result<(), ApexPlanError> {
    if cancel.is_some_and(|c| c.is_cancelled()) {
        Err(ApexPlanError::Cancelled)
    } else {
        Ok(())
    }
}

/// Try [`ConstructionRung::Direct`] first, then Setback, then Open.
/// Records the successful rung (and any fallback reason) in
/// [`PlanTelemetry`].  Polls `cancel` at each rung boundary.
fn plan_with_rung_ladder(
    brep: &BRep,
    apex: VertexHandle,
    boundary: &BlendBoundaryRef,
    fillet_specs: &[FilletSpec],
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &popkern_core::tolerance::Tolerance,
    cancel: Option<&CancellationToken>,
) -> Result<ApexPlan, ApexPlanError> {
    // --- Rung 1: Direct ---
    check_cancelled(cancel)?;
    match compute_per_fillet_terminations(brep, apex, boundary, fillet_surfaces, tolerance) {
        Ok(terminations) => {
            let mut telemetry = PlanTelemetry::new();
            let mut cap_kind =
                decide_cap_kind(&terminations, fillet_specs, tolerance, &mut telemetry)?;
            let mut target_topology = build_target_topology(&cap_kind, &terminations);
            if let CapKind::Sphere { radius, .. } = &cap_kind
                && brep_has_all_fillet_cylinder_faces(brep, fillet_surfaces)
                && let Some(max_gap) = sphere_anchor_support_gap(brep, &target_topology, *radius)
            {
                if popkern_core::trace::trace_enabled("sphere-fit") {
                    eprintln!(
                        "[spherefit] apex {apex:?}: Sphere DEMOTED to Gregory — anchor gap {max_gap:.4e} > support {:.4e}",
                        radius * 0.9
                    );
                }
                telemetry.set_value("sphere_anchor_support_gap", max_gap);
                telemetry.cap_kind_reason = format!(
                    "{}; Sphere target anchor gap {:.3e} exceeds relocation support {:.3e} — Gregory cap chosen for constructive topology compatibility",
                    telemetry.cap_kind_reason,
                    max_gap,
                    radius * 0.9
                );
                cap_kind = gregory_cap_kind(&terminations, tolerance);
                target_topology = build_target_topology(&cap_kind, &terminations);
            }
            if popkern_core::trace::trace_enabled("sphere-fit") {
                eprintln!(
                    "[spherefit] apex {apex:?}: Direct rung OK, sphere={}",
                    cap_kind.is_sphere()
                );
            }
            telemetry.rung = ConstructionRung::Direct;
            Ok(ApexPlan {
                apex_vertex: apex,
                cap_kind,
                target_topology,
                telemetry,
            })
        }
        Err(direct_err) => {
            if popkern_core::trace::trace_enabled("sphere-fit") {
                eprintln!("[spherefit] apex {apex:?}: Direct rung DECLINED: {direct_err:?}");
            }
            // --- Rung 2: Setback ---
            //
            // Setback strategy: tolerate individual per-fillet
            // sphere-fit failures.  Successful fillets contribute
            // apex_vertices/edges; failed fillets are skipped.
            // Requires ≥3 successful fillets.
            check_cancelled(cancel)?;
            if let Ok(plan) = plan_setback_rung(
                brep,
                apex,
                boundary,
                fillet_specs,
                fillet_surfaces,
                tolerance,
                &direct_err,
            ) {
                return Ok(plan);
            }

            // --- Rung 3: Open ---
            //
            // No cap construction.  Plan has empty apex topology;
            // builder leaves apex as raw rolling-ball terminations.
            check_cancelled(cancel)?;
            Ok(plan_open_rung(apex, fillet_surfaces.len(), &direct_err))
        }
    }
}

fn brep_has_all_fillet_cylinder_faces(
    brep: &BRep,
    fillet_surfaces: &[FilletSurfaceResult],
) -> bool {
    fillet_surfaces
        .iter()
        .all(|fillet| matching_cylinder_face_exists(brep, fillet))
}

fn matching_cylinder_face_exists(brep: &BRep, fillet: &FilletSurfaceResult) -> bool {
    let SurfaceGeometry::Cylinder {
        origin: f_origin,
        axis: f_axis,
        radius: f_radius,
    } = &fillet.surface
    else {
        return false;
    };
    let Some(f_axis_unit) = f_axis.normalized() else {
        return false;
    };
    brep.faces.values().any(|face| {
        let Some(surface) = face.surface.as_ref() else {
            return false;
        };
        let SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } = &surface.geometry
        else {
            return false;
        };
        if (*radius - *f_radius).abs() > f_radius.abs() * 1e-4 + 1e-9 {
            return false;
        }
        let Some(axis_unit) = axis.normalized() else {
            return false;
        };
        if axis_unit.dot(&f_axis_unit).abs() < 1.0 - 1e-3 {
            return false;
        }
        let d = popkern_core::geometry::Vec3::new(
            origin.x - f_origin.x,
            origin.y - f_origin.y,
            origin.z - f_origin.z,
        );
        let along = d.dot(&f_axis_unit);
        let perp = popkern_core::geometry::Vec3::new(
            d.x - along * f_axis_unit.x,
            d.y - along * f_axis_unit.y,
            d.z - along * f_axis_unit.z,
        );
        perp.length() <= 1e-3
    })
}

fn sphere_anchor_support_gap(
    brep: &BRep,
    target_topology: &TargetTopology,
    radius: f64,
) -> Option<f64> {
    if radius <= 0.0 {
        return Some(f64::INFINITY);
    }
    let support_radius = radius * 0.9;
    let mut worst = 0.0_f64;
    for target in &target_topology.apex_vertices {
        let nearest = brep
            .vertices
            .values()
            .map(|vertex| vertex.position.distance_to(target))
            .fold(f64::INFINITY, f64::min);
        if nearest > worst {
            worst = nearest;
        }
    }
    (worst > support_radius).then_some(worst)
}

// Type alias to keep the rung-ladder signature concise without
// re-importing BlendBoundary at the top of plan.rs (it's only
// referenced inside this private fn).
type BlendBoundaryRef = crate::blend::boundary::BlendBoundary;
