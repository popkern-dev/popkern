use popkern_core::geometry::Vec3;
use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};
use popkern_core::geometry::surface::SurfaceEval;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::mapping::geometry::GeometryCache;
use crate::step::parser::types::EntityId;
use crate::step::schema::geometry::{
    ConicalSurface, CylindricalSurface, SeamCurve, SphericalSurface, StepCircle, StepEllipse,
    StepPcurve, SurfaceCurve, ToroidalSurface,
};
use crate::step::schema::{as_entity_ref, as_list, param_at};

/// A resolved STEP PCURVE with the entity provenance needed by topology
/// diagnostics and association repair.
#[derive(Debug, Clone)]
pub(crate) struct ResolvedStepPcurve {
    /// STEP `PCURVE` entity id.
    pub pcurve_id: EntityId,
    /// STEP 2D basis curve referenced by the pcurve's definitional representation.
    pub reference_curve_id: EntityId,
    /// Kernel pcurve geometry in the owning surface's parameter frame.
    pub pcurve: Pcurve,
}

impl GeometryCache<'_> {
    /// Resolve every pcurve associated with a surface curve on `surface_id`.
    ///
    /// STEP permits multiple associated pcurves. Import callers that know the
    /// owning half-edge should evaluate all returned candidates against the
    /// recovered edge endpoints rather than accepting the first match.
    ///
    /// A pcurve whose basis surface is a different entity from
    /// `surface_id` is accepted when the two surfaces agree along it
    /// within `linear_tolerance` (see [`Self::basis_agrees_along`]).
    pub(crate) fn resolve_surface_curve_pcurves(
        &mut self,
        curve_id: EntityId,
        surface_id: EntityId,
        linear_tolerance: f64,
        log: &mut DiagnosticLog,
    ) -> Vec<ResolvedStepPcurve> {
        let Some(entity) = self.graph.get(curve_id) else {
            return Vec::new();
        };
        let associated_geometry = match entity.type_name.as_str() {
            "SURFACE_CURVE" => {
                SurfaceCurve::from_params(&entity.params).map(|curve| curve.associated_geometry)
            }
            "SEAM_CURVE" => {
                SeamCurve::from_params(&entity.params).map(|curve| curve.associated_geometry)
            }
            _ => return Vec::new(),
        }
        .map_err(|e| log.error(Phase::MapGeometry, Some(curve_id), format!("{e}")))
        .ok();

        let Some(associated_geometry) = associated_geometry else {
            return Vec::new();
        };

        let mut resolved = Vec::new();
        let mut nearest_foreign: Option<(f64, EntityId)> = None;
        for pcurve_id in associated_geometry {
            match self.resolve_pcurve_on_surface(pcurve_id, surface_id, linear_tolerance, log) {
                Ok(pcurve) => resolved.push(pcurve),
                Err(Some(gap)) if nearest_foreign.is_none_or(|(g, _)| gap < g) => {
                    nearest_foreign = Some((gap, pcurve_id));
                }
                Err(_) => {}
            }
        }
        // A surface curve carries one pcurve per adjacent face, so a
        // foreign-basis pcurve that disagrees is normally the other
        // face's. Only when none fits this face is it a defect.
        if resolved.is_empty()
            && let Some((gap, pcurve_id)) = nearest_foreign
        {
            log.warn(
                Phase::MapGeometry,
                Some(pcurve_id),
                format!(
                    "no PCURVE of #{curve_id} fits face surface #{surface_id}: nearest \
                     foreign-basis pcurve disagrees by {gap:.3e} > {linear_tolerance:.3e}"
                ),
            );
        }
        resolved
    }

    fn resolve_pcurve_on_surface(
        &mut self,
        pcurve_id: EntityId,
        surface_id: EntityId,
        linear_tolerance: f64,
        log: &mut DiagnosticLog,
    ) -> Result<ResolvedStepPcurve, Option<f64>> {
        let entity = self.graph.get(pcurve_id).ok_or(None)?;
        if entity.type_name != "PCURVE" {
            return Err(None);
        }
        let pcurve = StepPcurve::from_params(&entity.params)
            .map_err(|e| log.error(Phase::MapGeometry, Some(pcurve_id), format!("{e}")))
            .map_err(|()| None)?;
        let basis_surface = pcurve.basis_surface;
        let curve_2d_id = self
            .definitional_representation_curve(pcurve.reference_curve, log)
            .ok_or(None)?;
        let pcurve = match self.resolve_parametric_conic(curve_2d_id, log) {
            Some(conic) => conic,
            None => self.resolve_pcurve_geometry(curve_2d_id, log)?,
        };
        // The pcurve's coordinates live in its basis surface's chart.
        let transform = self
            .surface_chart_transform(basis_surface, log)
            .unwrap_or_default();
        let pcurve = apply_chart_transform(pcurve, transform);
        if basis_surface != surface_id {
            let gap = self
                .basis_gap_along(&pcurve, (basis_surface, surface_id), log)
                .ok_or(None)?;
            if gap > linear_tolerance {
                return Err(Some(gap));
            }
        }
        Ok(ResolvedStepPcurve {
            pcurve_id,
            reference_curve_id: curve_2d_id,
            pcurve,
        })
    }

    /// A 2-D CIRCLE or ELLIPSE pcurve as an exact elliptical pcurve,
    /// parameterized like the STEP conic (by angle from the placement's
    /// x direction); `None` for other entity types.
    fn resolve_parametric_conic(
        &mut self,
        curve_2d_id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<Pcurve> {
        let entity = self.graph.get(curve_2d_id)?;
        let (position, a, b) = match entity.type_name.as_str() {
            "CIRCLE" => {
                let c = StepCircle::from_params(&entity.params).ok()?;
                (c.position, c.radius, c.radius)
            }
            "ELLIPSE" => {
                let e = StepEllipse::from_params(&entity.params).ok()?;
                (e.position, e.semi_axis_1, e.semi_axis_2)
            }
            _ => return None,
        };
        let (center, _normal, x_axis) = self.resolve_placement(position, log)?;
        let tau = std::f64::consts::TAU;
        Some(Pcurve::elliptical(
            (center.x, center.y),
            a,
            b,
            x_axis.y.atan2(x_axis.x),
            0.0,
            tau,
            0.0,
            tau,
        ))
    }

    /// A 2-D LINE or B-spline pcurve.
    fn resolve_pcurve_geometry(
        &mut self,
        curve_2d_id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Result<Pcurve, Option<f64>> {
        let curve = self.resolve_curve(curve_2d_id, log).ok_or(None)?;
        match curve.geometry {
            CurveGeometry::Nurbs(nurbs) => Ok(Pcurve::from_nurbs(nurbs)),
            CurveGeometry::Line { origin, direction } => Ok(Pcurve::linear(
                (origin.x, origin.y),
                (origin.x + direction.x, origin.y + direction.y),
                curve.t_start,
                curve.t_end,
            )),
            other => {
                log.warn(
                    Phase::MapGeometry,
                    Some(curve_2d_id),
                    format!("unsupported pcurve basis geometry: {other:?}"),
                );
                Err(None)
            }
        }
    }

    /// How far a pcurve authored on its own copy of a surface (a
    /// separate basis-surface entity, as OpenCASCADE writes every
    /// PCURVE) is from being valid on the face surface: the largest
    /// distance between the two surfaces, resolved into popkern charts,
    /// at the pcurve's uv samples. This is agreement in the
    /// parameterization the pcurve relies on, not merely as point sets.
    /// `None` when either surface does not resolve.
    fn basis_gap_along(
        &mut self,
        pcurve: &Pcurve,
        (basis_surface, face_surface): (EntityId, EntityId),
        log: &mut DiagnosticLog,
    ) -> Option<f64> {
        let basis = self.resolve_surface(basis_surface, true, log)?;
        let face = self.resolve_surface(face_surface, true, log)?;
        let mut gap = 0.0_f64;
        for k in 0..=BASIS_AGREEMENT_SAMPLES {
            let s = k as f64 / BASIS_AGREEMENT_SAMPLES as f64;
            let uv = pcurve.point_at(pcurve.t_start + s * (pcurve.t_end - pcurve.t_start));
            let a = basis.geometry.point_at(uv.x, uv.y);
            let b = face.geometry.point_at(uv.x, uv.y);
            gap = gap.max(a.distance_to(&b));
        }
        Some(gap)
    }

    fn definitional_representation_curve(
        &self,
        representation_id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<EntityId> {
        let entity = self.graph.get(representation_id)?;
        if entity.type_name != "DEFINITIONAL_REPRESENTATION" {
            log.warn(
                Phase::MapGeometry,
                Some(representation_id),
                "PCURVE reference is not a DEFINITIONAL_REPRESENTATION",
            );
            return None;
        }
        let items = param_at(&entity.params, 1, "DEFINITIONAL_REPRESENTATION")
            .and_then(|param| as_list(param, "DEFINITIONAL_REPRESENTATION"))
            .map_err(|e| log.error(Phase::MapGeometry, Some(representation_id), format!("{e}")))
            .ok()?;
        items
            .first()
            .and_then(|param| as_entity_ref(param, "DEFINITIONAL_REPRESENTATION").ok())
    }

    /// Chart transform mapping STEP pcurve coordinates to popkern's
    /// chart for the resolved surface: a u offset for analytic
    /// surfaces whose placement ref_direction differs from the
    /// kernel chart's u origin (`orthonormal_basis(axis).0`), plus
    /// the latitude→colatitude flip for spheres (STEP v is measured
    /// from the equator, popkern's from the +axis pole:
    /// v_popkern = π/2 − v_step).
    fn surface_chart_transform(
        &mut self,
        surface_id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<ChartTransform> {
        let entity = self.graph.get(surface_id)?;
        let (placement, latitude_flip) = match entity.type_name.as_str() {
            "CYLINDRICAL_SURFACE" => (
                CylindricalSurface::from_params(&entity.params)
                    .ok()?
                    .position,
                false,
            ),
            "CONICAL_SURFACE" => (
                ConicalSurface::from_params(&entity.params).ok()?.position,
                false,
            ),
            "TOROIDAL_SURFACE" => (
                ToroidalSurface::from_params(&entity.params).ok()?.position,
                false,
            ),
            "SPHERICAL_SURFACE" => (
                SphericalSurface::from_params(&entity.params).ok()?.position,
                true,
            ),
            _ => return None,
        };
        let (_origin, axis, ref_direction) = self.resolve_placement(placement, log)?;
        Some(ChartTransform {
            u_offset: analytic_u_frame_offset(axis, ref_direction).unwrap_or(0.0),
            latitude_flip,
        })
    }
}

/// Uniform pcurve samples (plus one) compared between a pcurve's own
/// basis surface and the face surface.
const BASIS_AGREEMENT_SAMPLES: usize = 32;

/// STEP→popkern pcurve coordinate transform for one basis surface.
#[derive(Debug, Clone, Copy, Default)]
struct ChartTransform {
    u_offset: f64,
    latitude_flip: bool,
}

fn apply_chart_transform(pcurve: Pcurve, transform: ChartTransform) -> Pcurve {
    let pcurve = translate_pcurve_u(pcurve, transform.u_offset);
    if !transform.latitude_flip {
        return pcurve;
    }
    flip_pcurve_latitude(pcurve)
}

/// v → π/2 − v (STEP sphere latitude to popkern colatitude). Import
/// produces Linear, NURBS and (from 2-D conics) Elliptical pcurves (see
/// `resolve_pcurve_on_surface`); all three are closed under this
/// reflection. An ellipse reflects exactly by negating its rotation and
/// sweep angles: `c + R(θ)(a cos α, b sin α)` maps to
/// `c' + R(−θ)(a cos(−α), b sin(−α))`.
fn flip_pcurve_latitude(mut pcurve: Pcurve) -> Pcurve {
    let half_pi = std::f64::consts::FRAC_PI_2;
    match &mut pcurve.geometry {
        PcurveGeometry::Linear { start, end } => {
            start.1 = half_pi - start.1;
            end.1 = half_pi - end.1;
        }
        PcurveGeometry::Nurbs(nurbs) => {
            for point in &mut nurbs.control_points {
                point.y = half_pi - point.y;
            }
        }
        PcurveGeometry::Elliptical {
            center,
            rotation,
            start_angle,
            end_angle,
            ..
        } => {
            center.1 = half_pi - center.1;
            *rotation = -*rotation;
            *start_angle = -*start_angle;
            *end_angle = -*end_angle;
        }
        // Not produced by import; leave unchanged rather than guessing
        // a convention.
        PcurveGeometry::Sinusoidal { .. } | PcurveGeometry::CylinderArcViaRationalQuad { .. } => {}
    }
    pcurve
}

fn analytic_u_frame_offset(axis: Vec3, step_ref_direction: Vec3) -> Option<f64> {
    let axis = axis.normalized()?;
    let step_u = (step_ref_direction - axis * step_ref_direction.dot(&axis)).normalized()?;
    let (kernel_u, kernel_v) = orthonormal_basis(axis);
    Some(step_u.dot(&kernel_v).atan2(step_u.dot(&kernel_u)))
}

fn translate_pcurve_u(mut pcurve: Pcurve, offset: f64) -> Pcurve {
    if offset.abs() <= 1e-14 {
        return pcurve;
    }
    match &mut pcurve.geometry {
        PcurveGeometry::Linear { start, end } => {
            start.0 += offset;
            end.0 += offset;
        }
        PcurveGeometry::Elliptical { center, .. } => {
            center.0 += offset;
        }
        PcurveGeometry::Sinusoidal { u_offset, .. } => {
            *u_offset += offset;
        }
        PcurveGeometry::Nurbs(nurbs) => {
            for point in &mut nurbs.control_points {
                point.x += offset;
            }
        }
        PcurveGeometry::CylinderArcViaRationalQuad {
            theta_a, theta_b, ..
        } => {
            *theta_a += offset;
            *theta_b += offset;
        }
    }
    pcurve
}

#[cfg(test)]
#[path = "pcurves_tests.rs"]
mod tests;
