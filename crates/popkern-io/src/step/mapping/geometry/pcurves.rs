use popkern_core::geometry::Vec3;
use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::mapping::geometry::GeometryCache;
use crate::step::parser::types::EntityId;
use crate::step::schema::geometry::{
    ConicalSurface, CylindricalSurface, SeamCurve, SphericalSurface, StepPcurve, SurfaceCurve,
    ToroidalSurface,
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
    pub(crate) fn resolve_surface_curve_pcurves(
        &mut self,
        curve_id: EntityId,
        surface_id: EntityId,
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

        associated_geometry
            .into_iter()
            .filter_map(|pcurve_id| self.resolve_pcurve_on_surface(pcurve_id, surface_id, log))
            .collect()
    }

    fn resolve_pcurve_on_surface(
        &mut self,
        pcurve_id: EntityId,
        surface_id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<ResolvedStepPcurve> {
        let entity = self.graph.get(pcurve_id)?;
        if entity.type_name != "PCURVE" {
            return None;
        }
        let pcurve = StepPcurve::from_params(&entity.params)
            .map_err(|e| log.error(Phase::MapGeometry, Some(pcurve_id), format!("{e}")))
            .ok()?;
        if pcurve.basis_surface != surface_id {
            return None;
        }
        let curve_2d_id = self.definitional_representation_curve(pcurve.reference_curve, log)?;
        let curve = self.resolve_curve(curve_2d_id, log)?;
        let pcurve = match curve.geometry {
            CurveGeometry::Nurbs(nurbs) => Some(Pcurve::from_nurbs(nurbs)),
            CurveGeometry::Line { origin, direction } => Some(Pcurve::linear(
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
                None
            }
        }?;
        let transform = self
            .surface_chart_transform(surface_id, log)
            .unwrap_or_default();
        Some(ResolvedStepPcurve {
            pcurve_id,
            reference_curve_id: curve_2d_id,
            pcurve: apply_chart_transform(pcurve, transform),
        })
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
/// only produces Linear and NURBS pcurves (see
/// `resolve_pcurve_on_surface`); both are affine-closed under this
/// map.
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
        // Not produced by import for sphere bases; leave unchanged
        // rather than guessing a convention.
        PcurveGeometry::Elliptical { .. }
        | PcurveGeometry::Sinusoidal { .. }
        | PcurveGeometry::CylinderArcViaRationalQuad { .. } => {}
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
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::nurbs_curve::NurbsCurve;

    use crate::step::document::{OpaqueBag, StepDocument, StepProductNode, UnitContext};
    use crate::step::export::StepExportOptions;
    use crate::step::export::{ExportContext, ref_, str_param};
    use crate::step::parser::entities::parse_step_file;
    use crate::step::parser::graph::EntityGraph;
    use crate::step::serializer::serialize_part21;

    #[test]
    fn serialized_rational_pcurve_roundtrip_preserves_weights() {
        let mut export = ExportContext::new();
        let surface_id = export.emit("PLANE", vec![str_param(""), ref_(999)]);
        let pcurve_id = export.export_pcurve(
            &NurbsCurve {
                degree: 2,
                knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
                control_points: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(0.5, 0.5, 0.0),
                    Point3::new(1.0, 0.0, 0.0),
                ],
                weights: vec![1.0, 0.5, 1.0],
            },
            surface_id,
        );
        let bytes = serialize_part21(
            &export.entities,
            &StepDocument {
                products: Vec::<StepProductNode>::new(),
                roots: Vec::new(),
                units: UnitContext::default(),
                opaque_bag: OpaqueBag::default(),
            },
            &StepExportOptions::default(),
        )
        .expect("serialized pcurve STEP should write");
        let parsed = parse_step_file(&bytes).expect("serialized pcurve STEP should parse");
        let mut log = DiagnosticLog::new();
        let graph = EntityGraph::build(parsed.entities, &mut log);
        let mut geometry = GeometryCache::new(&graph);

        let resolved = geometry
            .resolve_pcurve_on_surface(pcurve_id, surface_id, &mut log)
            .expect("pcurve should resolve on its exported surface");

        assert_eq!(resolved.pcurve_id, pcurve_id);
        assert_ne!(resolved.reference_curve_id, resolved.pcurve_id);

        match resolved.pcurve.geometry {
            popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(nurbs) => {
                assert_eq!(nurbs.weights, vec![1.0, 0.5, 1.0]);
            }
            other => panic!("expected NURBS pcurve, got {other:?}"),
        }
    }
}
