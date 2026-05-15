use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::pcurve::Pcurve;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::mapping::geometry::GeometryCache;
use crate::step::parser::types::EntityId;
use crate::step::schema::geometry::{SeamCurve, StepPcurve, SurfaceCurve};
use crate::step::schema::{as_entity_ref, as_list, param_at};

impl GeometryCache<'_> {
    /// Resolve the pcurve associated with a surface curve on `surface_id`.
    pub fn resolve_surface_curve_pcurve(
        &mut self,
        curve_id: EntityId,
        surface_id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<Pcurve> {
        self.resolve_surface_curve_pcurves(curve_id, surface_id, log)
            .into_iter()
            .next()
    }

    /// Resolve every pcurve associated with a surface curve on `surface_id`.
    ///
    /// STEP permits multiple associated pcurves. Import callers that know the
    /// owning half-edge should evaluate all returned candidates against the
    /// recovered edge endpoints rather than accepting the first match.
    pub fn resolve_surface_curve_pcurves(
        &mut self,
        curve_id: EntityId,
        surface_id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Vec<Pcurve> {
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
    ) -> Option<Pcurve> {
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
        match curve.geometry {
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
        }
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
}
