//! Complex-form rational B-spline curve resolution.
//!
//! Split from `mod.rs` (structural line ratchet). The standard
//! Part-21 complex instance distributes spline data over the
//! B_SPLINE_CURVE / B_SPLINE_CURVE_WITH_KNOTS partials with weights
//! on RATIONAL_B_SPLINE_CURVE.

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::nurbs_curve::NurbsCurve;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::parser::types::EntityId;
use crate::step::schema::as_real_list;

use super::{GeometryCache, expand_knots};

impl<'a> GeometryCache<'a> {
    pub(super) fn resolve_bspline_curve_complex(
        &mut self,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<CurveOnEdge> {
        let entity = self.graph.get(id)?.clone();
        let part = |name: &str| {
            entity
                .complex_parts
                .iter()
                .find(|(n, _)| n == name)
                .map(|(_, p)| p.as_slice())
        };
        let bc = part("B_SPLINE_CURVE")?;
        let bk = part("B_SPLINE_CURVE_WITH_KNOTS")?;
        let degree = match bc.first() {
            Some(crate::step::parser::types::Parameter::Int(d)) => *d as usize,
            _ => {
                log.error(
                    Phase::MapGeometry,
                    Some(id),
                    "B_SPLINE_CURVE partial: missing degree".to_string(),
                );
                return None;
            }
        };
        let cp_ids = crate::step::schema::as_ref_list(bc.get(1)?, "B_SPLINE_CURVE control points")
            .map_err(|e| log.error(Phase::MapGeometry, Some(id), format!("{e}")))
            .ok()?;
        let mut control_points = Vec::with_capacity(cp_ids.len());
        for cp_id in cp_ids {
            control_points.push(self.resolve_point(cp_id, log)?);
        }
        let mults = crate::step::schema::as_int_list(bk.first()?, "knot multiplicities")
            .map_err(|e| log.error(Phase::MapGeometry, Some(id), format!("{e}")))
            .ok()?;
        let knot_vals = as_real_list(bk.get(1)?, "knot values")
            .map_err(|e| log.error(Phase::MapGeometry, Some(id), format!("{e}")))
            .ok()?;
        let knots = expand_knots(&knot_vals, &mults);
        let mut weights = resolve_rational_bspline_curve_weights(&entity, log)
            .unwrap_or_else(|| vec![1.0; control_points.len()]);
        if weights.len() != control_points.len() {
            log.warn(
                Phase::MapGeometry,
                Some(id),
                format!(
                    "weight count {} != control point count {}; resetting to 1",
                    weights.len(),
                    control_points.len()
                ),
            );
            weights = vec![1.0; control_points.len()];
        }
        let t_start = knots.first().copied().unwrap_or(0.0);
        let t_end = knots.last().copied().unwrap_or(1.0);
        Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(NurbsCurve {
                degree,
                knots,
                control_points,
                weights,
            }),
            t_start,
            t_end,
        })
    }
}

pub(super) fn resolve_rational_bspline_curve_weights(
    entity: &crate::step::parser::types::EntityInstance,
    log: &mut DiagnosticLog,
) -> Option<Vec<f64>> {
    entity
        .complex_parts
        .iter()
        .find(|(name, _)| name == "RATIONAL_B_SPLINE_CURVE")
        .and_then(|(_, params)| params.first())
        .and_then(|weights| {
            as_real_list(weights, "RATIONAL_B_SPLINE_CURVE")
                .map_err(|e| {
                    log.warn(Phase::MapGeometry, Some(entity.id), format!("{e}"));
                    e
                })
                .ok()
        })
}
