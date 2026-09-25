//! Export pcurve and surface-curve entities to STEP.
//!
//! Creates PCURVE, DEFINITIONAL_REPRESENTATION, and SURFACE_CURVE entities
//! that associate 2-D parametric trim curves with face surfaces, enabling
//! OpenCASCADE-based viewers (f3d, FreeCAD) to properly trim faces.

use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use crate::step::export::{ExportContext, enum_, list, ref_, str_param};
use crate::step::parser::types::{EntityId, Parameter};

impl ExportContext {
    /// Lazily emit the shared 2-D parametric representation context.
    ///
    /// STEP pcurves live in a `DEFINITIONAL_REPRESENTATION` that references a
    /// `( PARAMETRIC_REPRESENTATION_CONTEXT() REPRESENTATION_CONTEXT(...) )`
    /// complex entity.  We emit it once and cache the ID.
    pub(crate) fn ensure_param_repr_context(&mut self) -> EntityId {
        if let Some(id) = self.param_repr_context_id {
            return id;
        }
        let id = self.alloc_id();
        self.entities
            .push(crate::step::parser::types::EntityInstance {
                id,
                // Alphabetical complex-instance order per Part 21;
                // readers (OpenCASCADE) require the context to carry
                // GEOMETRIC_REPRESENTATION_CONTEXT — a bare
                // parametric context is "illegal type" for
                // DEFINITIONAL_REPRESENTATION.context_of_items.
                type_name: "GEOMETRIC_REPRESENTATION_CONTEXT".into(),
                params: vec![crate::step::export::int(2)],
                complex_parts: vec![
                    ("PARAMETRIC_REPRESENTATION_CONTEXT".into(), vec![]),
                    (
                        "REPRESENTATION_CONTEXT".into(),
                        vec![str_param("2D SPACE"), str_param("")],
                    ),
                ],
            });
        self.param_repr_context_id = Some(id);
        id
    }

    /// Export a pcurve NURBS as a PCURVE entity on the given surface.
    ///
    /// Emits:
    /// 1. A 2-D B_SPLINE_CURVE_WITH_KNOTS (control points have z = 0)
    /// 2. A DEFINITIONAL_REPRESENTATION wrapping the 2-D curve
    /// 3. A PCURVE referencing the surface and definitional representation
    pub(crate) fn export_pcurve(&mut self, pcurve: &NurbsCurve, surface_id: EntityId) -> EntityId {
        // The pcurve NURBS uses Point3(u, v, 0) convention — export as-is
        // since STEP AP203/214 uses 3-D CARTESIAN_POINT for pcurve control
        // points (z coordinate is 0).
        let curve_2d_id = self.export_nurbs_curve(pcurve);

        // Parametric representation context (shared, emitted once)
        let ctx_id = self.ensure_param_repr_context();

        // DEFINITIONAL_REPRESENTATION('', (#curve_2d), #context)
        let def_repr_id = self.emit(
            "DEFINITIONAL_REPRESENTATION",
            vec![str_param(""), list(vec![ref_(curve_2d_id)]), ref_(ctx_id)],
        );

        // PCURVE('', #surface, #def_repr)
        self.emit(
            "PCURVE",
            vec![str_param(""), ref_(surface_id), ref_(def_repr_id)],
        )
    }

    /// Create a SURFACE_CURVE wrapping a 3-D curve with associated pcurves.
    ///
    /// Returns a SURFACE_CURVE entity ID that can be used in place of the
    /// raw 3-D curve ID in an EDGE_CURVE.
    pub(crate) fn export_surface_curve(
        &mut self,
        curve_3d_id: EntityId,
        pcurve_ids: &[EntityId],
    ) -> EntityId {
        let pcurve_refs: Vec<Parameter> = pcurve_ids.iter().map(|&id| ref_(id)).collect();
        self.emit(
            "SURFACE_CURVE",
            vec![
                str_param(""),
                ref_(curve_3d_id),
                list(pcurve_refs),
                enum_("PCURVE_S1"),
            ],
        )
    }

    /// Create a SEAM_CURVE wrapping a 3-D curve with associated pcurves.
    ///
    /// SEAM_CURVE has the same structure as SURFACE_CURVE but indicates that
    /// both half-edges belong to the same face (e.g., a cylinder or torus seam).
    /// OpenCASCADE handles seam topology correctly with this entity type.
    pub(crate) fn export_seam_curve(
        &mut self,
        curve_3d_id: EntityId,
        pcurve_ids: &[EntityId],
    ) -> EntityId {
        let pcurve_refs: Vec<Parameter> = pcurve_ids.iter().map(|&id| ref_(id)).collect();
        self.emit(
            "SEAM_CURVE",
            vec![
                str_param(""),
                ref_(curve_3d_id),
                list(pcurve_refs),
                enum_("PCURVE_S1"),
            ],
        )
    }
}

/// Fraction of the model's linear tolerance an exported analytic
/// pcurve may deviate from its source, in model space: the written
/// trim must be indistinguishable from the kernel's at the model's
/// own accuracy (PK-550: the former fixed 17-sample cubic opened real
/// seam gaps of 7e-6 to 1.2e-2 on re-import).
const ANALYTIC_PCURVE_EXPORT_FRACTION: f64 = 1e-3;
/// Pcurve samples used to estimate the surface's uv metric scale.
const METRIC_SCALE_SAMPLES: usize = 9;

/// The NURBS pcurve written for `pcurve` on `surface`: exact for NURBS
/// and linear forms; analytic arcs and sinusoids are fitted
/// co-parameterized to within `ANALYTIC_PCURVE_EXPORT_FRACTION` of
/// `linear_tolerance`, converted to uv through the surface's metric
/// scale along the pcurve (1 on a unit-axis plane, the radius for a
/// cylinder's angular u).
pub(crate) fn export_pcurve_nurbs(
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
    linear_tolerance: f64,
) -> NurbsCurve {
    let uv_tolerance =
        ANALYTIC_PCURVE_EXPORT_FRACTION * linear_tolerance / uv_metric_scale(pcurve, surface);
    pcurve.to_nurbs_within(uv_tolerance).0
}

/// Largest `|S_u|` or `|S_v|` along the pcurve: a uv displacement `d`
/// moves the surface point by at most about `scale · d`. Falls back to
/// 1 for a degenerate or non-finite metric.
fn uv_metric_scale(pcurve: &Pcurve, surface: &SurfaceGeometry) -> f64 {
    let mut scale = 0.0_f64;
    for k in 0..METRIC_SCALE_SAMPLES {
        let s = k as f64 / (METRIC_SCALE_SAMPLES - 1) as f64;
        let uv = pcurve.point_at(pcurve.t_start + s * (pcurve.t_end - pcurve.t_start));
        let su = surface.partial_u(uv.x, uv.y).length();
        let sv = surface.partial_v(uv.x, uv.y).length();
        scale = scale.max(su).max(sv);
    }
    if scale.is_finite() && scale > 0.0 {
        scale
    } else {
        1.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;

    #[test]
    fn export_pcurve_creates_three_entities() {
        let mut ctx = ExportContext::new();

        // Create a dummy surface
        let surface_id = ctx.emit("PLANE", vec![str_param(""), ref_(1)]);

        let pcurve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        };

        let pcurve_id = ctx.export_pcurve(&pcurve, surface_id);

        // Should have created: PLANE, CARTESIAN_POINTs, B_SPLINE_CURVE_WITH_KNOTS,
        // PARAMETRIC_REPRESENTATION_CONTEXT, DEFINITIONAL_REPRESENTATION, PCURVE
        let entity = ctx.entities.iter().find(|e| e.id == pcurve_id);
        assert!(entity.is_some());
        assert_eq!(
            entity.as_ref().map(|e| e.type_name.as_str()),
            Some("PCURVE")
        );

        // Verify DEFINITIONAL_REPRESENTATION exists
        let def_repr = ctx
            .entities
            .iter()
            .find(|e| e.type_name == "DEFINITIONAL_REPRESENTATION");
        assert!(def_repr.is_some());

        // Verify the 2D context exists (a complex instance led by
        // GEOMETRIC_REPRESENTATION_CONTEXT with a
        // PARAMETRIC_REPRESENTATION_CONTEXT component).
        let ctx_ent = ctx
            .entities
            .iter()
            .find(|e| e.type_name == "GEOMETRIC_REPRESENTATION_CONTEXT");
        assert!(ctx_ent.is_some());
        assert!(
            ctx_ent
                .unwrap()
                .complex_parts
                .iter()
                .any(|(t, _)| t == "PARAMETRIC_REPRESENTATION_CONTEXT")
        );
    }

    #[test]
    fn parametric_context_is_shared() {
        let mut ctx = ExportContext::new();
        let id1 = ctx.ensure_param_repr_context();
        let id2 = ctx.ensure_param_repr_context();
        assert_eq!(id1, id2, "Context should be cached and reused");
    }

    #[test]
    fn export_surface_curve_wraps_pcurves() {
        let mut ctx = ExportContext::new();
        let curve_3d = ctx.emit("LINE", vec![str_param("")]);
        let pcurve1 = ctx.emit("PCURVE", vec![str_param("")]);
        let pcurve2 = ctx.emit("PCURVE", vec![str_param("")]);

        let sc_id = ctx.export_surface_curve(curve_3d, &[pcurve1, pcurve2]);

        let entity = ctx.entities.iter().find(|e| e.id == sc_id);
        assert!(entity.is_some());
        assert_eq!(
            entity.as_ref().map(|e| e.type_name.as_str()),
            Some("SURFACE_CURVE")
        );

        // Verify it references the 3D curve and pcurves
        let params = &entity.as_ref().map(|e| &e.params);
        assert!(params.is_some());
        let params = params.as_ref().map(|p| p.len());
        assert_eq!(params, Some(4));
    }

    #[test]
    fn export_seam_curve_wraps_pcurves() {
        let mut ctx = ExportContext::new();
        let curve_3d = ctx.emit("LINE", vec![str_param("")]);
        let pcurve1 = ctx.emit("PCURVE", vec![str_param("")]);
        let pcurve2 = ctx.emit("PCURVE", vec![str_param("")]);

        let sc_id = ctx.export_seam_curve(curve_3d, &[pcurve1, pcurve2]);

        let entity = ctx.entities.iter().find(|e| e.id == sc_id);
        assert!(entity.is_some());
        assert_eq!(
            entity.as_ref().map(|e| e.type_name.as_str()),
            Some("SEAM_CURVE")
        );

        // Verify same 4-parameter structure as SURFACE_CURVE
        let params = &entity.as_ref().map(|e| &e.params);
        assert!(params.is_some());
        let params = params.as_ref().map(|p| p.len());
        assert_eq!(params, Some(4));
    }
}
