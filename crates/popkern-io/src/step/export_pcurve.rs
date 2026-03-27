//! Export pcurve and surface-curve entities to STEP.
//!
//! Creates PCURVE, DEFINITIONAL_REPRESENTATION, and SURFACE_CURVE entities
//! that associate 2-D parametric trim curves with face surfaces, enabling
//! OpenCASCADE-based viewers (f3d, FreeCAD) to properly trim faces.

use popkern_core::geometry::nurbs_curve::NurbsCurve;

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
                type_name: "PARAMETRIC_REPRESENTATION_CONTEXT".into(),
                params: vec![],
                complex_parts: vec![(
                    "REPRESENTATION_CONTEXT".into(),
                    vec![str_param("2D SPACE"), str_param("")],
                )],
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

        // Verify parametric context exists
        let ctx_ent = ctx
            .entities
            .iter()
            .find(|e| e.type_name == "PARAMETRIC_REPRESENTATION_CONTEXT");
        assert!(ctx_ent.is_some());
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
