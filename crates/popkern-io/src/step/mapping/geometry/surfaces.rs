//! Surface resolution helpers for [`GeometryCache`].

use popkern_core::geometry::Vec3;
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_math::PortableMath;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::parser::types::EntityId;
use crate::step::schema::geometry::{
    BSplineSurfaceWithKnots, ConicalSurface, CylindricalSurface, SphericalSurface, StepPlane,
    ToroidalSurface,
};
use crate::step::schema::{as_list, as_real_list};

use super::{GeometryCache, expand_knots};

impl<'a> GeometryCache<'a> {
    pub(super) fn resolve_plane(
        &mut self,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<SurfaceGeometry> {
        let entity = self.graph.get(id)?;
        let plane = match StepPlane::from_params(&entity.params) {
            Ok(p) => p,
            Err(e) => {
                log.error(Phase::MapGeometry, Some(id), format!("{e}"));
                return None;
            }
        };
        let (origin, z_axis, x_axis) = self.resolve_placement(plane.position, log)?;
        let y_axis = z_axis.cross(&x_axis);
        let v_axis = y_axis.normalized().unwrap_or(Vec3::new(0.0, 1.0, 0.0));
        let u_axis = x_axis.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
        Some(SurfaceGeometry::Plane {
            origin,
            normal: z_axis,
            u_axis,
            v_axis,
        })
    }

    pub(super) fn resolve_cylindrical_surface(
        &mut self,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<SurfaceGeometry> {
        let entity = self.graph.get(id)?;
        let cyl = match CylindricalSurface::from_params(&entity.params) {
            Ok(c) => c,
            Err(e) => {
                log.error(Phase::MapGeometry, Some(id), format!("{e}"));
                return None;
            }
        };
        let (origin, axis, _x) = self.resolve_placement(cyl.position, log)?;
        Some(SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius: cyl.radius,
        })
    }

    pub(super) fn resolve_conical_surface(
        &mut self,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<SurfaceGeometry> {
        let entity = self.graph.get(id)?;
        let cone = match ConicalSurface::from_params(&entity.params) {
            Ok(c) => c,
            Err(e) => {
                log.error(Phase::MapGeometry, Some(id), format!("{e}"));
                return None;
            }
        };
        let (origin, axis, _x) = self.resolve_placement(cone.position, log)?;
        let half_angle = cone.semi_angle * self.radians_per_angle;
        if !(half_angle > 0.0 && half_angle < std::f64::consts::FRAC_PI_2) {
            log.error(
                Phase::MapGeometry,
                Some(id),
                format!("cone semi-angle {half_angle} rad is outside (0, π/2)"),
            );
            return None;
        }
        // STEP places a cone where its radius is R; popkern's apex is
        // R / tan α back along the axis (measured: drill-point cones
        // with R = 4e-5 sat 2.4e-5 off their own edges, and chamfer or
        // countersink cones with a real R were displaced by it).
        Some(SurfaceGeometry::Cone {
            apex: origin - axis * cone_apex_distance(cone.radius, half_angle),
            axis,
            half_angle,
        })
    }

    pub(super) fn resolve_spherical_surface(
        &mut self,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<SurfaceGeometry> {
        let entity = self.graph.get(id)?;
        let sph = match SphericalSurface::from_params(&entity.params) {
            Ok(s) => s,
            Err(e) => {
                log.error(Phase::MapGeometry, Some(id), format!("{e}"));
                return None;
            }
        };
        // The placement's axis is the chart pole: dropping it put a
        // tilted sphere's poles on or inside its trimmed faces.
        let (origin, axis, _x) = self.resolve_placement(sph.position, log)?;
        Some(SurfaceGeometry::Sphere {
            center: origin,
            radius: sph.radius,
            axis,
        })
    }

    pub(super) fn resolve_toroidal_surface(
        &mut self,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<SurfaceGeometry> {
        let entity = self.graph.get(id)?;
        let tor = match ToroidalSurface::from_params(&entity.params) {
            Ok(t) => t,
            Err(e) => {
                log.error(Phase::MapGeometry, Some(id), format!("{e}"));
                return None;
            }
        };
        let (origin, axis, _x) = self.resolve_placement(tor.position, log)?;
        Some(SurfaceGeometry::Torus {
            center: origin,
            axis,
            major_radius: tor.major_radius,
            minor_radius: tor.minor_radius,
        })
    }

    pub(super) fn resolve_bspline_surface(
        &mut self,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<SurfaceGeometry> {
        let params = self.graph.get(id)?.params.clone();
        self.resolve_bspline_surface_from(id, &params, log)
    }

    /// Resolve the standard Part-21 complex form for rational B-spline
    /// surfaces: leading partial `BOUNDED_SURFACE`, spline data
    /// distributed over the `B_SPLINE_SURFACE` /
    /// `B_SPLINE_SURFACE_WITH_KNOTS` partials, weights on
    /// `RATIONAL_B_SPLINE_SURFACE` (the form the exporter writes and
    /// OCC itself emits — PK-155).
    pub(super) fn resolve_bspline_surface_complex(
        &mut self,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<SurfaceGeometry> {
        let entity = self.graph.get(id)?;
        let part = |name: &str| {
            entity
                .complex_parts
                .iter()
                .find(|(n, _)| n == name)
                .map(|(_, p)| p.clone())
        };
        let bs = part("B_SPLINE_SURFACE")?;
        let bk = part("B_SPLINE_SURFACE_WITH_KNOTS")?;
        if bs.len() < 7 || bk.len() < 5 {
            log.error(
                Phase::MapGeometry,
                Some(id),
                format!(
                    "complex B-spline surface: expected 7 B_SPLINE_SURFACE and 5 \
                     B_SPLINE_SURFACE_WITH_KNOTS parameters, got {} and {}",
                    bs.len(),
                    bk.len()
                ),
            );
            return None;
        }
        // Reassemble the monolithic parameter layout the schema
        // parser expects: name, degrees, control points, form, flags,
        // multiplicities, knots, knot spec.
        let mut params = Vec::with_capacity(13);
        params.push(crate::step::parser::types::Parameter::String(String::new()));
        params.extend(bs.iter().take(7).cloned());
        params.extend(bk.iter().take(5).cloned());
        self.resolve_bspline_surface_from(id, &params, log)
    }

    pub(super) fn resolve_bspline_surface_from(
        &mut self,
        id: EntityId,
        params: &[crate::step::parser::types::Parameter],
        log: &mut DiagnosticLog,
    ) -> Option<SurfaceGeometry> {
        let entity = self.graph.get(id)?;
        let bs = match BSplineSurfaceWithKnots::from_params(params) {
            Ok(b) => b,
            Err(e) => {
                log.error(Phase::MapGeometry, Some(id), format!("{e}"));
                return None;
            }
        };

        let num_u = bs.control_points.len();
        let num_v = if num_u > 0 {
            bs.control_points[0].len()
        } else {
            0
        };

        let mut control_points = Vec::with_capacity(num_u * num_v);
        for row in &bs.control_points {
            for &cp_id in row {
                control_points.push(self.resolve_point(cp_id, log)?);
            }
        }

        let knots_u = expand_knots(&bs.u_knots, &bs.u_multiplicities);
        let knots_v = expand_knots(&bs.v_knots, &bs.v_multiplicities);

        let mut weights = resolve_rational_bspline_surface_weights(entity, log)
            .or_else(|| bs.weights.map(|w| w.iter().flatten().copied().collect()))
            .unwrap_or_else(|| vec![1.0; control_points.len()]);
        if weights.len() != control_points.len() {
            log.warn(
                Phase::MapGeometry,
                Some(id),
                format!(
                    "rational B-spline surface weight count {} does not match control point count {}; using unit weights",
                    weights.len(),
                    control_points.len()
                ),
            );
            weights = vec![1.0; control_points.len()];
        }

        Some(SurfaceGeometry::Nurbs(NurbsSurface {
            degree_u: bs.u_degree as usize,
            degree_v: bs.v_degree as usize,
            knots_u,
            knots_v,
            control_points,
            num_u,
            num_v,
            weights,
        }))
    }
}

fn resolve_rational_bspline_surface_weights(
    entity: &crate::step::parser::types::EntityInstance,
    log: &mut DiagnosticLog,
) -> Option<Vec<f64>> {
    entity
        .complex_parts
        .iter()
        .find(|(name, _)| name == "RATIONAL_B_SPLINE_SURFACE")
        .and_then(|(_, params)| params.first())
        .and_then(|weight_grid| {
            let rows = as_list(weight_grid, "RATIONAL_B_SPLINE_SURFACE")
                .map_err(|e| {
                    log.warn(Phase::MapGeometry, Some(entity.id), format!("{e}"));
                    e
                })
                .ok()?;
            let mut weights = Vec::new();
            for row in rows {
                weights.extend(
                    as_real_list(row, "RATIONAL_B_SPLINE_SURFACE")
                        .map_err(|e| {
                            log.warn(Phase::MapGeometry, Some(entity.id), format!("{e}"));
                            e
                        })
                        .ok()?,
                );
            }
            Some(weights)
        })
}

/// The distance along a STEP cone's axis from its placement (where its
/// radius is `radius`) back to the apex.
pub(super) fn cone_apex_distance(radius: f64, half_angle: f64) -> f64 {
    radius / half_angle.ptan()
}
