//! Surface resolution helpers for [`GeometryCache`].

use popkern_core::geometry::Vec3;
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::surface::SurfaceGeometry;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::parser::types::EntityId;
use crate::step::schema::geometry::{
    BSplineSurfaceWithKnots, ConicalSurface, CylindricalSurface, SphericalSurface, StepPlane,
    ToroidalSurface,
};

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
        Some(SurfaceGeometry::Cone {
            apex: origin,
            axis,
            half_angle: cone.semi_angle,
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
        let (origin, _axis, _x) = self.resolve_placement(sph.position, log)?;
        Some(SurfaceGeometry::Sphere {
            center: origin,
            radius: sph.radius,
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
        let entity = self.graph.get(id)?;
        let bs = match BSplineSurfaceWithKnots::from_params(&entity.params) {
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

        let weights = match &bs.weights {
            Some(w) => w.iter().flatten().copied().collect(),
            None => vec![1.0; control_points.len()],
        };

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
