//! Cone constructor.

use std::f64::consts::TAU;

use crate::error::CoreError;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

impl BRep {
    /// Build a cone with 3-lune topology.
    ///
    /// The cone has its apex at (0, 0, height) and its base circle centered
    /// at the origin with the given `radius`. The two vertices are the apex
    /// and the base center (0, 0, 0).
    ///
    /// # Topology
    ///
    /// V=2, E=3, F=3, L=3, S=1.
    ///
    /// # Geometry
    ///
    /// All three faces are assigned `SurfaceGeometry::Cone` with the apex,
    /// downward axis, and half-angle derived from `radius` and `height`.
    /// Edges have line curves (seam lines from base center to apex).
    /// Each half-edge carries a pcurve (degree-1 NURBS line in UV space)
    /// representing a generator line at constant u, running from v=0 (apex)
    /// to v=height (base) or vice versa.
    ///
    /// # Errors
    ///
    /// Returns [`CoreError`] if `radius` or `height` is not positive,
    /// or if any Euler operation or validation fails.
    pub fn cone(radius: f64, height: f64, tolerance: Tolerance) -> Result<Self, CoreError> {
        if radius <= 0.0 {
            return Err(CoreError::PreconditionFailed(
                "cone radius must be positive".into(),
            ));
        }
        if height <= 0.0 {
            return Err(CoreError::PreconditionFailed(
                "cone height must be positive".into(),
            ));
        }

        let mut brep = BRep::new(tolerance);
        let base = Point3::new(0.0, 0.0, 0.0);
        let apex = Point3::new(0.0, 0.0, height);

        // Build 3-lune topology (identical structure to sphere).
        let mvfs = brep.mvfs(base)?;
        let mev1 = brep.mev_at_vertex(mvfs.vertex, mvfs.loop_, apex)?;
        let mef1 = brep.mef(mev1.half_edge, mev1.half_edge_twin)?;
        let he_in_new_loop = mef1.half_edge_twin;
        let next_he = brep.half_edges[he_in_new_loop].next;
        let _mef2 = brep.mef(he_in_new_loop, next_he)?;

        // Attach line curves to all edges.
        brep.attach_line_curves();

        // Attach cone surface geometry to every face.
        // The axis points from apex toward the base (downward along -Z).
        // The half-angle is the angle between the axis and the surface generator.
        let half_angle = (radius / height).atan();
        let cone_geom = SurfaceGeometry::Cone {
            apex,
            axis: Vec3::new(0.0, 0.0, -1.0),
            half_angle,
        };
        let face_handles: Vec<_> = brep.faces.keys().collect();
        for fh in &face_handles {
            brep.faces[*fh].surface = Some(SurfaceOnFace {
                geometry: cone_geom.clone(),
                orientation: true,
            });
        }

        // Attach pcurves to all half-edges.
        // Cone parameterization with axis=(0,0,-1), apex=(0,0,h):
        //   point = apex + axis*v + perp1*(v*tan(ha)*cos(u)) + perp2*(v*tan(ha)*sin(u))
        //   => z = h - v, so v = h - z.  At apex (z=h): v=0. At base (z=0): v=h.
        // Each lune face has 2 half-edges. For face i:
        //   base->apex half-edge: u=u_left, from (u_left, h) to (u_left, 0)
        //   apex->base half-edge: u=u_right, from (u_right, 0) to (u_right, h)
        let base_vertex = mvfs.vertex;
        for (face_idx, &fh) in face_handles.iter().enumerate() {
            let u_left = face_idx as f64 * TAU / 3.0;
            let u_right = (face_idx as f64 + 1.0) * TAU / 3.0;

            let outer_loop = brep.faces[fh].outer_loop;
            let loop_data = brep
                .loops
                .get(outer_loop)
                .ok_or_else(|| CoreError::InvalidHandle("face outer loop invalid".into()))?;
            let start_he = loop_data.half_edge;

            if !brep.half_edges.contains_key(start_he) {
                continue;
            }

            // Walk the loop (exactly 2 half-edges per lune face).
            let mut he_handles = vec![start_he];
            let mut current = brep.half_edges[start_he].next;
            while current != start_he {
                he_handles.push(current);
                current = brep.half_edges[current].next;
            }

            for &heh in &he_handles {
                let origin = brep.half_edges[heh].origin;
                let is_from_base = origin == base_vertex;

                let pcurve = if is_from_base {
                    // Base->Apex: (u_left, h) -> (u_left, 0)
                    NurbsCurve {
                        degree: 1,
                        knots: vec![0.0, 0.0, 1.0, 1.0],
                        control_points: vec![
                            Point3::new(u_left, height, 0.0),
                            Point3::new(u_left, 0.0, 0.0),
                        ],
                        weights: vec![1.0, 1.0],
                    }
                } else {
                    // Apex->Base: (u_right, 0) -> (u_right, h)
                    NurbsCurve {
                        degree: 1,
                        knots: vec![0.0, 0.0, 1.0, 1.0],
                        control_points: vec![
                            Point3::new(u_right, 0.0, 0.0),
                            Point3::new(u_right, height, 0.0),
                        ],
                        weights: vec![1.0, 1.0],
                    }
                };

                brep.half_edges[heh].pcurve = Some(pcurve);
            }
        }

        // Validate the result.
        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::surface::SurfaceGeometry;
    use crate::tolerance::Tolerance;

    #[test]
    fn cone_topology() {
        let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
        assert_eq!(cone.vertices.len(), 2);
        assert_eq!(cone.edges.len(), 3);
        assert_eq!(cone.faces.len(), 3);
    }

    #[test]
    fn cone_passes_validation() {
        let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
        assert!(cone.validate().is_ok(), "Cone should pass validation");
    }

    #[test]
    fn cone_has_edge_geometry() {
        let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
        for edge in cone.edges.values() {
            assert!(
                edge.curve.is_some(),
                "Every edge should have curve geometry"
            );
        }
    }

    #[test]
    fn cone_vertex_positions() {
        let cone = BRep::cone(1.0, 3.0, Tolerance::default()).unwrap();
        let positions: Vec<_> = cone.vertices.values().map(|v| v.position).collect();
        let has_base = positions
            .iter()
            .any(|p| p.z.abs() < 1e-10 && p.x.abs() < 1e-10 && p.y.abs() < 1e-10);
        let has_apex = positions
            .iter()
            .any(|p| (p.z - 3.0).abs() < 1e-10 && p.x.abs() < 1e-10 && p.y.abs() < 1e-10);
        assert!(has_base, "Should have base vertex");
        assert!(has_apex, "Should have apex vertex");
    }

    #[test]
    fn cone_invalid_params() {
        assert!(BRep::cone(0.0, 1.0, Tolerance::default()).is_err());
        assert!(BRep::cone(-1.0, 1.0, Tolerance::default()).is_err());
        assert!(BRep::cone(1.0, 0.0, Tolerance::default()).is_err());
        assert!(BRep::cone(1.0, -1.0, Tolerance::default()).is_err());
    }

    #[test]
    fn cone_all_faces_have_surface_geometry() {
        let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
        for face in cone.faces.values() {
            assert!(
                face.surface.is_some(),
                "every cone face should have surface geometry"
            );
        }
    }

    #[test]
    fn cone_faces_have_cone_surface_geometry() {
        let cone = BRep::cone(1.0, 2.0, Tolerance::default()).unwrap();
        let expected_half_angle = (1.0_f64 / 2.0_f64).atan();
        for face in cone.faces.values() {
            let surface = face
                .surface
                .as_ref()
                .expect("cone face should have surface geometry");
            match &surface.geometry {
                SurfaceGeometry::Cone {
                    apex,
                    axis,
                    half_angle,
                } => {
                    assert!(
                        (apex.z - 2.0).abs() < 1e-10,
                        "cone apex z should be at height"
                    );
                    assert!(
                        apex.x.abs() < 1e-10 && apex.y.abs() < 1e-10,
                        "cone apex should be on z-axis"
                    );
                    assert!(
                        (axis.z + 1.0).abs() < 1e-10
                            && axis.x.abs() < 1e-10
                            && axis.y.abs() < 1e-10,
                        "cone axis should point in -Z direction"
                    );
                    assert!(
                        (*half_angle - expected_half_angle).abs() < 1e-10,
                        "cone half-angle should match atan(radius/height)"
                    );
                }
                other => panic!(
                    "expected SurfaceGeometry::Cone, got {:?}",
                    std::mem::discriminant(other)
                ),
            }
        }
    }

    #[test]
    fn cone_surface_geometry_matches_params() {
        // Test with different radius/height to ensure half_angle is computed correctly.
        let radius = 3.0_f64;
        let height = 4.0_f64;
        let cone = BRep::cone(radius, height, Tolerance::default()).unwrap();
        let expected_half_angle = (radius / height).atan();
        for face in cone.faces.values() {
            if let Some(ref sof) = face.surface {
                if let SurfaceGeometry::Cone { half_angle, .. } = sof.geometry {
                    assert!(
                        (half_angle - expected_half_angle).abs() < 1e-10,
                        "half_angle mismatch: got {half_angle}, expected {expected_half_angle}"
                    );
                }
            }
        }
    }
}
