//! Export popkern geometry types to STEP entity instances.
//!
//! Converts [`Point3`], [`Vec3`], [`CurveOnEdge`], and [`SurfaceOnFace`] into
//! the corresponding STEP geometric entities (CARTESIAN_POINT, DIRECTION,
//! LINE, PLANE, CYLINDRICAL_SURFACE, etc.).

use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveOnEdge, NurbsCurve, Point3, SurfaceOnFace, Vec3};

/// True when the STEP plane placement written by `export_surface`
/// (ẑ = u_axis × v_axis, chart-true) OPPOSES the stored popkern
/// normal — the ADVANCED_FACE same_sense flag must flip to keep the
/// face's material side unchanged.
pub(crate) fn step_plane_sense_flip(surface: &SurfaceGeometry) -> bool {
    if let SurfaceGeometry::Plane {
        normal,
        u_axis,
        v_axis,
        ..
    } = surface
    {
        u_axis.cross(v_axis).dot(normal) < 0.0
    } else {
        false
    }
}

/// Convert a pcurve (already in NURBS form) from popkern's chart
/// convention to the STEP parameterization of the written surface.
/// Spheres differ: popkern v is COLATITUDE from the +Z pole,
/// STEP SPHERICAL_SURFACE v is LATITUDE from the equator
/// (v_step = π/2 − v). Plane/cylinder/torus charts match the
/// written placements exactly (plane via the chart-true ẑ above;
/// cylinder/torus via ref_direction = the kernel chart's u origin).
pub(crate) fn pcurve_to_step_chart(surface: &SurfaceGeometry, mut nurbs: NurbsCurve) -> NurbsCurve {
    if matches!(surface, SurfaceGeometry::Sphere { .. }) {
        for p in &mut nurbs.control_points {
            p.y = std::f64::consts::FRAC_PI_2 - p.y;
        }
    }
    nurbs
}

use crate::step::export::{ExportContext, PointKey, list, real, ref_, str_param};
use crate::step::export_nurbs::compute_perpendicular;
use crate::step::parser::types::EntityId;

impl ExportContext {
    /// Export a 3D point, deduplicating by exact float bits.
    pub(crate) fn export_point(&mut self, p: &Point3) -> EntityId {
        let key = PointKey::from_point(p);
        if let Some(&id) = self.point_cache.get(&key) {
            return id;
        }
        let id = self.emit(
            "CARTESIAN_POINT",
            vec![str_param(""), list(vec![real(p.x), real(p.y), real(p.z)])],
        );
        self.point_cache.insert(key, id);
        id
    }

    /// Export a 3D direction vector, deduplicating by exact float bits.
    pub(crate) fn export_direction(&mut self, v: &Vec3) -> EntityId {
        let key = PointKey::from_vec(v);
        if let Some(&id) = self.direction_cache.get(&key) {
            return id;
        }
        let id = self.emit(
            "DIRECTION",
            vec![str_param(""), list(vec![real(v.x), real(v.y), real(v.z)])],
        );
        self.direction_cache.insert(key, id);
        id
    }

    /// Export an AXIS2_PLACEMENT_3D from origin, Z axis, and X axis.
    pub(crate) fn export_placement(
        &mut self,
        origin: &Point3,
        z_axis: &Vec3,
        x_axis: &Vec3,
    ) -> EntityId {
        let pt = self.export_point(origin);
        let z = self.export_direction(z_axis);
        let x = self.export_direction(x_axis);
        self.emit(
            "AXIS2_PLACEMENT_3D",
            vec![str_param(""), ref_(pt), ref_(z), ref_(x)],
        )
    }

    /// Export a curve geometry to STEP entities.
    pub(crate) fn export_curve(&mut self, curve: &CurveOnEdge) -> EntityId {
        match &curve.geometry {
            CurveGeometry::Line { origin, direction } => {
                let pt = self.export_point(origin);
                let dir = self.export_direction(direction);
                let mag = direction.length();
                let vec_id = self.emit("VECTOR", vec![str_param(""), ref_(dir), real(mag)]);
                self.emit("LINE", vec![str_param(""), ref_(pt), ref_(vec_id)])
            }
            CurveGeometry::Circle {
                center,
                normal,
                radius,
            } => {
                let x_axis = compute_perpendicular(normal);
                let placement = self.export_placement(center, normal, &x_axis);
                self.emit(
                    "CIRCLE",
                    vec![str_param(""), ref_(placement), real(*radius)],
                )
            }
            CurveGeometry::Arc {
                center,
                normal,
                radius,
                ..
            } => {
                // Export as a full CIRCLE; trimming is handled by edge
                // parameter bounds in the topology layer.
                let x_axis = compute_perpendicular(normal);
                let placement = self.export_placement(center, normal, &x_axis);
                self.emit(
                    "CIRCLE",
                    vec![str_param(""), ref_(placement), real(*radius)],
                )
            }
            CurveGeometry::Ellipse { .. }
            | CurveGeometry::Parabola { .. }
            | CurveGeometry::Hyperbola { .. } => {
                // STEP has native ELLIPSE / PARABOLA / HYPERBOLA
                // entities but our exporter currently routes all
                // through the NURBS path. Parabolas are exactly
                // representable as a degree-2 Bézier; ellipses and
                // hyperbolas use a sample-and-fit / standard
                // 9-CP rational approximation.
                use popkern_core::geometry::CurveEval;
                let range = curve.geometry.parameter_range();
                let nurbs = curve.geometry.to_nurbs_in_range(range);
                self.export_nurbs_curve(&nurbs)
            }
            CurveGeometry::Nurbs(_) => {
                let nurbs = curve
                    .geometry
                    .to_nurbs_in_range((curve.t_start, curve.t_end));
                self.export_nurbs_curve(&nurbs)
            }
        }
    }

    /// Export a surface geometry to STEP entities.
    pub(crate) fn export_surface(&mut self, surface: &SurfaceOnFace) -> EntityId {
        match &surface.geometry {
            SurfaceGeometry::Plane {
                origin,
                normal,
                u_axis,
                v_axis,
            } => {
                // STEP readers derive the chart's v direction as
                // ẑ × x̂ — the stored v_axis is never written. A
                // left-handed popkern chart (v_axis = −(normal ×
                // u_axis)) exported with the stored normal as ẑ
                // REFLECTS every pcurve trim on the face (measured:
                // plate faces displaced by their full extent in
                // OpenCASCADE). Write ẑ = u × v so the derived chart
                // matches the pcurves exactly; the ADVANCED_FACE
                // same_sense flag carries the orientation difference
                // (see `step_plane_sense_flip`).
                let chart_normal = u_axis.cross(v_axis).normalized().unwrap_or_else(|| *normal);
                let placement = self.export_placement(origin, &chart_normal, u_axis);
                self.emit("PLANE", vec![str_param(""), ref_(placement)])
            }
            SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } => {
                let x_axis = compute_perpendicular(axis);
                let placement = self.export_placement(origin, axis, &x_axis);
                self.emit(
                    "CYLINDRICAL_SURFACE",
                    vec![str_param(""), ref_(placement), real(*radius)],
                )
            }
            SurfaceGeometry::Cone {
                apex,
                axis,
                half_angle,
            } => {
                let x_axis = compute_perpendicular(axis);
                let placement = self.export_placement(apex, axis, &x_axis);
                self.emit(
                    "CONICAL_SURFACE",
                    vec![
                        str_param(""),
                        ref_(placement),
                        real(0.0), // radius at apex
                        real(*half_angle),
                    ],
                )
            }
            SurfaceGeometry::Sphere { center, radius } => {
                let z = Vec3::new(0.0, 0.0, 1.0);
                let x = Vec3::new(1.0, 0.0, 0.0);
                let placement = self.export_placement(center, &z, &x);
                self.emit(
                    "SPHERICAL_SURFACE",
                    vec![str_param(""), ref_(placement), real(*radius)],
                )
            }
            SurfaceGeometry::Torus {
                center,
                axis,
                major_radius,
                minor_radius,
            } => {
                let x_axis = compute_perpendicular(axis);
                let placement = self.export_placement(center, axis, &x_axis);
                self.emit(
                    "TOROIDAL_SURFACE",
                    vec![
                        str_param(""),
                        ref_(placement),
                        real(*major_radius),
                        real(*minor_radius),
                    ],
                )
            }
            SurfaceGeometry::Nurbs(nurbs) => self.export_nurbs_surface(nurbs),
            // STEP AP203/242 has no direct Gregory entity. Convert
            // to a NURBS fit (via `to_nurbs()`) and export as
            // B_SPLINE_SURFACE. Loses the rational-inner-CP exact
            // form but preserves the geometry to tessellation
            // fidelity — the standard interchange strategy.
            SurfaceGeometry::Gregory(_) => {
                let nurbs = surface.geometry.to_nurbs();
                self.export_nurbs_surface(&nurbs)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::curve::CurveOnEdge;

    #[test]
    fn plane_sense_flip_detects_left_handed_chart() {
        let right = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        assert!(!step_plane_sense_flip(&right));
        let left = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(0.0, 1.0, 0.0),
            v_axis: Vec3::new(1.0, 0.0, 0.0),
        };
        assert!(step_plane_sense_flip(&left));
    }

    /// The written plane placement must satisfy the STEP reader's
    /// derivation ŷ = ẑ × x̂ == v_axis even for left-handed popkern
    /// charts (v_axis = −(normal × u_axis)) — otherwise every pcurve
    /// trim on the face is REFLECTED in OpenCASCADE (measured on the
    /// flagship plate: trims displaced by the full face extent).
    #[test]
    fn left_handed_plane_exports_chart_true_placement() {
        let mut ctx = ExportContext::new();
        let left = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(0.0, 1.0, 0.0),
            v_axis: Vec3::new(1.0, 0.0, 0.0),
        };
        ctx.export_surface(&SurfaceOnFace {
            geometry: left,
            orientation: true,
        });
        // ẑ' = u × v = Y × X = −Z; readers derive ŷ = ẑ' × x̂ = X = v_axis.
        let dirs: Vec<String> = ctx
            .entities
            .iter()
            .filter(|e| e.type_name == "DIRECTION")
            .map(|e| format!("{:?}", e.params))
            .collect();
        assert!(
            dirs.iter().any(|d| d.contains("-1")),
            "left-handed chart must export ẑ = u×v = −normal, got {dirs:?}"
        );
    }

    #[test]
    fn sphere_pcurve_converts_colatitude_to_latitude() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let nurbs = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.3, 0.0, 0.0), Point3::new(0.3, 1.2, 0.0)],
            weights: vec![1.0, 1.0],
        };
        let out = pcurve_to_step_chart(&sphere, nurbs);
        // popkern colatitude 0 (north pole) → STEP latitude π/2.
        assert!((out.control_points[0].y - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
        assert!((out.control_points[1].y - (std::f64::consts::FRAC_PI_2 - 1.2)).abs() < 1e-12);
        // u untouched.
        assert!((out.control_points[0].x - 0.3).abs() < 1e-12);
    }

    #[test]
    fn export_point_deduplication() {
        let mut ctx = ExportContext::new();
        let p = Point3::new(1.0, 2.0, 3.0);
        let id1 = ctx.export_point(&p);
        let id2 = ctx.export_point(&p);
        assert_eq!(id1, id2, "Same point should produce same entity ID");
        assert_eq!(ctx.entities.len(), 1, "Should only emit one entity");
    }

    #[test]
    fn export_point_different() {
        let mut ctx = ExportContext::new();
        let id1 = ctx.export_point(&Point3::new(1.0, 2.0, 3.0));
        let id2 = ctx.export_point(&Point3::new(4.0, 5.0, 6.0));
        assert_ne!(id1, id2);
        assert_eq!(ctx.entities.len(), 2);
    }

    #[test]
    fn export_direction() {
        let mut ctx = ExportContext::new();
        let v = Vec3::new(0.0, 0.0, 1.0);
        let id = ctx.export_direction(&v);
        assert_eq!(ctx.entities[0].type_name, "DIRECTION");
        let id2 = ctx.export_direction(&v);
        assert_eq!(id, id2);
    }

    #[test]
    fn export_placement() {
        let mut ctx = ExportContext::new();
        let origin = Point3::new(0.0, 0.0, 0.0);
        let z = Vec3::new(0.0, 0.0, 1.0);
        let x = Vec3::new(1.0, 0.0, 0.0);
        let id = ctx.export_placement(&origin, &z, &x);
        let entity = ctx.entities.iter().find(|e| e.id == id);
        assert!(entity.is_some());
        assert_eq!(entity.unwrap().type_name, "AXIS2_PLACEMENT_3D");
    }

    #[test]
    fn export_line_curve() {
        let mut ctx = ExportContext::new();
        let curve = CurveOnEdge {
            geometry: CurveGeometry::Line {
                origin: Point3::new(0.0, 0.0, 0.0),
                direction: Vec3::new(1.0, 0.0, 0.0),
            },
            t_start: 0.0,
            t_end: 1.0,
        };
        let id = ctx.export_curve(&curve);
        let entity = ctx.entities.iter().find(|e| e.id == id);
        assert!(entity.is_some());
        assert_eq!(entity.unwrap().type_name, "LINE");
    }

    #[test]
    fn export_circle_curve() {
        let mut ctx = ExportContext::new();
        let curve = CurveOnEdge {
            geometry: CurveGeometry::Circle {
                center: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                radius: 5.0,
            },
            t_start: 0.0,
            t_end: std::f64::consts::TAU,
        };
        let id = ctx.export_curve(&curve);
        let entity = ctx.entities.iter().find(|e| e.id == id);
        assert!(entity.is_some());
        assert_eq!(entity.unwrap().type_name, "CIRCLE");
    }

    #[test]
    fn export_plane_surface() {
        let mut ctx = ExportContext::new();
        let surface = SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
            },
            orientation: true,
        };
        let id = ctx.export_surface(&surface);
        let entity = ctx.entities.iter().find(|e| e.id == id);
        assert!(entity.is_some());
        assert_eq!(entity.unwrap().type_name, "PLANE");
    }

    #[test]
    fn export_cylinder_surface() {
        let mut ctx = ExportContext::new();
        let surface = SurfaceOnFace {
            geometry: SurfaceGeometry::Cylinder {
                origin: Point3::new(0.0, 0.0, 0.0),
                axis: Vec3::new(0.0, 0.0, 1.0),
                radius: 2.0,
            },
            orientation: true,
        };
        let id = ctx.export_surface(&surface);
        let entity = ctx.entities.iter().find(|e| e.id == id);
        assert!(entity.is_some());
        assert_eq!(entity.unwrap().type_name, "CYLINDRICAL_SURFACE");
    }

    #[test]
    fn export_sphere_surface() {
        let mut ctx = ExportContext::new();
        let surface = SurfaceOnFace {
            geometry: SurfaceGeometry::Sphere {
                center: Point3::new(1.0, 2.0, 3.0),
                radius: 4.0,
            },
            orientation: true,
        };
        let id = ctx.export_surface(&surface);
        let entity = ctx.entities.iter().find(|e| e.id == id);
        assert!(entity.is_some());
        assert_eq!(entity.unwrap().type_name, "SPHERICAL_SURFACE");
    }
}
