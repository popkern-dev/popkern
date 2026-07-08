//! Geometry mapping: STEP geometric entities to popkern-core geometry types.
//!
//! [`GeometryCache`] lazily resolves and caches frequently-referenced entities
//! (points, directions, axis placements) while resolving curves and surfaces
//! on demand since they are typically referenced only once per edge/face.
//!
//! Surface resolution is in the [`surfaces`] submodule.

mod pcurves;
mod surfaces;
#[cfg(test)]
mod tests_resolve;

pub(crate) use pcurves::ResolvedStepPcurve;

use std::collections::HashMap;

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceOnFace;
use popkern_core::geometry::{Point3, Vec3};

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::parser::graph::EntityGraph;
use crate::step::parser::types::EntityId;
use crate::step::schema::as_real_list;
use crate::step::schema::geometry::{
    Axis2Placement3d, BSplineCurveWithKnots, CartesianPoint, Direction, StepCircle, StepLine,
    StepVector, SurfaceCurve,
};

/// Lazy-resolving cache for STEP geometry entities.
///
/// Points, directions, and axis placements are cached because they are
/// referenced many times. Curves and surfaces are resolved on demand.
pub struct GeometryCache<'a> {
    graph: &'a EntityGraph,
    points: HashMap<EntityId, Point3>,
    directions: HashMap<EntityId, Vec3>,
    placements: HashMap<EntityId, (Point3, Vec3, Vec3)>,
}

impl<'a> GeometryCache<'a> {
    /// Create a new cache backed by the given entity graph.
    pub fn new(graph: &'a EntityGraph) -> Self {
        Self {
            graph,
            points: HashMap::new(),
            directions: HashMap::new(),
            placements: HashMap::new(),
        }
    }

    /// Resolve a `CARTESIAN_POINT` entity to a [`Point3`].
    pub fn resolve_point(&mut self, id: EntityId, log: &mut DiagnosticLog) -> Option<Point3> {
        if let Some(&pt) = self.points.get(&id) {
            return Some(pt);
        }
        let pt = self.resolve_point_inner(id, log)?;
        self.points.insert(id, pt);
        Some(pt)
    }

    /// Resolve a `DIRECTION` entity to a normalised [`Vec3`].
    pub fn resolve_direction(&mut self, id: EntityId, log: &mut DiagnosticLog) -> Option<Vec3> {
        if let Some(&dir) = self.directions.get(&id) {
            return Some(dir);
        }
        let dir = self.resolve_direction_inner(id, log)?;
        self.directions.insert(id, dir);
        Some(dir)
    }

    /// Resolve an `AXIS2_PLACEMENT_3D` to (origin, z_axis, x_axis).
    pub fn resolve_placement(
        &mut self,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<(Point3, Vec3, Vec3)> {
        if let Some(&p) = self.placements.get(&id) {
            return Some(p);
        }
        let p = self.resolve_placement_inner(id, log)?;
        self.placements.insert(id, p);
        Some(p)
    }

    /// Resolve a curve entity (LINE, CIRCLE, B_SPLINE_CURVE_WITH_KNOTS, SURFACE_CURVE)
    /// to a [`CurveOnEdge`].
    pub fn resolve_curve(&mut self, id: EntityId, log: &mut DiagnosticLog) -> Option<CurveOnEdge> {
        let entity = self.graph.get(id)?;
        match entity.type_name.as_str() {
            "LINE" => self.resolve_line(id, log),
            "CIRCLE" => self.resolve_circle(id, log),
            "B_SPLINE_CURVE_WITH_KNOTS" => self.resolve_bspline_curve(id, log),
            // Standard Part-21 complex form for rational B-splines:
            // the leading partial is BOUNDED_CURVE, the spline data
            // lives in the B_SPLINE_CURVE / B_SPLINE_CURVE_WITH_KNOTS
            // partials.
            "BOUNDED_CURVE"
                if entity
                    .complex_parts
                    .iter()
                    .any(|(n, _)| n == "B_SPLINE_CURVE_WITH_KNOTS") =>
            {
                self.resolve_bspline_curve_complex(id, log)
            }
            "SURFACE_CURVE" | "SEAM_CURVE" => self.resolve_surface_curve(id, log),
            other => {
                log.warn(
                    Phase::MapGeometry,
                    Some(id),
                    format!("unsupported curve type: {other}"),
                );
                None
            }
        }
    }

    /// Resolve a surface entity to a [`SurfaceOnFace`].
    pub fn resolve_surface(
        &mut self,
        id: EntityId,
        same_sense: bool,
        log: &mut DiagnosticLog,
    ) -> Option<SurfaceOnFace> {
        let entity = self.graph.get(id)?;
        let geometry = match entity.type_name.as_str() {
            "PLANE" => self.resolve_plane(id, log)?,
            "CYLINDRICAL_SURFACE" => self.resolve_cylindrical_surface(id, log)?,
            "CONICAL_SURFACE" => self.resolve_conical_surface(id, log)?,
            "SPHERICAL_SURFACE" => self.resolve_spherical_surface(id, log)?,
            "TOROIDAL_SURFACE" => self.resolve_toroidal_surface(id, log)?,
            "B_SPLINE_SURFACE_WITH_KNOTS" => self.resolve_bspline_surface(id, log)?,
            other => {
                log.warn(
                    Phase::MapGeometry,
                    Some(id),
                    format!("unsupported surface type: {other}"),
                );
                return None;
            }
        };
        Some(SurfaceOnFace {
            geometry,
            orientation: same_sense,
        })
    }
}

// ---------------------------------------------------------------------------
// Private resolution helpers
// ---------------------------------------------------------------------------

impl<'a> GeometryCache<'a> {
    fn resolve_point_inner(&self, id: EntityId, log: &mut DiagnosticLog) -> Option<Point3> {
        let entity = self.graph.get(id)?;
        match CartesianPoint::from_params(&entity.params) {
            Ok(cp) => Some(Point3::new(cp.coords[0], cp.coords[1], cp.coords[2])),
            Err(e) => {
                log.error(Phase::MapGeometry, Some(id), format!("{e}"));
                None
            }
        }
    }

    fn resolve_direction_inner(&self, id: EntityId, log: &mut DiagnosticLog) -> Option<Vec3> {
        let entity = self.graph.get(id)?;
        match Direction::from_params(&entity.params) {
            Ok(d) => {
                let v = Vec3::new(d.ratios[0], d.ratios[1], d.ratios[2]);
                v.normalized().or_else(|| {
                    log.warn(
                        Phase::MapGeometry,
                        Some(id),
                        "zero-length direction, using (0,0,1)",
                    );
                    Some(Vec3::new(0.0, 0.0, 1.0))
                })
            }
            Err(e) => {
                log.error(Phase::MapGeometry, Some(id), format!("{e}"));
                None
            }
        }
    }

    fn resolve_placement_inner(
        &mut self,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<(Point3, Vec3, Vec3)> {
        let entity = self.graph.get(id)?;
        let a2p = match Axis2Placement3d::from_params(&entity.params) {
            Ok(a) => a,
            Err(e) => {
                log.error(Phase::MapGeometry, Some(id), format!("{e}"));
                return None;
            }
        };

        let origin = self.resolve_point(a2p.location, log)?;
        let z_axis = match a2p.axis {
            Some(axis_id) => self.resolve_direction(axis_id, log)?,
            None => Vec3::new(0.0, 0.0, 1.0),
        };
        let x_axis = match a2p.ref_direction {
            Some(ref_id) => self.resolve_direction(ref_id, log)?,
            None => Vec3::new(1.0, 0.0, 0.0),
        };
        Some((origin, z_axis, x_axis))
    }

    fn resolve_line(&mut self, id: EntityId, log: &mut DiagnosticLog) -> Option<CurveOnEdge> {
        let entity = self.graph.get(id)?;
        let line = match StepLine::from_params(&entity.params) {
            Ok(l) => l,
            Err(e) => {
                log.error(Phase::MapGeometry, Some(id), format!("{e}"));
                return None;
            }
        };
        let origin = self.resolve_point(line.point, log)?;
        let direction = self.resolve_vector_direction(line.direction, log)?;
        Some(CurveOnEdge {
            geometry: CurveGeometry::Line { origin, direction },
            t_start: 0.0,
            t_end: 1.0,
        })
    }

    fn resolve_circle(&mut self, id: EntityId, log: &mut DiagnosticLog) -> Option<CurveOnEdge> {
        let entity = self.graph.get(id)?;
        let circle = match StepCircle::from_params(&entity.params) {
            Ok(c) => c,
            Err(e) => {
                log.error(Phase::MapGeometry, Some(id), format!("{e}"));
                return None;
            }
        };
        let (center, normal, _x_axis) = self.resolve_placement(circle.position, log)?;
        Some(CurveOnEdge {
            geometry: CurveGeometry::Circle {
                center,
                normal,
                radius: circle.radius,
            },
            t_start: 0.0,
            t_end: std::f64::consts::TAU,
        })
    }

    fn resolve_bspline_curve(
        &mut self,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<CurveOnEdge> {
        let entity = self.graph.get(id)?;
        let bs = match BSplineCurveWithKnots::from_params(&entity.params) {
            Ok(b) => b,
            Err(e) => {
                log.error(Phase::MapGeometry, Some(id), format!("{e}"));
                return None;
            }
        };

        let mut control_points = Vec::with_capacity(bs.control_points.len());
        for &cp_id in &bs.control_points {
            control_points.push(self.resolve_point(cp_id, log)?);
        }

        // Expand knots from (value, multiplicity) pairs
        let knots = expand_knots(&bs.knots, &bs.knot_multiplicities);

        let mut weights = resolve_rational_bspline_curve_weights(entity, log)
            .or_else(|| bs.weights.clone())
            .unwrap_or_else(|| vec![1.0; control_points.len()]);
        if weights.len() != control_points.len() {
            log.warn(
                Phase::MapGeometry,
                Some(id),
                format!(
                    "rational B-spline curve weight count {} does not match control point count {}; using unit weights",
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
                degree: bs.degree as usize,
                knots,
                control_points,
                weights,
            }),
            t_start,
            t_end,
        })
    }

    fn resolve_surface_curve(
        &mut self,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<CurveOnEdge> {
        let entity = self.graph.get(id)?;
        let sc = match SurfaceCurve::from_params(&entity.params) {
            Ok(s) => s,
            Err(e) => {
                log.error(Phase::MapGeometry, Some(id), format!("{e}"));
                return None;
            }
        };
        // Delegate to the 3D curve
        self.resolve_curve(sc.curve_3d, log)
    }

    /// Resolve a VECTOR entity to get its direction (normalised).
    fn resolve_vector_direction(&mut self, id: EntityId, log: &mut DiagnosticLog) -> Option<Vec3> {
        let entity = self.graph.get(id)?;
        match StepVector::from_params(&entity.params) {
            Ok(v) => self
                .resolve_direction(v.orientation, log)
                .map(|direction| direction * v.magnitude),
            Err(e) => {
                log.error(Phase::MapGeometry, Some(id), format!("{e}"));
                None
            }
        }
    }
}

impl<'a> GeometryCache<'a> {
    /// Resolve the standard Part-21 complex form of a (rational)
    /// B-spline curve: degree + control points from the
    /// B_SPLINE_CURVE partial, knots from B_SPLINE_CURVE_WITH_KNOTS,
    /// weights from RATIONAL_B_SPLINE_CURVE (defaulting to 1s).
    fn resolve_bspline_curve_complex(
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

fn resolve_rational_bspline_curve_weights(
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

/// Expand (knot_values, multiplicities) into a flat knot vector.
fn expand_knots(values: &[f64], multiplicities: &[i64]) -> Vec<f64> {
    let mut knots = Vec::new();
    for (val, &mult) in values.iter().zip(multiplicities.iter()) {
        for _ in 0..mult {
            knots.push(*val);
        }
    }
    knots
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step::parser::types::{EntityInstance, Parameter};
    use popkern_core::geometry::surface::SurfaceGeometry;

    fn ent(id: EntityId, ty: &str, p: Vec<Parameter>) -> EntityInstance {
        EntityInstance {
            id,
            type_name: ty.into(),
            params: p,
            complex_parts: vec![],
        }
    }
    fn s(v: &str) -> Parameter {
        Parameter::String(v.into())
    }
    fn r(v: f64) -> Parameter {
        Parameter::Real(v)
    }
    fn er(id: u64) -> Parameter {
        Parameter::EntityRef(id)
    }
    fn l(items: Vec<Parameter>) -> Parameter {
        Parameter::List(items)
    }

    fn graph(entities: Vec<EntityInstance>) -> EntityGraph {
        EntityGraph::build(entities, &mut DiagnosticLog::new())
    }

    /// Build the standard origin-Z-X placement (point #1, dir #2 = Z, dir #3 = X, a2p #4).
    fn placement_entities(px: f64, py: f64, pz: f64) -> Vec<EntityInstance> {
        vec![
            ent(
                1,
                "CARTESIAN_POINT",
                vec![s(""), l(vec![r(px), r(py), r(pz)])],
            ),
            ent(2, "DIRECTION", vec![s(""), l(vec![r(0.0), r(0.0), r(1.0)])]),
            ent(3, "DIRECTION", vec![s(""), l(vec![r(1.0), r(0.0), r(0.0)])]),
            ent(4, "AXIS2_PLACEMENT_3D", vec![s(""), er(1), er(2), er(3)]),
        ]
    }

    #[test]
    fn resolve_cartesian_point() {
        let g = graph(vec![ent(
            1,
            "CARTESIAN_POINT",
            vec![s("p"), l(vec![r(1.0), r(2.0), r(3.0)])],
        )]);
        let mut c = GeometryCache::new(&g);
        let mut log = DiagnosticLog::new();
        let pt = c.resolve_point(1, &mut log).map(|p| (p.x, p.y, p.z));
        assert_eq!(pt, Some((1.0, 2.0, 3.0)));
        assert!(c.resolve_point(1, &mut log).is_some()); // cached
    }

    #[test]
    fn resolve_direction_normalized() {
        let g = graph(vec![ent(
            1,
            "DIRECTION",
            vec![s("d"), l(vec![r(0.0), r(0.0), r(2.0)])],
        )]);
        let mut c = GeometryCache::new(&g);
        let mut log = DiagnosticLog::new();
        let d = c.resolve_direction(1, &mut log);
        assert!(d.is_some());
        let v = d.unwrap();
        assert!(v.x.abs() < 1e-10 && v.y.abs() < 1e-10 && (v.z - 1.0).abs() < 1e-10);
    }

    #[test]
    fn resolve_axis2_placement() {
        let ents = placement_entities(10.0, 20.0, 30.0);
        // Override point to (10, 20, 30)
        let g = graph(ents);
        let mut c = GeometryCache::new(&g);
        let mut log = DiagnosticLog::new();
        let (origin, z, x) = c.resolve_placement(4, &mut log).unwrap();
        assert!((origin.x - 10.0).abs() < 1e-10);
        assert!((z.z - 1.0).abs() < 1e-10);
        assert!((x.x - 1.0).abs() < 1e-10);
    }

    #[test]
    fn resolve_plane_surface() {
        let mut ents = placement_entities(0.0, 0.0, 0.0);
        ents.push(ent(5, "PLANE", vec![s(""), er(4)]));
        let g = graph(ents);
        let mut c = GeometryCache::new(&g);
        let mut log = DiagnosticLog::new();
        let sof = c.resolve_surface(5, true, &mut log).unwrap();
        assert!(sof.orientation);
        assert!(matches!(sof.geometry, SurfaceGeometry::Plane { .. }));
    }

    #[test]
    fn expand_knots_works() {
        assert_eq!(
            expand_knots(&[0.0, 1.0], &[3, 3]),
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
        );
    }
}
