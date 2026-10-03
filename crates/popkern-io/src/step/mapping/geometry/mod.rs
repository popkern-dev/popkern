//! Geometry mapping: STEP geometric entities to popkern-core geometry types.
//!
//! [`GeometryCache`] lazily resolves and caches frequently-referenced entities
//! (points, directions, axis placements) while resolving curves and surfaces
//! on demand since they are typically referenced only once per edge/face.
//!
//! Surface resolution is in the [`surfaces`] submodule.

mod bspline_implicit_knots;
mod curves_complex;
mod pcurve_chart;
mod pcurves;
mod placement_axes;
mod surfaces;
mod surfaces_swept;
use curves_complex::resolve_rational_bspline_curve_weights;
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
use crate::step::schema::geometry::{
    Axis2Placement2d, Axis2Placement3d, BSplineCurveWithKnots, CartesianPoint, Direction,
    StepCircle, StepLine, StepVector, SurfaceCurve,
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
    /// Radians per file plane-angle unit (π/180 for a degree file).
    radians_per_angle: f64,
    /// Box of the file's 3-D points (bounds extrusion patches).
    model_box: Option<(Point3, Point3)>,
    /// How far a vertex may sit from its edge curve's end and still be
    /// that end (the file's declared precision, PK-639).
    trim_allowance: f64,
}

impl<'a> GeometryCache<'a> {
    /// Create a new cache backed by the given entity graph.
    pub fn new(graph: &'a EntityGraph) -> Self {
        Self {
            graph,
            points: HashMap::new(),
            directions: HashMap::new(),
            placements: HashMap::new(),
            radians_per_angle: 1.0,
            model_box: None,
            trim_allowance: 0.0,
        }
    }

    /// Accept vertices within `allowance` of their edge curve's ends
    /// (the file's declared precision; trim recovery only).
    pub(crate) fn set_trim_allowance(&mut self, allowance: f64) {
        self.trim_allowance = allowance;
    }

    /// The trim-recovery allowance.
    pub(crate) fn trim_allowance(&self) -> f64 {
        self.trim_allowance
    }

    /// Read angle-valued attributes (a cone's semi-angle) in a file
    /// plane-angle unit of `radians_per_angle` radians.
    pub(crate) fn set_radians_per_angle(&mut self, radians_per_angle: f64) {
        self.radians_per_angle = radians_per_angle;
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
            "ELLIPSE" => self.resolve_ellipse(id, log),
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
            // INTERSECTION_CURVE is a SURFACE_CURVE subtype with the same
            // attributes (its associated geometry names the two surfaces
            // it lies on); the 3-D curve is the edge geometry either way.
            t if bspline_implicit_knots::IMPLICIT_CURVES.contains(&t) => {
                self.resolve_implicit_knot_curve(id, log)
            }
            "BOUNDED_CURVE"
                if entity.complex_parts.iter().any(|(n, _)| {
                    bspline_implicit_knots::IMPLICIT_CURVES.contains(&n.as_str())
                }) =>
            {
                self.resolve_implicit_knot_curve(id, log)
            }
            "SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE" => {
                self.resolve_surface_curve(id, log)
            }
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
            "SURFACE_OF_REVOLUTION" => self.resolve_surface_of_revolution(id, log)?,
            t if bspline_implicit_knots::IMPLICIT_SURFACES.contains(&t) => {
                self.resolve_implicit_knot_surface(id, log)?
            }
            "BOUNDED_SURFACE"
                if entity.complex_parts.iter().any(|(n, _)| {
                    bspline_implicit_knots::IMPLICIT_SURFACES.contains(&n.as_str())
                }) =>
            {
                self.resolve_implicit_knot_surface(id, log)?
            }
            "SURFACE_OF_LINEAR_EXTRUSION" => self.resolve_surface_of_linear_extrusion(id, log)?,
            // Standard Part-21 complex form for rational B-spline
            // surfaces (mirrors the BOUNDED_CURVE arm above).
            "BOUNDED_SURFACE"
                if entity
                    .complex_parts
                    .iter()
                    .any(|(n, _)| n == "B_SPLINE_SURFACE_WITH_KNOTS") =>
            {
                self.resolve_bspline_surface_complex(id, log)?
            }
            other => {
                log.warn(
                    Phase::MapGeometry,
                    Some(id),
                    format!("unsupported surface type: {other}"),
                );
                return None;
            }
        };
        // The kernel's orientation flag says the outer loop runs
        // counter-clockwise in the surface CHART (the winding-honesty
        // contract); STEP's `same_sense` says it runs counter-clockwise
        // about STEP's surface normal. They agree where the chart normal
        // ∂u × ∂v is STEP's normal — every elementary surface except
        // the sphere, whose colatitude chart normal points INWARD.
        // Without the flip, every imported sphere face whose loop the
        // audit cannot measure (a cap around a pole) kept a backwards
        // flag (measured: jack_in_the_box's ball knobs).
        let chart_opposes_step = matches!(
            geometry,
            popkern_core::geometry::surface::SurfaceGeometry::Sphere { .. }
        );
        Some(SurfaceOnFace {
            geometry,
            orientation: same_sense != chart_opposes_step,
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
        if entity.type_name == "AXIS2_PLACEMENT_2D" {
            // A 2-D placement (a pcurve conic's frame) is the z = 0
            // plane's frame: z is +Z, x the given ref direction.
            let a2p = match Axis2Placement2d::from_params(&entity.params) {
                Ok(a) => a,
                Err(e) => {
                    log.error(Phase::MapGeometry, Some(id), format!("{e}"));
                    return None;
                }
            };
            let origin = self.resolve_point(a2p.location, log)?;
            let x_axis = match a2p.ref_direction {
                Some(ref_id) => self.resolve_direction(ref_id, log)?,
                None => Vec3::new(1.0, 0.0, 0.0),
            };
            return Some((origin, Vec3::new(0.0, 0.0, 1.0), x_axis));
        }
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
            None => placement_axes::default_ref_direction(z_axis),
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
#[path = "tests_cache.rs"]
mod tests;
