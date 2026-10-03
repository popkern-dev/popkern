//! Swept STEP surfaces and the ellipse curve (PK-636), mapped exactly
//! onto kernel forms:
//!
//! - `SURFACE_OF_REVOLUTION(name, swept_curve, AXIS1_PLACEMENT)` — the
//!   profile's NURBS revolved a full turn (a rational NURBS surface);
//! - `SURFACE_OF_LINEAR_EXTRUSION(name, swept_curve, VECTOR)` — the
//!   profile's NURBS swept along the vector, degree 1 in v;
//! - `ELLIPSE(name, AXIS2_PLACEMENT_3D, a, b)` — the kernel ellipse.
//!
//! STEP leaves an extrusion's v unbounded. The NURBS spans
//! `|v| ≤ D / |V|`, D the diagonal of the box holding every 3-D point
//! of the file and the profile's control hull: every point of the
//! solid lies in that box, so the face's whole trim is inside the patch
//! (the patch is exact over any v range; the bound only keeps it
//! finite). Measured on the CAx-IF corpus: 24 revolved B-spline faces
//! (io1-tc) and 9 extruded ellipses (jack_in_the_box, moon_buggy) were
//! dropped as unsupported and their solids meshed open.
//!
//! The kernel NURBS parameterizations differ from STEP's (a rational
//! arc is not linear in angle), so a file pcurve on these faces fails
//! the fidelity check and is synthesized by projection, which the
//! import log records.

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::mapping::geometry::GeometryCache;
use crate::step::parser::types::{EntityId, Parameter};
use crate::step::schema::geometry::{Axis1Placement, StepEllipse, StepVector};

impl GeometryCache<'_> {
    /// `ELLIPSE` as a full-turn kernel ellipse in its placement frame.
    pub(super) fn resolve_ellipse(
        &mut self,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<CurveOnEdge> {
        let entity = self.graph.get(id)?;
        let ellipse = StepEllipse::from_params(&entity.params)
            .map_err(|e| log.error(Phase::MapGeometry, Some(id), format!("{e}")))
            .ok()?;
        let (center, z, x) = self.resolve_placement(ellipse.position, log)?;
        let x = (x - z * x.dot(&z)).normalized()?;
        let y = z.cross(&x);
        Some(CurveOnEdge {
            geometry: CurveGeometry::Ellipse {
                center,
                u_basis: x * ellipse.semi_axis_1,
                v_basis: y * ellipse.semi_axis_2,
                start_angle: 0.0,
                end_angle: std::f64::consts::TAU,
            },
            t_start: 0.0,
            t_end: std::f64::consts::TAU,
        })
    }

    /// `SURFACE_OF_REVOLUTION`: the profile revolved a full turn.
    pub(super) fn resolve_surface_of_revolution(
        &mut self,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<SurfaceGeometry> {
        let entity = self.graph.get(id)?;
        let (Some(Parameter::EntityRef(curve_id)), Some(Parameter::EntityRef(axis_id))) =
            (entity.params.get(1), entity.params.get(2))
        else {
            log.error(
                Phase::MapGeometry,
                Some(id),
                "SURFACE_OF_REVOLUTION needs a swept curve and an axis",
            );
            return None;
        };
        let profile = self.bounded_profile(*curve_id, id, log)?;
        let axis_entity = self.graph.get(*axis_id)?;
        let axis = Axis1Placement::from_params(&axis_entity.params)
            .map_err(|e| log.error(Phase::MapGeometry, Some(*axis_id), format!("{e}")))
            .ok()?;
        let origin = self.resolve_point(axis.location, log)?;
        let direction = match axis.axis {
            Some(d) => self.resolve_direction(d, log)?,
            None => Vec3::new(0.0, 0.0, 1.0),
        };
        popkern_core::geometry::revolution::revolution_surface(
            &profile,
            &origin,
            &direction,
            std::f64::consts::TAU,
        )
        .map(SurfaceGeometry::Nurbs)
        .map_err(|e| {
            log.error(
                Phase::MapGeometry,
                Some(id),
                format!("surface of revolution: {e}"),
            )
        })
        .ok()
    }

    /// `SURFACE_OF_LINEAR_EXTRUSION`: the profile swept along the
    /// vector, degree 1 in v over `|v| ≤ D / |V|` (module doc).
    pub(super) fn resolve_surface_of_linear_extrusion(
        &mut self,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<SurfaceGeometry> {
        let entity = self.graph.get(id)?;
        let (Some(Parameter::EntityRef(curve_id)), Some(Parameter::EntityRef(vector_id))) =
            (entity.params.get(1), entity.params.get(2))
        else {
            log.error(
                Phase::MapGeometry,
                Some(id),
                "SURFACE_OF_LINEAR_EXTRUSION needs a swept curve and a vector",
            );
            return None;
        };
        let profile = self.bounded_profile(*curve_id, id, log)?;
        let vector_entity = self.graph.get(*vector_id)?;
        let vector = StepVector::from_params(&vector_entity.params)
            .map_err(|e| log.error(Phase::MapGeometry, Some(*vector_id), format!("{e}")))
            .ok()?;
        let sweep = self.resolve_direction(vector.orientation, log)? * vector.magnitude;
        if sweep.length() <= 0.0 || !sweep.length().is_finite() {
            log.error(
                Phase::MapGeometry,
                Some(id),
                "extrusion vector has no length",
            );
            return None;
        }
        let (mut lo, mut hi) = self.model_box();
        for p in &profile.control_points {
            lo = Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        let reach = 1.1 * lo.distance_to(&hi).max(1e-9) / sweep.length();
        let n = profile.control_points.len();
        let mut control_points = Vec::with_capacity(2 * n);
        let mut weights = Vec::with_capacity(2 * n);
        for (p, w) in profile.control_points.iter().zip(&profile.weights) {
            control_points.push(*p - sweep * reach);
            control_points.push(*p + sweep * reach);
            weights.extend([*w, *w]);
        }
        Some(SurfaceGeometry::Nurbs(NurbsSurface {
            degree_u: profile.degree,
            degree_v: 1,
            knots_u: profile.knots.clone(),
            knots_v: vec![-reach, -reach, reach, reach],
            control_points,
            num_u: n,
            num_v: 2,
            weights,
        }))
    }

    /// The swept curve as a bounded NURBS; an unbounded line has no
    /// finite profile and is declined.
    fn bounded_profile(
        &mut self,
        curve_id: EntityId,
        surface_id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<popkern_core::geometry::nurbs_curve::NurbsCurve> {
        let curve = self.resolve_curve(curve_id, log)?;
        if matches!(curve.geometry, CurveGeometry::Line { .. }) {
            log.warn(
                Phase::MapGeometry,
                Some(surface_id),
                "swept surface over an unbounded LINE profile is not supported",
            );
            return None;
        }
        Some(
            curve
                .geometry
                .to_nurbs_in_range((curve.t_start, curve.t_end)),
        )
    }

    /// The box of every 3-D `CARTESIAN_POINT` in the file (cached).
    fn model_box(&mut self) -> (Point3, Point3) {
        if let Some(b) = self.model_box {
            return b;
        }
        let mut lo = Point3::new(f64::MAX, f64::MAX, f64::MAX);
        let mut hi = Point3::new(f64::MIN, f64::MIN, f64::MIN);
        for (_, entity) in self.graph.iter() {
            if entity.type_name != "CARTESIAN_POINT" {
                continue;
            }
            let Some(Parameter::List(coords)) = entity.params.get(1) else {
                continue;
            };
            let c: Vec<f64> = coords
                .iter()
                .filter_map(|p| match p {
                    Parameter::Real(v) => Some(*v),
                    Parameter::Int(v) => Some(*v as f64),
                    _ => None,
                })
                .collect();
            if c.len() == 3 && c.iter().all(|v| v.is_finite()) {
                lo = Point3::new(lo.x.min(c[0]), lo.y.min(c[1]), lo.z.min(c[2]));
                hi = Point3::new(hi.x.max(c[0]), hi.y.max(c[1]), hi.z.max(c[2]));
            }
        }
        if lo.x > hi.x {
            lo = Point3::new(0.0, 0.0, 0.0);
            hi = lo;
        }
        self.model_box = Some((lo, hi));
        (lo, hi)
    }
}
