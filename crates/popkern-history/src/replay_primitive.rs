//! Replay of [`Feature::Primitive`](crate::feature::Feature::Primitive)
//! (PK-1083): a primitive built in its local frame from an exact
//! profile, the same sweep kernel the scripts' profiles use (so planes,
//! cylinders, cones, spheres and tori stay exact), named by role, then
//! placed by its frame.
//!
//! Local frame: axis +z, reference +x. A box stands on a corner at the
//! origin (or is centred on it); a cylinder and cone stand on the
//! origin; a sphere and torus are centred on it. A sector (`angle`
//! below 360) of a cylinder, cone, sphere or torus turns its section
//! from the reference direction, its two cut faces meeting on the axis
//! (or, for a torus, apart).
//!
//! Names (under the feature's operation): the bottom face `Cap 0`, the
//! top `Cap 1` (a sphere's latitude caps too), a sector's start and
//! end faces `Cap 2` and `Cap 3`, every other face `Lateral i` in
//! order of angle about the axis, then height. An edge is an
//! `Intersection` and a vertex a `Junction` of the faces they bound
//! (by name), indexed in order along the axis and round it.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_math::PortableMath;
use popkern_sketch::plane::SketchPlane;
use popkern_sketch::profile::{AnalyticProfile, PathSegment, ProfileShape};

use crate::error::ReplayError;
use crate::failure::FeatureFailure;
use crate::feature_primitive::{PrimitiveDesc, PrimitiveFrame, PrimitiveShape};
use crate::parameter::{ParameterId, ParameterTable};
use crate::replay_side_effects::ExecuteSideEffects;
use crate::resolve::resolve_param;

#[path = "replay_primitive_polyhedra.rs"]
mod replay_primitive_polyhedra;

fn invalid(shape: &str, parameter: &str, value: f64, reason: &str) -> ReplayError {
    ReplayError::FeatureFailed {
        label: shape.into(),
        error: Box::new(FeatureFailure::InvalidParameter {
            parameter: parameter.into(),
            value: value.to_string(),
            reason: reason.into(),
        }),
    }
}

fn failed(shape: &str, what: &str, e: impl std::fmt::Display) -> ReplayError {
    ReplayError::FeatureFailed {
        label: shape.into(),
        error: FeatureFailure::sweep(shape, &format!("{what}: {e}")).into(),
    }
}

/// Build, name and place the primitive; the working body becomes it.
pub(super) fn execute(
    brep: &mut BRep,
    desc: &PrimitiveDesc,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
    table: &ParameterTable,
) -> Result<ExecuteSideEffects, ReplayError> {
    let shape = desc.shape.spelling();
    let value = |id: ParameterId| resolve_param(table, id);
    let positive = |role: &str, id: ParameterId| -> Result<f64, ReplayError> {
        let v = value(id)?;
        if v.is_finite() && v > 0.0 {
            Ok(v)
        } else {
            Err(invalid(shape, role, v, "must be positive"))
        }
    };
    // A sector's angle; a whole turn (360) is no sector.
    let sweep = match desc.shape.angle() {
        Some(id) => {
            let a = value(id)?;
            if !(a > 0.0 && a <= std::f64::consts::TAU + 1e-12) {
                return Err(invalid(
                    shape,
                    "angle",
                    a.to_degrees(),
                    "needs 0 < angle <= 360",
                ));
            }
            (a < std::f64::consts::TAU - 1e-12).then_some(a)
        }
        None => None,
    };
    let mut local = match &desc.shape {
        PrimitiveShape::Box { x, y, z, centred } => {
            let (x, y, z) = (positive("x", *x)?, positive("y", *y)?, positive("z", *z)?);
            let (x0, y0, z0) = if *centred {
                (-x / 2.0, -y / 2.0, -z / 2.0)
            } else {
                (0.0, 0.0, 0.0)
            };
            let profile = AnalyticProfile::new(ProfileShape::Polygon {
                points: vec![[x0, y0], [x0 + x, y0], [x0 + x, y0 + y], [x0, y0 + y]],
            });
            extruded(shape, &profile, z0, z, *tolerance)?
        }
        PrimitiveShape::Cylinder { radius, height, .. } => {
            let r = positive("radius", *radius)?;
            let h = positive("height", *height)?;
            if sweep.is_some() {
                let segments = vec![
                    PathSegment::Line { to: [r, 0.0] },
                    PathSegment::Line { to: [r, h] },
                    PathSegment::Line { to: [0.0, h] },
                ];
                revolved(shape, &section(segments, [0.0, 0.0]), sweep, *tolerance)?
            } else {
                let profile = AnalyticProfile::new(ProfileShape::Circle { radius: r });
                extruded(shape, &profile, 0.0, h, *tolerance)?
            }
        }
        PrimitiveShape::Cone { r1, r2, height, .. } => {
            let (a, b) = (value(*r1)?, value(*r2)?);
            let h = positive("height", *height)?;
            for (role, r) in [("r1", a), ("r2", b)] {
                if !(r.is_finite() && r >= 0.0) {
                    return Err(invalid(shape, role, r, "must not be negative"));
                }
            }
            if a <= 0.0 && b <= 0.0 {
                return Err(invalid(shape, "r1", a, "r1 and r2 cannot both be 0"));
            }
            let mut segments = Vec::new();
            if a > 0.0 {
                segments.push(PathSegment::Line { to: [a, 0.0] });
            }
            segments.push(PathSegment::Line { to: [b, h] });
            if b > 0.0 {
                segments.push(PathSegment::Line { to: [0.0, h] });
            }
            revolved(shape, &section(segments, [0.0, 0.0]), sweep, *tolerance)?
        }
        PrimitiveShape::Sphere {
            radius, latitudes, ..
        } => {
            let r = positive("radius", *radius)?;
            let half_pi = std::f64::consts::FRAC_PI_2;
            let (lo, hi) = match latitudes {
                Some([lo, hi]) => (value(*lo)?, value(*hi)?),
                None => (-half_pi, half_pi),
            };
            if !(lo >= -half_pi - 1e-12 && hi <= half_pi + 1e-12 && lo < hi) {
                return Err(invalid(
                    shape,
                    "latitudes",
                    lo.to_degrees(),
                    "needs -90 <= lo < hi <= 90 degrees",
                ));
            }
            let on = |a: f64| [r * a.pcos(), r * a.psin()];
            let pole = |a: f64| (a.abs() - half_pi).abs() <= 1e-12;
            let start = [0.0, r * lo.psin()];
            let mut segments = Vec::new();
            if !pole(lo) {
                segments.push(PathSegment::Line { to: on(lo) });
            }
            segments.push(PathSegment::Arc {
                via: on((lo + hi) / 2.0),
                to: if pole(hi) { [0.0, r] } else { on(hi) },
            });
            if !pole(hi) {
                segments.push(PathSegment::Line {
                    to: [0.0, r * hi.psin()],
                });
            }
            revolved(shape, &section(segments, start), sweep, *tolerance)?
        }
        PrimitiveShape::Prism {
            sides,
            radius,
            height,
            flats,
        } => {
            let r = positive("radius", *radius)?;
            let h = positive("height", *height)?;
            replay_primitive_polyhedra::prism(shape, *sides, r, h, *flats, *tolerance)?
        }
        PrimitiveShape::Pyramid {
            sides,
            r1,
            r2,
            height,
            flats,
        } => {
            let (a, b) = (value(*r1)?, value(*r2)?);
            for (role, r) in [("r1", a), ("r2", b)] {
                if !(r.is_finite() && r >= 0.0) {
                    return Err(invalid(shape, role, r, "must not be negative"));
                }
            }
            let h = positive("height", *height)?;
            replay_primitive_polyhedra::pyramid(shape, *sides, (a, b), h, *flats, *tolerance)?
        }
        PrimitiveShape::Wedge { x, y, z, top_x } => {
            let size = [positive("x", *x)?, positive("y", *y)?, positive("z", *z)?];
            let t = value(*top_x)?;
            if !(t.is_finite() && t >= 0.0) {
                return Err(invalid(shape, "top_x", t, "must not be negative"));
            }
            replay_primitive_polyhedra::wedge(shape, size, t, *tolerance)?
        }
        PrimitiveShape::Torus { major, minor, .. } => {
            let big = positive("major", *major)?;
            let small = positive("minor", *minor)?;
            if small >= big {
                return Err(invalid(
                    shape,
                    "minor",
                    small,
                    "must be less than major (a spindle or horn torus is not this primitive)",
                ));
            }
            let profile = AnalyticProfile {
                shape: ProfileShape::Circle { radius: small },
                center: [big, 0.0],
            };
            revolved(shape, &profile, sweep, *tolerance)?
        }
    };
    super::replay_primitive_names::name_by_role(&mut local, tracking, sweep.is_some());
    *brep = placed(shape, &local, &desc.frame)?;
    Ok(ExecuteSideEffects::empty())
}

/// A closed path of the xz section from `start` (x radial, z axial).
fn section(segments: Vec<PathSegment<f64>>, start: [f64; 2]) -> AnalyticProfile<f64> {
    AnalyticProfile::new(ProfileShape::Path { start, segments })
}

/// `profile` on the xy plane swept along +z from `z0` by `h`.
fn extruded(
    shape: &str,
    profile: &AnalyticProfile<f64>,
    z0: f64,
    h: f64,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let wire = profile
        .placed_wire(&SketchPlane::xy())
        .map_err(|e| failed(shape, "profile", e))?;
    // The sweep carries the profile from the path's start: start it at
    // the wire's sampled centroid, as the scripts' exact extrude does.
    let (t0, t1) = wire.parameter_range();
    let n = 64;
    let mut sum = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..n {
        sum = sum
            + (wire.point_at(t0 + (t1 - t0) * i as f64 / n as f64) - Point3::new(0.0, 0.0, 0.0));
    }
    let mean = Point3::new(0.0, 0.0, 0.0) + sum * (1.0 / n as f64);
    let lift = Vec3::new(0.0, 0.0, 1.0);
    let wire = NurbsCurve {
        control_points: wire.control_points.iter().map(|p| *p + lift * z0).collect(),
        ..wire
    };
    let (a, b) = (mean + lift * z0, mean + lift * (z0 + h));
    let path = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    };
    popkern_sweep::sweep(
        &popkern_sweep::SweepProfile::Wire(wire),
        &path,
        &popkern_sweep::SweepOptions::default(),
        tolerance,
        None,
    )
    .map(|r| r.brep)
    .map_err(|e| failed(shape, "sweep", e))
}

/// `profile` on the xz plane turned about +z (whole turn by default).
fn revolved(
    shape: &str,
    profile: &AnalyticProfile<f64>,
    angle: Option<f64>,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let wire = profile
        .placed_wire(&SketchPlane::xz())
        .map_err(|e| failed(shape, "profile", e))?;
    popkern_sweep::revolve(
        &popkern_sweep::SweepProfile::Wire(wire),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        angle.unwrap_or(std::f64::consts::TAU),
        tolerance,
        None,
    )
    .map(|r| r.brep)
    .map_err(|e| failed(shape, "revolve", e))
}

/// The local body placed by `frame`.
fn placed(shape: &str, local: &BRep, frame: &PrimitiveFrame) -> Result<BRep, ReplayError> {
    if frame.is_world() {
        return Ok(local.clone());
    }
    let v = |a: [f64; 3]| Vec3::new(a[0], a[1], a[2]);
    let axis = v(frame.axis);
    let reference = v(frame.reference);
    let reference = reference - axis * (reference.dot(&axis) / axis.dot(&axis).max(1e-300));
    let placement = popkern_core::geometry::transform::Transform3::from_axis_placement(
        Point3::new(frame.origin[0], frame.origin[1], frame.origin[2]),
        axis,
        reference,
    )
    .map_err(|_| {
        invalid(
            shape,
            "axis",
            axis.length(),
            "the axis must be non-zero and the reference not along it",
        )
    })?;
    local
        .transform(&placement)
        .map_err(|e| failed(shape, "placement", e))
}
