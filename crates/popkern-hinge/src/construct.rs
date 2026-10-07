//! Construction helpers: a profile compiled to an exact planar sheet
//! and extruded (lines become planes, arcs cylinders), booleans and
//! placements with the hinge stage named, and the validation every
//! served body passes.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sketch::{AnalyticProfile, ProfileShape, SketchPlane};

use crate::HingeError;

/// A sketch plane through `origin` spanned by `u` and `v` (unit,
/// orthogonal); its normal is `u × v`.
pub(crate) fn plane(origin: Point3, u: Vec3, v: Vec3) -> SketchPlane {
    SketchPlane {
        origin,
        u_axis: u,
        v_axis: v,
        normal: u.cross(&v),
    }
}

/// Extrude `shape` (in `plane`'s coordinates) along the plane normal
/// by `depth`, into a validated solid.
pub(crate) fn extrude(
    shape: ProfileShape<f64>,
    plane: SketchPlane,
    depth: f64,
    stage: &str,
    tolerance: Tolerance,
) -> Result<BRep, HingeError> {
    let normal = plane.normal;
    let sheet = AnalyticProfile {
        shape,
        center: [0.0, 0.0],
    }
    .to_brep(plane, tolerance, None)
    .map_err(|e| construction(stage, format!("profile compile failed: {e}")))?;
    let face = sheet
        .regions
        .first()
        .map(|r| r.face)
        .ok_or_else(|| construction(stage, "the profile produced no face".into()))?;
    let mut brep = sheet.brep;
    brep.extrude(face, normal, depth, None)
        .map_err(|e| construction(stage, format!("extrude failed: {e}")))?;
    validated(brep, stage)
}

/// A cylinder of `radius` on the Z axis from `z0` to `z1`.
pub(crate) fn z_cylinder(
    radius: f64,
    z0: f64,
    z1: f64,
    stage: &str,
    tolerance: Tolerance,
) -> Result<BRep, HingeError> {
    extrude(
        ProfileShape::Circle { radius },
        plane(
            Point3::new(0.0, 0.0, z0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ),
        z1 - z0,
        stage,
        tolerance,
    )
}

fn boolean(
    op: popkern_bool::BoolOperation,
    body: &BRep,
    tool: &BRep,
    stage: &str,
    tolerance: Tolerance,
) -> Result<BRep, HingeError> {
    let result = popkern_bool::boolean(op, body.clone(), tool.clone(), tolerance, None).map_err(
        |source| HingeError::Boolean {
            stage: stage.to_string(),
            source: Box::new(source),
        },
    )?;
    validated(result, stage)
}

/// `body − tool`.
pub(crate) fn subtract(
    body: &BRep,
    tool: &BRep,
    stage: &str,
    tolerance: Tolerance,
) -> Result<BRep, HingeError> {
    boolean(
        popkern_bool::BoolOperation::Subtraction,
        body,
        tool,
        stage,
        tolerance,
    )
}

/// `body ∪ other`.
pub(crate) fn union(
    body: &BRep,
    other: &BRep,
    stage: &str,
    tolerance: Tolerance,
) -> Result<BRep, HingeError> {
    boolean(
        popkern_bool::BoolOperation::Union,
        body,
        other,
        stage,
        tolerance,
    )
}

/// `body` moved by `t`.
pub(crate) fn placed(body: &BRep, t: &Transform3, stage: &str) -> Result<BRep, HingeError> {
    let moved = body
        .transform(t)
        .map_err(|e| construction(stage, format!("placement failed: {e}")))?;
    validated(moved, stage)
}

/// A served body: non-empty, structurally and geometrically valid,
/// every edge curve on its faces.
pub(crate) fn validated(brep: BRep, stage: &str) -> Result<BRep, HingeError> {
    if brep.faces.is_empty() {
        return Err(construction(
            stage,
            "produced an empty body (0 faces)".into(),
        ));
    }
    let invalid = |source| HingeError::Validation {
        stage: stage.to_string(),
        source,
    };
    brep.validate().map_err(invalid)?;
    brep.validate_curve_surface().map_err(invalid)?;
    Ok(brep)
}

pub(crate) fn construction(stage: &str, detail: String) -> HingeError {
    HingeError::Construction {
        stage: stage.to_string(),
        detail,
    }
}
