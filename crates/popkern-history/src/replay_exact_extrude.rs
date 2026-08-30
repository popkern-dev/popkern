//! Exact-wire extrude for analytic profiles (design D2, spec §3.1).
//!
//! `popkern-part` extrudes an analytic profile by sweeping its EXACT
//! closed wire (`AnalyticProfile::local_wire`) along a line — one
//! periodic cylinder face for a circle. The tree's `BRep::extrude`
//! sweeps every edge of the compiled sketch FACE, whose full-circle
//! boundary the sketch builder splits into two arcs (self-loop edges
//! are avoided there), so the same document replayed to two
//! half-cylinders: a fillet recorded by the direct build's selection
//! then failed on the co-cylindrical seam lines (measured on the boss
//! plate, PK-382). Recording and direct building must be one code
//! path (parity fixtures, design §13), so a blind extrude of an
//! analytic sketch face takes the very arithmetic `popkern-part`
//! uses; every other extrude keeps `BRep::extrude`.
//!
//! The exact path names the result like a snapshot: every face, edge
//! and vertex `Generated/…` in slotmap order under the extrude's
//! operation (the sketch sheet is replaced wholesale).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, OperationId};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sketch::SketchPlane;

use crate::ReplayError;
use crate::feature::{ExtrudeEnd, ExtrudeSecondDir, Feature, SketchPlaneDesc, Vec3Desc};
use crate::feature_sketch::SketchPayload;
use crate::parameter::{ParameterId, ParameterTable};

/// The exact profile behind an extrude's face.
#[derive(Debug, Clone)]
pub(crate) struct ExactProfileWire {
    /// The unplaced wire lifted onto the plane (as `popkern-part` sweeps it).
    pub wire: NurbsCurve,
    /// Placement on the plane.
    pub center: [f64; 2],
    /// The sketch plane.
    pub plane: SketchPlane,
}

/// Find the analytic profile that produced `face`, from the features
/// in scope (`(feature, operation)` pairs).
pub(crate) fn exact_profile_for<'a>(
    features: impl IntoIterator<Item = (&'a Feature, OperationId)>,
    face: &EntityName,
    brep: &BRep,
    param_table: &ParameterTable,
) -> Option<ExactProfileWire> {
    let (feature, _) = features.into_iter().find(|(_, op)| *op == face.operation)?;
    let Feature::Sketch {
        plane,
        sketch_data: SketchPayload::Analytic(profile),
    } = feature
    else {
        return None;
    };
    let profile = crate::feature_sketch::resolve_profile(profile, param_table).ok()?;
    let plane = crate::resolve::resolve_sketch_plane(plane, brep, param_table).ok()?;
    let wire = profile.local_wire(&plane).ok()?;
    let _ = SketchPlaneDesc::XY;
    Some(ExactProfileWire {
        wire,
        center: profile.center,
        plane,
    })
}

/// Whether this extrude takes the exact path: blind, no draft, no
/// second direction, direction absent or along the plane normal.
pub(crate) fn takes_exact_path(
    direction: Option<&Vec3Desc>,
    end_condition: &ExtrudeEnd,
    draft_angle: Option<&ParameterId>,
    second_direction: Option<&ExtrudeSecondDir>,
    plane: &SketchPlane,
) -> bool {
    if draft_angle.is_some() || second_direction.is_some() {
        return false;
    }
    if !matches!(end_condition, ExtrudeEnd::Blind { .. }) {
        return false;
    }
    match direction {
        None => true,
        Some(Vec3Desc::Literal(v)) => {
            let d = Vec3::new(v[0], v[1], v[2]);
            d.normalized()
                .map(|d| d.dot(&plane.normal).abs() > 1.0 - 1e-9)
                .unwrap_or(false)
        }
        Some(_) => false,
    }
}

/// Sweep the exact wire by `distance` along the plane normal (or its
/// negation for a negative distance), replacing `brep`, and name
/// every entity under this operation.
pub(crate) fn execute_extrude_exact(
    brep: &mut BRep,
    exact: &ExactProfileWire,
    distance: f64,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
) -> Result<(), ReplayError> {
    let failed = |error: String| ReplayError::FeatureFailed {
        label: "extrude".into(),
        error: error.into(),
    };
    // The same centroid arithmetic as `popkern_part::Profile::extrude`:
    // the mean of 64 wire samples, offset by the placement.
    let (t0, t1) = exact.wire.parameter_range();
    let n = 64;
    let mut acc = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..n {
        let t = t0 + (t1 - t0) * (i as f64) / (n as f64);
        let p = exact.wire.point_at(t);
        acc = acc + (p - Point3::new(0.0, 0.0, 0.0));
    }
    let mean = Point3::new(0.0, 0.0, 0.0) + acc * (1.0 / n as f64);
    let placement = exact.plane.u_axis * exact.center[0] + exact.plane.v_axis * exact.center[1];
    let start = mean + placement;
    let (a, b) = if distance >= 0.0 {
        (start, start + exact.plane.normal * distance)
    } else {
        (start + exact.plane.normal * distance, start)
    };
    let path = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    };
    let result = popkern_sweep::sweep(
        &popkern_sweep::SweepProfile::Wire(exact.wire.clone()),
        &path,
        &popkern_sweep::SweepOptions::default(),
        *tolerance,
        None,
    )
    .map_err(|e| failed(format!("exact profile sweep: {e}")))?;
    *brep = result.brep;
    crate::feature_snapshot::own_all_entities(brep, tracking);
    Ok(())
}
