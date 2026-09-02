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
use crate::feature::{ExtrudeSecondDir, Feature, SketchPlaneDesc, Vec3Desc};
use crate::feature_sketch::SketchPayload;
use crate::parameter::{ParameterId, ParameterTable};

/// The exact profile behind an extrude's face.
#[derive(Debug, Clone)]
pub(crate) struct ExactProfileWire {
    /// The unplaced wire lifted onto the plane (as `popkern-part` sweeps it).
    pub wire: NurbsCurve,
    /// The same wire in the plane's own 2-D frame (`z = 0`): what the
    /// sweep kernel's wire contract takes on any plane but XY.
    pub flat: NurbsCurve,
    /// Placement on the plane.
    pub center: [f64; 2],
    /// The sketch plane.
    pub plane: SketchPlane,
    /// The profile is a single full circle — the class whose generic
    /// extrude splits the lateral (PK-387) and whose drafted extrude
    /// has an exact arm (PK-395: draft the one periodic lateral).
    pub circular: bool,
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
    let Feature::Sketch { plane, sketch_data } = feature else {
        return None;
    };
    let profile = match sketch_data {
        SketchPayload::Analytic(profile) => {
            crate::feature_sketch::resolve_profile(profile, param_table).ok()?
        }
        SketchPayload::Constraint(sketch_json) => single_full_circle_profile(sketch_json)?,
    };
    let circular = matches!(profile.shape, popkern_sketch::ProfileShape::Circle { .. });
    let plane = crate::resolve::resolve_sketch_plane(plane, brep, param_table).ok()?;
    let wire = profile.local_wire(&plane).ok()?;
    let flat = profile.flat_wire().ok()?;
    let _ = SketchPlaneDesc::XY;
    Some(ExactProfileWire {
        wire,
        flat,
        center: profile.center,
        plane,
        circular,
    })
}

/// Whether this extrude takes the exact path: blind, no draft, no
/// second direction, direction absent or along the plane normal.
/// A solved constraint sketch whose only non-construction curve is a
/// full circle, as the exact circle profile. The sketch builder splits
/// such a circle into two arcs (its edge graph avoids self-loop edges),
/// so `BRep::extrude` of the compiled face makes two co-cylindrical
/// laterals; the exact wire sweep makes one periodic cylinder, the same
/// body `popkern-part` serves.
fn single_full_circle_profile(sketch_json: &str) -> Option<popkern_sketch::AnalyticProfile<f64>> {
    let mut sketch: popkern_constraint::Sketch = serde_json::from_str(sketch_json).ok()?;
    sketch.solve().ok()?;
    let curves = sketch.non_construction_curves();
    let [popkern_constraint::ElementHandle::Arc(arc)] = curves.as_slice() else {
        return None;
    };
    let (cx, cy, r, start_angle, end_angle) = sketch.arc_params(*arc)?;
    let full = ((end_angle - start_angle).abs() - std::f64::consts::TAU).abs() < 1.0e-9;
    if !full || !(r.is_finite() && r > 0.0) {
        return None;
    }
    Some(popkern_sketch::AnalyticProfile {
        shape: popkern_sketch::ProfileShape::Circle { radius: r },
        center: [cx, cy],
    })
}

/// Whether the exact wire sweep can serve this extrude: no draft on
/// either side (the sweep has no draft arm — that class stays on the
/// generic extrude and is recorded as a fallback), and every direction
/// along the sketch normal. The end condition is resolved to a distance
/// by the caller, so every end kind qualifies.
pub(crate) fn takes_exact_path(
    direction: Option<&Vec3Desc>,
    draft_angle: Option<&ParameterId>,
    second_direction: Option<&ExtrudeSecondDir>,
    circular: bool,
    plane: &SketchPlane,
) -> bool {
    if draft_angle.is_some() && !(circular && second_direction.is_none()) {
        // The wire sweep has no draft arm for general profiles; a full
        // circle's single periodic lateral is drafted exactly by the
        // same operator the generic path uses (PK-395).
        return false;
    }
    if let Some(second) = second_direction
        && (second.draft_angle.is_some() || !along_normal(Some(&second.direction), plane))
    {
        return false;
    }
    along_normal(direction, plane)
}

/// Whether `plane` is the world XY frame at some height: `u = x`,
/// `v = y`, `normal = z` (to 1e-12), so a wire lifted onto it is
/// already in the sweep's path-frame coordinates.
fn world_xy_frame(plane: &SketchPlane) -> bool {
    let close = |a: Vec3, b: Vec3| (a - b).length() < 1.0e-12;
    close(plane.u_axis, Vec3::new(1.0, 0.0, 0.0))
        && close(plane.v_axis, Vec3::new(0.0, 1.0, 0.0))
        && close(plane.normal, Vec3::new(0.0, 0.0, 1.0))
}

/// `None` (the face normal) or a literal parallel to the sketch normal.
fn along_normal(direction: Option<&Vec3Desc>, plane: &SketchPlane) -> bool {
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

/// Signed extent along the sketch normal for a resolved distance in
/// `direction` (`None` = the normal itself).
pub(crate) fn signed_extent(
    direction: Option<&Vec3Desc>,
    plane: &SketchPlane,
    distance: f64,
) -> f64 {
    match direction {
        Some(Vec3Desc::Literal(v)) if Vec3::new(v[0], v[1], v[2]).dot(&plane.normal) < 0.0 => {
            -distance
        }
        _ => distance,
    }
}

/// Serve the extrude through the exact periodic wire sweep when
/// [`takes_exact_path`] admits it. Every end condition is resolved to a
/// distance along the sketch normal (up-to-surface, offset, mid-plane,
/// through-all, …) and a second direction extends the same sweep the
/// other way: one lateral per smooth profile span by construction.
/// Returns `Ok(false)` when the generic face extrude must run instead.
#[allow(clippy::too_many_arguments)]
pub(crate) fn try_execute_exact(
    brep: &mut BRep,
    exact: &ExactProfileWire,
    face: &EntityName,
    direction: Option<&Vec3Desc>,
    end_condition: &crate::feature::ExtrudeEnd,
    draft_angle: Option<&ParameterId>,
    second_direction: Option<&ExtrudeSecondDir>,
    host_brep: Option<&BRep>,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
    param_table: &ParameterTable,
) -> Result<bool, ReplayError> {
    if !takes_exact_path(
        direction,
        draft_angle,
        second_direction,
        exact.circular,
        &exact.plane,
    ) {
        return Ok(false);
    }
    let draft = match draft_angle {
        Some(param) => {
            let angle = crate::resolve::resolve_param(param_table, *param)?;
            (angle.abs() > 1e-10).then_some(angle)
        }
        None => None,
    };
    let face_handle = crate::resolve::resolve_face(face, brep)?;
    let centroid = crate::replay::replay_fns::face_centroid(brep, face_handle);
    let dir = if signed_extent(direction, &exact.plane, 1.0) > 0.0 {
        exact.plane.normal
    } else {
        -exact.plane.normal
    };
    let dist = crate::resolve::resolve_extrude_distance(
        end_condition,
        centroid,
        dir,
        host_brep,
        brep,
        param_table,
        tolerance,
    )?;
    let forward = signed_extent(direction, &exact.plane, dist);
    // Mid-plane: the resolved distance is the half extent, both ways.
    let (mut lo, mut hi) = if matches!(end_condition, crate::feature::ExtrudeEnd::MidPlane { .. }) {
        (-forward.abs(), forward.abs())
    } else {
        (forward.min(0.0), forward.max(0.0))
    };
    if let Some(second) = second_direction {
        let dir2 = crate::resolve::resolve_direction(&second.direction, brep)?;
        let dist2 = crate::resolve::resolve_extrude_distance(
            &second.end_condition,
            centroid,
            dir2,
            host_brep,
            brep,
            param_table,
            tolerance,
        )?;
        let back = signed_extent(Some(&second.direction), &exact.plane, dist2);
        lo = lo.min(back);
        hi = hi.max(back);
    }
    execute_extrude_exact(
        brep,
        exact,
        face.operation,
        (lo, hi),
        draft,
        tracking,
        tolerance,
    )?;
    Ok(true)
}

pub(crate) fn execute_extrude_exact(
    brep: &mut BRep,
    exact: &ExactProfileWire,
    sketch_op: OperationId,
    span: (f64, f64),
    draft: Option<f64>,
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
    let (lo, hi) = (span.0.min(span.1), span.0.max(span.1));
    if !(hi.is_finite() && lo.is_finite()) || hi - lo <= 0.0 {
        return Err(failed(format!(
            "exact profile sweep: empty extent [{lo}, {hi}] along the sketch normal"
        )));
    }
    let line_path = |a: Point3, b: Point3| NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    };
    let sweep = |wire: &NurbsCurve, path: &NurbsCurve| {
        popkern_sweep::sweep(
            &popkern_sweep::SweepProfile::Wire(wire.clone()),
            path,
            &popkern_sweep::SweepOptions::default(),
            *tolerance,
            None,
        )
        .map_err(|e| failed(format!("exact profile sweep: {e}")))
    };
    // The sweep kernel's wire contract is the PATH FRAME's 2-D plane
    // (profile x → frame normal, y → binormal; z dropped). On the
    // world XY frame the lifted wire IS that plane — the very
    // arithmetic `popkern_part::Profile::extrude` performs, kept
    // bit-identical (parity fixtures). On any other plane the wire is
    // swept in the plane's own frame along local z and the body is
    // placed by the plane's rigid frame afterwards; lifting first
    // collapsed every span along the plane's v axis (PK-416).
    if world_xy_frame(&exact.plane) {
        let (a, b) = (
            start + exact.plane.normal * lo,
            start + exact.plane.normal * hi,
        );
        *brep = sweep(&exact.wire, &line_path(a, b))?.brep;
    } else {
        let (ft0, ft1) = exact.flat.parameter_range();
        let mut flat_acc = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..n {
            let t = ft0 + (ft1 - ft0) * (i as f64) / (n as f64);
            flat_acc = flat_acc + (exact.flat.point_at(t) - Point3::new(0.0, 0.0, 0.0));
        }
        let flat_mean = flat_acc * (1.0 / n as f64);
        let (cx, cy) = (flat_mean.x + exact.center[0], flat_mean.y + exact.center[1]);
        let local = sweep(
            &exact.flat,
            &line_path(Point3::new(cx, cy, lo), Point3::new(cx, cy, hi)),
        )?
        .brep;
        let placement = popkern_core::geometry::transform::Transform3::from_axis_placement(
            exact.plane.origin,
            exact.plane.normal,
            exact.plane.u_axis,
        )
        .map_err(|e| failed(format!("sketch plane frame: {e}")))?;
        *brep = local
            .transform(&placement)
            .map_err(|e| failed(format!("placing the swept body on its sketch plane: {e}")))?;
    }
    // Draft: the swept body has ONE periodic lateral; drafting it with
    // the SAME operator the generic path applies to its split laterals
    // yields one exact Cone (popkern-draft's cylindrical surface arm),
    // where the generic route serves two co-conical halves (PK-395).
    if let Some(angle) = draft {
        let laterals: Vec<popkern_core::topology::FaceHandle> = brep
            .faces
            .iter()
            .filter(|(_, f)| {
                matches!(
                    f.surface.as_ref().map(|s| &s.geometry),
                    Some(popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. })
                )
            })
            .map(|(h, _)| h)
            .collect();
        let pull = exact.plane.normal * if hi >= 0.0 { 1.0 } else { -1.0 };
        let specs: Vec<popkern_draft::DraftSpec> = laterals
            .iter()
            .map(|&face| popkern_draft::DraftSpec {
                face,
                angle,
                pull_direction: pull,
                neutral: popkern_draft::NeutralReference::Plane {
                    origin: start,
                    normal: pull,
                },
                angle_variation: popkern_draft::AngleVariation::Constant,
            })
            .collect();
        let drafted = popkern_draft::draft_faces(
            brep,
            &specs,
            &popkern_draft::DraftOptions::default(),
            *tolerance,
            None,
        )
        .map_err(|e| failed(format!("draft of the periodic lateral: {e}")))?;
        *brep = drafted.brep;
    }
    // The wire sweep places the profile at `start` (its sampled
    // centroid plus the placement); the name frame measures heights
    // from there along the sketch normal.
    let placed_wire = {
        let shift = start - Point3::new(0.0, 0.0, 0.0) - (mean - Point3::new(0.0, 0.0, 0.0));
        NurbsCurve {
            control_points: exact
                .wire
                .control_points
                .iter()
                .map(|p| *p + shift)
                .collect(),
            ..exact.wire.clone()
        }
    };
    crate::replay_exact_extrude_names::name_exact_extrude(
        brep,
        tracking,
        sketch_op,
        &placed_wire,
        start,
        exact.plane.normal,
        (lo, hi),
    );
    Ok(())
}
