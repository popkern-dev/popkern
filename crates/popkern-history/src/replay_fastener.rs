//! Replay arms for the fastener features (PK-460):
//! [`Feature::Fastener`](crate::feature::Feature::Fastener) bodies
//! plus the `TappedHole` / `ClearanceHole` cuts, executed through
//! [`popkern_fastener`].

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fastener::builders;
use popkern_fastener::{ClearanceHoleOptions, HoleDepth, TapOptions, ThreadSpec};

use crate::error::ReplayError;
use crate::feature::AxisDesc;
use crate::feature_fastener::{ClearanceFitDesc, FastenerKind, ThreadDesc};
use crate::parameter::{ParameterId, ParameterTable};
use crate::replay::replay_fns::resolve_param;
use crate::replay_side_effects::ExecuteSideEffects;
use crate::resolve::resolve_axis;

fn thread_spec(
    thread: &ThreadDesc,
    param_table: &ParameterTable,
) -> Result<ThreadSpec, ReplayError> {
    let mut spec = ThreadSpec::custom(
        resolve_param(param_table, thread.major_diameter)?,
        resolve_param(param_table, thread.pitch)?,
    )
    .with_starts(thread.starts)
    .with_millimetre(millimetre_of(param_table, thread));
    if thread.left_handed {
        spec = spec.left_handed();
    }
    Ok(spec)
}

/// Builder length units per millimetre for this thread's document
/// (PK-462). Parameters are stored in SI. A major diameter bound in a
/// real length unit (mm, cm, in) is a true SI length, so a millimetre
/// is 0.001; one bound in metres follows the default-unit convention
/// that fastener sizes are millimetre MAGNITUDES (an "M6" records 6.0),
/// so a millimetre is 1.0. The ISO size tables are looked up in
/// millimetres through this scale.
fn millimetre_of(param_table: &ParameterTable, thread: &ThreadDesc) -> f64 {
    use crate::unit::Unit;
    match param_table.get(thread.major_diameter).map(|p| p.unit) {
        Some(Unit::Millimeters | Unit::Centimeters | Unit::Inches) => 1.0e-3,
        _ => 1.0,
    }
}

fn fastener_failed(label: &'static str, e: &popkern_fastener::FastenerError) -> ReplayError {
    ReplayError::FeatureFailed {
        label: label.into(),
        error: crate::failure::FeatureFailure::other(e.to_string()).into(),
    }
}

fn optional_param(
    table: &ParameterTable,
    id: Option<ParameterId>,
) -> Result<Option<f64>, ReplayError> {
    id.map(|id| resolve_param(table, id)).transpose()
}

/// Dispatch the three fastener feature arms (single entry so the
/// main replay match stays one arm).
pub(super) fn execute_any(
    feature: &crate::feature::Feature,
    brep: &mut BRep,
    tolerance: &Tolerance,
    param_table: &ParameterTable,
) -> Result<ExecuteSideEffects, ReplayError> {
    use crate::feature::Feature;
    match feature {
        Feature::Fastener {
            kind,
            thread,
            length,
            thread_length,
            clearance,
            tip_chamfer,
            corner_chamfer,
            mouth_chamfer,
        } => execute_fastener(
            brep,
            *kind,
            thread,
            *length,
            *thread_length,
            *clearance,
            FastenerChamfers {
                tip: *tip_chamfer,
                corner: *corner_chamfer,
                mouth: *mouth_chamfer,
            },
            tolerance,
            param_table,
        ),
        Feature::TappedHole {
            axis,
            thread,
            depth,
            clearance,
            entry_chamfer,
        } => execute_tapped_hole(
            brep,
            axis,
            thread,
            *depth,
            *clearance,
            *entry_chamfer,
            tolerance,
            param_table,
        ),
        Feature::ClearanceHole {
            axis,
            thread,
            fit,
            depth,
            countersink,
            counterbore,
            entry_chamfer,
        } => execute_clearance_hole(
            brep,
            axis,
            thread,
            fit,
            *depth,
            *countersink,
            *counterbore,
            *entry_chamfer,
            tolerance,
            param_table,
        ),
        _ => Err(ReplayError::FeatureFailed {
            label: "fastener".into(),
            error: crate::failure::FeatureFailure::other(
                "not a fastener feature (dispatch invariant)".to_string(),
            )
            .into(),
        }),
    }
}

/// The `Feature::Fastener` arm: build a fastener body in the
/// insertion frame (new body — the incoming brep is replaced).
/// The chamfer toggles a fastener feature records: tip (screws,
/// rods), hex corners and thread mouths (nuts).
#[derive(Clone, Copy)]
struct FastenerChamfers {
    tip: bool,
    corner: bool,
    mouth: bool,
}

#[allow(clippy::too_many_arguments)]
fn execute_fastener(
    brep: &mut BRep,
    kind: FastenerKind,
    thread: &ThreadDesc,
    length: ParameterId,
    thread_length: Option<ParameterId>,
    clearance: Option<ParameterId>,
    chamfers: FastenerChamfers,
    tolerance: &Tolerance,
    param_table: &ParameterTable,
) -> Result<ExecuteSideEffects, ReplayError> {
    let tip_chamfer = chamfers.tip;
    let spec = thread_spec(thread, param_table)?;
    let length = resolve_param(param_table, length)?;
    let thread_length = optional_param(param_table, thread_length)?;
    let clearance = optional_param(param_table, clearance)?.unwrap_or(0.0);
    let label = kind.spelling();
    let built = match kind {
        FastenerKind::ThreadedRod => builders::threaded_rod(
            &builders::RodParams {
                spec,
                length,
                tip_chamfer,
            },
            *tolerance,
        ),
        FastenerKind::HexBolt => builders::hex_bolt(
            &builders::HexBoltParams {
                spec,
                length,
                thread_length,
                across_flats: None,
                head_height: None,
                head_chamfer: true,
                tip_chamfer,
            },
            *tolerance,
        ),
        FastenerKind::SocketCapScrew => builders::socket_cap_screw(
            &builders::SocketCapParams {
                spec,
                length,
                thread_length,
                tip_chamfer,
            },
            *tolerance,
        ),
        FastenerKind::CountersunkScrew => builders::countersunk_screw(
            &builders::CountersunkParams {
                spec,
                length,
                thread_length,
                tip_chamfer,
            },
            *tolerance,
        ),
        FastenerKind::HexNut => builders::hex_nut(
            &builders::NutParams {
                spec,
                thickness: (length > 0.0).then_some(length),
                across_flats: None,
                clearance,
                corner_chamfer: chamfers.corner,
                mouth_chamfer: chamfers.mouth,
            },
            *tolerance,
        ),
    }
    .map_err(|e| fastener_failed(label, &e))?;
    *brep = built;
    Ok(ExecuteSideEffects::empty())
}

/// The `Feature::TappedHole` arm: cut a threaded hole into the body.
#[allow(clippy::too_many_arguments)]
fn execute_tapped_hole(
    brep: &mut BRep,
    axis: &AxisDesc,
    thread: &ThreadDesc,
    depth: Option<ParameterId>,
    clearance: ParameterId,
    entry_chamfer: bool,
    tolerance: &Tolerance,
    param_table: &ParameterTable,
) -> Result<ExecuteSideEffects, ReplayError> {
    let spec = thread_spec(thread, param_table)?;
    let (origin, direction) = resolve_axis(axis, brep)?;
    let depth = match optional_param(param_table, depth)? {
        Some(d) => HoleDepth::Blind(d),
        None => HoleDepth::Through,
    };
    let options = TapOptions {
        clearance: resolve_param(param_table, clearance)?,
        entry_chamfer,
    };
    let result = popkern_fastener::holes::tapped_hole(
        brep, &spec, origin, direction, depth, &options, *tolerance,
    )
    .map_err(|e| fastener_failed("tapped_hole", &e))?;
    *brep = result;
    Ok(ExecuteSideEffects::empty())
}

/// The `Feature::ClearanceHole` arm: cut a clearance hole into the
/// body.
#[allow(clippy::too_many_arguments)]
fn execute_clearance_hole(
    brep: &mut BRep,
    axis: &AxisDesc,
    thread: &ThreadDesc,
    fit: &ClearanceFitDesc,
    depth: Option<ParameterId>,
    countersink: Option<ParameterId>,
    counterbore: Option<ParameterId>,
    entry_chamfer: bool,
    tolerance: &Tolerance,
    param_table: &ParameterTable,
) -> Result<ExecuteSideEffects, ReplayError> {
    let spec = thread_spec(thread, param_table)?;
    let (origin, direction) = resolve_axis(axis, brep)?;
    let depth = match optional_param(param_table, depth)? {
        Some(d) => HoleDepth::Blind(d),
        None => HoleDepth::Through,
    };
    let fit = match fit {
        ClearanceFitDesc::Close => popkern_fastener::ClearanceFit::Close,
        ClearanceFitDesc::Normal => popkern_fastener::ClearanceFit::Normal,
        ClearanceFitDesc::Loose => popkern_fastener::ClearanceFit::Loose,
        ClearanceFitDesc::Custom(id) => {
            popkern_fastener::ClearanceFit::Custom(resolve_param(param_table, *id)?)
        }
    };
    let options = ClearanceHoleOptions {
        fit,
        countersink: optional_param(param_table, countersink)?,
        counterbore: optional_param(param_table, counterbore)?,
        entry_chamfer,
    };
    let result = popkern_fastener::holes::clearance_hole(
        brep, &spec, origin, direction, depth, &options, *tolerance,
    )
    .map_err(|e| fastener_failed("clearance_hole", &e))?;
    *brep = result;
    Ok(ExecuteSideEffects::empty())
}
