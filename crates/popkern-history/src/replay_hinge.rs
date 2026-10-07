//! Replay arms for the hinge features (PK-693):
//! [`Feature::Hinge`](crate::feature::Feature::Hinge) components and
//! [`Feature::LivingHinge`](crate::feature::Feature::LivingHinge)
//! bodies, built through [`popkern_hinge`]. Both are new bodies: the
//! incoming brep is replaced.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_hinge::{
    HingeComponent, HingeHoles, HingeKind, HingeParams, LivingHingeParams, LivingHingeStyle,
};

use crate::error::ReplayError;
use crate::feature::Feature;
use crate::feature_hinge::{
    HingeComponentDesc, HingeDesc, HingeKindDesc, LivingHingeDesc, LivingHingeStyleDesc,
    SeamHingeDesc,
};
use crate::parameter::{ParameterId, ParameterTable};
use crate::replay::replay_fns::resolve_param;
use crate::replay_side_effects::ExecuteSideEffects;

fn failed(label: &str, detail: String) -> ReplayError {
    ReplayError::FeatureFailed {
        label: label.into(),
        error: crate::failure::FeatureFailure::other(detail).into(),
    }
}

/// A count parameter rounded to a whole number, like pattern counts.
fn count(table: &ParameterTable, id: ParameterId, role: &str) -> Result<u32, ReplayError> {
    let v = resolve_param(table, id)?.round();
    if !(v.is_finite() && v >= 0.0 && v <= f64::from(u32::MAX)) {
        return Err(failed("hinge", format!("{role} = {v} is not a count")));
    }
    Ok(v as u32)
}

fn hinge_params(spec: &HingeDesc, table: &ParameterTable) -> Result<HingeParams, ReplayError> {
    let p = |id| resolve_param(table, id);
    let holes = match &spec.holes {
        Some(h) => Some(HingeHoles {
            per_leaf: count(table, h.per_leaf, "holes")?,
            thread: super::replay_fastener::thread_spec(&h.thread, table)?,
            countersink: h.countersink,
        }),
        None => None,
    };
    Ok(HingeParams {
        length: p(spec.length)?,
        leaf_width: p(spec.leaf_width)?,
        thickness: p(spec.thickness)?,
        knuckles: count(table, spec.knuckles, "knuckles")?,
        pin_diameter: p(spec.pin_diameter)?,
        knuckle_diameter: spec.knuckle_diameter.map(p).transpose()?,
        clearance: p(spec.clearance)?,
        gap: p(spec.gap)?,
        // Angle parameters are stored in radians.
        angle_deg: p(spec.angle)?.to_degrees(),
        holes,
    })
}

/// A seam hinge's parameters as the plate-hinge builder takes them (no
/// leaf width or open angle: the bodies are the leaves).
fn seam_params(spec: &SeamHingeDesc, table: &ParameterTable) -> Result<HingeParams, ReplayError> {
    let p = |id| resolve_param(table, id);
    Ok(HingeParams {
        length: p(spec.length)?,
        leaf_width: 0.0,
        thickness: p(spec.wall)?,
        knuckles: count(table, spec.knuckles, "knuckles")?,
        pin_diameter: p(spec.pin_diameter)?,
        knuckle_diameter: spec.knuckle_diameter.map(p).transpose()?,
        clearance: p(spec.clearance)?,
        gap: p(spec.gap)?,
        angle_deg: 180.0,
        holes: None,
    })
}

fn living_params(
    style: LivingHingeStyleDesc,
    spec: &LivingHingeDesc,
    table: &ParameterTable,
) -> Result<LivingHingeParams, ReplayError> {
    let p = |id| resolve_param(table, id);
    let style_param = |id: Option<ParameterId>, role: &str| match id {
        Some(id) => p(id),
        None => Err(failed(
            "living_hinge",
            format!("the {role} parameter is missing for this style"),
        )),
    };
    let (style, web, root, pitch, width, bridge) = match style {
        LivingHingeStyleDesc::Straight => (
            LivingHingeStyle::Straight,
            style_param(spec.web, "web")?,
            style_param(spec.root_radius, "root_radius")?,
            0.0,
            0.0,
            0.0,
        ),
        LivingHingeStyleDesc::Slotted => (
            LivingHingeStyle::Slotted,
            0.0,
            0.0,
            style_param(spec.slot_pitch, "slot_pitch")?,
            style_param(spec.slot_width, "slot_width")?,
            style_param(spec.bridge, "bridge")?,
        ),
    };
    Ok(LivingHingeParams {
        style,
        length: p(spec.length)?,
        span: p(spec.span)?,
        panel_width: p(spec.panel_width)?,
        thickness: p(spec.thickness)?,
        web,
        root_radius: root,
        slot_pitch: pitch,
        slot_width: width,
        bridge,
    })
}

/// Dispatch the two hinge feature arms.
pub(super) fn execute_any(
    feature: &Feature,
    brep: &mut BRep,
    tolerance: &Tolerance,
    param_table: &ParameterTable,
) -> Result<ExecuteSideEffects, ReplayError> {
    let built = match feature {
        Feature::Hinge {
            kind,
            component,
            spec,
            ..
        } => {
            let kind = match kind {
                HingeKindDesc::Butt => HingeKind::Butt,
                HingeKindDesc::PrintInPlace => HingeKind::PrintInPlace,
            };
            let part = match component {
                HingeComponentDesc::LeafA => HingeComponent::LeafA,
                HingeComponentDesc::LeafB => HingeComponent::LeafB,
                HingeComponentDesc::Pin => HingeComponent::Pin,
            };
            let params = hinge_params(spec, param_table)?;
            popkern_hinge::hinge_component(kind, part, &params, *tolerance)
                .map_err(|e| failed(&format!("hinge {}", component.spelling()), e.to_string()))?
        }
        Feature::SeamHinge {
            kind,
            component,
            seam,
            spec,
            ..
        } => {
            let kind = match kind {
                HingeKindDesc::Butt => HingeKind::Butt,
                HingeKindDesc::PrintInPlace => HingeKind::PrintInPlace,
            };
            let part = match component {
                HingeComponentDesc::LeafA => HingeComponent::LeafA,
                HingeComponentDesc::LeafB => HingeComponent::LeafB,
                HingeComponentDesc::Pin => HingeComponent::Pin,
            };
            let params = seam_params(spec, param_table)?;
            let v = |a: [f64; 3]| popkern_core::geometry::Vec3::new(a[0], a[1], a[2]);
            let seam = popkern_hinge::Seam {
                origin: popkern_core::geometry::Point3::new(
                    seam.origin[0],
                    seam.origin[1],
                    seam.origin[2],
                ),
                axis: v(seam.axis),
                outward: v(seam.outward),
                up: v(seam.up),
            };
            let host = (part != HingeComponent::Pin).then_some(&*brep);
            popkern_hinge::seam_hinge_component(kind, part, &seam, &params, host, *tolerance)
                .map_err(|e| {
                    failed(
                        &format!("seam hinge {}", component.spelling()),
                        e.to_string(),
                    )
                })?
        }
        Feature::LivingHinge { style, spec } => {
            let params = living_params(*style, spec, param_table)?;
            popkern_hinge::living_hinge(&params, *tolerance)
                .map_err(|e| failed("living_hinge", e.to_string()))?
        }
        _ => {
            return Err(failed(
                "hinge",
                "not a hinge feature (dispatch invariant)".to_string(),
            ));
        }
    };
    *brep = built;
    Ok(ExecuteSideEffects::empty())
}
