//! Round-trip emit for the hinge features (PK-693). The owning
//! component emits `m.hinge(kind, length, #{ only: …, … }).<component>`;
//! every other component emits `m.hinge_component("<owner>", "<component>")`
//! so it shares the owner's parameters, as recorded. Options still at
//! their default — a dimension following the length, a default literal —
//! are left out: the call re-derives exactly that default.

use popkern_history::FeatureTree;
use popkern_history::expression::{BinOp, Expression};
use popkern_history::feature::Feature;
use popkern_history::feature_hinge::{
    HingeDesc, HingeKindDesc, LivingHingeDesc, LivingHingeStyleDesc,
};
use popkern_history::parameter::ParameterId;

use super::{Statement, label_option, options};
use crate::emit::EmitError;
use crate::emit_spelling::{quoted, scalar};

/// Whether `id` is still the default that follows `of × factor`.
fn follows(tree: &FeatureTree, id: ParameterId, of: ParameterId, factor: f64) -> bool {
    matches!(
        tree.parameters().get(id).and_then(|p| p.expression.as_ref()),
        Some(Expression::BinOp(a, BinOp::Mul, b))
            if matches!(**a, Expression::Param(x) if x == of)
                && matches!(**b, Expression::Literal(f) if f.to_bits() == factor.to_bits())
    )
}

/// Whether `id` is an undriven literal equal to `value` (display units).
fn is_literal(tree: &FeatureTree, id: ParameterId, value: f64) -> bool {
    tree.parameters()
        .get(id)
        .is_some_and(|p| p.expression.is_none() && (p.unit.from_si(p.value) - value).abs() < 1e-12)
}

fn push_unless(
    pairs: &mut Vec<(&'static str, String)>,
    key: &'static str,
    default: bool,
    tree: &FeatureTree,
    id: ParameterId,
    context: &str,
) -> Result<(), EmitError> {
    if !default {
        pairs.push((key, scalar(tree, id, context)?));
    }
    Ok(())
}

fn hinge_pairs(
    tree: &FeatureTree,
    kind: HingeKindDesc,
    s: &HingeDesc,
    context: &str,
) -> Result<Vec<(&'static str, String)>, EmitError> {
    let print = kind == HingeKindDesc::PrintInPlace;
    let l = s.length;
    let mut pairs = Vec::new();
    push_unless(
        &mut pairs,
        "leaf_width",
        follows(tree, s.leaf_width, l, 0.25),
        tree,
        s.leaf_width,
        context,
    )?;
    push_unless(
        &mut pairs,
        "thickness",
        follows(tree, s.thickness, l, 1.0 / 30.0),
        tree,
        s.thickness,
        context,
    )?;
    push_unless(
        &mut pairs,
        "knuckles",
        is_literal(tree, s.knuckles, 5.0),
        tree,
        s.knuckles,
        context,
    )?;
    push_unless(
        &mut pairs,
        "pin",
        follows(tree, s.pin_diameter, l, 1.0 / 15.0),
        tree,
        s.pin_diameter,
        context,
    )?;
    if let Some(id) = s.knuckle_diameter {
        pairs.push(("knuckle_diameter", scalar(tree, id, context)?));
    }
    let clearance_default = if print {
        follows(tree, s.clearance, l, 1.0 / 200.0)
    } else {
        is_literal(tree, s.clearance, 0.0)
    };
    push_unless(
        &mut pairs,
        "clearance",
        clearance_default,
        tree,
        s.clearance,
        context,
    )?;
    let gap_factor = if print { 1.0 / 200.0 } else { 1.0 / 400.0 };
    push_unless(
        &mut pairs,
        "gap",
        follows(tree, s.gap, l, gap_factor),
        tree,
        s.gap,
        context,
    )?;
    push_unless(
        &mut pairs,
        "angle",
        is_literal(tree, s.angle, 180.0),
        tree,
        s.angle,
        context,
    )?;
    if let Some(h) = &s.holes {
        pairs.push((
            "screw",
            format!(
                "#{{ major: {}, pitch: {} }}",
                scalar(tree, h.thread.major_diameter, context)?,
                scalar(tree, h.thread.pitch, context)?
            ),
        ));
        push_unless(
            &mut pairs,
            "holes",
            is_literal(tree, h.per_leaf, 2.0),
            tree,
            h.per_leaf,
            context,
        )?;
        if !h.countersink {
            pairs.push(("countersink", "false".to_string()));
        }
    }
    Ok(pairs)
}

/// The statement for one hinge component.
pub(super) fn hinge_statement(
    tree: &FeatureTree,
    context: &str,
    label: &str,
    feature: &Feature,
) -> Result<Statement, EmitError> {
    if matches!(feature, Feature::LivingHinge { .. }) {
        return living_hinge_statement(tree, context, label, feature);
    }
    let Feature::Hinge {
        kind,
        component,
        spec,
        owner,
    } = feature
    else {
        return Err(EmitError::not_applicable(
            label,
            "hinge",
            "not a hinge feature",
        ));
    };
    let call = if *owner {
        let mut pairs = hinge_pairs(tree, *kind, spec, context)?;
        pairs.push(("only", quoted(component.spelling())));
        pairs.push(label_option(label));
        let kind_text = match kind {
            HingeKindDesc::Butt => "butt",
            HingeKindDesc::PrintInPlace => "print_in_place",
        };
        format!(
            "m.hinge({}, {}{}).{}",
            quoted(kind_text),
            scalar(tree, spec.length, context)?,
            options(pairs),
            component.spelling()
        )
    } else {
        let owner_label = owner_label(tree, spec).ok_or_else(|| {
            EmitError::not_applicable(label, "hinge", "the hinge's owning component is gone")
        })?;
        format!(
            "m.hinge_component({}, {}{})",
            quoted(&owner_label),
            quoted(component.spelling()),
            options(vec![label_option(label)])
        )
    };
    Ok(Statement {
        call,
        starts_body: true,
    })
}

/// The label of the feature that owns `spec`'s parameters.
fn owner_label(tree: &FeatureTree, spec: &HingeDesc) -> Option<String> {
    tree.timeline().iter().find_map(|&id| {
        let node = tree.feature(id)?;
        match &node.feature {
            Feature::Hinge {
                spec: s,
                owner: true,
                ..
            } if s == spec => Some(node.label.clone()),
            _ => None,
        }
    })
}

/// The `m.living_hinge(..)` statement.
fn living_hinge_statement(
    tree: &FeatureTree,
    context: &str,
    label: &str,
    feature: &Feature,
) -> Result<Statement, EmitError> {
    let Feature::LivingHinge { style, spec } = feature else {
        return Err(EmitError::not_applicable(
            label,
            "living_hinge",
            "not a living hinge",
        ));
    };
    let s: &LivingHingeDesc = spec;
    let mut pairs = Vec::new();
    push_unless(
        &mut pairs,
        "panel_width",
        follows(tree, s.panel_width, s.span, 2.0),
        tree,
        s.panel_width,
        context,
    )?;
    push_unless(
        &mut pairs,
        "thickness",
        follows(tree, s.thickness, s.span, 0.2),
        tree,
        s.thickness,
        context,
    )?;
    let optional = [
        ("web", s.web, s.thickness, 0.2),
        ("root_radius", s.root_radius, s.thickness, 0.4),
        ("slot_pitch", s.slot_pitch, s.span, 0.125),
        ("slot_width", s.slot_width, s.span, 1.0 / 24.0),
        ("bridge", s.bridge, s.length, 0.1),
    ];
    for (key, id, of, factor) in optional {
        if let Some(id) = id {
            push_unless(
                &mut pairs,
                key,
                follows(tree, id, of, factor),
                tree,
                id,
                context,
            )?;
        }
    }
    pairs.push(label_option(label));
    let style_text = match style {
        LivingHingeStyleDesc::Straight => "straight",
        LivingHingeStyleDesc::Slotted => "slotted",
    };
    Ok(Statement {
        call: format!(
            "m.living_hinge({}, {}, {}{})",
            quoted(style_text),
            scalar(tree, s.span, context)?,
            scalar(tree, s.length, context)?,
            options(pairs)
        ),
        starts_body: true,
    })
}

fn vec3_text(v: [f64; 3]) -> String {
    format!(
        "[{}, {}, {}]",
        crate::emit_spelling::float(v[0]),
        crate::emit_spelling::float(v[1]),
        crate::emit_spelling::float(v[2])
    )
}

/// The label of the seam hinge feature owning `spec`'s parameters.
fn seam_owner_label(
    tree: &FeatureTree,
    spec: &popkern_history::feature_hinge::SeamHingeDesc,
) -> Option<String> {
    tree.timeline().iter().find_map(|&id| {
        let node = tree.feature(id)?;
        match &node.feature {
            Feature::SeamHinge {
                spec: s,
                owner: true,
                ..
            } if s == spec => Some(node.label.clone()),
            _ => None,
        }
    })
}

/// A seam hinge component: body A's knuckles chain `.seam_hinge(…)` on
/// its body with the seam spelled out; body B's chain
/// `.seam_hinge_component(owner, "leaf_b")`; the pin is the model-level
/// `m.seam_hinge_component(owner, "pin")` (`Err(statement)` marks it).
pub(super) fn seam_hinge_call(
    tree: &FeatureTree,
    context: &str,
    label: &str,
    feature: &Feature,
) -> Result<Result<String, Statement>, EmitError> {
    use popkern_history::feature_hinge::HingeComponentDesc;
    let Feature::SeamHinge {
        kind,
        component,
        seam,
        spec,
        owner,
    } = feature
    else {
        return Err(EmitError::not_applicable(
            label,
            "seam_hinge",
            "not a seam hinge",
        ));
    };
    if *owner {
        let print = *kind == HingeKindDesc::PrintInPlace;
        let l = spec.length;
        let factor = if print { 1.0 / 200.0 } else { 1.0 / 400.0 };
        let mut pairs = vec![
            ("outward", vec3_text(seam.outward)),
            ("up", vec3_text(seam.up)),
        ];
        push_unless(
            &mut pairs,
            "knuckles",
            is_literal(tree, spec.knuckles, 5.0),
            tree,
            spec.knuckles,
            context,
        )?;
        push_unless(
            &mut pairs,
            "pin",
            follows(tree, spec.pin_diameter, l, 1.0 / 15.0),
            tree,
            spec.pin_diameter,
            context,
        )?;
        push_unless(
            &mut pairs,
            "wall",
            follows(tree, spec.wall, l, 1.0 / 30.0),
            tree,
            spec.wall,
            context,
        )?;
        if let Some(id) = spec.knuckle_diameter {
            pairs.push(("knuckle_diameter", scalar(tree, id, context)?));
        }
        push_unless(
            &mut pairs,
            "clearance",
            follows(tree, spec.clearance, l, factor),
            tree,
            spec.clearance,
            context,
        )?;
        push_unless(
            &mut pairs,
            "gap",
            follows(tree, spec.gap, l, factor),
            tree,
            spec.gap,
            context,
        )?;
        if *component != HingeComponentDesc::LeafA {
            pairs.push(("component", quoted(component.spelling())));
        }
        pairs.push(label_option(label));
        let kind_text = if print { "print_in_place" } else { "butt" };
        let call = format!(
            "seam_hinge({}, axis_at({}, {}), {}{})",
            quoted(kind_text),
            vec3_text(seam.origin),
            vec3_text(seam.axis),
            scalar(tree, l, context)?,
            options(pairs)
        );
        // A pin that owns its parameters is a new body: model level.
        return Ok(if *component == HingeComponentDesc::Pin {
            Err(Statement {
                call: format!("m.{call}"),
                starts_body: true,
            })
        } else {
            Ok(format!(".{call}"))
        });
    }
    let owner_label = seam_owner_label(tree, spec).ok_or_else(|| {
        EmitError::not_applicable(label, "seam_hinge", "the hinge's owning component is gone")
    })?;
    let tail = format!(
        "seam_hinge_component({}, {}{})",
        quoted(&owner_label),
        quoted(component.spelling()),
        options(vec![label_option(label)])
    );
    if *component == HingeComponentDesc::Pin {
        Ok(Err(Statement {
            call: format!("m.{tail}"),
            starts_body: true,
        }))
    } else {
        Ok(Ok(format!(".{tail}")))
    }
}
