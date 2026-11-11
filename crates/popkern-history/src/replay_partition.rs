//! Replay of the two-operand boolean and the partition's features
//! (design: operations on the cell complex, 5.4): a split keeping the
//! pieces a rule selects, an N-ary boolean, an imprint. Each consumes
//! its tools and names its result after the target and the tools.

use popkern_bool::{BoolError, BoolOptions};
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::ReplayError;
use crate::failure::FeatureFailure;
use crate::feature::{Feature, PieceRule};
use crate::replay_report::FeatureEvidence;
use crate::replay_side_effects::ExecuteSideEffects;
use crate::working_bodies::WorkingBodies;

/// Replay a two-operand boolean `feature` on `brep`.
pub(crate) fn execute_boolean(
    brep: &mut BRep,
    feature: &Feature,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
    body_breps: &WorkingBodies,
) -> Result<ExecuteSideEffects, ReplayError> {
    let Feature::Boolean {
        operation,
        tool_body,
        non_manifold,
    } = feature
    else {
        return Err(failed_with("boolean", "not a boolean feature".into()));
    };
    let tool_brep = body_breps.get(tool_body).cloned().ok_or_else(|| {
        failed_with(
            "boolean",
            format!("tool body {tool_body:?} not found").into(),
        )
    })?;
    // Drain stale registry entries so the firings taken after the call
    // belong to this boolean.
    let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
    let _ = popkern_bool::sub_tolerance::take_sub_tolerance_events();
    let options = BoolOptions {
        non_manifold: *non_manifold,
    };
    let bodies = popkern_bool::boolean_bodies(
        crate::replay::replay_fns::to_bool_operation(*operation),
        brep.clone(),
        tool_brep,
        *tolerance,
        Some(tracking),
        &options,
    )
    .map_err(|e| {
        let op = format!("{operation:?}").to_lowercase();
        failed_with("boolean", FeatureFailure::boolean(&op, &e).into())
    })?;
    let mut evidence = FeatureEvidence::empty();
    evidence.record_boolean_fallbacks(popkern_bool::merge_fallback::take_merge_fallback_events());
    evidence.record_sub_tolerance(popkern_bool::sub_tolerance::take_sub_tolerance_events());
    Ok(ExecuteSideEffects {
        consumed_bodies: vec![*tool_body],
        evidence,
        written_bodies: Vec::new(),
        served_bodies: bodies,
    })
}

/// Replay `feature` (a split, N-ary boolean or imprint) on `brep`.
pub(crate) fn execute(
    brep: &mut BRep,
    feature: &Feature,
    tracking: &mut TrackingContext,
    tolerance: &Tolerance,
    body_breps: &WorkingBodies,
) -> Result<ExecuteSideEffects, ReplayError> {
    let (label, tools) = match feature {
        Feature::Split { tools, .. } => ("split", tools.clone()),
        Feature::NaryBoolean { tools, .. } => ("nary_boolean", tools.clone()),
        Feature::Imprint { tool } => ("imprint", vec![*tool]),
        _ => return Err(failed_with("partition", "not a partition feature".into())),
    };
    let tool_breps: Vec<BRep> = tools
        .iter()
        .map(|id| {
            body_breps
                .get(id)
                .cloned()
                .ok_or_else(|| failed_with(label, format!("tool body {id:?} not found").into()))
        })
        .collect::<Result<_, _>>()?;
    let refs: Vec<&BRep> = tool_breps.iter().collect();
    // Drain stale registry entries so the firings taken after the call
    // belong to this feature.
    let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
    let _ = popkern_bool::sub_tolerance::take_sub_tolerance_events();
    let mut classes = Vec::new();
    let (mut result, touching_pieces) = match feature {
        Feature::Split { keep, .. } => (split(brep, &refs, keep, tolerance)?, false),
        Feature::NaryBoolean {
            expr, non_manifold, ..
        } => {
            let mut bodies = vec![&*brep];
            bodies.extend_from_slice(&refs);
            let options = BoolOptions {
                non_manifold: *non_manifold,
            };
            let served = popkern_bool::nary_boolean_bodies(&bodies, expr, *tolerance, &options)
                .map_err(|e| failed(label, &e))?;
            all_solids(label, served)?
        }
        _ => {
            let imprinted = popkern_bool::ops::imprint::imprint(brep, &tool_breps[0], tolerance)
                .map_err(|e| failed(label, &e))?;
            classes = imprinted.classes;
            (imprinted.imprinted_target, false)
        }
    };
    popkern_bool::name_from_operands(&mut result, brep, &refs, tolerance, tracking);
    let mut evidence = FeatureEvidence::empty();
    if matches!(feature, Feature::Imprint { .. }) {
        // Naming keeps face handles: each class is keyed by the name
        // its face now carries.
        let classes = classes
            .into_iter()
            .filter_map(|(face, class)| {
                let id = result.faces.get(face)?.name?;
                Some((result.name_pool.get(id)?.clone(), class))
            })
            .collect();
        evidence.imprint = Some(crate::replay_report::ImprintRecord { classes });
    }
    evidence.record_boolean_fallbacks(popkern_bool::merge_fallback::take_merge_fallback_events());
    // A gap at or under the tolerance the selection closed (PK-1162).
    evidence.record_sub_tolerance(popkern_bool::sub_tolerance::take_sub_tolerance_events());
    // Named as one, then served as its pieces: a split serves a body
    // per kept piece, an N-ary boolean a body per touching piece
    // (PK-1120); an imprint changes its body in place.
    let served_bodies = match feature {
        Feature::Imprint { .. } => {
            *brep = result;
            Vec::new()
        }
        Feature::Split { .. } => {
            popkern_bool::pieces::piece_bodies(&result).map_err(|e| failed(label, &e))?
        }
        _ if touching_pieces => result
            .piece_bodies()
            .map_err(|e| failed_with(label, format!("{e}").into()))?,
        _ => vec![result],
    };
    Ok(ExecuteSideEffects {
        consumed_bodies: tools,
        evidence,
        written_bodies: Vec::new(),
        served_bodies,
    })
}

/// `target` cut by `tools`, the pieces `keep` selects as one body of
/// several solids (named as one, then served as a body per piece).
fn split(
    target: &BRep,
    tools: &[&BRep],
    keep: &PieceRule,
    tolerance: &Tolerance,
) -> Result<BRep, ReplayError> {
    let cut = popkern_bool::split(target, tools, *tolerance).map_err(|e| failed("split", &e))?;
    let bodies: Vec<&BRep> = cut.pieces.iter().map(|p| &p.body).collect();
    let inside: Vec<popkern_bool::BodySet> = cut.pieces.iter().map(|p| p.inside_tools).collect();
    let membership = Some((inside.as_slice(), tools.len()));
    let kept = crate::feature_pieces::select_pieces("split", keep, &bodies, membership, tolerance)
        .map_err(|e| failed_with("split", e))?;
    let mut pieces = cut
        .pieces
        .into_iter()
        .enumerate()
        .filter(|(i, _)| kept.contains(i));
    let Some((_, first)) = pieces.next() else {
        return Err(failed_with("split", "no piece kept".into()));
    };
    pieces.try_fold(first.body, |body, (_, piece)| {
        popkern_bool::concatenate_solids(body, piece.body).map_err(|e| failed("split", &e))
    })
}

/// Bodies served for one feature (split manifold pieces that touch)
/// as one body of several solids, to be named as one; `true` when
/// there were several, to be served again as their pieces.
fn all_solids(label: &str, bodies: Vec<BRep>) -> Result<(BRep, bool), ReplayError> {
    let several = bodies.len() > 1;
    let mut bodies = bodies.into_iter();
    let Some(first) = bodies.next() else {
        return Err(failed_with(label, "no body served".into()));
    };
    let body = bodies.try_fold(first, |body, piece| {
        popkern_bool::concatenate_solids(body, piece).map_err(|e| failed(label, &e))
    })?;
    Ok((body, several))
}

fn failed(label: &str, error: &BoolError) -> ReplayError {
    failed_with(label, FeatureFailure::boolean(label, error).into())
}

fn failed_with(label: &str, error: Box<FeatureFailure>) -> ReplayError {
    ReplayError::FeatureFailed {
        label: label.into(),
        error,
    }
}
