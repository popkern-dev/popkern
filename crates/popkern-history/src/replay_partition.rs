//! Replay of the partition's features (design: operations on the cell
//! complex, 5.4): a split keeping the pieces a rule selects, an N-ary
//! boolean, an imprint. Each consumes its tools, as a boolean consumes
//! its tool, and names its result after the target and the tools.

use popkern_bool::{BoolError, BoolOptions, PointClassification};
use popkern_core::geometry::Point3;
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::ReplayError;
use crate::failure::FeatureFailure;
use crate::feature::{Feature, PieceRule};
use crate::replay_report::FeatureEvidence;
use crate::replay_side_effects::ExecuteSideEffects;
use crate::working_bodies::WorkingBodies;

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
    let mut result = match feature {
        Feature::Split { keep, .. } => split(brep, &refs, keep, tolerance)?,
        Feature::NaryBoolean {
            expr, non_manifold, ..
        } => {
            let mut bodies = vec![&*brep];
            bodies.extend_from_slice(&refs);
            let options = BoolOptions {
                non_manifold: *non_manifold,
            };
            popkern_bool::nary_boolean(&bodies, expr, *tolerance, &options)
                .map_err(|e| failed(label, &e))?
        }
        _ => {
            popkern_bool::ops::imprint::imprint(brep, &tool_breps[0], tolerance)
                .map_err(|e| failed(label, &e))?
                .imprinted_target
        }
    };
    popkern_bool::name_from_operands(&mut result, brep, &refs, tolerance, tracking);
    let mut evidence = FeatureEvidence::empty();
    evidence.record_boolean_fallbacks(popkern_bool::merge_fallback::take_merge_fallback_events());
    *brep = result;
    Ok(ExecuteSideEffects {
        consumed_bodies: tools,
        evidence,
    })
}

/// `target` cut by `tools`, the pieces `keep` selects as one body of
/// several solids.
fn split(
    target: &BRep,
    tools: &[&BRep],
    keep: &PieceRule,
    tolerance: &Tolerance,
) -> Result<BRep, ReplayError> {
    let cut = popkern_bool::split(target, tools, *tolerance).map_err(|e| failed("split", &e))?;
    let count = cut.pieces.len();
    let kept: Vec<usize> = match keep {
        PieceRule::All => (0..count).collect(),
        PieceRule::InsideTools(expr) => (0..count)
            .filter(|&i| expr.holds(&cut.pieces[i].inside_tools, tools.len()))
            .collect(),
        PieceRule::Containing([x, y, z]) => {
            let p = Point3::new(*x, *y, *z);
            (0..count)
                .filter(|&i| {
                    let body = &cut.pieces[i].body;
                    let boxes = popkern_bool::adaptive_precision::body_face_aabbs(body);
                    popkern_bool::adaptive_precision::robust_point_in_solid(
                        &p, body, &boxes, tolerance,
                    ) == PointClassification::Inside
                })
                .collect()
        }
        PieceRule::Largest(n) => {
            let mut by_volume: Vec<(usize, f64)> = (0..count)
                .map(|i| {
                    let v =
                        popkern_measure::body_volume(&cut.pieces[i].body, tolerance).unwrap_or(0.0);
                    (i, v)
                })
                .collect();
            // Largest first; equal volumes keep the partition's order.
            by_volume.sort_by(|a, b| b.1.total_cmp(&a.1));
            let mut kept: Vec<usize> = by_volume.into_iter().take(*n).map(|(i, _)| i).collect();
            kept.sort_unstable();
            kept
        }
    };
    if kept.is_empty() {
        return Err(ReplayError::FeatureFailed {
            label: "split".into(),
            error: FeatureFailure::PieceRuleMatchedNothing {
                feature: "split".into(),
                rule: format!("{keep:?}"),
                pieces: count,
            }
            .into(),
        });
    }
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

fn failed(label: &str, error: &BoolError) -> ReplayError {
    failed_with(label, FeatureFailure::boolean(label, error).into())
}

fn failed_with(label: &str, error: Box<FeatureFailure>) -> ReplayError {
    ReplayError::FeatureFailed {
        label: label.into(),
        error,
    }
}
