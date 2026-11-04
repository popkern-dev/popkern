//! Replay of the piece features (design: operations on the cell
//! complex, 5.4): `KeepPieces` keeps the solids of the current body a
//! rule selects; `SeparatePieces` moves the solids of a source body a
//! rule selects into the feature's new body, the source keeping the
//! rest. Pieces are taken whole (`popkern_bool::pieces`): every kept
//! face, edge and vertex keeps its name, so features after this one
//! resolve the same entities.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::ReplayError;
use crate::failure::FeatureFailure;
use crate::feature::{Feature, PieceRule};
use crate::feature_pieces::select_pieces;
use crate::replay_report::FeatureEvidence;
use crate::replay_side_effects::ExecuteSideEffects;
use crate::working_bodies::WorkingBodies;

/// Replay a `KeepPieces` or `SeparatePieces` `feature` on `brep`.
pub(crate) fn execute(
    brep: &mut BRep,
    feature: &Feature,
    tolerance: &Tolerance,
    body_breps: &WorkingBodies,
) -> Result<ExecuteSideEffects, ReplayError> {
    match feature {
        Feature::KeepPieces { rule } => {
            let (kept, _) = partition(brep, "keep_pieces", rule, tolerance)?;
            *brep = kept;
            Ok(ExecuteSideEffects::empty())
        }
        Feature::SeparatePieces { source, rule } => {
            let label = "separate_pieces";
            let from = body_breps
                .get(source)
                .ok_or_else(|| failed(label, format!("source body {source:?} not found").into()))?;
            let (moved, rest) = partition(from, label, rule, tolerance)?;
            if rest.solids.is_empty() {
                return Err(failed(
                    label,
                    FeatureFailure::InvalidParameter {
                        parameter: "rule".into(),
                        value: format!("{rule:?}"),
                        reason: "the rule selects every piece: the source body would be empty"
                            .into(),
                    }
                    .into(),
                ));
            }
            *brep = moved;
            Ok(ExecuteSideEffects {
                consumed_bodies: Vec::new(),
                evidence: FeatureEvidence::empty(),
                written_bodies: vec![(*source, rest)],
            })
        }
        _ => Err(failed("pieces", "not a piece feature".into())),
    }
}

/// `body` split into the pieces `rule` selects and the others.
fn partition(
    body: &BRep,
    label: &str,
    rule: &PieceRule,
    tolerance: &Tolerance,
) -> Result<(BRep, BRep), ReplayError> {
    let handles = popkern_bool::pieces::piece_handles(body);
    let bodies = popkern_bool::pieces::piece_bodies(body).map_err(|e| boolean(label, &e))?;
    let refs: Vec<&BRep> = bodies.iter().collect();
    let kept = select_pieces(label, rule, &refs, None, tolerance).map_err(|e| failed(label, e))?;
    let (selected, others): (Vec<_>, Vec<_>) = handles
        .iter()
        .enumerate()
        .partition(|(i, _)| kept.contains(i));
    let take = |set: Vec<(usize, &popkern_core::topology::SolidHandle)>| {
        let handles: Vec<_> = set.into_iter().map(|(_, h)| *h).collect();
        popkern_bool::pieces::retain_pieces(body, &handles).map_err(|e| boolean(label, &e))
    };
    Ok((take(selected)?, take(others)?))
}

fn boolean(label: &str, error: &popkern_bool::BoolError) -> ReplayError {
    failed(label, FeatureFailure::boolean(label, error).into())
}

fn failed(label: &str, error: Box<FeatureFailure>) -> ReplayError {
    ReplayError::FeatureFailed {
        label: label.into(),
        error,
    }
}
