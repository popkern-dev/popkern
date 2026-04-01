//! Top-level pattern execution engine.
//!
//! Orchestrates the full pattern lifecycle: compute transforms, replay or
//! copy-transform each instance, boolean-merge into the host body, and
//! collect per-instance outcomes.

use std::collections::HashMap;

use popkern_core::history::TrackingContext;
use popkern_core::naming::OperationId;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::auto_merge;
use crate::error::ReplayError;
use crate::feature::{Feature, OperationMode};
use crate::node::{BodyId, FeatureId, FeatureNode};
use crate::parameter::ParameterTable;
use crate::pattern::{InstanceOutcome, PatternResult};
use crate::pattern_replay;
use crate::pattern_transform;

/// Aggregated context for pattern execution, bundling the many parameters
/// needed by [`execute_pattern`] into a single struct to satisfy clippy's
/// `too_many_arguments` lint.
pub(crate) struct PatternContext<'a> {
    /// The pattern feature definition (cloned from the tree).
    pub pattern_feature: &'a Feature,
    /// The seed feature definitions (cloned from the tree).
    pub seed_features: &'a [Feature],
    /// OperationIds of the seed features.
    pub seed_op_ids: &'a [OperationId],
    /// Existing instance OperationId assignments.
    pub instance_op_ids: &'a mut HashMap<(u32, u32), OperationId>,
    /// Mutable map of all body BReps.
    pub body_breps: &'a mut HashMap<BodyId, BRep>,
    /// Parameter table for value resolution.
    pub param_table: &'a ParameterTable,
    /// Geometric tolerance.
    pub tolerance: &'a Tolerance,
    /// Counter for allocating new OperationIds.
    pub next_op_id: &'a mut u64,
    /// OperationId of the pattern feature itself.
    pub pattern_op_id: OperationId,
    /// Body that the seed features built their tool body on.
    pub seed_body: BodyId,
    /// All feature nodes (for host body detection).
    pub features: &'a slotmap::SlotMap<FeatureId, FeatureNode>,
    /// FeatureIds of the seed features.
    pub seed_feature_ids: &'a [FeatureId],
}

/// Execute a pattern feature, producing geometry for all instances.
///
/// This is the top-level function called from `tree_replay.rs`.  It computes
/// transforms, replays or copy-transforms each instance, booleans the result
/// into the host body, and returns a [`PatternResult`] summarising outcomes.
///
/// # Errors
///
/// Returns [`ReplayError`] if transform computation fails or all instances
/// fail fatally.
#[allow(clippy::too_many_lines)]
pub(crate) fn execute_pattern(ctx: &mut PatternContext<'_>) -> Result<PatternResult, ReplayError> {
    // 1. Determine host body BRep for transform computation.
    let host_brep =
        ctx.body_breps
            .get(&ctx.seed_body)
            .ok_or_else(|| ReplayError::FeatureFailed {
                label: "pattern".into(),
                error: format!("seed body {:?} has no BRep", ctx.seed_body),
            })?;

    // 2. Compute per-instance transforms.
    let transforms = pattern_transform::compute_pattern_transforms(
        ctx.pattern_feature,
        host_brep,
        ctx.param_table,
    )?;

    // 3. Extract pattern flags.
    let (geometry_pattern, instance_overrides) = extract_pattern_flags(ctx.pattern_feature);

    // 4. Determine boolean operation from seed features' mode.
    let bool_op = seed_bool_operation(ctx.seed_features);

    // 5. Detect host body for the boolean merge.
    let host_body = if !ctx.seed_feature_ids.is_empty() {
        auto_merge::detect_host_body(ctx.seed_feature_ids[0], ctx.features, ctx.body_breps)
    } else {
        None
    };

    // 6. Snapshot the tool body BEFORE instance 0 consumes it.
    let tool_snapshot = ctx.body_breps.get(&ctx.seed_body).cloned();

    // 7. Allocate/reuse OperationIds for all instances.
    let num_seeds = ctx.seed_op_ids.len();
    for inst_idx in 0..transforms.len() as u32 {
        for seed_idx in 0..num_seeds as u32 {
            let key = (inst_idx, seed_idx);
            ctx.instance_op_ids.entry(key).or_insert_with(|| {
                let id = OperationId(*ctx.next_op_id);
                *ctx.next_op_id += 1;
                id
            });
        }
    }

    // 8. Execute instances.
    let mut outcomes: Vec<(usize, InstanceOutcome)> = Vec::with_capacity(transforms.len());

    for (inst_idx, transform) in transforms.iter().enumerate() {
        // Check suppression.
        if instance_overrides
            .get(&inst_idx)
            .is_some_and(|ov| ov.suppressed)
        {
            outcomes.push((inst_idx, InstanceOutcome::Skipped));
            continue;
        }

        if inst_idx == 0 {
            // Instance 0: the seed features already built the tool body.
            // Boolean it into the host body.
            match boolean_into_host(
                ctx.seed_body,
                host_body,
                bool_op,
                ctx.body_breps,
                ctx.pattern_op_id,
                ctx.tolerance,
            ) {
                Ok(()) => outcomes.push((inst_idx, InstanceOutcome::Replayed)),
                Err(e) => outcomes.push((
                    inst_idx,
                    InstanceOutcome::Failed {
                        error: e.to_string(),
                    },
                )),
            }
            continue;
        }

        // Instances 1..N: try feature-replay, then geometry fallback.
        let host_target = host_body.unwrap_or(ctx.seed_body);
        let host_backup = ctx.body_breps.get(&host_target).cloned();
        let mut success = false;

        if !geometry_pattern {
            match try_feature_replay(
                inst_idx,
                ctx.seed_features,
                ctx.seed_op_ids,
                ctx.instance_op_ids,
                ctx.body_breps,
                ctx.param_table,
                ctx.tolerance,
                host_target,
                bool_op,
                ctx.pattern_op_id,
            ) {
                Ok(()) => {
                    outcomes.push((inst_idx, InstanceOutcome::Replayed));
                    success = true;
                }
                Err(_) => {
                    if let Some(backup) = &host_backup {
                        ctx.body_breps.insert(host_target, backup.clone());
                    }
                }
            }
        }

        if !success {
            match try_geometry_fallback(
                &tool_snapshot,
                transform,
                ctx.body_breps,
                host_target,
                bool_op,
                ctx.pattern_op_id,
                ctx.tolerance,
            ) {
                Ok(()) => {
                    outcomes.push((inst_idx, InstanceOutcome::GeometryFallback));
                }
                Err(e) => {
                    if let Some(backup) = &host_backup {
                        ctx.body_breps.insert(host_target, backup.clone());
                    }
                    outcomes.push((
                        inst_idx,
                        InstanceOutcome::Failed {
                            error: e.to_string(),
                        },
                    ));
                }
            }
        }
    }

    Ok(build_pattern_result(outcomes))
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Extract `geometry_pattern` flag and `instance_overrides` from a pattern
/// feature variant.
fn extract_pattern_flags(
    feature: &Feature,
) -> (
    bool,
    &std::collections::HashMap<usize, crate::pattern::InstanceOverride>,
) {
    static EMPTY: std::sync::LazyLock<
        std::collections::HashMap<usize, crate::pattern::InstanceOverride>,
    > = std::sync::LazyLock::new(std::collections::HashMap::new);

    match feature {
        Feature::LinearPattern {
            geometry_pattern,
            instance_overrides,
            ..
        }
        | Feature::CircularPattern {
            geometry_pattern,
            instance_overrides,
            ..
        }
        | Feature::CurvePattern {
            geometry_pattern,
            instance_overrides,
            ..
        } => (*geometry_pattern, instance_overrides),
        Feature::Mirror { .. } => (true, &EMPTY),
        _ => (false, &EMPTY),
    }
}

/// Determine the boolean operation from the seed features' modes.
///
/// Returns `Union` if any seed is `Boss`, `Subtraction` if any is `Cut`,
/// and `Union` as default.
fn seed_bool_operation(seed_features: &[Feature]) -> popkern_bool::BoolOperation {
    for feat in seed_features {
        match auto_merge::feature_mode(feat) {
            OperationMode::Cut => return popkern_bool::BoolOperation::Subtraction,
            OperationMode::Boss => return popkern_bool::BoolOperation::Union,
            OperationMode::NewBody => {}
        }
    }
    popkern_bool::BoolOperation::Union
}

/// Boolean-merge the tool body into the host body.
fn boolean_into_host(
    tool_body: BodyId,
    host_body: Option<BodyId>,
    bool_op: popkern_bool::BoolOperation,
    body_breps: &mut HashMap<BodyId, BRep>,
    op_id: OperationId,
    tolerance: &Tolerance,
) -> Result<(), ReplayError> {
    let host = match host_body {
        Some(h) if h != tool_body => h,
        _ => return Ok(()),
    };
    let mut ctx = TrackingContext::new(op_id);
    auto_merge::execute_auto_merge(
        match bool_op {
            popkern_bool::BoolOperation::Union => OperationMode::Boss,
            popkern_bool::BoolOperation::Subtraction => OperationMode::Cut,
            popkern_bool::BoolOperation::Intersection => OperationMode::Boss,
        },
        tool_body,
        host,
        body_breps,
        &mut ctx,
        tolerance,
    )
}

/// Try the feature-replay path for a single instance.
#[allow(clippy::too_many_arguments)]
fn try_feature_replay(
    inst_idx: usize,
    seed_features: &[Feature],
    seed_op_ids: &[OperationId],
    instance_op_ids: &HashMap<(u32, u32), OperationId>,
    body_breps: &mut HashMap<BodyId, BRep>,
    param_table: &ParameterTable,
    tolerance: &Tolerance,
    host_body: BodyId,
    bool_op: popkern_bool::BoolOperation,
    pattern_op_id: OperationId,
) -> Result<(), ReplayError> {
    let mut op_remap = HashMap::new();
    for (seed_idx, seed_op) in seed_op_ids.iter().enumerate() {
        let key = (inst_idx as u32, seed_idx as u32);
        if let Some(&inst_op) = instance_op_ids.get(&key) {
            op_remap.insert(*seed_op, inst_op);
        }
    }

    let remapped: Vec<(Feature, OperationId)> =
        seed_features
            .iter()
            .enumerate()
            .map(|(seed_idx, feat)| {
                let cloned = pattern_replay::clone_seed_feature(feat, &op_remap);
                let key = (inst_idx as u32, seed_idx as u32);
                let op = instance_op_ids.get(&key).copied().ok_or_else(|| {
                    ReplayError::FeatureFailed {
                        label: String::new(),
                        error: format!(
                            "missing instance OperationId for instance {inst_idx}, seed {seed_idx}"
                        ),
                    }
                })?;
                Ok((cloned, op))
            })
            .collect::<Result<Vec<_>, ReplayError>>()?;

    let tool_brep =
        pattern_replay::replay_seed_chain(&remapped, tolerance, param_table, body_breps)?;

    let host_brep = body_breps
        .remove(&host_body)
        .ok_or_else(|| ReplayError::FeatureFailed {
            label: "pattern".into(),
            error: format!("host body {host_body:?} has no BRep"),
        })?;

    let mut ctx = TrackingContext::new(pattern_op_id);
    let result = popkern_bool::boolean(bool_op, host_brep, tool_brep, *tolerance, Some(&mut ctx))
        .map_err(|e| ReplayError::FeatureFailed {
        label: "pattern".into(),
        error: e.to_string(),
    })?;

    body_breps.insert(host_body, result);
    Ok(())
}

/// Try the geometry fallback path for a single instance.
///
/// Clones the tool body snapshot, transforms it, and booleans it into the
/// host body.
fn try_geometry_fallback(
    tool_snapshot: &Option<BRep>,
    transform: &popkern_core::geometry::transform::Transform3,
    body_breps: &mut HashMap<BodyId, BRep>,
    host_body: BodyId,
    bool_op: popkern_bool::BoolOperation,
    pattern_op_id: OperationId,
    tolerance: &Tolerance,
) -> Result<(), ReplayError> {
    let snapshot = tool_snapshot
        .as_ref()
        .ok_or_else(|| ReplayError::FeatureFailed {
            label: "pattern".into(),
            error: "no tool body snapshot for geometry fallback".into(),
        })?;

    let tool_brep = snapshot
        .transform(transform)
        .map_err(|e| ReplayError::FeatureFailed {
            label: "pattern".into(),
            error: format!("transform failed: {e}"),
        })?;

    let host_brep = body_breps
        .remove(&host_body)
        .ok_or_else(|| ReplayError::FeatureFailed {
            label: "pattern".into(),
            error: format!("host body {host_body:?} has no BRep for geometry fallback"),
        })?;

    let mut ctx = TrackingContext::new(pattern_op_id);
    let result = popkern_bool::boolean(bool_op, host_brep, tool_brep, *tolerance, Some(&mut ctx))
        .map_err(|e| ReplayError::FeatureFailed {
        label: "pattern".into(),
        error: e.to_string(),
    })?;

    body_breps.insert(host_body, result);
    Ok(())
}

/// Build a [`PatternResult`] from per-instance outcomes.
fn build_pattern_result(outcomes: Vec<(usize, InstanceOutcome)>) -> PatternResult {
    let mut failed: Vec<(usize, String)> = Vec::new();
    let mut succeeded: Vec<(usize, InstanceOutcome)> = Vec::new();

    for (idx, outcome) in outcomes {
        match &outcome {
            InstanceOutcome::Failed { error } => failed.push((idx, error.clone())),
            _ => succeeded.push((idx, outcome)),
        }
    }

    if failed.is_empty() {
        PatternResult::AllSucceeded {
            count: succeeded.len(),
        }
    } else if succeeded.is_empty() {
        PatternResult::AllFailed { errors: failed }
    } else {
        PatternResult::PartialSuccess { succeeded, failed }
    }
}
