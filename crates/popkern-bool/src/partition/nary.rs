//! A boolean of N bodies by set expression, with the recorded pairwise
//! fallback (design 7.3, I.3).
//!
//! The partition answers in one pass. When it declines (a stage of its
//! complex cannot be built), the expression is a boolean formula and
//! two-operand booleans evaluate it: the structural nodes fold directly
//! (a body, ∧, ∨, −, all of them, any of them); the counting nodes
//! (exactly one, odd, at least k) by their disjunctive form over the
//! membership subsets, for at most [`FALLBACK_COUNTING_LIMIT`] bodies.
//! A fallback is recorded (`PartitionFallback` shadow line), never
//! silent; one that itself fails returns the partition's decline.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::{BodySet, SetExpr, partition};
use crate::error::{BoolError, BoolOperation};
use crate::options::BoolOptions;

/// The most bodies a counting expression's fallback expands over (its
/// disjunctive form has up to 2ⁿ terms).
pub const FALLBACK_COUNTING_LIMIT: usize = 6;

/// The body `expr` keeps of `bodies`: the partition's selection, or the
/// pairwise fallback when the partition declines.
///
/// # Errors
///
/// The partition's errors that are answers rather than declines
/// (non-manifold result, unsupported body kind, too many bodies, a bad
/// expression, cancellation) as they come; a decline when the fallback
/// fails too or cannot expand the expression.
pub fn nary_boolean(
    bodies: &[&BRep],
    expr: &SetExpr,
    tolerance: Tolerance,
    options: &BoolOptions,
) -> Result<BRep, BoolError> {
    let pieces = nary_pieces(bodies, expr, tolerance, options);
    crate::split_serve::serve_one(pieces, options.non_manifold, &tolerance)
}

/// [`nary_boolean`] as a list of bodies: one body, except under
/// `SplitToManifold` a result whose manifold pieces touch, each piece
/// its own valid body (as [`crate::boolean_bodies`]).
///
/// # Errors
///
/// As [`nary_boolean`], never [`BoolError::SplitIntoBodies`].
pub fn nary_boolean_bodies(
    bodies: &[&BRep],
    expr: &SetExpr,
    tolerance: Tolerance,
    options: &BoolOptions,
) -> Result<Vec<BRep>, BoolError> {
    let pieces = nary_pieces(bodies, expr, tolerance, options);
    crate::split_serve::serve_bodies(pieces, options.non_manifold, &tolerance)
}

/// [`nary_boolean`] as the pipeline answers it: a non-manifold result
/// as manifold pieces in one body under the split and general policies.
fn nary_pieces(
    bodies: &[&BRep],
    expr: &SetExpr,
    tolerance: Tolerance,
    options: &BoolOptions,
) -> Result<BRep, BoolError> {
    let pieces = general_operands(bodies, options)?;
    let bodies: Vec<&BRep> = match &pieces {
        Some(pieces) => pieces.iter().collect(),
        None => bodies.to_vec(),
    };
    let served = partition(&bodies, tolerance).and_then(|p| p.select_pieces(expr, options));
    or_fold(served, &bodies, expr, tolerance, options)
}

/// Under the general policy, the bodies as manifold pieces (a general
/// body taken apart); the entry point glues the result. `None`
/// otherwise.
fn general_operands(
    bodies: &[&BRep],
    options: &BoolOptions,
) -> Result<Option<Vec<BRep>>, BoolError> {
    if options.non_manifold != crate::options::NonManifoldPolicy::General {
        return Ok(None);
    }
    let pieces = bodies
        .iter()
        .enumerate()
        .map(|(i, b)| crate::general::operand_pieces((*b).clone(), i))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(Some(pieces))
}

/// The partition's answer, or on a decline the fold of `expr` by
/// two-operand booleans (the seam [`nary_boolean`] answers through).
pub(super) fn or_fold(
    served: Result<BRep, BoolError>,
    bodies: &[&BRep],
    expr: &SetExpr,
    tolerance: Tolerance,
    options: &BoolOptions,
) -> Result<BRep, BoolError> {
    let decline = match decline_of(served) {
        Ok(body) => return Ok(body),
        Err(Ok(answer)) => return Err(answer),
        Err(Err(decline)) => decline,
    };
    match fold(bodies, expr, tolerance, options) {
        Some(Ok(body)) => {
            crate::complex_shadow::append_log_line(&format!(
                "PartitionFallback {} bodies {expr:?}: {decline}",
                bodies.len()
            ));
            Ok(body)
        }
        _ => Err(decline),
    }
}

/// [`nary_boolean`] by a predicate on each region's bodies (the
/// code-only form); its fallback is the disjunctive form, for at most
/// [`FALLBACK_COUNTING_LIMIT`] bodies.
///
/// # Errors
///
/// As [`nary_boolean`].
pub fn nary_boolean_where(
    bodies: &[&BRep],
    keep: impl Fn(&BodySet) -> bool,
    tolerance: Tolerance,
    options: &BoolOptions,
) -> Result<BRep, BoolError> {
    let pieces = nary_where_pieces(bodies, keep, tolerance, options);
    crate::split_serve::serve_one(pieces, options.non_manifold, &tolerance)
}

/// [`nary_boolean_where`] as a list of bodies (as
/// [`nary_boolean_bodies`]).
///
/// # Errors
///
/// As [`nary_boolean_where`], never [`BoolError::SplitIntoBodies`].
pub fn nary_boolean_where_bodies(
    bodies: &[&BRep],
    keep: impl Fn(&BodySet) -> bool,
    tolerance: Tolerance,
    options: &BoolOptions,
) -> Result<Vec<BRep>, BoolError> {
    let pieces = nary_where_pieces(bodies, keep, tolerance, options);
    crate::split_serve::serve_bodies(pieces, options.non_manifold, &tolerance)
}

/// [`nary_boolean_where`] as the pipeline answers it (see
/// [`nary_pieces`]).
fn nary_where_pieces(
    bodies: &[&BRep],
    keep: impl Fn(&BodySet) -> bool,
    tolerance: Tolerance,
    options: &BoolOptions,
) -> Result<BRep, BoolError> {
    let pieces = general_operands(bodies, options)?;
    let bodies: Vec<&BRep> = match &pieces {
        Some(pieces) => pieces.iter().collect(),
        None => bodies.to_vec(),
    };
    nary_boolean_where_under(&bodies, keep, tolerance, options)
}

fn nary_boolean_where_under(
    bodies: &[&BRep],
    keep: impl Fn(&BodySet) -> bool,
    tolerance: Tolerance,
    options: &BoolOptions,
) -> Result<BRep, BoolError> {
    let decline = match decline_of(
        partition(bodies, tolerance).and_then(|p| p.select_where_pieces(&keep, options)),
    ) {
        Ok(body) => return Ok(body),
        Err(Ok(answer)) => return Err(answer),
        Err(Err(decline)) => decline,
    };
    if bodies.len() > FALLBACK_COUNTING_LIMIT {
        return Err(decline);
    }
    match disjunctive(bodies, &keep, tolerance, options) {
        Ok(body) => {
            crate::complex_shadow::append_log_line(&format!(
                "PartitionFallback {} bodies (predicate): {decline}",
                bodies.len()
            ));
            Ok(body)
        }
        Err(_) => Err(decline),
    }
}

/// The partition's body; else its error, split into an answer (a
/// non-manifold result, an unsupported body kind, too many bodies, a
/// bad expression, cancellation: `Ok`) or a decline the fallback may
/// rescue (`Err`).
fn decline_of(served: Result<BRep, BoolError>) -> Result<BRep, Result<BoolError, BoolError>> {
    served.map_err(|e| match e {
        BoolError::NonManifoldResult { .. }
        | BoolError::UnsupportedBodyKind { .. }
        | BoolError::TooManyBodies { .. }
        | BoolError::DegenerateInput(_)
        | BoolError::Cancelled => Ok(e),
        e => Err(e),
    })
}

/// `a op b` with an empty body as the identity it is.
fn combine(
    op: BoolOperation,
    a: BRep,
    b: BRep,
    tolerance: Tolerance,
    options: &BoolOptions,
) -> Result<BRep, BoolError> {
    let empty = |x: &BRep| x.faces.is_empty();
    match (op, empty(&a), empty(&b)) {
        (BoolOperation::Union, true, _) => Ok(b),
        (BoolOperation::Union, _, true) => Ok(a),
        (BoolOperation::Intersection, true, _) | (BoolOperation::Subtraction, true, _) => Ok(a),
        (BoolOperation::Intersection, _, true) => Ok(b),
        (BoolOperation::Subtraction, _, true) => Ok(a),
        _ => crate::bool_api_core::boolean_under(op, a, b, tolerance, None, options.non_manifold),
    }
}

/// `expr` by two-operand booleans; `None` when it cannot be expanded.
fn fold(
    bodies: &[&BRep],
    expr: &SetExpr,
    tolerance: Tolerance,
    options: &BoolOptions,
) -> Option<Result<BRep, BoolError>> {
    let n = bodies.len();
    let all = |op| {
        let mut acc = BRep::new(tolerance);
        for (k, b) in bodies.iter().enumerate() {
            acc = if k == 0 {
                (*b).clone()
            } else {
                combine(op, acc, (*b).clone(), tolerance, options)?
            };
        }
        Ok(acc)
    };
    let pair = |op, x: &SetExpr, y: &SetExpr| -> Option<Result<BRep, BoolError>> {
        let a = fold(bodies, x, tolerance, options)?;
        let b = fold(bodies, y, tolerance, options)?;
        Some(a.and_then(|a| b.and_then(|b| combine(op, a, b, tolerance, options))))
    };
    match expr {
        SetExpr::Body(i) => bodies.get(*i).map(|b| Ok((*b).clone())),
        SetExpr::Any => Some(all(BoolOperation::Union)),
        SetExpr::All => Some(all(BoolOperation::Intersection)),
        SetExpr::And(x, y) => pair(BoolOperation::Intersection, x, y),
        SetExpr::Or(x, y) => pair(BoolOperation::Union, x, y),
        SetExpr::Minus(x, y) => pair(BoolOperation::Subtraction, x, y),
        counting => {
            if n > FALLBACK_COUNTING_LIMIT {
                return None;
            }
            Some(disjunctive(
                bodies,
                |m| counting.holds(m, n),
                tolerance,
                options,
            ))
        }
    }
}

/// ∪ over the membership subsets `keep` accepts of (∩ the subset − ∪
/// the rest).
fn disjunctive(
    bodies: &[&BRep],
    keep: impl Fn(&BodySet) -> bool,
    tolerance: Tolerance,
    options: &BoolOptions,
) -> Result<BRep, BoolError> {
    let n = bodies.len();
    let mut acc = BRep::new(tolerance);
    for mask in 1u64..(1u64 << n) {
        let inside: Vec<usize> = (0..n).filter(|&i| mask & (1 << i) != 0).collect();
        if !keep(&BodySet::of(&inside)) {
            continue;
        }
        let mut term = bodies[inside[0]].clone();
        for &i in &inside[1..] {
            term = combine(
                BoolOperation::Intersection,
                term,
                bodies[i].clone(),
                tolerance,
                options,
            )?;
        }
        for i in (0..n).filter(|i| mask & (1 << i) == 0) {
            term = combine(
                BoolOperation::Subtraction,
                term,
                bodies[i].clone(),
                tolerance,
                options,
            )?;
        }
        acc = combine(BoolOperation::Union, acc, term, tolerance, options)?;
    }
    Ok(acc)
}
