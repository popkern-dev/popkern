//! The non-manifold corpus under `NonManifoldPolicy::SplitToManifold`
//! (PK-964 Stage 6, PK-1115): every row that declines by default serves
//! as separate manifold bodies that touch — each `validate_all` clean,
//! together the exact volume of the touching bodies, under every
//! placement — and the single-body boolean declines it as several
//! bodies. Rows that serve by default serve the same one body under the
//! split policy.

mod nonmanifold_support;

use nonmanifold_support::{Expect, build, cases, placements, two_cut_nm_edge};
use popkern_bool::complex_shadow::ComplexRescueOutcome;
use popkern_bool::{
    BoolError, BoolOperation, BoolOptions, NonManifoldPolicy, boolean_bodies, boolean_with_options,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

const SPLIT: BoolOptions = BoolOptions {
    non_manifold: NonManifoldPolicy::SplitToManifold,
};

fn volume(b: &BRep) -> f64 {
    popkern_measure::body_volume(b, &Tolerance::default()).expect("volume")
}

/// A split serve: several bodies, each valid, together the volume of
/// `want`.
fn judge_split(outcome: &Result<Vec<BRep>, BoolError>, want: f64) -> Result<(), String> {
    let bodies = outcome.as_ref().map_err(|e| format!("declined: {e}"))?;
    if bodies.len() < 2 {
        return Err(format!("{} body served", bodies.len()));
    }
    let mut v = 0.0;
    for body in bodies {
        if let Some(first) = body.validate_all().first() {
            return Err(format!("validate_all: {first}"));
        }
        v += volume(body);
    }
    if (v - want).abs() > 1e-5 * want {
        return Err(format!("volume {v}, expected {want}"));
    }
    Ok(())
}

#[test]
fn nonmanifold_rows_serve_split_as_separate_bodies() {
    let mut rows = Vec::new();
    let mut failures = 0;
    let mut events = Vec::new();
    for (pname, place, scale3) in placements() {
        let tol = Tolerance::default();
        for case in cases() {
            let a = build(case.a, &place, tol);
            let b = build(case.b, &place, tol);
            let (va, vb) = (volume(&a), volume(&b));
            let single = boolean_with_options(case.op, a.clone(), b.clone(), tol, None, &SPLIT);
            let outcome = boolean_bodies(case.op, a, b, tol, None, &SPLIT);
            let verdict = match case.expect {
                // Touching operands: the union is both, side by side, as
                // separate bodies the single-body boolean cannot serve.
                Expect::NonManifold => match single {
                    Err(BoolError::SplitIntoBodies { .. }) => judge_split(&outcome, va + vb),
                    other => Err(format!("single body: {:?}", other.map(|b| b.faces.len()))),
                },
                Expect::Valid { volume: v, .. } => match outcome.as_deref() {
                    Ok([body])
                        if body.validate_all().is_empty()
                            && (volume(body) - v * scale3).abs() <= 1e-5 * v * scale3 =>
                    {
                        Ok(())
                    }
                    Ok(bodies) => Err(format!(
                        "served {} bodies, volume {}",
                        bodies.len(),
                        bodies.iter().map(volume).sum::<f64>()
                    )),
                    Err(e) => Err(format!("declined: {e}")),
                },
                Expect::Empty => match outcome.as_deref() {
                    Ok([body]) if body.faces.is_empty() => Ok(()),
                    other => Err(format!("{:?}", other.map(|b| b.len()))),
                },
            };
            if verdict.is_err() {
                failures += 1;
            }
            rows.push(format!("{pname:>8} {:<28} {verdict:?}", case.name));
        }
        // The corpus rows' rescue events, before the two-cut row (whose
        // split declines by design).
        events.extend(popkern_bool::complex_shadow::take_complex_rescue_events());
        // Two overshooting cuts: one solid touching itself along an edge.
        // It has no separate bodies to split into, so the split policy
        // declines it; the general policy serves it as one solid region
        // whose edge is a ring of four.
        let block = || {
            build(
                nonmanifold_support::Shape::Cube((0.0, 0.0, 0.0), 2.0),
                &place,
                tol,
            )
        };
        let first_tool = || {
            build(
                nonmanifold_support::Shape::Cube((1.0, -0.5, -0.5), 1.5),
                &place,
                tol,
            )
        };
        let second_tool = || {
            build(
                nonmanifold_support::Shape::Cube((-0.5, 1.0, -0.5), 1.5),
                &place,
                tol,
            )
        };
        let general = BoolOptions {
            non_manifold: NonManifoldPolicy::General,
        };
        let verdict = boolean_with_options(
            BoolOperation::Subtraction,
            block(),
            first_tool(),
            tol,
            None,
            &SPLIT,
        )
        .map_err(|e| format!("first cut: {e}"))
        .and_then(|s1| {
            match boolean_bodies(
                BoolOperation::Subtraction,
                s1.clone(),
                second_tool(),
                tol,
                None,
                &SPLIT,
            ) {
                Err(BoolError::NonManifoldResult { .. }) => {}
                other => {
                    return Err(format!("split: {:?}", other.map(|b| b.len())));
                }
            }
            let body = boolean_with_options(
                BoolOperation::Subtraction,
                s1,
                second_tool(),
                tol,
                None,
                &general,
            )
            .map_err(|e| format!("general: {e}"))?;
            if let Some(issue) = body.validate_all().first() {
                return Err(format!("general invalid: {issue}"));
            }
            let v = volume(&body);
            if body.solid_regions().len() != 1 || (v - 6.0 * scale3).abs() > 1e-5 * 6.0 * scale3 {
                return Err(format!(
                    "general: {} regions, volume {v}",
                    body.solid_regions().len()
                ));
            }
            Ok(())
        });
        let _ = two_cut_nm_edge;
        let _ = popkern_bool::complex_shadow::take_complex_rescue_events();
        if verdict.is_err() {
            failures += 1;
        }
        rows.push(format!("{pname:>8} {:<28} {verdict:?}", "two_cut_nm_edge"));
    }
    assert_eq!(failures, 0, "split corpus:\n{}", rows.join("\n"));
    // Every split serves from the complex's own pairing (PK-987), never
    // from graph-merge's coincident copies: no split was rejected, and
    // each non-manifold verdict was followed by a served split.
    let nm = events
        .iter()
        .filter(|e| matches!(e.outcome, ComplexRescueOutcome::NonManifold))
        .count();
    let served = events
        .iter()
        .filter(|e| matches!(e.outcome, ComplexRescueOutcome::Served))
        .count();
    assert!(nm > 0, "the corpus reaches the split");
    assert!(
        events
            .iter()
            .all(|e| !matches!(e.outcome, ComplexRescueOutcome::SplitRejected(_))),
        "{events:?}"
    );
    assert!(served >= nm, "{events:?}");
}

/// Negative control: the split bodies of two cubes on an edge, held in
/// one body, fail `validate_all` (the coincident copies).
#[test]
fn touching_split_bodies_in_one_body_fail_validate_all() {
    let tol = Tolerance::default();
    let place = popkern_core::geometry::Transform3::identity();
    let a = build(
        nonmanifold_support::Shape::Cube((0.0, 0.0, 0.0), 1.0),
        &place,
        tol,
    );
    let b = build(
        nonmanifold_support::Shape::Cube((1.0, 1.0, 0.0), 1.0),
        &place,
        tol,
    );
    let bodies = boolean_bodies(BoolOperation::Union, a, b, tol, None, &SPLIT).expect("splits");
    assert_eq!(bodies.len(), 2);
    let mut body = BRep::new(tol);
    for one in &bodies {
        assert!(one.validate_all().is_empty());
        body.append_body(one);
    }
    assert!(!body.validate_all().is_empty());
}

/// The corpus under `NonManifoldPolicy::General` (PK-1115 to PK-1117):
/// every row that touches itself serves one general body that validates
/// with the volume of the touching bodies, its contacts all held by
/// topology (rings, shared vertices, vertex loops
/// on faces, edges split where a vertex touches them, contact curves
/// imprinted across faces). The split pieces glue from their geometry
/// alone (`BRep::glue_touching`, the split bodies appended) into a
/// valid body; and the general body splits back into valid manifold
/// bodies of the same volume (`BRep::split_to_manifold`, the imprinted
/// edges kept).
#[test]
fn nonmanifold_rows_serve_general_bodies_with_nothing_recorded() {
    const GENERAL: BoolOptions = BoolOptions {
        non_manifold: NonManifoldPolicy::General,
    };
    let mut rows = Vec::new();
    let mut failures = 0;
    for (pname, place, _) in placements() {
        let tol = Tolerance::default();
        for case in cases() {
            if !matches!(case.expect, Expect::NonManifold) {
                continue;
            }
            let a = build(case.a, &place, tol);
            let b = build(case.b, &place, tol);
            let want = volume(&a) + volume(&b);
            let touch = boolean_bodies(case.op, a.clone(), b.clone(), tol, None, &SPLIT)
                .map_err(|e| format!("split declined: {e}"))
                .and_then(|bodies| {
                    let mut pieces = BRep::new(tol);
                    for one in &bodies {
                        pieces.append_body(one);
                    }
                    let done = pieces
                        .glue_touching()
                        .map_err(|e| format!("touch glue: {e}"))?;
                    match pieces.validate_all().first() {
                        Some(issue) => Err(format!("touch-glued invalid: {issue}")),
                        None => Ok(done),
                    }
                });
            let verdict = touch.and_then(|_| {
                let body = boolean_with_options(case.op, a, b, tol, None, &GENERAL)
                    .map_err(|e| format!("declined: {e}"))?;
                let issues = body.validate_all();
                if !issues.is_empty() {
                    Err(format!("invalid: {}", issues[0]))
                } else if (volume(&body) - want).abs() > 1e-5 * want {
                    Err(format!("volume {} vs {want}", volume(&body)))
                } else {
                    // And back into manifold bodies (imprints kept).
                    let back = body
                        .split_to_manifold()
                        .map_err(|e| format!("split back: {e}"))?;
                    let total: f64 = back.iter().map(volume).sum();
                    match back
                        .iter()
                        .find_map(|b| b.validate_all().into_iter().next())
                    {
                        Some(issue) => Err(format!("split back invalid: {issue}")),
                        None if (total - want).abs() > 1e-5 * want => {
                            Err(format!("split back volume {total}"))
                        }
                        None => Ok(()),
                    }
                }
            });
            if verdict.is_err() {
                failures += 1;
            }
            rows.push(format!("{pname:>8} {:<28} {verdict:?}", case.name));
        }
    }
    assert_eq!(failures, 0, "general corpus:\n{}", rows.join("\n"));
}
