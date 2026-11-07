//! The non-manifold corpus under `NonManifoldPolicy::SplitToManifold`
//! (PK-964 Stage 6): every row that declines by default serves as
//! manifold pieces meeting at recorded contacts — `validate_all` clean
//! with the contacts recorded, the exact volume of the touching bodies,
//! under every placement. Rows that serve by default serve the same
//! under the split policy.

mod nonmanifold_support;

use nonmanifold_support::{Expect, build, cases, placements, two_cut_nm_edge};
use popkern_bool::complex_shadow::ComplexRescueOutcome;
use popkern_bool::{
    BoolError, BoolOperation, BoolOptions, NonManifoldPolicy, boolean_with_options,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

const SPLIT: BoolOptions = BoolOptions {
    non_manifold: NonManifoldPolicy::SplitToManifold,
};

fn volume(b: &BRep) -> f64 {
    popkern_measure::body_volume(b, &Tolerance::default()).expect("volume")
}

/// A split serve: valid with its contacts, contacts present, and the
/// volume of `want`.
fn judge_split(outcome: &Result<BRep, BoolError>, want: f64) -> Result<(), String> {
    let body = outcome.as_ref().map_err(|e| format!("declined: {e}"))?;
    if body.topology_contacts.is_empty() {
        return Err("no contact recorded".into());
    }
    let issues = body.validate_all();
    if let Some(first) = issues.first() {
        return Err(format!("validate_all: {first}"));
    }
    let v = volume(body);
    if (v - want).abs() > 1e-5 * want {
        return Err(format!("volume {v}, expected {want}"));
    }
    Ok(())
}

#[test]
fn nonmanifold_rows_serve_split_with_recorded_contacts() {
    let mut rows = Vec::new();
    let mut failures = 0;
    for (pname, place, scale3) in placements() {
        let tol = Tolerance::default();
        for case in cases() {
            let a = build(case.a, &place, tol);
            let b = build(case.b, &place, tol);
            let (va, vb) = (volume(&a), volume(&b));
            let outcome = boolean_with_options(case.op, a, b, tol, None, &SPLIT);
            let verdict = match case.expect {
                // Touching operands: the union is both, side by side.
                Expect::NonManifold => judge_split(&outcome, va + vb),
                Expect::Valid { volume: v, .. } => match &outcome {
                    Ok(body)
                        if body.validate_all().is_empty()
                            && (volume(body) - v * scale3).abs() <= 1e-5 * v * scale3 =>
                    {
                        Ok(())
                    }
                    Ok(body) => Err(format!("served volume {}", volume(body))),
                    Err(e) => Err(format!("declined: {e}")),
                },
                Expect::Empty => match &outcome {
                    Ok(body) if body.faces.is_empty() => Ok(()),
                    other => Err(format!("{:?}", other.as_ref().map(|b| b.faces.len()))),
                },
            };
            if verdict.is_err() {
                failures += 1;
            }
            rows.push(format!("{pname:>8} {:<28} {verdict:?}", case.name));
        }
        // Two overshooting cuts: one body whose regions meet at an edge.
        let first = boolean_with_options(
            BoolOperation::Subtraction,
            build(
                nonmanifold_support::Shape::Cube((0.0, 0.0, 0.0), 2.0),
                &place,
                tol,
            ),
            build(
                nonmanifold_support::Shape::Cube((1.0, -0.5, -0.5), 1.5),
                &place,
                tol,
            ),
            tol,
            None,
            &SPLIT,
        );
        let verdict = first.map_err(|e| format!("first cut: {e}")).and_then(|s1| {
            let cut = boolean_with_options(
                BoolOperation::Subtraction,
                s1,
                build(
                    nonmanifold_support::Shape::Cube((-0.5, 1.0, -0.5), 1.5),
                    &place,
                    tol,
                ),
                tol,
                None,
                &SPLIT,
            );
            judge_split(&cut, 6.0 * scale3)
        });
        let _ = two_cut_nm_edge;
        if verdict.is_err() {
            failures += 1;
        }
        rows.push(format!("{pname:>8} {:<28} {verdict:?}", "two_cut_nm_edge"));
    }
    assert_eq!(failures, 0, "split corpus:\n{}", rows.join("\n"));
    // Every split serves from the complex's own pairing (PK-987), never
    // from graph-merge's coincident copies: no split was rejected, and
    // each non-manifold verdict was followed by a served split.
    let events = popkern_bool::complex_shadow::take_complex_rescue_events();
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

/// Negative control: dropping the recorded contacts from a split serve
/// makes `validate_all` reject the coincident copies again.
#[test]
fn an_unrecorded_split_fails_validate_all() {
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
    let mut body =
        boolean_with_options(BoolOperation::Union, a, b, tol, None, &SPLIT).expect("splits");
    assert!(body.validate_all().is_empty());
    body.topology_contacts.clear();
    assert!(!body.validate_all().is_empty());
}

/// Rows whose general body still records contacts: a contact line
/// across a face's interior has no topology until it is imprinted
/// (PK-1117), and vertex contacts of some placements the vertex glue
/// misses there too. Each still serves a valid general body.
const STILL_RECORDED: &[(&str, &str)] = &[
    ("identity", "edge_on_face_union"),
    ("rotated", "edge_touch_union"),
    ("rotated", "line_contact_union"),
    ("rotated", "edge_on_face_union"),
    ("mm", "edge_on_face_union"),
    ("far", "line_contact_union"),
    ("far", "edge_on_face_union"),
];

/// The corpus under `NonManifoldPolicy::General` (PK-1115, PK-1116):
/// every row that touches itself serves one general body that validates
/// with the volume of the touching bodies and nothing recorded, its
/// contacts all held by topology (rings, shared vertices, vertex loops
/// on faces, edges split where a vertex touches them), except the rows
/// [`STILL_RECORDED`] names; a listed row that records nothing fails, so
/// the list shrinks as PK-1117 lands.
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
            let listed = STILL_RECORDED.contains(&(pname, case.name));
            let verdict = match boolean_with_options(case.op, a, b, tol, None, &GENERAL) {
                Ok(body) => {
                    let issues = body.validate_all();
                    let recorded = body.topology_contacts.len();
                    if !issues.is_empty() {
                        Err(format!("invalid: {}", issues[0]))
                    } else if (volume(&body) - want).abs() > 1e-5 * want {
                        Err(format!("volume {} vs {want}", volume(&body)))
                    } else if listed && recorded == 0 {
                        Err("records nothing now: drop it from STILL_RECORDED".into())
                    } else if !listed && recorded > 0 {
                        Err(format!("{recorded} recorded"))
                    } else {
                        Ok(())
                    }
                }
                Err(e) => Err(format!("declined: {e}")),
            };
            if verdict.is_err() {
                failures += 1;
            }
            rows.push(format!("{pname:>8} {:<28} {verdict:?}", case.name));
        }
    }
    assert_eq!(failures, 0, "general corpus:\n{}", rows.join("\n"));
}
