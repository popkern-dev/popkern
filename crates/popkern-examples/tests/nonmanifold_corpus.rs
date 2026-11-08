//! Non-manifold and contact boolean corpus (PK-967).
//!
//! Every configuration where operands touch, lie flush, or leave a
//! non-manifold result, each under rigid motions, scales and a far
//! origin. A row's expectation is one of:
//!
//! - `Valid { shells, volume }`: a served body that passes
//!   `validate_all`, with the expected shell count and volume;
//! - `NonManifold`: the result is genuinely non-manifold (bodies
//!   touching along an edge or at a point), so the boolean declines
//!   with `BoolError::NonManifoldResult` naming the contacts (decision
//!   D1: decline by default);
//! - `Empty`: a regularized empty result (a bare body).
//!
//! Volumes are exact for the box rows and checked by identity
//! (vol A∪B + vol A∩B = vol A + vol B) everywhere. A case whose
//! outcome differs from its expectation fails with the whole table
//! printed, so a regression shows its neighbours.

mod nonmanifold_support;

use nonmanifold_support::{Expect, Shape, build, cases, placements, two_cut_nm_edge};
use popkern_bool::{BoolError, BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// What a boolean outcome is, judged against the oracles.
fn judge(outcome: &Result<BRep, BoolError>, expect: Expect, scale3: f64) -> Result<(), String> {
    match (outcome, expect) {
        (Err(BoolError::NonManifoldResult { contacts }), Expect::NonManifold)
            if !contacts.is_empty() =>
        {
            Ok(())
        }
        (Ok(body), Expect::Empty) if body.faces.is_empty() && body.back_shell_count() == 0 => {
            Ok(())
        }
        (Ok(body), Expect::Valid { shells, volume }) => {
            let issues = body.validate_all();
            if !issues.is_empty() {
                return Err(format!("validate_all: {}", issues[0]));
            }
            if body.back_shell_count() != shells {
                return Err(format!(
                    "{} shells, expected {shells}",
                    body.back_shell_count()
                ));
            }
            let v = popkern_measure::body_volume(body, &Tolerance::default())
                .map_err(|e| format!("volume: {e}"))?;
            let want = volume * scale3;
            if (v - want).abs() > 1.0e-6 * want.max(1e-300) * 10.0 {
                return Err(format!("volume {v}, expected {want}"));
            }
            Ok(())
        }
        (Ok(body), e) => Err(format!(
            "served V{} E{} F{} shells {} (validate_all {}), expected {e:?}",
            body.vertices.len(),
            body.edges.len(),
            body.faces.len(),
            body.back_shell_count(),
            body.validate_all().len()
        )),
        (Err(err), e) => Err(format!("declined: {err}; expected {e:?}")),
    }
}

#[test]
fn nonmanifold_corpus() {
    let mut rows = Vec::new();
    let mut failures = 0;
    for (pname, place, scale3) in placements() {
        let tol = Tolerance::default();
        for case in cases() {
            let a = build(case.a, &place, tol);
            let b = build(case.b, &place, tol);
            let outcome = boolean(case.op, a, b, tol, None);
            let verdict = judge(&outcome, case.expect, scale3);
            if verdict.is_err() {
                failures += 1;
            }
            rows.push(format!("{pname:>8} {:<28} {:?}", case.name, verdict));
        }
        let nm = two_cut_nm_edge(&place, tol);
        let verdict = judge(&nm, Expect::NonManifold, scale3);
        if verdict.is_err() {
            failures += 1;
        }
        rows.push(format!(
            "{pname:>8} {:<28} {:?}",
            "two_cut_nm_edge", verdict
        ));
    }
    assert_eq!(failures, 0, "non-manifold corpus:\n{}", rows.join("\n"));
}

/// The inclusion-exclusion identity holds wherever ∪ and ∩ both serve.
#[test]
fn union_and_intersection_volumes_add_up() {
    let tol = Tolerance::default();
    let id = Transform3::identity();
    let vol = |b: &BRep| popkern_measure::body_volume(b, &tol).expect("volume");
    let pairs = [
        (
            Shape::Cube((0.0, 0.0, 0.0), 1.0),
            Shape::Cube((0.5, 0.5, 0.5), 1.0),
        ),
        (
            Shape::Cube((0.0, 0.0, 0.0), 2.0),
            Shape::Cube((1.0, 0.0, 0.0), 1.0),
        ),
        (
            Shape::Cube((0.0, 0.0, 0.0), 1.0),
            Shape::Cube((1.0, 0.5, 0.5), 1.0),
        ),
    ];
    for (a, b) in pairs {
        let (va, vb) = (vol(&build(a, &id, tol)), vol(&build(b, &id, tol)));
        let u = boolean(
            BoolOperation::Union,
            build(a, &id, tol),
            build(b, &id, tol),
            tol,
            None,
        )
        .expect("union");
        let i = boolean(
            BoolOperation::Intersection,
            build(a, &id, tol),
            build(b, &id, tol),
            tol,
            None,
        )
        .expect("intersection");
        let vi = if i.faces.is_empty() { 0.0 } else { vol(&i) };
        assert!((vol(&u) + vi - va - vb).abs() < 1e-9, "{a:?} {b:?}");
    }
}
