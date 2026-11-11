//! Slots closed by N-ary booleans are reported as a two-operand union
//! reports them (PK-1162, a follow-up of PK-1131).
//!
//! Three blocks in a row, [0, 10]³, then each next one a gap g further
//! along x, at identity and turned. A selection keeping both sides of a
//! gap at or under the linear tolerance fuses them across it and
//! records a `SlotClosed` per gap; one keeping a single side leaves a
//! face where the gap was and records nothing; a gap over the tolerance
//! is never closed. Before, only `boolean()`'s union reported a slot:
//! `nary_boolean` and `Partition::select` closed the same gaps
//! silently.

use popkern_bool::BoolOptions;
use popkern_bool::partition::{SetExpr, nary_boolean, nary_boolean_where, partition};
use popkern_bool::sub_tolerance::{SubToleranceKind, take_sub_tolerance_events};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn block(x0: f64, tol: Tolerance, place: &Transform3) -> BRep {
    let rect = |z: f64| {
        vec![
            Point3::new(x0, 0.0, z),
            Point3::new(x0 + 10.0, 0.0, z),
            Point3::new(x0 + 10.0, 10.0, z),
            Point3::new(x0, 10.0, z),
        ]
    };
    BRep::polygon_frustum(&rect(0.0), &rect(10.0), tol)
        .and_then(|b| b.transform(place))
        .expect("block")
}

/// The widths of the slots recorded since the last drain.
fn slots() -> Vec<f64> {
    take_sub_tolerance_events()
        .into_iter()
        .filter(|e| e.kind == SubToleranceKind::SlotClosed)
        .map(|e| e.measured_width)
        .collect()
}

#[test]
fn an_nary_selection_reports_each_gap_it_closes() {
    let tol = Tolerance::default();
    let turned =
        Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.7, 0.3), 0.6)
            .expect("axis");
    let options = BoolOptions::default();
    let mut problems = Vec::new();
    for (place_name, place) in [("identity", Transform3::identity()), ("turned", turned)] {
        for gap in [0.0, 5e-7, 1e-6, 1e-5] {
            let blocks: Vec<BRep> = (0..3)
                .map(|k| block(k as f64 * (10.0 + gap), tol, &place))
                .collect();
            let refs: Vec<&BRep> = blocks.iter().collect();
            let closes = gap > 0.0 && tol.within_linear(gap);
            // Each expression with the gaps it keeps both sides of:
            // all three keep both (A|B, B|C), A or B one, A alone or
            // all but B none.
            let cases = [
                ("any", SetExpr::Any, 2),
                ("a or b", SetExpr::body(0).or(SetExpr::body(1)), 1),
                ("a", SetExpr::body(0), 0),
                ("any minus b", SetExpr::Any.minus(SetExpr::body(1)), 0),
            ];
            for (expr_name, expr, kept) in cases {
                let name = format!("{place_name} gap {gap:e} {expr_name}");
                let _ = slots();
                match nary_boolean(&refs, &expr, tol, &options) {
                    Ok(body) => {
                        if let Some(first) = body.validate_all().first() {
                            problems.push(format!("{name}: {first}"));
                        }
                    }
                    Err(e) => problems.push(format!("{name}: declined {e}")),
                }
                let want = if closes { kept } else { 0 };
                let got = slots();
                if got.len() != want || got.iter().any(|w| (w - gap).abs() > 1e-6 * gap) {
                    problems.push(format!("{name}: slots {got:?}, want {want} of {gap:e}"));
                }
            }
            // The predicate form and a partition selected directly.
            let _ = slots();
            let any = nary_boolean_where(&refs, |m| !m.is_empty(), tol, &options);
            let want = if closes { 2 } else { 0 };
            let got = slots();
            if any.is_err() || got.len() != want {
                problems.push(format!(
                    "{place_name} gap {gap:e} where: {:?}, slots {got:?}",
                    any.err()
                ));
            }
            let selected = partition(&refs, tol).and_then(|p| p.select(&SetExpr::Any, &options));
            let got = slots();
            if selected.is_err() || got.len() != want {
                problems.push(format!(
                    "{place_name} gap {gap:e} select: {:?}, slots {got:?}",
                    selected.err()
                ));
            }
        }
    }
    assert!(problems.is_empty(), "{problems:#?}");
}
