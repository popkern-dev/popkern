//! Template cells (design §4): every curated tangency corpus case under
//! random placements and self-turns beyond the corpus's fixed set,
//! judged against the case's own closed form.

use popkern_bool::{BoolError, BoolOptions, NonManifoldPolicy, boolean_with_options};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[path = "../tangency_corpus_support/mod.rs"]
mod support;

use super::rng::{Rng, stream_id};
use support::{Case, Expect, apply_turn, cases};

fn judge(
    outcome: &Result<BRep, BoolError>,
    expect: &Expect,
    place: &Transform3,
    scale: f64,
    tol: &Tolerance,
) -> Result<(), String> {
    let near = |p: &Point3, q: &Point3| p.distance_to(q) <= 1e-5 * scale.max(1.0);
    match (outcome, expect) {
        (Ok(body), Expect::Body { volume, area }) => {
            if body.faces.is_empty() {
                return Err("served an empty body".into());
            }
            if let Some(first) = body.validate_all().first() {
                return Err(format!("validate_all: {first}"));
            }
            let v = popkern_measure::body_volume(body, tol).map_err(|e| format!("volume: {e}"))?;
            let a =
                popkern_measure::body_surface_area(body, tol).map_err(|e| format!("area: {e}"))?;
            let (wv, wa) = (volume * scale.powi(3), area * scale * scale);
            let bv = super::judge::MEASURE_REL * wv + tol.linear() * wa;
            let ba = super::judge::MEASURE_REL * wa + 4.0 * tol.linear() * wa.sqrt();
            if (v - wv).abs() > bv || (a - wa).abs() > ba {
                return Err(format!("volume {v} want {wv}, area {a} want {wa}"));
            }
            Ok(())
        }
        (Ok(body), Expect::Empty) if body.faces.is_empty() => Ok(()),
        (Err(BoolError::NonManifoldResult { contacts }), Expect::NonManifold(at)) => {
            let want: Vec<Point3> = at.iter().map(|p| place.apply_point(p)).collect();
            let missing = want
                .iter()
                .filter(|w| !contacts.iter().any(|c| near(&c.location, w)))
                .count();
            let extra = contacts
                .iter()
                .filter(|c| !want.iter().any(|w| near(&c.location, w)))
                .count();
            if missing == 0 && extra == 0 {
                Ok(())
            } else {
                Err(format!("contacts: {missing} missing, {extra} unexpected"))
            }
        }
        (Ok(body), e) => Err(format!("served F{}, expected {e:?}", body.faces.len())),
        (Err(err), e) => Err(format!("declined `{err}`, expected {e:?}")),
    }
}

/// Run every template `random` times under random placements; the
/// failures as ratchet lines (`template-<case>[policy]/template/-`).
pub fn run_templates(seed: u64, random: usize) -> (Vec<(String, String)>, usize, usize) {
    let tol = Tolerance::default();
    let mut failures = Vec::new();
    let (mut runs, mut passed) = (0, 0);
    for case in cases() {
        let rows: Vec<(NonManifoldPolicy, Expect, &str)> =
            std::iter::once((NonManifoldPolicy::Decline, case.decline.clone(), "Decline"))
                .chain(
                    case.general
                        .clone()
                        .map(|g| (NonManifoldPolicy::General, g, "General")),
                )
                .collect();
        for (policy, expect, label) in rows {
            for k in 0..random as u64 {
                let mut rng = Rng::new(
                    seed,
                    stream_id(&format!("template-{}-{label}", case.name)),
                    k,
                );
                let (place, scale) = placement(&mut rng);
                let built = build(&case, &mut rng);
                let (a, b) = (
                    built.0.transform(&place).expect("place"),
                    built.1.transform(&place).expect("place"),
                );
                runs += 1;
                let outcome = boolean_with_options(
                    case.op,
                    a,
                    b,
                    tol,
                    None,
                    &BoolOptions {
                        non_manifold: policy,
                    },
                );
                match judge(&outcome, &expect, &place, scale, &tol) {
                    Ok(()) => passed += 1,
                    Err(e) => failures.push((
                        format!("template-{}[{label}]/template/-", case.name),
                        format!("{e} [seed {seed} k {k}]"),
                    )),
                }
            }
        }
    }
    (failures, runs, passed)
}

fn placement(rng: &mut Rng) -> (Transform3, f64) {
    let scale = [1.0, 1.0, 1e-3, 1e3][rng.below(4)];
    let axis = rng.direction();
    let turn = Transform3::rotation_around_axis(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(axis[0], axis[1], axis[2]),
        rng.range(0.0, std::f64::consts::TAU),
    )
    .expect("axis");
    let shift = Transform3::translation(
        rng.range(-50.0, 50.0),
        rng.range(-50.0, 50.0),
        rng.range(-50.0, 50.0),
    );
    (
        Transform3::scale_uniform(scale).compose(&shift.compose(&turn)),
        scale,
    )
}

fn build(case: &Case, rng: &mut Rng) -> (BRep, BRep) {
    let pair = (case.build)(Tolerance::default());
    if case.turns.is_empty() || rng.unit() < 0.3 {
        return pair;
    }
    // The case's own self-turns, as given: some are discrete symmetries.
    let turn = case.turns[rng.below(case.turns.len())];
    apply_turn(pair, &turn)
}
