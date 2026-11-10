//! Tangency corpus (PK-1124): operands that touch, are tangent, or meet
//! in curves that cross, under scene placements (turns, a far
//! translation, millimetre and kilometre scales), self-turns that
//! move each operand's seams and poles off and onto its contacts, and
//! for a union or common the operands swapped.
//!
//! Every run is judged against the case's closed form for its policy:
//! a served body must pass `validate_all` and measure its volume and
//! area; a non-manifold result must decline naming exactly its contact
//! points; an empty result must be empty. Placement luck is the class
//! these cases expose (the inscribed ball's intersection served EMPTY
//! with its poles on the cube's faces and right when turned), so a row
//! passes only at every placement.
//!
//! A ratchet keeps the table honest: a `Gate` row must pass everywhere,
//! and a `Target` row (the exit gate of a later stage of the tangency
//! foundations plan) must not pass everywhere unnoticed: when it does,
//! the test fails asking for it to be promoted to a gate.

mod tangency_corpus_support;

use popkern_bool::{BoolError, BoolOptions, NonManifoldPolicy, boolean_with_options};
use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use tangency_corpus_support::{Case, Expect, SelfTurn, Status, apply_turn, cases, placements};

/// Relative volume and area error a served body may carry.
const MEASURE_REL: f64 = 1e-6;

/// One run's verdict: `Ok(worst relative measure error)` or why it failed.
fn judge(
    outcome: &Result<BRep, BoolError>,
    expect: &Expect,
    place: &popkern_core::geometry::Transform3,
    scale: f64,
    tol: &Tolerance,
) -> Result<f64, String> {
    match (outcome, expect) {
        (Ok(body), Expect::Body { volume, area }) => {
            if body.faces.is_empty() {
                return Err("served an empty body".into());
            }
            let issues = body.validate_all();
            if let Some(first) = issues.first() {
                return Err(format!("validate_all: {first}"));
            }
            let measured = body;
            let v =
                popkern_measure::body_volume(measured, tol).map_err(|e| format!("volume: {e}"))?;
            let a = popkern_measure::body_surface_area(measured, tol)
                .map_err(|e| format!("area: {e}"))?;
            let (want_v, want_a) = (volume * scale.powi(3), area * scale * scale);
            let (ev, ea) = ((v - want_v).abs() / want_v, (a - want_a).abs() / want_a);
            // The boundary is known to the tolerance, so the volume is
            // known to about tolerance × area beyond the measure's own
            // precision (a pole cap 1e-4 high: 4%).
            let bound_v = MEASURE_REL * want_v + tol.linear() * want_a;
            // And the area to about tolerance × perimeter (a square's
            // for that area).
            let bound_a = MEASURE_REL * want_a + 4.0 * tol.linear() * want_a.sqrt();
            if (v - want_v).abs() > bound_v || (a - want_a).abs() > bound_a {
                return Err(format!(
                    "volume {v} (want {want_v}, rel {ev:.1e}; mesh reads {}), area {a} (want {want_a}, rel {ea:.1e})",
                    mesh_volume(measured).map_or("?".into(), |m| format!("{m}"))
                ));
            }
            Ok(ev.max(ea))
        }
        (Ok(body), Expect::Empty) if body.faces.is_empty() => Ok(0.0),
        (Err(BoolError::NonManifoldResult { contacts }), Expect::NonManifold(at)) => {
            let want: Vec<Point3> = at.iter().map(|p| place.apply_point(p)).collect();
            let near = |p: &Point3, q: &Point3| p.distance_to(q) <= 1e-5 * scale.max(1.0);
            let missing: Vec<_> = want
                .iter()
                .filter(|w| !contacts.iter().any(|c| near(&c.location, w)))
                .collect();
            let extra: Vec<_> = contacts
                .iter()
                .filter(|c| !want.iter().any(|w| near(&c.location, w)))
                .map(|c| c.location)
                .collect();
            if missing.is_empty() && extra.is_empty() {
                Ok(0.0)
            } else {
                Err(format!(
                    "contacts: missing {missing:?}, unexpected {extra:?}"
                ))
            }
        }
        (Ok(body), e) => Err(format!(
            "served V{} E{} F{}, expected {e:?}",
            body.vertices.len(),
            body.edges.len(),
            body.faces.len()
        )),
        (Err(err), e) => Err(format!("declined `{err}`, expected {e:?}")),
    }
}

/// Seconds a run may take before it counts as failed: a hang, not a
/// slow run (a kilometre-scale run under full-suite load took over 60).
const RUN_LIMIT_S: u64 = 300;

/// One boolean and its verdict, in a thread with a deadline: a run that
/// does not finish fails its row (the thread is left to the process's
/// end), so one stuck case cannot hide the rest of the table.
#[allow(clippy::too_many_arguments)]
fn run_with_watchdog(
    op: popkern_bool::BoolOperation,
    a: BRep,
    b: BRep,
    tol: Tolerance,
    options: BoolOptions,
    expect: Expect,
    place: popkern_core::geometry::Transform3,
    scale: f64,
) -> Result<f64, String> {
    let (send, receive) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let outcome = boolean_with_options(op, a, b, tol, None, &options);
        let _ = send.send(judge(&outcome, &expect, &place, scale, &tol));
    });
    receive
        .recv_timeout(std::time::Duration::from_secs(RUN_LIMIT_S))
        .unwrap_or_else(|_| Err(format!("did not finish in {RUN_LIMIT_S} s")))
}

/// The volume of a fine tessellation, independent of the measure's
/// quadrature: a disagreement attributes a wrong volume to the measure,
/// an agreement to the body.
fn mesh_volume(body: &BRep) -> Option<f64> {
    // A chord of 1e-4 of the body's size: fine enough to attribute a
    // wrong volume (the chord's error is under 1e-7 of it), coarse
    // enough to stay a fraction of a second.
    let size = popkern_measure::body_bounding_box(body)
        .map(|b| b.min.distance_to(&b.max))
        .unwrap_or(1.0);
    let options = popkern_mesh::TessellationOptions {
        chord_height: 1e-4 * size.max(1e-9),
        ..Default::default()
    };
    let mesh = popkern_mesh::tessellate(body, &options).ok()?;
    let p = |i: usize| mesh.vertices[i].position;
    Some(
        mesh.triangles
            .iter()
            .map(|t| {
                let (a, b, c) = (p(t[0]), p(t[1]), p(t[2]));
                (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                    + a.z * (b.x * c.y - b.y * c.x))
                    / 6.0
            })
            .sum::<f64>()
            .abs(),
    )
}

struct RowResult {
    name: String,
    status: Status,
    runs: usize,
    passed: usize,
    worst: f64,
    first_failure: Option<String>,
}

fn run_case(case: &Case, policy: NonManifoldPolicy, expect: &Expect) -> RowResult {
    let mut row = RowResult {
        name: format!("{} [{policy:?}]", case.name),
        status: if policy == NonManifoldPolicy::General {
            case.general_status
        } else {
            case.status
        },
        runs: 0,
        passed: 0,
        worst: 0.0,
        first_failure: None,
    };
    let base = Tolerance::default();
    // Each self-turn, and for a symmetric operation the operands
    // swapped (the plan's 8.2 probe: the answer cannot depend on which
    // operand is first).
    let symmetric = matches!(
        case.op,
        popkern_bool::BoolOperation::Union | popkern_bool::BoolOperation::Intersection
    );
    let turns: Vec<(Option<&SelfTurn>, bool)> = std::iter::once((None, false))
        .chain(case.turns.iter().map(|t| (Some(t), false)))
        .chain(symmetric.then_some((None, true)))
        .collect();
    let only_place = std::env::var("TANGENCY_CORPUS_PLACE").ok();
    for (pname, place, scale) in placements() {
        if only_place.as_deref().is_some_and(|o| o != pname) {
            continue;
        }
        let tol = base.with_linear(base.linear() * scale.max(1.0));
        for &(turn, swapped) in &turns {
            let pair = (case.build)(base);
            let pair = match turn {
                Some(t) => apply_turn(pair, t),
                None => pair,
            };
            let pair = if swapped { (pair.1, pair.0) } else { pair };
            let a = pair.0.transform(&place).expect("place a");
            let b = pair.1.transform(&place).expect("place b");
            let options = BoolOptions {
                non_manifold: policy,
            };
            let started = std::time::Instant::now();
            let verdict =
                run_with_watchdog(case.op, a, b, tol, options, expect.clone(), place, scale);
            if std::env::var_os("TANGENCY_CORPUS_TRACE").is_some() {
                eprintln!(
                    "  {} {pname} {turn:?} swapped {swapped}: {:.0} ms {}",
                    row.name,
                    started.elapsed().as_secs_f64() * 1e3,
                    verdict.as_ref().map_or_else(|e| e.as_str(), |_| "ok")
                );
            }
            row.runs += 1;
            match verdict {
                Ok(err) => {
                    row.passed += 1;
                    row.worst = row.worst.max(err);
                }
                Err(why) => {
                    if row.first_failure.is_none() {
                        row.first_failure =
                            Some(format!("{pname} turn {turn:?} swapped {swapped}: {why}"));
                    }
                }
            }
        }
    }
    row
}

#[test]
fn tangency_corpus() {
    let mut rows = Vec::new();
    // `TANGENCY_CORPUS_CASE=<name>` runs one case: for reproducing a
    // row, not a gate (the ratchet then reports only that case).
    let only = std::env::var("TANGENCY_CORPUS_CASE").ok();
    let jobs: Vec<(Case, NonManifoldPolicy, Expect)> = cases()
        .into_iter()
        .filter(|c| only.as_deref().is_none_or(|o| o == c.name))
        .flat_map(|c| {
            let mut out = vec![(NonManifoldPolicy::Decline, c.decline.clone())];
            if let Some(g) = &c.general {
                out.push((NonManifoldPolicy::General, g.clone()));
            }
            out.into_iter()
                .map(|(p, e)| (c.clone(), p, e))
                .collect::<Vec<_>>()
        })
        .collect();
    // The rows run side by side (each boolean is independent), at most
    // eight at once, and report in the table's order.
    let width = std::thread::available_parallelism().map_or(4, |n| n.get().clamp(1, 8));
    for chunk in jobs.chunks(width) {
        let done: Vec<RowResult> = std::thread::scope(|s| {
            let handles: Vec<_> = chunk
                .iter()
                .map(|(c, p, e)| s.spawn(move || run_case(c, *p, e)))
                .collect();
            handles
                .into_iter()
                .map(|h| h.join().expect("row thread"))
                .collect()
        });
        rows.extend(done);
    }
    let mut problems = Vec::new();
    for r in &rows {
        let all = r.passed == r.runs;
        eprintln!(
            "{:<44} {:>10?} {:>3}/{:<3} worst {:.1e} {}",
            r.name,
            r.status,
            r.passed,
            r.runs,
            r.worst,
            r.first_failure.as_deref().unwrap_or("")
        );
        match (r.status, all) {
            (Status::Gate, false) => problems.push(format!(
                "gate {} failed: {}",
                r.name,
                r.first_failure.as_deref().unwrap_or("")
            )),
            (Status::Target(stage), true) => problems.push(format!(
                "target {} ({stage}) now passes everywhere: promote it to a gate",
                r.name
            )),
            _ => {}
        }
    }
    assert!(problems.is_empty(), "{}", problems.join("\n"));
}
