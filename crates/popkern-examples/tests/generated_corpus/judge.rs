//! Running a scene and judging it (design §6.2, §6.3).

use popkern_bool::{
    BoolError, BoolOperation, BoolOptions, NonManifoldPolicy, boolean_with_options,
};
use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};
use serde::Serialize;

use super::scene::{Config, Placed, Scene, place};

/// The measure's own relative precision (as the tangency corpus): the
/// design's 1e-9 for analytic-bounded results is recorded per run as the
/// observed error, and tightened once the inventory shows it is met.
pub const MEASURE_REL: f64 = 1e-6;

/// The operations a scene runs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize)]
pub enum Op {
    Union,
    Common,
    AminusB,
    BminusA,
}

pub const OPS: [Op; 4] = [Op::Union, Op::Common, Op::AminusB, Op::BminusA];

/// What an operation should give.
#[derive(Debug, Clone, PartialEq)]
pub enum Want {
    /// A served solid of this volume (`None`: any volume, identities
    /// decide).
    Body(Option<f64>),
    /// Empty.
    Empty,
    /// A non-manifold decline naming the contact.
    NonManifold,
}

/// One oracle failure.
#[derive(Debug, Clone, Serialize)]
pub struct Failure {
    pub oracle: &'static str,
    pub op: Option<Op>,
    pub detail: String,
}

fn band(offset: f64) -> i8 {
    if offset.abs() <= 1.0 {
        0
    } else if offset > 0.0 {
        1
    } else {
        -1
    }
}

/// The expected outcome of `op` under `Decline` (design §6.3).
pub fn want(scene: &Scene, op: Op, va: f64, vb: f64) -> Want {
    use Op::*;
    let b = band(scene.offset);
    // A curve tangency answers as a point touch of its orientation.
    let config = match scene.cell.config {
        Config::K if scene.aligned => Config::I,
        Config::K => Config::E,
        c => c,
    };
    match (config, b) {
        (Config::D, _) | (Config::E, 1) => match op {
            Union => Want::Body(Some(va + vb)),
            Common => Want::Empty,
            AminusB => Want::Body(Some(va)),
            BminusA => Want::Body(Some(vb)),
        },
        (Config::E, 0) => match op {
            Union => Want::NonManifold,
            Common => Want::Empty,
            AminusB => Want::Body(Some(va)),
            BminusA => Want::Body(Some(vb)),
        },
        (Config::I, 0) => match op {
            Union => Want::Body(Some(va)),
            Common => Want::Body(Some(vb)),
            AminusB => Want::NonManifold,
            BminusA => Want::Empty,
        },
        (Config::N, _) | (Config::I, 1) => match op {
            Union => Want::Body(Some(va)),
            Common => Want::Body(Some(vb)),
            AminusB => Want::Body(Some(va - vb)),
            BminusA => Want::Empty,
        },
        // A coincident patch: a pure symmetry gives A back; a prism
        // (box, cylinder) shifted along its axis by a fraction f gives
        // the closed forms; other shifts the identities decide.
        (Config::O, _) if scene.shift == 0.0 => match op {
            Union | Common => Want::Body(Some(va)),
            AminusB | BminusA => Want::Empty,
        },
        (Config::O, _)
            if matches!(
                scene.cell.a,
                super::operand::Kind::P | super::operand::Kind::Y
            ) =>
        {
            let f = scene.shift;
            match op {
                Union => Want::Body(Some(va * (1.0 + f))),
                Common => Want::Body(Some(va * (1.0 - f))),
                AminusB | BminusA => Want::Body(Some(va * f)),
            }
        }
        // A tangent crossing at its touch: the gap and overlap sectors
        // are wedges meeting at the contact. With normals opposed the
        // union's outside and the common's inside pinch there; with
        // them aligned (B inside A there) both differences do.
        (Config::V, 0) => {
            let pinched = if scene.aligned {
                matches!(op, AminusB | BminusA)
            } else {
                matches!(op, Union | Common)
            };
            if pinched {
                Want::NonManifold
            } else {
                Want::Body(None)
            }
        }
        // Overlapping near misses, tangent crossings away from the
        // touch, transversal crossings: valid bodies, volumes decided
        // by the identities.
        _ => Want::Body(None),
    }
}

pub(super) fn op_of(op: Op) -> (BoolOperation, bool) {
    match op {
        Op::Union => (BoolOperation::Union, false),
        Op::Common => (BoolOperation::Intersection, false),
        Op::AminusB => (BoolOperation::Subtraction, false),
        Op::BminusA => (BoolOperation::Subtraction, true),
    }
}

/// The outcome of one boolean.
pub struct Outcome {
    pub op: Op,
    pub result: Result<BRep, BoolError>,
    pub volume: Option<f64>,
}

pub fn run(
    placed: &Placed,
    op: Op,
    policy: NonManifoldPolicy,
    tol: Tolerance,
) -> Result<BRep, BoolError> {
    let (bop, swap) = op_of(op);
    let (x, y) = if swap {
        (placed.b.clone(), placed.a.clone())
    } else {
        (placed.a.clone(), placed.b.clone())
    };
    boolean_with_options(
        bop,
        x,
        y,
        tol,
        None,
        &BoolOptions {
            non_manifold: policy,
        },
    )
}

pub(super) fn volume(body: &BRep, tol: &Tolerance) -> Option<f64> {
    if body.faces.is_empty() {
        return Some(0.0);
    }
    popkern_measure::body_volume(body, tol).ok()
}

/// The scene's verdict: every oracle failure (empty: passed).
#[derive(Debug, Clone, Serialize)]
pub struct Verdict {
    pub failures: Vec<Failure>,
    /// The worst relative volume error against a known volume.
    pub worst_rel: f64,
}

pub(super) fn fail(out: &mut Vec<Failure>, oracle: &'static str, op: Option<Op>, detail: String) {
    out.push(Failure { oracle, op, detail });
}

/// Run and judge one scene (no OCC oracle).
pub fn judge(scene: &Scene) -> Verdict {
    judge_with(scene, false)
}

/// Run and judge one scene; `occ` adds the OCC oracle on the union. A
/// panic anywhere in the kernel is the scene's failure (oracle `panic`,
/// the message as its detail), so one scene never loses its pair's report.
pub fn judge_with(scene: &Scene, occ: bool) -> Verdict {
    caught(|| judge_scene(scene, occ))
}

/// `judge` under `catch_unwind`: a panic becomes a `panic` failure.
pub fn caught(judge: impl FnOnce() -> Verdict) -> Verdict {
    match std::panic::catch_unwind(std::panic::AssertUnwindSafe(judge)) {
        Ok(v) => v,
        Err(payload) => {
            let msg = payload
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| payload.downcast_ref::<&str>().map(|s| s.to_string()))
                .unwrap_or_else(|| "non-string panic payload".to_string());
            Verdict {
                failures: vec![Failure {
                    oracle: "panic",
                    op: None,
                    detail: msg,
                }],
                worst_rel: 0.0,
            }
        }
    }
}

fn judge_scene(scene: &Scene, occ: bool) -> Verdict {
    let tol = Tolerance::default();
    let placed = place(scene);
    if let Some(third) = &scene.third {
        let sim = popkern_core::geometry::Transform3 {
            matrix: scene.scene,
        };
        let (c, _) = super::three::third_body(scene.cell.a, third, &sim);
        let s2 = scene.scale * scene.scale;
        let area = (placed.a_op.area * (1.0 + third.scale * third.scale)
            + placed.b_op.area * scene.b_scale * scene.b_scale)
            * s2;
        return super::three::judge_three(&placed.a, &placed.b, &c, area);
    }
    let s3 = scene.scale.powi(3);
    let s2 = scene.scale * scene.scale;
    let bs3 = scene.b_scale.powi(3);
    let va = placed.a_op.volume * s3;
    let vb = placed.b_op.volume * bs3 * s3;
    let area = (placed.a_op.area + placed.b_op.area * scene.b_scale * scene.b_scale) * s2;
    let bound = |v: f64| MEASURE_REL * v.abs().max(va.min(vb)) + 4.0 * tol.linear() * area;
    let mut failures = Vec::new();
    let mut worst = 0.0f64;
    let mut volumes = std::collections::HashMap::new();
    let contact = scene.contact.map(|c| Point3::new(c[0], c[1], c[2]));
    let mut union_sig = None;
    let sig = |r: &Result<BRep, BoolError>| match r {
        Ok(b) => format!(
            "V{} E{} F{} {:?}",
            b.vertices.len(),
            b.edges.len(),
            b.faces.len(),
            volume(b, &tol).map(f64::to_bits)
        ),
        Err(e) => format!("{e}"),
    };
    for op in super::judge::OPS {
        let expect = want(scene, op, va, vb);
        let result = run(&placed, op, NonManifoldPolicy::Decline, tol);
        if op == Op::Union {
            union_sig = Some(sig(&result));
        }
        match (&result, &expect) {
            (Ok(body), Want::Empty) => {
                if !body.faces.is_empty() {
                    fail(
                        &mut failures,
                        "policy",
                        Some(op),
                        format!("served F{} where empty is expected", body.faces.len()),
                    );
                } else {
                    volumes.insert(op, 0.0);
                }
            }
            (Ok(body), Want::Body(v)) => {
                if body.faces.is_empty() {
                    fail(
                        &mut failures,
                        "policy",
                        Some(op),
                        "served empty where a body is expected".into(),
                    );
                    continue;
                }
                let issues = body.validate_all();
                if let Some(first) = issues.first() {
                    fail(
                        &mut failures,
                        "validity",
                        Some(op),
                        format!("validate_all: {first}"),
                    );
                }
                if !matches!(body.body_kind(), BodyKind::Solid) {
                    fail(
                        &mut failures,
                        "validity",
                        Some(op),
                        format!("kind {:?}", body.body_kind()),
                    );
                }
                match volume(body, &tol) {
                    Some(got) => {
                        volumes.insert(op, got);
                        if let Some(want) = v {
                            let err = (got - want).abs();
                            worst = worst.max(err / want.abs().max(1e-300));
                            if err > bound(*want) {
                                fail(
                                    &mut failures,
                                    "closed-form",
                                    Some(op),
                                    format!("volume {got} want {want}"),
                                );
                            }
                        }
                        // Two volumes (design §6.2 oracle 8).
                        // An inscribed tessellation loses up to area × chord.
                        if let Some((mv, chord)) = super::judge_policies::mesh_volume(body) {
                            let mb = 1.5 * area * chord + bound(got);
                            if (mv - got).abs() > mb {
                                fail(
                                    &mut failures,
                                    "two-volumes",
                                    Some(op),
                                    format!("B-rep {got} mesh {mv}"),
                                );
                            }
                        }
                    }
                    None => fail(&mut failures, "validity", Some(op), "volume failed".into()),
                }
            }
            (Err(BoolError::NonManifoldResult { contacts }), Want::NonManifold) => {
                let reach = 1e-5 * scene.scale.max(1.0);
                let curve: Vec<Point3> = scene
                    .curve
                    .iter()
                    .map(|p| Point3::new(p[0], p[1], p[2]))
                    .collect();
                let on_contact = |q: &Point3| {
                    if curve.is_empty() {
                        contact.is_some_and(|c| q.distance_to(&c) <= reach)
                    } else {
                        // A circle's chords sag below it by r(1 − cos(π/16)).
                        let sag = curve[0].distance_to(&curve[1]) * 0.1;
                        super::geometry::polyline_distance(*q, &curve, scene.curve_closed)
                            <= reach + sag
                    }
                };
                if contact.is_some() && !contacts.iter().any(|k| on_contact(&k.location)) {
                    fail(
                        &mut failures,
                        "policy",
                        Some(op),
                        format!(
                            "contacts {:?} miss {contact:?}",
                            contacts.iter().map(|k| k.location).collect::<Vec<_>>()
                        ),
                    );
                }
                if let Some(v) = super::judge_policies::policies(
                    &placed,
                    op,
                    va,
                    vb,
                    scene,
                    &mut failures,
                    bound,
                ) {
                    volumes.insert(op, v);
                }
            }
            (Ok(body), Want::NonManifold) => fail(
                &mut failures,
                "policy",
                Some(op),
                format!(
                    "served F{} where a non-manifold decline is expected",
                    body.faces.len()
                ),
            ),
            (Err(e), w) => fail(
                &mut failures,
                "policy",
                Some(op),
                format!("declined `{e}`, expected {w:?}"),
            ),
        }
    }
    // OCC (oracle 9): the union and A measured by OpenCASCADE, as a ratio.
    if occ
        && let Some(&u) = volumes.get(&Op::Union)
        && u > 0.0
        && let Ok(Ok(union)) =
            std::panic::catch_unwind(|| run(&placed, Op::Union, NonManifoldPolicy::Decline, tol))
    {
        let tag = format!("{}-{}-{}", scene.cell.name(), scene.seed, scene.index);
        match (
            super::occ::measure(&union, &format!("{tag}-u")),
            super::occ::measure(&placed.a, &format!("{tag}-a")),
        ) {
            (Ok(ou), Ok(oa)) if oa > 0.0 => {
                let (ours, theirs) = (u / va, ou / oa);
                if (ours - theirs).abs() > 1e-6 * ours.abs().max(1.0) + bound(u) / va {
                    fail(
                        &mut failures,
                        "occ",
                        Some(Op::Union),
                        format!("V(A∪B)/V(A): popkern {ours} OCC {theirs}"),
                    );
                }
            }
            (Err(e), _) | (_, Err(e)) => fail(&mut failures, "occ", Some(Op::Union), e),
            _ => {}
        }
    }
    // Set identities (oracle 2).
    let g = |op: Op| volumes.get(&op).copied();
    if let (Some(u), Some(c)) = (g(Op::Union), g(Op::Common))
        && (u + c - va - vb).abs() > bound(va + vb)
    {
        fail(
            &mut failures,
            "identity",
            None,
            format!("V(A∪B)+V(A∩B) = {} vs V(A)+V(B) = {}", u + c, va + vb),
        );
    }
    if let (Some(d), Some(c)) = (g(Op::AminusB), g(Op::Common))
        && (d - (va - c)).abs() > bound(va)
    {
        fail(
            &mut failures,
            "identity",
            None,
            format!("V(A−B) = {d} vs V(A)−V(A∩B) = {}", va - c),
        );
    }
    if let (Some(d), Some(c)) = (g(Op::BminusA), g(Op::Common))
        && (d - (vb - c)).abs() > bound(vb)
    {
        fail(
            &mut failures,
            "identity",
            None,
            format!("V(B−A) = {d} vs V(B)−V(A∩B) = {}", vb - c),
        );
    }
    // Determinism (oracle 7): the union again, bit for bit.
    let again = sig(&run(&placed, Op::Union, NonManifoldPolicy::Decline, tol));
    if let Some(first) = union_sig
        && again != first
    {
        fail(
            &mut failures,
            "determinism",
            Some(Op::Union),
            format!("{first} vs {again}"),
        );
    }
    Verdict {
        failures,
        worst_rel: worst,
    }
}
