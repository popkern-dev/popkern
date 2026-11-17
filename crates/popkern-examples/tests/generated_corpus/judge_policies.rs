//! The policy oracle and the mesh volume (split from `judge.rs` for the
//! file-size limit).

use popkern_bool::{BoolError, BoolOptions, NonManifoldPolicy, boolean_bodies};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::geometry::mesh_signed_volume;
use super::judge::{Failure, Op, fail, op_of, run, volume};
use super::scene::{Config, Placed, Scene};

/// Under `General` one general body of the expected volume; under
/// `SplitToManifold` the pieces, their volumes summing to it. Returns the
/// general body's volume (for the identities) when it served.
pub(super) fn policies(
    placed: &Placed,
    op: Op,
    va: f64,
    vb: f64,
    scene: &Scene,
    failures: &mut Vec<Failure>,
    bound: impl Fn(f64) -> f64,
) -> Option<f64> {
    let tol = Tolerance::default();
    let config = match scene.cell.config {
        Config::K if scene.aligned => Config::I,
        Config::K => Config::E,
        c => c,
    };
    let total = match (config, op) {
        (Config::E, Op::Union) => Some(va + vb),
        (Config::I, Op::AminusB) => Some(va - vb),
        _ => None,
    };
    let mut served = None;
    match run(placed, op, NonManifoldPolicy::General, tol) {
        Ok(body) => {
            let issues = body.validate_all();
            if let Some(first) = issues.first() {
                fail(
                    failures,
                    "policy-general",
                    Some(op),
                    format!("validate_all: {first}"),
                );
            }
            match (volume(&body, &tol), total) {
                (Some(v), Some(t)) if (v - t).abs() > bound(t) => fail(
                    failures,
                    "policy-general",
                    Some(op),
                    format!("volume {v} want {t}"),
                ),
                (None, _) => fail(failures, "policy-general", Some(op), "volume failed".into()),
                (v, _) => served = v,
            }
        }
        Err(e) => fail(
            failures,
            "policy-general",
            Some(op),
            format!("declined `{e}`"),
        ),
    }
    let total = total.or(served)?;
    let (bop, swap) = op_of(op);
    let (x, y) = if swap {
        (placed.b.clone(), placed.a.clone())
    } else {
        (placed.a.clone(), placed.b.clone())
    };
    match boolean_bodies(
        bop,
        x,
        y,
        tol,
        None,
        &BoolOptions {
            non_manifold: NonManifoldPolicy::SplitToManifold,
        },
    ) {
        Ok(bodies) => {
            let sum: f64 = bodies.iter().filter_map(|b| volume(b, &tol)).sum();
            if bodies.iter().any(|b| !b.validate_all().is_empty()) {
                fail(
                    failures,
                    "policy-split",
                    Some(op),
                    "a piece fails validate_all".into(),
                );
            }
            if (sum - total).abs() > bound(total) {
                fail(
                    failures,
                    "policy-split",
                    Some(op),
                    format!("pieces sum {sum} want {total}"),
                );
            }
        }
        Err(BoolError::SplitIntoBodies { .. }) if (config, op) == (Config::I, Op::AminusB) => {
            // A solid touching itself has no separate bodies.
        }
        Err(BoolError::SplitIntoBodies { .. }) if scene.cell.config == Config::V => {}
        Err(e) => fail(
            failures,
            "policy-split",
            Some(op),
            format!("declined `{e}`"),
        ),
    }
    served
}

pub(super) fn mesh_volume(body: &BRep) -> Option<(f64, f64)> {
    // The size from a coarse tessellation (`body_bounding_box` reads
    // vertices only, PK-1258).
    // The vertex box has the right order of magnitude for a coarse chord.
    let rough = popkern_measure::body_bounding_box(body)
        .map(|b| b.min.distance_to(&b.max))
        .unwrap_or(1.0)
        .max(1e-9);
    let coarse = popkern_mesh::tessellate(
        body,
        &popkern_mesh::TessellationOptions {
            chord_height: 1e-2 * rough,
            ..Default::default()
        },
    )
    .ok()?;
    let (lo, hi) = super::distance::mesh_box(&coarse);
    let chord = 5e-4 * lo.distance_to(&hi).max(1e-9);
    let options = popkern_mesh::TessellationOptions {
        chord_height: chord,
        ..Default::default()
    };
    let mesh = popkern_mesh::tessellate(body, &options).ok()?;
    Some((mesh_signed_volume(&mesh).abs(), chord))
}
