//! Three-body scenes (configuration X3, design §5): A and B crossing, C a
//! second operand of A's kind placed over a point inside both, so all
//! three overlap. The partition's N-ary results are judged against the
//! two-body folds (the identity question of plan W4).

use popkern_bool::{
    BoolError, BoolOperation, BoolOptions, NonManifoldPolicy, SetExpr, boolean_with_options,
    nary_boolean,
};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};
use serde::{Deserialize, Serialize};

use super::geometry::inside;
use super::judge::{Failure, MEASURE_REL, Verdict};
use super::operand::operand;
use super::rng::Rng;

/// The third body of a three-body scene.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Third {
    pub variant: u8,
    pub scale: f64,
    /// Its rigid placement (row-major 4×4), before the scene similarity.
    pub place: [[f64; 4]; 4],
}

/// C's placement for A (in its frame) and B (placed): a point inside both,
/// C centred there, turned randomly. `None` when A and B share no probe.
pub fn place_third(
    a_mesh: &popkern_mesh::TriMesh,
    b_mesh: &popkern_mesh::TriMesh,
    a_kind: super::operand::Kind,
    a_variant: u8,
    rng: &mut Rng,
) -> Option<Third> {
    let common: Vec<Point3> = b_mesh
        .vertices
        .iter()
        .step_by((b_mesh.vertices.len() / 200).max(1))
        .map(|v| v.position)
        .filter(|p| inside(a_mesh, *p))
        .collect();
    let at = *common.get(rng.below(common.len().max(1)))?;
    let variant = (a_variant + 1) % 3;
    let c = operand(a_kind, variant);
    let scale = 0.6;
    let axis = rng.direction();
    let turn = Transform3::rotation_around_axis(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(axis[0], axis[1], axis[2]),
        rng.range(0.0, std::f64::consts::TAU),
    )
    .ok()?;
    let (lo, hi) = super::distance::mesh_box(&c.mesh);
    let centre = Point3::new(
        0.5 * (lo.x + hi.x),
        0.5 * (lo.y + hi.y),
        0.5 * (lo.z + hi.z),
    );
    let to_origin =
        Transform3::translation(-centre.x * scale, -centre.y * scale, -centre.z * scale);
    let place = Transform3::translation(at.x, at.y, at.z).compose(&turn.compose(&to_origin));
    Some(Third {
        variant,
        scale,
        place: place.matrix,
    })
}

/// The third body placed in the scene's final frame.
pub fn third_body(
    kind: super::operand::Kind,
    third: &Third,
    sim: &Transform3,
) -> (BRep, popkern_mesh::TriMesh) {
    let c = operand(kind, third.variant);
    let t = sim
        .compose(&Transform3 {
            matrix: third.place,
        })
        .compose(&Transform3::scale_uniform(third.scale));
    let brep = c.brep.transform(&t).expect("place C");
    let mut mesh = c.mesh.clone();
    for v in &mut mesh.vertices {
        v.position = t.apply_point(&v.position);
    }
    (brep, mesh)
}

fn volume(b: &BRep, tol: &Tolerance) -> Option<f64> {
    if b.faces.is_empty() {
        return Some(0.0);
    }
    popkern_measure::body_volume(b, tol).ok()
}

/// Judge the N-ary results against their pairwise folds.
pub fn judge_three(a: &BRep, b: &BRep, c: &BRep, area: f64) -> Verdict {
    let tol = Tolerance::default();
    let opts = BoolOptions {
        non_manifold: NonManifoldPolicy::Decline,
    };
    let mut failures = Vec::new();
    let mut fail = |oracle: &'static str, detail: String| {
        failures.push(Failure {
            oracle,
            op: None,
            detail,
        })
    };
    let fold =
        |op: BoolOperation, x: BRep, y: BRep| boolean_with_options(op, x, y, tol, None, &opts);
    let cases: [(&str, SetExpr, Result<BRep, BoolError>); 3] = [
        (
            "union",
            SetExpr::body(0).or(SetExpr::body(1)).or(SetExpr::body(2)),
            fold(BoolOperation::Union, a.clone(), b.clone())
                .and_then(|ab| fold(BoolOperation::Union, ab, c.clone())),
        ),
        (
            "common",
            SetExpr::body(0).and(SetExpr::body(1)).and(SetExpr::body(2)),
            fold(BoolOperation::Intersection, a.clone(), b.clone())
                .and_then(|ab| fold(BoolOperation::Intersection, ab, c.clone())),
        ),
        (
            "union-minus",
            SetExpr::body(0)
                .or(SetExpr::body(1))
                .minus(SetExpr::body(2)),
            fold(BoolOperation::Union, a.clone(), b.clone())
                .and_then(|ab| fold(BoolOperation::Subtraction, ab, c.clone())),
        ),
    ];
    for (name, expr, folded) in cases {
        let nary = nary_boolean(&[a, b, c], &expr, tol, &opts);
        let check = |r: &Result<BRep, BoolError>,
                     path: &str,
                     fail: &mut dyn FnMut(&'static str, String)|
         -> Option<f64> {
            match r {
                Ok(body) => {
                    if let Some(first) = body.validate_all().first() {
                        fail(
                            "three-validity",
                            format!("{name} {path}: validate_all: {first}"),
                        );
                    }
                    if !body.faces.is_empty() && !matches!(body.body_kind(), BodyKind::Solid) {
                        fail(
                            "three-validity",
                            format!("{name} {path}: kind {:?}", body.body_kind()),
                        );
                    }
                    volume(body, &tol)
                }
                Err(e) => {
                    fail("three-declined", format!("{name} {path}: `{e}`"));
                    None
                }
            }
        };
        let vn = check(&nary, "n-ary", &mut fail);
        let vf = check(&folded, "fold", &mut fail);
        if let (Some(x), Some(y)) = (vn, vf) {
            let bound = MEASURE_REL * x.abs().max(y.abs()) + 4.0 * tol.linear() * area;
            if (x - y).abs() > bound {
                fail("three-fold", format!("{name}: n-ary {x} vs fold {y}"));
            }
        }
        // Determinism of the n-ary path.
        let again = nary_boolean(
            &[a, b, c],
            &match name {
                "union" => SetExpr::body(0).or(SetExpr::body(1)).or(SetExpr::body(2)),
                "common" => SetExpr::body(0).and(SetExpr::body(1)).and(SetExpr::body(2)),
                _ => SetExpr::body(0)
                    .or(SetExpr::body(1))
                    .minus(SetExpr::body(2)),
            },
            tol,
            &opts,
        );
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
        if sig(&nary) != sig(&again) {
            fail(
                "determinism",
                format!("{name} n-ary: {} vs {}", sig(&nary), sig(&again)),
            );
        }
    }
    Verdict {
        failures,
        worst_rel: 0.0,
    }
}
