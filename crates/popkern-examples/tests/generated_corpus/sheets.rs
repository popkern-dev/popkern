//! Sheet scenes (kind H, design §3): a planar quad and a bicubic NURBS
//! dome against every solid kind, crossing, tangent or disjoint, judged
//! by the sheet booleans' own identities.

use popkern_bool::{BoolOperation, SheetBooleanResult, boolean_with_sheets};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use serde::{Deserialize, Serialize};

use super::geometry::{Steer, inside, op_point};
use super::judge::{Failure, MEASURE_REL, Verdict};
use super::operand::{Kind, operand};
use super::rng::{Rng, stream_id};

/// The two sheet operands.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SheetKind {
    Quad,
    Dome,
}

/// The sheet configurations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SheetConfig {
    Cross,
    Tangent,
    Apart,
}

/// A sheet scene: plain data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SheetScene {
    pub sheet: SheetKind,
    pub solid: Kind,
    pub config: SheetConfig,
    pub variant: u8,
    pub index: u64,
    /// The sheet's placement (row-major 4×4) in the solid's frame.
    pub place: [[f64; 4]; 4],
    /// The quad's plane distance from a ball's centre, for the disc area
    /// closed form (negative: no closed form).
    pub sphere_d: f64,
}

/// The sheet, centred at the origin in the xy plane, half-width `h`.
fn sheet(kind: SheetKind, h: f64, tol: Tolerance) -> BRep {
    match kind {
        SheetKind::Quad => BRep::sheet_quad(-h, h, -h, h, 0.0, tol).expect("quad"),
        SheetKind::Dome => {
            let src = format!(
                r#"let m = model();
let z = [[0.0, 0.1, 0.12, 0.0], [0.08, 0.3, 0.25, 0.05], [0.05, 0.28, 0.35, 0.1], [0.0, 0.06, 0.1, 0.02]];
let net = [];
for i in 0..4 {{ let row = []; for j in 0..4 {{ row.push([{w} * i / 3.0 - {h}, {w} * j / 3.0 - {h}, z[i][j] * {h}]); }} net.push(row); }}
m.sheet_nurbs(net, #{{ degree: 3 }});
m"#,
                w = 2.0 * h,
                h = h
            );
            let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
            let run = host
                .run(&src, "dome_sheet")
                .unwrap_or_else(|e| panic!("dome sheet: {e}"));
            run.model.expect("model").to_parts_all().expect("parts")[0]
                .brep()
                .clone()
        }
    }
}

/// Draw a sheet scene: `None` when the draw is rejected.
pub fn draw_sheet(
    sheet_kind: SheetKind,
    solid: Kind,
    config: SheetConfig,
    seed: u64,
    index: u64,
) -> Option<SheetScene> {
    let mut rng = Rng::new(
        seed,
        stream_id(&format!("{sheet_kind:?}x{solid:?}-{config:?}")),
        index,
    );
    let variant = rng.below(3) as u8;
    let op = operand(solid, variant);
    let (place, sphere_d) = match config {
        SheetConfig::Cross => {
            let (lo, hi) = super::distance::mesh_box(&op.mesh);
            let mut at = None;
            for _ in 0..64 {
                let p = Point3::new(
                    rng.range(lo.x, hi.x),
                    rng.range(lo.y, hi.y),
                    rng.range(lo.z, hi.z),
                );
                if inside(&op.mesh, p) {
                    at = Some(p);
                    break;
                }
            }
            let at = at?;
            let d = rng.direction();
            let n = Vec3::new(d[0], d[1], d[2]);
            let turn = rotation_onto(Vec3::new(0.0, 0.0, 1.0), n);
            let place = Transform3::translation(at.x, at.y, at.z).compose(&turn);
            let sphere_d = if solid == Kind::S && sheet_kind == SheetKind::Quad {
                (at - Point3::new(0.0, 0.0, 0.0)).dot(&n).abs()
            } else {
                -1.0
            };
            (place, sphere_d)
        }
        SheetConfig::Tangent => {
            if sheet_kind != SheetKind::Quad {
                return None;
            }
            let sp = op_point(&op, Steer::Free, &mut rng)?;
            // The plane touches only where the solid is convex there.
            let (lo, hi) = super::geometry::tangent_eigen(&sp.t, sp.n);
            if !(hi < -1e-3 && lo < -1e-3) {
                return None;
            }
            let turn = rotation_onto(Vec3::new(0.0, 0.0, 1.0), sp.n);
            (
                Transform3::translation(sp.p.x, sp.p.y, sp.p.z).compose(&turn),
                -1.0,
            )
        }
        SheetConfig::Apart => {
            let d = rng.direction();
            let far = 2.0 * op.size;
            (
                Transform3::translation(d[0] * far, d[1] * far, d[2] * far),
                -1.0,
            )
        }
    };
    Some(SheetScene {
        sheet: sheet_kind,
        solid,
        config,
        variant,
        index,
        place: place.matrix,
        sphere_d,
    })
}

fn rotation_onto(from: Vec3, to: Vec3) -> Transform3 {
    let axis = from.cross(&to);
    let c = from.dot(&to).clamp(-1.0, 1.0);
    if axis.length() < 1e-12 {
        if c > 0.0 {
            return Transform3::identity();
        }
        return Transform3::rotation_x(std::f64::consts::PI);
    }
    Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), axis, c.acos()).expect("axis")
}

fn area(b: &BRep, tol: &Tolerance) -> f64 {
    if b.faces.is_empty() {
        return 0.0;
    }
    popkern_measure::body_surface_area(b, tol).unwrap_or(f64::NAN)
}

/// Judge a sheet scene.
pub fn judge_sheet(scene: &SheetScene) -> Verdict {
    super::judge::caught(|| judge_sheet_scene(scene))
}

fn judge_sheet_scene(scene: &SheetScene) -> Verdict {
    let tol = Tolerance::default();
    let op = operand(scene.solid, scene.variant);
    let h = 1.5 * op.size;
    let s = sheet(scene.sheet, h, tol)
        .transform(&Transform3 {
            matrix: scene.place,
        })
        .expect("place sheet");
    let solid = op.brep.clone();
    let sheet_area = area(&s, &tol);
    let mut failures = Vec::new();
    let mut fail = |oracle: &'static str, detail: String| {
        failures.push(Failure {
            oracle,
            op: None,
            detail,
        })
    };
    let bound = 4.0 * tol.linear() * 8.0 * h + MEASURE_REL * sheet_area;
    let sheet_result = |op: BoolOperation,
                        fail: &mut dyn FnMut(&'static str, String),
                        name: &str|
     -> Option<f64> {
        match boolean_with_sheets(op, &s, &solid, tol) {
            Ok(SheetBooleanResult::OpenShell(b)) | Ok(SheetBooleanResult::Solid(b)) => {
                if let Some(first) = b.validate_all().first() {
                    fail("sheet-validity", format!("{name}: validate_all: {first}"));
                }
                Some(area(&b, &tol))
            }
            Ok(SheetBooleanResult::Pieces(ps)) => Some(ps.iter().map(|b| area(b, &tol)).sum()),
            Ok(other) => {
                fail(
                    "sheet-shape",
                    format!("{name}: {:?}", std::mem::discriminant(&other)),
                );
                None
            }
            Err(e) => {
                fail("sheet-declined", format!("{name}: `{e}`"));
                None
            }
        }
    };
    let inside_area = sheet_result(BoolOperation::Intersection, &mut fail, "sheet ∩ solid");
    let outside_area = sheet_result(BoolOperation::Subtraction, &mut fail, "sheet − solid");
    if let (Some(i), Some(o)) = (inside_area, outside_area)
        && (i + o - sheet_area).abs() > bound
    {
        fail(
            "sheet-identity",
            format!("area in {i} + out {o} = {} vs sheet {sheet_area}", i + o),
        );
    }
    if let Some(i) = inside_area {
        match scene.config {
            SheetConfig::Tangent | SheetConfig::Apart if i > bound => fail(
                "sheet-closed-form",
                format!("a touching or apart sheet has {i} inside"),
            ),
            SheetConfig::Cross if scene.sphere_d >= 0.0 => {
                // A plane through a ball cuts a disc of radius √(r² − d²).
                let r = op
                    .brep
                    .faces
                    .values()
                    .find_map(|f| match f.surface.as_ref().map(|s| &s.geometry) {
                        Some(popkern_core::geometry::surface::SurfaceGeometry::Sphere {
                            radius,
                            ..
                        }) => Some(*radius),
                        _ => None,
                    })
                    .unwrap_or(f64::NAN);
                let want = std::f64::consts::PI * (r * r - scene.sphere_d * scene.sphere_d);
                if (i - want).abs() > bound {
                    fail("sheet-closed-form", format!("disc area {i} want {want}"));
                }
            }
            _ => {}
        }
    }
    // The solid split by the sheet: its pieces hold its volume.
    match boolean_with_sheets(BoolOperation::Subtraction, &solid, &s, tol) {
        Ok(SheetBooleanResult::Pieces(ps)) => {
            let sum: f64 = ps
                .iter()
                .filter_map(|b| popkern_measure::body_volume(b, &tol).ok())
                .sum();
            if (sum - op.volume).abs() > MEASURE_REL * op.volume + 4.0 * tol.linear() * op.area {
                fail(
                    "sheet-split",
                    format!("pieces {sum} vs solid {}", op.volume),
                );
            }
            for p in &ps {
                if let Some(first) = p.validate_all().first() {
                    fail("sheet-validity", format!("solid − sheet piece: {first}"));
                }
            }
        }
        Ok(SheetBooleanResult::Solid(b)) => {
            let v = popkern_measure::body_volume(&b, &tol).unwrap_or(f64::NAN);
            if (v - op.volume).abs() > MEASURE_REL * op.volume + 4.0 * tol.linear() * op.area {
                fail("sheet-split", format!("solid {v} vs {}", op.volume));
            }
        }
        Ok(other) => fail(
            "sheet-shape",
            format!("solid − sheet: {:?}", std::mem::discriminant(&other)),
        ),
        Err(e) => fail("sheet-declined", format!("solid − sheet: `{e}`")),
    }
    Verdict {
        failures,
        worst_rel: 0.0,
    }
}
