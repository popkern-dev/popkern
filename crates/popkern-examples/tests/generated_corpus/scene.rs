//! Scenes (design §4–§5): a cell, a seed and an index give a descriptor;
//! the descriptor gives the two operands placed exactly.

use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use serde::{Deserialize, Serialize};

use super::geometry::{Steer, inside, op_point, rotate_tensor, tangent_eigen};
use super::operand::{Kind, Operand, operand};
use super::rng::{Rng, stream_id};

/// The configurations (design §4).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Config {
    /// Disjoint.
    D,
    /// External point touch.
    E,
    /// Internal point touch.
    I,
    /// Tangent crossing point (relative curvature indefinite).
    V,
    /// Transversal crossing.
    X,
    /// Containment.
    N,
    /// Tangent along a curve (a template: a ruling or a circle).
    K,
    /// Three bodies: A and B crossing, C (A's kind) overlapping both.
    X3,
    /// Coincident patch: B is A's own operand moved by a symmetry (the
    /// surfaces coincide, the charts differ) or shifted along an axis
    /// (faces coincide in part).
    O,
}

/// A (pair, configuration) cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct Cell {
    pub a: Kind,
    pub b: Kind,
    pub config: Config,
}

impl Cell {
    pub fn name(&self) -> String {
        format!("{:?}x{:?}-{:?}", self.a, self.b, self.config)
    }
}

/// The near-miss offsets, in tolerances (positive apart, negative
/// overlapping); 0 is the exact touch.
pub const OFFSETS: [f64; 7] = [0.0, 0.5, -0.5, 2.0, -2.0, 10.0, -10.0];

/// One ladder rung or random draw's knobs (design §5).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Rung {
    /// Near-miss offset in tolerances.
    pub offset: f64,
    /// Scene scale.
    pub scale: f64,
    /// Where A's contact sits on its chart.
    pub steer_a: Steer,
    /// Where B's contact sits on its chart.
    pub steer_b: Steer,
}

impl Rung {
    pub const FREE: Rung = Rung {
        offset: 0.0,
        scale: 1.0,
        steer_a: Steer::Free,
        steer_b: Steer::Free,
    };
}

/// A scene descriptor: plain data, replayable from its JSON.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Scene {
    pub cell: Cell,
    pub seed: u64,
    pub index: u64,
    pub variant_a: u8,
    pub variant_b: u8,
    /// B's uniform scale before placement.
    pub b_scale: f64,
    /// B's rigid placement against A (row-major 4×4).
    pub place_b: [[f64; 4]; 4],
    /// The scene's similarity applied to both (row-major 4×4).
    pub scene: [[f64; 4]; 4],
    /// The scene's scale factor (for the expected volumes).
    pub scale: f64,
    /// The near-miss offset in tolerances (touching configurations).
    pub offset: f64,
    /// The rung that drew it.
    pub rung: Rung,
    /// Touching configurations: the outward normals at the contact are
    /// aligned (B inside A there) rather than opposed.
    pub aligned: bool,
    /// Coincident patch: the shift along A's axis as a fraction of A's
    /// extent there (0: a pure symmetry, the surfaces coincide whole).
    pub shift: f64,
    /// The contact point in the final frame (touching configurations).
    pub contact: Option<[f64; 3]>,
    /// A curve contact's samples in the final frame (K).
    pub curve: Vec<[f64; 3]>,
    /// Whether the curve closes (a circle).
    pub curve_closed: bool,
    /// The third body (X3).
    pub third: Option<super::three::Third>,
}

fn rotation_onto(from: Vec3, to: Vec3) -> Transform3 {
    let axis = from.cross(&to);
    let c = from.dot(&to).clamp(-1.0, 1.0);
    if axis.length() < 1e-12 {
        if c > 0.0 {
            return Transform3::identity();
        }
        let perp = if from.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        let axis = from.cross(&perp);
        return Transform3::rotation_around_axis(
            Point3::new(0.0, 0.0, 0.0),
            axis,
            std::f64::consts::PI,
        )
        .expect("axis");
    }
    Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), axis, c.acos()).expect("axis")
}

fn r3(t: &Transform3) -> [[f64; 3]; 3] {
    let m = t.matrix;
    [
        [m[0][0], m[0][1], m[0][2]],
        [m[1][0], m[1][1], m[1][2]],
        [m[2][0], m[2][1], m[2][2]],
    ]
}

fn transform_mesh(mesh: &popkern_mesh::TriMesh, t: &Transform3) -> popkern_mesh::TriMesh {
    let mut out = mesh.clone();
    for v in &mut out.vertices {
        v.position = t.apply_point(&v.position);
    }
    out
}

/// The scene's operands, placed: `(A, B, A mesh, B mesh)`.
pub struct Placed {
    pub a: BRep,
    pub b: BRep,
    pub a_op: std::sync::Arc<Operand>,
    pub b_op: std::sync::Arc<Operand>,
    pub a_mesh: popkern_mesh::TriMesh,
    pub b_mesh: popkern_mesh::TriMesh,
}

fn matrix(t: &Transform3) -> [[f64; 4]; 4] {
    t.matrix
}

fn from_matrix(m: [[f64; 4]; 4]) -> Transform3 {
    Transform3 { matrix: m }
}

/// Build the placed operands of a descriptor.
pub fn place(scene: &Scene) -> Placed {
    let a_op = operand(scene.cell.a, scene.variant_a);
    let b_op = operand(scene.cell.b, scene.variant_b);
    let s = from_matrix(scene.scene);
    let pb = from_matrix(scene.place_b).compose(&Transform3::scale_uniform(scene.b_scale));
    let b_total = s.compose(&pb);
    let a = a_op.brep.transform(&s).expect("place A");
    let b = b_op.brep.transform(&b_total).expect("place B");
    Placed {
        a_mesh: transform_mesh(&a_op.mesh, &s),
        b_mesh: transform_mesh(&b_op.mesh, &b_total),
        a,
        b,
        a_op,
        b_op,
    }
}

/// Why a draw was rejected (counted per cell, design §4).
#[derive(Debug, Clone, PartialEq)]
pub enum Reject {
    NoSample,
    Classification,
    SecondApproach(f64),
    NoApproach(f64),
    NotOverlapping,
    NotContained,
    NotDisjoint,
    /// A touch whose B lies on the wrong side of A away from the contact
    /// (an external touch inside A, an internal one outside it).
    WrongSide,
}

/// Draw scene `index` of `cell` under `seed`: `Err` when the draw is
/// rejected by its own checks (design §6.1).
pub fn draw(cell: Cell, seed: u64, index: u64, rung: Rung) -> Result<Scene, Reject> {
    let (offset, scale) = (rung.offset, rung.scale);
    let tol = Tolerance::default().linear();
    let mut rng = Rng::new(seed, stream_id(&cell.name()), index);
    let (va, mut vb) = (rng.below(3) as u8, rng.below(3) as u8);
    if cell.config == Config::O {
        if cell.a != cell.b {
            return Err(Reject::Classification);
        }
        vb = va;
    }
    let (a_op, b_op) = (operand(cell.a, va), operand(cell.b, vb));
    // The scene similarity: a random turn, a translation, the scale.
    let axis = rng.direction();
    let turn = Transform3::rotation_around_axis(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(axis[0], axis[1], axis[2]),
        rng.range(0.0, std::f64::consts::TAU),
    )
    .expect("axis");
    let shift = Transform3::translation(
        rng.range(-3.0, 3.0),
        rng.range(-3.0, 3.0),
        rng.range(-3.0, 3.0),
    );
    let sim = Transform3::scale_uniform(scale).compose(&shift.compose(&turn));
    // Offsets act in the final frame: a δ of `offset` tolerances.
    let delta = offset * tol / scale;
    let mut b_scale = 1.0;
    let mut aligned = false;
    let mut shift = 0.0;
    let mut curve: Vec<Point3> = Vec::new();
    let mut curve_closed = false;
    let (place_b, contact) = match cell.config {
        Config::E | Config::V | Config::I => {
            let pa = op_point(&a_op, rung.steer_a, &mut rng).ok_or(Reject::NoSample)?;
            if cell.config == Config::I {
                // B fits inside A at the contact: shrink it.
                b_scale = rng.range(0.15, 0.35) * a_op.size / b_op.size;
            }
            let pb = op_point(&b_op, rung.steer_b, &mut rng).ok_or(Reject::NoSample)?;
            // V crosses at its touch: from a saddle against anything, or
            // with normals aligned against mixed curvature (Viviani), so
            // both orientations are drawn.
            let internal = match cell.config {
                Config::I => true,
                Config::V => rng.unit() < 0.5,
                _ => false,
            };
            aligned = internal;
            let target = if internal { pa.n } else { pa.n * -1.0 };
            let r1 = rotation_onto(pb.n, target);
            let spin = Transform3::rotation_around_axis(
                Point3::new(0.0, 0.0, 0.0),
                pa.n,
                rng.range(0.0, std::f64::consts::TAU),
            )
            .expect("axis");
            let rot = spin.compose(&r1);
            let q = rot.apply_point(&Point3::new(
                pb.p.x * b_scale,
                pb.p.y * b_scale,
                pb.p.z * b_scale,
            ));
            let along = if internal { -delta } else { delta };
            let goal = pa.p + pa.n * along;
            let place =
                Transform3::translation(goal.x - q.x, goal.y - q.y, goal.z - q.z).compose(&rot);
            // Relative curvature in A's frame (design §4).
            let tb = rotate_tensor(&r3(&rot), &pb.t);
            let tb = tb.map(|row| row.map(|x| x / b_scale));
            let mut d = [[0.0; 3]; 3];
            for i in 0..3 {
                for j in 0..3 {
                    d[i][j] = if internal {
                        pa.t[i][j] - tb[i][j]
                    } else {
                        -(pa.t[i][j] + tb[i][j])
                    };
                }
            }
            let (lo, hi) = tangent_eigen(&d, pa.n);
            let scale_k = lo.abs().max(hi.abs()).max(1e-12);
            let floor = 0.05 * scale_k + 1e-3 / a_op.size;
            let definite = lo > floor;
            let indefinite = lo < -floor && hi > floor;
            let ok = match cell.config {
                Config::V => indefinite,
                _ => definite,
            };
            if !ok {
                return Err(Reject::Classification);
            }
            (place, Some(sim.apply_point(&pa.p)))
        }
        Config::X | Config::X3 => {
            let pa = op_point(&a_op, rung.steer_a, &mut rng).ok_or(Reject::NoSample)?;
            let axis = rng.direction();
            let rot = Transform3::rotation_around_axis(
                Point3::new(0.0, 0.0, 0.0),
                Vec3::new(axis[0], axis[1], axis[2]),
                rng.range(0.0, std::f64::consts::TAU),
            )
            .expect("axis");
            let depth = rng.range(-0.3, 0.3) * b_op.size;
            let goal = pa.p + pa.n * depth;
            (
                Transform3::translation(goal.x, goal.y, goal.z).compose(&rot),
                None,
            )
        }
        Config::K => {
            let t = super::k_templates::place(&a_op, &b_op, delta, &mut rng)
                .ok_or(Reject::Classification)?;
            b_scale = t.b_scale;
            aligned = t.aligned;
            curve_closed = t.curve.len() == 16;
            curve = t.curve.iter().map(|p| sim.apply_point(p)).collect();
            let c = curve[curve.len() / 2];
            (t.place, Some(c))
        }
        Config::O => {
            use super::operand::Kind as K;
            let z = Vec3::new(0.0, 0.0, 1.0);
            let origin = Point3::new(0.0, 0.0, 0.0);
            let turn_z =
                |angle: f64| Transform3::rotation_around_axis(origin, z, angle).expect("axis");
            let (lo, hi) = super::distance::mesh_box(&a_op.mesh);
            match cell.a {
                K::P | K::W | K::L | K::G => {
                    // Shift along x: the faces square to it stay coplanar.
                    shift = rng.range(0.2, 0.8);
                    (
                        Transform3::translation(shift * (hi.x - lo.x), 0.0, 0.0),
                        None,
                    )
                }
                K::Y => {
                    shift = rng.range(0.2, 0.8);
                    (
                        Transform3::translation(0.0, 0.0, shift * (hi.z - lo.z)),
                        None,
                    )
                }
                K::S => {
                    let d = rng.direction();
                    (
                        Transform3::rotation_around_axis(
                            origin,
                            Vec3::new(d[0], d[1], d[2]),
                            rng.range(0.1, std::f64::consts::TAU - 0.1),
                        )
                        .expect("axis"),
                        None,
                    )
                }
                // An ellipsoid whose x and y scales differ keeps its shape
                // only under a half turn.
                K::B if va == 1 => (turn_z(std::f64::consts::PI), None),
                _ => (turn_z(rng.range(0.1, std::f64::consts::TAU - 0.1)), None),
            }
        }
        Config::D => {
            let dir = rng.direction();
            let far = 1.5 * (a_op.size + b_op.size);
            (
                Transform3::translation(dir[0] * far, dir[1] * far, dir[2] * far),
                None,
            )
        }
        Config::N => {
            b_scale = rng.range(0.1, 0.25) * a_op.size / b_op.size;
            let (lo, hi) = super::distance::mesh_box(&a_op.mesh);
            let mut found = None;
            for _ in 0..64 {
                let p = Point3::new(
                    rng.range(lo.x, hi.x),
                    rng.range(lo.y, hi.y),
                    rng.range(lo.z, hi.z),
                );
                if inside(&a_op.mesh, p) {
                    found = Some(p);
                    break;
                }
            }
            let p = found.ok_or(Reject::NoSample)?;
            let axis = rng.direction();
            let rot = Transform3::rotation_around_axis(
                Point3::new(0.0, 0.0, 0.0),
                Vec3::new(axis[0], axis[1], axis[2]),
                rng.range(0.0, std::f64::consts::TAU),
            )
            .expect("axis");
            (Transform3::translation(p.x, p.y, p.z).compose(&rot), None)
        }
    };
    let scene = Scene {
        cell,
        seed,
        index,
        variant_a: va,
        variant_b: vb,
        b_scale,
        place_b: matrix(&place_b),
        scene: matrix(&sim),
        scale,
        offset,
        rung,
        aligned,
        shift,
        contact: contact.map(|p| [p.x, p.y, p.z]),
        curve: curve.iter().map(|p| [p.x, p.y, p.z]).collect(),
        curve_closed,
        third: None,
    };
    let scene = if cell.config == Config::X3 {
        let placed = place(&scene);
        // C over a point inside both, in A's frame (before the similarity).
        let inv = sim.inverse().expect("similarity");
        let mut a_mesh = placed.a_mesh.clone();
        let mut b_mesh = placed.b_mesh.clone();
        for v in a_mesh.vertices.iter_mut().chain(b_mesh.vertices.iter_mut()) {
            v.position = inv.apply_point(&v.position);
        }
        let third = super::three::place_third(&a_mesh, &b_mesh, cell.a, va, &mut rng);
        Scene {
            third: Some(third.ok_or(Reject::NotOverlapping)?),
            ..scene
        }
    } else {
        scene
    };
    super::scene_check::check(&scene)?;
    Ok(scene)
}
