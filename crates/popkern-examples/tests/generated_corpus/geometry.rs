//! The generator's own geometry: surface samples with outward normals and
//! second fundamental forms, and tessellation distances. Independent of
//! the kernel's classifiers: it reads surface evaluation and meshes only.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::BRep;
use popkern_mesh::TriMesh;

pub use super::distance::{
    Grid, closest_approach, closest_approach_except, inside, polyline_distance,
};
use super::rng::Rng;

/// The tessellation's signed divergence volume.
pub fn mesh_signed_volume(mesh: &TriMesh) -> f64 {
    let p = |i: usize| mesh.vertices[i].position;
    mesh.triangles
        .iter()
        .map(|t| {
            let (a, b, c) = (p(t[0]), p(t[1]), p(t[2]));
            (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                + a.z * (b.x * c.y - b.y * c.x))
                / 6.0
        })
        .sum()
}

/// A point on a surface with its out-of-material normal and its second
/// fundamental form as a symmetric tensor on the tangent plane (`II(x,
/// x) = xᵀ T x` for tangent `x`, taken against the outward normal).
#[derive(Debug, Clone, Copy)]
pub struct SurfacePoint {
    pub p: Point3,
    pub n: Vec3,
    pub t: [[f64; 3]; 3],
    /// The face's chart parameters.
    pub uv: (f64, f64),
}

fn v(p: Vec3) -> [f64; 3] {
    [p.x, p.y, p.z]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// The outward normal and the II tensor at `(u, v)` given an outward hint.
/// Near a singular chart row the derivatives are read a little inside,
/// toward `(u_in, v_in)`; the point itself is evaluated exactly.
pub fn surface_point(
    s: &SurfaceGeometry,
    (u, w): (f64, f64),
    toward: (f64, f64),
    hint: Vec3,
) -> Option<SurfacePoint> {
    let p = s.point_at(u, w);
    let mut at = (u, w);
    let mut frame = None;
    for k in 0..8 {
        let (su, sv) = s.first_partials(at.0, at.1);
        let n = su.cross(&sv);
        if n.length() > 1e-10 * (su.length() * sv.length()).max(1e-30) && su.length() > 1e-12 {
            frame = Some((su, sv, n.normalized()?));
            break;
        }
        let f = 1e-7 * 10f64.powi(k);
        at = (u + (toward.0 - u) * f, w + (toward.1 - w) * f);
    }
    let (su, sv, mut n) = frame?;
    if n.dot(&hint) < 0.0 {
        n = n * -1.0;
    }
    let (suu, suv, svv) = (
        s.second_partial_uu(at.0, at.1),
        s.second_partial_uv(at.0, at.1),
        s.second_partial_vv(at.0, at.1),
    );
    let (e, f, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
    let (l, m, nn) = (suu.dot(&n), suv.dot(&n), svv.dot(&n));
    let det = e * g - f * f;
    if det.abs() < 1e-30 {
        return None;
    }
    // An orthonormal tangent frame and each frame vector's chart coords.
    let e1 = su.normalized()?;
    let e2 = n.cross(&e1);
    let coords = |d: Vec3| {
        let (r1, r2) = (d.dot(&su), d.dot(&sv));
        ((g * r1 - f * r2) / det, (e * r2 - f * r1) / det)
    };
    let (a1, a2) = (coords(e1), coords(e2));
    let ii =
        |x: (f64, f64), y: (f64, f64)| l * x.0 * y.0 + m * (x.0 * y.1 + x.1 * y.0) + nn * x.1 * y.1;
    let (k11, k12, k22) = (ii(a1, a1), ii(a1, a2), ii(a2, a2));
    let (b1, b2) = (v(e1), v(e2));
    let mut t = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            t[i][j] =
                k11 * b1[i] * b1[j] + k12 * (b1[i] * b2[j] + b2[i] * b1[j]) + k22 * b2[i] * b2[j];
        }
    }
    Some(SurfacePoint {
        p,
        n,
        t,
        uv: (u, w),
    })
}

/// A random point inside a relevant face: a random relevant triangle,
/// a random point of it, its chart parameters interpolated and the
/// surface evaluated there (kept only when it lies on the triangle within
/// the chord bound, so it is inside the face's trim).
pub fn sample_point(
    brep: &BRep,
    mesh: &TriMesh,
    relevant: &[usize],
    chord: f64,
    rng: &mut Rng,
) -> Option<SurfacePoint> {
    sample_point_weighted(brep, mesh, relevant, None, chord, rng)
}

/// A triangle of `relevant` drawn by area when `cumulative` is given.
fn pick(relevant: &[usize], cumulative: Option<&[f64]>, rng: &mut Rng) -> usize {
    match cumulative {
        Some(c) if c.len() == relevant.len() && c.last().is_some_and(|t| *t > 0.0) => {
            let x = rng.unit() * c[c.len() - 1];
            relevant[c.partition_point(|&a| a < x).min(relevant.len() - 1)]
        }
        _ => relevant[rng.below(relevant.len())],
    }
}

/// [`sample_point`] with area-weighted triangle choice.
pub fn sample_point_weighted(
    brep: &BRep,
    mesh: &TriMesh,
    relevant: &[usize],
    cumulative: Option<&[f64]>,
    chord: f64,
    rng: &mut Rng,
) -> Option<SurfacePoint> {
    for _ in 0..32 {
        let ti = pick(relevant, cumulative, rng);
        let tri = mesh.triangles[ti];
        let face = brep.faces.get(mesh.face_map[ti])?;
        let surface = &face.surface.as_ref()?.geometry;
        let (mut a, mut b) = (rng.unit(), rng.unit());
        if a + b > 1.0 {
            (a, b) = (1.0 - a, 1.0 - b);
        }
        let c = 1.0 - a - b;
        let vx = |i: usize| &mesh.vertices[tri[i]];
        let uv = (
            a * vx(0).uv.0 + b * vx(1).uv.0 + c * vx(2).uv.0,
            a * vx(0).uv.1 + b * vx(1).uv.1 + c * vx(2).uv.1,
        );
        let centre = (
            (vx(0).uv.0 + vx(1).uv.0 + vx(2).uv.0) / 3.0,
            (vx(0).uv.1 + vx(1).uv.1 + vx(2).uv.1) / 3.0,
        );
        let (p0, p1, p2) = (vx(0).position, vx(1).position, vx(2).position);
        let hint = (p1 - p0).cross(&(p2 - p0));
        let Some(sp) = surface_point(surface, uv, centre, hint) else {
            continue;
        };
        let q = Point3::new(
            a * p0.x + b * p1.x + c * p2.x,
            a * p0.y + b * p1.y + c * p2.y,
            a * p0.z + b * p1.z + c * p2.z,
        );
        if sp.p.distance_to(&q) <= 4.0 * chord + 1e-9 {
            return Some(sp);
        }
    }
    None
}

/// The 2×2 restriction of a tangent tensor to the plane with normal `n`,
/// and its eigenvalues (ascending).
pub fn tangent_eigen(t: &[[f64; 3]; 3], n: Vec3) -> (f64, f64) {
    let a = if n.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let e1 = n.cross(&a).normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let e2 = n.cross(&e1);
    let q = |x: Vec3, y: Vec3| {
        let (x, y) = (v(x), v(y));
        let mut s = 0.0;
        for i in 0..3 {
            for j in 0..3 {
                s += x[i] * t[i][j] * y[j];
            }
        }
        s
    };
    let (m11, m12, m22) = (q(e1, e1), q(e1, e2), q(e2, e2));
    let mean = 0.5 * (m11 + m22);
    let disc = (0.25 * (m11 - m22).powi(2) + m12 * m12).sqrt();
    (mean - disc, mean + disc)
}

/// `R T Rᵀ` for a rotation given as a 3×3 matrix.
pub fn rotate_tensor(r: &[[f64; 3]; 3], t: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut rt = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            rt[i][j] = (0..3).map(|k| r[i][k] * t[k][j]).sum();
        }
    }
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = (0..3).map(|k| rt[i][k] * r[j][k]).sum();
        }
    }
    out
}

/// `dot3` re-exported for the placement code.
pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    dot3(a, b)
}

/// Where a ladder rung puts the contact on an operand's chart.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum Steer {
    /// Anywhere (random).
    Free,
    /// At `d` (in space) from the chart's seam.
    Seam(f64),
    /// At `d` (in space) from a singular chart row (a pole).
    Pole(f64),
}

/// The surface's seam parameter (`u` at its range start) when the chart
/// closes along u, and its singular rows (v values whose row maps to one
/// point).
pub fn chart_events(s: &SurfaceGeometry, size: f64) -> (Option<f64>, Vec<f64>) {
    let ((u0, u1), (v0, v1)) = s.parameter_range();
    if !(u0.is_finite() && u1.is_finite()) {
        return (None, Vec::new());
    }
    // A cylinder's v range is unbounded: the seam needs only u's.
    let vm = if v0.is_finite() && v1.is_finite() {
        0.5 * (v0 + v1)
    } else {
        0.0
    };
    let closes = s.parameter_periodicity().0.is_some()
        || s.point_at(u0, vm).distance_to(&s.point_at(u1, vm)) <= 1e-9 * size;
    let poles = [v0, v1]
        .into_iter()
        .filter(|v| v.is_finite())
        .filter(|&v| {
            let (a, b) = (
                s.point_at(u0 + 0.13 * (u1 - u0), v),
                s.point_at(u0 + 0.61 * (u1 - u0), v),
            );
            a.distance_to(&b) <= 1e-9 * size
        })
        .collect();
    (closes.then_some(u0), poles)
}

/// A point on a relevant face placed by `steer`: the seam or pole rung
/// moves the parameter off the event by `d` in space (through the chart
/// speed), the other parameter random; kept only when it lies on the
/// face (within the chord of the face's own triangles).
pub fn steered_point(
    brep: &BRep,
    mesh: &TriMesh,
    relevant: &[usize],
    chord: f64,
    size: f64,
    steer: Steer,
    rng: &mut Rng,
) -> Option<SurfacePoint> {
    steered_point_weighted(brep, mesh, relevant, None, chord, size, steer, rng)
}

/// A point on an operand's relevant faces placed by `steer`, triangles
/// drawn by area.
pub fn op_point(op: &super::operand::Operand, steer: Steer, rng: &mut Rng) -> Option<SurfacePoint> {
    steered_point_weighted(
        &op.brep,
        &op.mesh,
        &op.relevant,
        Some(&op.relevant_area),
        op.chord,
        op.size,
        steer,
        rng,
    )
}

#[allow(clippy::too_many_arguments)]
fn steered_point_weighted(
    brep: &BRep,
    mesh: &TriMesh,
    relevant: &[usize],
    cumulative: Option<&[f64]>,
    chord: f64,
    size: f64,
    steer: Steer,
    rng: &mut Rng,
) -> Option<SurfacePoint> {
    if steer == Steer::Free {
        return sample_point_weighted(brep, mesh, relevant, cumulative, chord, rng);
    }
    for _ in 0..64 {
        let ti = pick(relevant, cumulative, rng);
        let face = brep.faces.get(mesh.face_map[ti])?;
        let surface = &face.surface.as_ref()?.geometry;
        let (seam, poles) = chart_events(surface, size);
        // The face's own parameter box, from its triangles.
        let mut lo = (f64::INFINITY, f64::INFINITY);
        let mut hi = (f64::NEG_INFINITY, f64::NEG_INFINITY);
        for (k, f) in mesh.face_map.iter().enumerate() {
            if *f == mesh.face_map[ti] {
                for &i in &mesh.triangles[k] {
                    let uv = mesh.vertices[i].uv;
                    lo = (lo.0.min(uv.0), lo.1.min(uv.1));
                    hi = (hi.0.max(uv.0), hi.1.max(uv.1));
                }
            }
        }
        let (u, v, toward) = match steer {
            Steer::Seam(d) => {
                let u0 = seam?;
                let v = rng.range(lo.1, hi.1);
                let (su, _) = surface.first_partials(u0, v);
                let du = d / su.length().max(1e-300);
                let u = if rng.unit() < 0.5 {
                    u0 + du
                } else {
                    surface.parameter_range().0.1 - du
                };
                (u, v, (u + 1e-3, v))
            }
            Steer::Pole(d) => {
                let &vp = poles.get(rng.below(poles.len().max(1)))?;
                let u = rng.range(lo.0, hi.0);
                // Space distance from the pole along the meridian.
                let ((_, _), (v0, v1)) = surface.parameter_range();
                let inward = if (vp - v0).abs() < (vp - v1).abs() {
                    1.0
                } else {
                    -1.0
                };
                let pole = surface.point_at(u, vp);
                let (mut a, mut b) = (0.0f64, 0.5 * (v1 - v0));
                for _ in 0..80 {
                    let m = 0.5 * (a + b);
                    if surface.point_at(u, vp + inward * m).distance_to(&pole) < d {
                        a = m;
                    } else {
                        b = m;
                    }
                }
                let v = vp + inward * 0.5 * (a + b);
                (u, v, (u, vp + inward * 0.5 * (v1 - v0)))
            }
            Steer::Free => unreachable!(),
        };
        let hint = {
            let t = mesh.triangles[ti];
            let (p0, p1, p2) = (
                mesh.vertices[t[0]].position,
                mesh.vertices[t[1]].position,
                mesh.vertices[t[2]].position,
            );
            (p1 - p0).cross(&(p2 - p0))
        };
        let Some(sp) = surface_point(surface, (u, v), toward, hint) else {
            continue;
        };
        // On this face: within a few chords of one of its triangles (a
        // vertex test fails a cylinder's strips, whose vertices are only at
        // the rims).
        let on_face = mesh.face_map.iter().enumerate().any(|(k, f)| {
            *f == mesh.face_map[ti] && {
                let t = mesh.triangles[k];
                let q = |i: usize| mesh.vertices[t[i]].position;
                super::distance::point_triangle_distance(sp.p, q(0), q(1), q(2))
                    <= 4.0 * chord + 1e-12
            }
        });
        if on_face {
            return Some(sp);
        }
    }
    None
}
