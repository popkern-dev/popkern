//! A cheap, sound-positive crossing witness for bi-nonic cap rings
//! (PK-547).
//!
//! The full gate ([`super::binonic_ring_crossing`]) clears clean rings in
//! microseconds through its separation certificate. A ring that really
//! crosses falls through to the surface-surface solver instead, which
//! took 27 s on the pyramid3_uniform welded-retry ring. That is too slow
//! to run per optimizer step, or to spend once per rejected candidate.
//!
//! The witness finds crossings cheaply and proves each one:
//! 1. Tessellate every patch on a coarse grid, and intersect the
//!    triangle edges of one patch with the triangles of the other
//!    (Möller–Trumbore, AABB-culled).
//! 2. Refine each hit that is not exempt by Gauss–Newton on
//!    `S_a(u, v) − S_b(s, t) = 0` to a true intersection point, within
//!    the gate's tolerance.
//! 3. Report it only when that point lies beyond the pair's exemption
//!    depth, under the same exemption rules as the gate.
//!
//! A witness is therefore a real crossing the gate would also decline.
//! `None` proves nothing: the coarse mesh can miss a crossing, and the
//! full gate stays the acceptance authority.

use super::super::RectangularGregoryBiNonic;
use super::ring_pair_crossing_certificate::distance_to_spoke;
use super::ring_pair_crossing_setup::{PairSetup, pair_setup, patch_aabbs, spoke_polylines};
use super::{BinonicRingCrossing, EXEMPTION_BAND_FACTOR};
use crate::geometry::{Point3, Vec3};

/// Grid points per patch side of the witness tessellation.
const WITNESS_GRID: usize = 13;
/// Gauss–Newton iterations of the witness refinement.
const REFINE_ITERATIONS: usize = 12;

/// The witness for a ring: the first proven crossing in scan order.
pub(super) fn ring_crossing_witness(
    ring: &[RectangularGregoryBiNonic],
    tolerance: f64,
) -> Option<BinonicRingCrossing> {
    let n = ring.len();
    if n < 2 {
        return None;
    }
    let band = EXEMPTION_BAND_FACTOR * tolerance;
    let aabbs = patch_aabbs(ring);
    let spokes = spoke_polylines(ring);
    let meshes: Vec<Mesh> = ring.iter().map(Mesh::of).collect();
    for a in 0..n {
        for b in (a + 1)..n {
            if let (Some(bb_a), Some(bb_b)) = (&aabbs[a], &aabbs[b])
                && !bb_a.overlaps(bb_b, band)
            {
                continue;
            }
            let setup = pair_setup(ring, &aabbs, &spokes, (a, b), band, false);
            let witness = pair_witness(ring, &meshes, (a, b), &setup, tolerance)
                .or_else(|| pair_witness(ring, &meshes, (b, a), &setup, tolerance));
            if let Some(mut crossing) = witness {
                // Report in the gate's (a < b) pair order.
                if crossing.patch_a != a {
                    std::mem::swap(&mut crossing.uv_a, &mut crossing.uv_b);
                    crossing.patch_a = a;
                    crossing.patch_b = b;
                }
                return Some(crossing);
            }
        }
    }
    None
}

/// A patch sampled on the witness grid.
struct Mesh {
    points: Vec<Point3>,
}

impl Mesh {
    fn of(patch: &RectangularGregoryBiNonic) -> Self {
        let g = WITNESS_GRID;
        let mut points = Vec::with_capacity(g * g);
        for i in 0..g {
            for j in 0..g {
                points.push(patch.point_at(param(i), param(j)));
            }
        }
        Self { points }
    }

    fn at(&self, i: usize, j: usize) -> Point3 {
        self.points[i * WITNESS_GRID + j]
    }

    /// The two triangles of every grid quad, with their `(u, v)` corners.
    fn triangles(&self) -> Vec<[(Point3, (f64, f64)); 3]> {
        let g = WITNESS_GRID;
        let corner = |i: usize, j: usize| (self.at(i, j), (param(i), param(j)));
        let mut out = Vec::with_capacity(2 * (g - 1) * (g - 1));
        for i in 0..g - 1 {
            for j in 0..g - 1 {
                out.push([corner(i, j), corner(i + 1, j), corner(i + 1, j + 1)]);
                out.push([corner(i, j), corner(i + 1, j + 1), corner(i, j + 1)]);
            }
        }
        out
    }

    /// Every triangle edge once: grid rows, columns and quad diagonals.
    fn segments(&self) -> Vec<[(Point3, (f64, f64)); 2]> {
        let g = WITNESS_GRID;
        let corner = |i: usize, j: usize| (self.at(i, j), (param(i), param(j)));
        let mut out = Vec::with_capacity(3 * g * g);
        for i in 0..g {
            for j in 0..g {
                if i + 1 < g {
                    out.push([corner(i, j), corner(i + 1, j)]);
                }
                if j + 1 < g {
                    out.push([corner(i, j), corner(i, j + 1)]);
                }
                if i + 1 < g && j + 1 < g {
                    out.push([corner(i, j), corner(i + 1, j + 1)]);
                }
            }
        }
        out
    }
}

fn param(i: usize) -> f64 {
    i as f64 / (WITNESS_GRID - 1) as f64
}

/// Segments of patch `pair.0` against triangles of patch `pair.1`.
fn pair_witness(
    ring: &[RectangularGregoryBiNonic],
    meshes: &[Mesh],
    (a, b): (usize, usize),
    setup: &PairSetup<'_>,
    tolerance: f64,
) -> Option<BinonicRingCrossing> {
    let triangles = meshes[b].triangles();
    let boxes: Vec<(Point3, Point3)> = triangles
        .iter()
        .map(|t| bounds(&[t[0].0, t[1].0, t[2].0]))
        .collect();
    for segment in meshes[a].segments() {
        let (lo, hi) = bounds(&[segment[0].0, segment[1].0]);
        for (triangle, (tlo, thi)) in triangles.iter().zip(boxes.iter()) {
            if hi.x < tlo.x || thi.x < lo.x || hi.y < tlo.y || thi.y < lo.y {
                continue;
            }
            if hi.z < tlo.z || thi.z < lo.z {
                continue;
            }
            let Some((s, beta, gamma)) = segment_triangle(segment[0].0, segment[1].0, triangle)
            else {
                continue;
            };
            let lerp = |p: (f64, f64), q: (f64, f64), f: f64| {
                (p.0 + f * (q.0 - p.0), p.1 + f * (q.1 - p.1))
            };
            let uv_a = lerp(segment[0].1, segment[1].1, s);
            let (t0, t1, t2) = (triangle[0].1, triangle[1].1, triangle[2].1);
            let alpha = 1.0 - beta - gamma;
            let uv_b = (
                alpha * t0.0 + beta * t1.0 + gamma * t2.0,
                alpha * t0.1 + beta * t1.1 + gamma * t2.1,
            );
            if let Some(crossing) = refine(ring, (a, b), uv_a, uv_b, setup, tolerance) {
                return Some(crossing);
            }
        }
    }
    None
}

/// Gauss–Newton on `F = S_a(u, v) − S_b(s, t)` (3 equations, 4
/// unknowns, minimum-norm steps), clamped to the unit square. A
/// converged point beyond the exemption depth is a proven crossing.
fn refine(
    ring: &[RectangularGregoryBiNonic],
    (a, b): (usize, usize),
    uv_a: (f64, f64),
    uv_b: (f64, f64),
    setup: &PairSetup<'_>,
    tolerance: f64,
) -> Option<BinonicRingCrossing> {
    let (pa, pb) = (&ring[a], &ring[b]);
    let mut x = [uv_a.0, uv_a.1, uv_b.0, uv_b.1];
    for _ in 0..REFINE_ITERATIONS {
        let f = pa.point_at(x[0], x[1]) - pb.point_at(x[2], x[3]);
        if f.length() <= tolerance {
            break;
        }
        let columns = [
            pa.partial_u(x[0], x[1]),
            pa.partial_v(x[0], x[1]),
            pb.partial_u(x[2], x[3]) * -1.0,
            pb.partial_v(x[2], x[3]) * -1.0,
        ];
        // Minimum-norm step: dx = Jᵀ (J Jᵀ)⁻¹ (−F).
        let jjt = gram(&columns);
        let y = solve3(jjt, f * -1.0)?;
        for (k, column) in columns.iter().enumerate() {
            x[k] = (x[k] + column.dot(&y)).clamp(0.0, 1.0);
        }
    }
    let point = pa.point_at(x[0], x[1]);
    if (point - pb.point_at(x[2], x[3])).length() > tolerance {
        return None;
    }
    let mut separation = point.distance_to(&setup.apex);
    for spoke in &setup.shared_spokes {
        separation = separation.min(distance_to_spoke(point, spoke.polyline, spoke.patch));
    }
    if separation <= setup.exempt_depth {
        return None;
    }
    Some(BinonicRingCrossing {
        patch_a: a,
        patch_b: b,
        relation: setup.relation,
        point,
        uv_a: (x[0], x[1]),
        uv_b: (x[2], x[3]),
        separation,
        exemption_band: setup.exempt_depth,
        coincident_region: false,
    })
}

fn bounds(points: &[Point3]) -> (Point3, Point3) {
    let mut lo = points[0];
    let mut hi = points[0];
    for p in &points[1..] {
        lo = Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
        hi = Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
    }
    (lo, hi)
}

/// Möller–Trumbore segment–triangle intersection: `(s, β, γ)` with the
/// hit at `p0 + s (p1 − p0)` and barycentrics `(1 − β − γ, β, γ)`.
fn segment_triangle(
    p0: Point3,
    p1: Point3,
    triangle: &[(Point3, (f64, f64)); 3],
) -> Option<(f64, f64, f64)> {
    let (v0, v1, v2) = (triangle[0].0, triangle[1].0, triangle[2].0);
    let direction = p1 - p0;
    let (e1, e2) = (v1 - v0, v2 - v0);
    let h = direction.cross(&e2);
    let det = e1.dot(&h);
    let scale = e1.length() * e2.length() * direction.length();
    if det.is_nan() || det.abs() <= 1e-12 * scale {
        return None;
    }
    let inv = 1.0 / det;
    let offset = p0 - v0;
    let beta = inv * offset.dot(&h);
    if !(0.0..=1.0).contains(&beta) {
        return None;
    }
    let q = offset.cross(&e1);
    let gamma = inv * direction.dot(&q);
    if gamma < 0.0 || beta + gamma > 1.0 {
        return None;
    }
    let s = inv * e2.dot(&q);
    (0.0..=1.0).contains(&s).then_some((s, beta, gamma))
}

/// `J Jᵀ` for the 3×4 Jacobian given by its four columns.
fn gram(columns: &[Vec3; 4]) -> [[f64; 3]; 3] {
    let mut m = [[0.0; 3]; 3];
    for c in columns {
        let v = [c.x, c.y, c.z];
        for (r, row) in m.iter_mut().enumerate() {
            for (k, cell) in row.iter_mut().enumerate() {
                *cell += v[r] * v[k];
            }
        }
    }
    m
}

/// Solve the 3×3 system `m y = r` by Cramer's rule (`None` if singular).
fn solve3(m: [[f64; 3]; 3], r: Vec3) -> Option<Vec3> {
    let det = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(m);
    let scale = m.iter().flatten().fold(0.0_f64, |acc, v| acc.max(v.abs()));
    if d.is_nan() || d.abs() <= 1e-14 * scale * scale * scale {
        return None;
    }
    let rhs = [r.x, r.y, r.z];
    let mut out = [0.0; 3];
    for (col, slot) in out.iter_mut().enumerate() {
        let mut replaced = m;
        for row in 0..3 {
            replaced[row][col] = rhs[row];
        }
        *slot = det(replaced) / d;
    }
    Some(Vec3::new(out[0], out[1], out[2]))
}
