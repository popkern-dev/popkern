//! Bent pairings for the triangles a face's affine and pair bounds leave
//! over the tolerance (PK-908; the bound itself is
//! `bent_triangle_chord_bound`, see `chord_bend` in popkern-core).
//!
//! A bent triangle pairs `L(y)` with `S(φ_T(y))`, `φ_T` the quadratic
//! Bézier map of its parameter triangle whose edge controls are `c_e`
//! (the midpoint for a straight edge). The certificate's two
//! directions stay covered:
//!
//! - Mesh to surface: every `L(y)` has the partner `S(φ_T(y))`, which
//!   lies on the face when `φ_T(T)` stays inside the face's parameter
//!   triangles. `φ_T(T)` lies in the hull of its control net, which is
//!   checked to lie in `T` and the triangles across its bent edges.
//! - Surface to mesh: an edge is bent only when both its triangles bend
//!   it with the same control, and is straight otherwise, so the
//!   piecewise map `Φ` (`φ_T` on bent triangles, the identity elsewhere)
//!   is continuous and the identity on the domain's boundary. It is
//!   then onto (degree one), and every surface point `S(Φ(y))` is
//!   paired with its triangle's `L(y)`; on an unbent triangle that is
//!   the triangle's own affine or pair bound.
//!
//! Any control gives a valid bound; the one chosen takes each edge's
//! midpoint to the surface point nearest the chord's midpoint, which
//! cancels the in-plane sag that the affine pairing charges.

use std::collections::HashMap;

use popkern_core::geometry::{Point3, SurfaceEval, SurfaceGeometry};

type Uv = (f64, f64);

/// One plain (single affine piece) triangle of the face.
#[derive(Debug, Clone, Copy)]
pub(super) struct BendTriangle {
    /// Welded vertex ids (equal for coincident vertices, so a seam's
    /// two sides share their edges).
    pub weld: [usize; 3],
    /// Unwrapped parameters.
    pub w: [Uv; 3],
    /// Emitted positions.
    pub p: [Point3; 3],
}

/// Rounds of rejection before the pass gives up on bending at all.
const MAX_ROUNDS: usize = 8;

/// Lower `bounds[k]` for the triangles over `tolerance` that a bent
/// pairing bounds more tightly. `bound(w, p, controls)` is the bent
/// bound, `None` where it does not apply or declines (the triangle then
/// keeps its unbent bound); `None` entries of `tris` (pole pieces)
/// never bend.
pub(super) fn bend_pass(
    surface: &SurfaceGeometry,
    tris: &[Option<BendTriangle>],
    bounds: &mut [f64],
    tolerance: f64,
    mut bound: impl FnMut([Uv; 3], [Point3; 3], [Option<Uv>; 3]) -> Option<f64>,
) {
    let mut users: HashMap<(usize, usize), Vec<(usize, usize)>> = HashMap::new();
    for (k, t) in tris.iter().enumerate() {
        let Some(t) = t else { continue };
        for e in 0..3 {
            users.entry(key(t, e)).or_default().push((k, e));
        }
    }
    // The midpoint displacement of every edge between two candidates.
    let mut bent: Vec<bool> = tris
        .iter()
        .zip(bounds.iter())
        .map(|(t, &b)| t.is_some() && b > tolerance)
        .collect();
    let mut offsets: HashMap<(usize, usize), Uv> = HashMap::new();
    for _ in 0..MAX_ROUNDS {
        let mut found = vec![None; tris.len()];
        let mut rejected = false;
        for k in 0..tris.len() {
            let Some(t) = tris[k].filter(|_| bent[k]) else {
                continue;
            };
            let neighbour = |e: usize| -> Option<usize> {
                match users.get(&key(&t, e)).map(Vec::as_slice) {
                    Some([(a, _), (b, _)]) => Some(if *a == k { *b } else { *a }),
                    _ => None,
                }
            };
            let mut controls = [None; 3];
            let mut across = [None; 3];
            for (e, c) in controls.iter_mut().enumerate() {
                let (a, b) = (t.w[e], t.w[(e + 1) % 3]);
                let m = (0.5 * (a.0 + b.0), 0.5 * (a.1 + b.1));
                let Some(j) = neighbour(e).filter(|&j| bent[j]) else {
                    continue;
                };
                let Some(tj) = tris[j] else { continue };
                let d = *offsets
                    .entry(key(&t, e))
                    .or_insert_with(|| foot(surface, &t, e).unwrap_or((0.0, 0.0)));
                *c = Some((m.0 + 2.0 * d.0, m.1 + 2.0 * d.1));
                across[e] = Some(aligned(&tj, &t, e));
            }
            let admissible = net_inside(t.w, controls, across);
            match admissible.then(|| bound(t.w, t.p, controls)).flatten() {
                Some(b) if b < bounds[k] => found[k] = Some(b),
                _ => {
                    bent[k] = false;
                    rejected = true;
                }
            }
        }
        if !rejected {
            for (b, f) in bounds.iter_mut().zip(found) {
                if let Some(f) = f {
                    *b = f;
                }
            }
            return;
        }
    }
}

/// The welded edge `e` (vertices `e`, `e + 1`) of `t`, unordered.
fn key(t: &BendTriangle, e: usize) -> (usize, usize) {
    let (a, b) = (t.weld[e], t.weld[(e + 1) % 3]);
    (a.min(b), a.max(b))
}

/// The displacement of edge `e`'s parameter midpoint that takes it to
/// the surface point nearest the chord's midpoint (Gauss–Newton); `None`
/// when it does not settle within a quarter of the edge.
fn foot(surface: &SurfaceGeometry, t: &BendTriangle, e: usize) -> Option<Uv> {
    let (i, j) = (e, (e + 1) % 3);
    let (a, b) = (t.w[i], t.w[j]);
    let m = (0.5 * (a.0 + b.0), 0.5 * (a.1 + b.1));
    let c = Point3::new(
        0.5 * (t.p[i].x + t.p[j].x),
        0.5 * (t.p[i].y + t.p[j].y),
        0.5 * (t.p[i].z + t.p[j].z),
    );
    let mut q = m;
    for _ in 0..8 {
        let r = surface.point_at(q.0, q.1) - c;
        let (su, sv) = (surface.partial_u(q.0, q.1), surface.partial_v(q.0, q.1));
        let (g, h, k) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
        let det = g * k - h * h;
        // NaN fails too.
        if det.abs().partial_cmp(&f64::MIN_POSITIVE) != Some(std::cmp::Ordering::Greater) {
            return None;
        }
        q.0 -= (k * su.dot(&r) - h * sv.dot(&r)) / det;
        q.1 -= (g * sv.dot(&r) - h * su.dot(&r)) / det;
    }
    let d = (q.0 - m.0, q.1 - m.1);
    let reach = 0.25 * ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt();
    (d.0.is_finite() && d.1.is_finite() && (d.0 * d.0 + d.1 * d.1).sqrt() <= reach).then_some(d)
}

/// The parameter triangle of `other` (across edge `e` of `t`) in `t`'s
/// unwrapped chart: shifted so their shared vertex coincides.
fn aligned(other: &BendTriangle, t: &BendTriangle, e: usize) -> [Uv; 3] {
    let shared = t.weld[e];
    let mine = t.w[e];
    let theirs = (0..3)
        .find(|&i| other.weld[i] == shared)
        .map_or(mine, |i| other.w[i]);
    let s = (mine.0 - theirs.0, mine.1 - theirs.1);
    other.w.map(|q| (q.0 + s.0, q.1 + s.1))
}

fn orient(a: Uv, b: Uv, c: Uv) -> f64 {
    let p = |q: Uv| robust::Coord { x: q.0, y: q.1 };
    robust::orient2d(p(a), p(b), p(c))
}

/// Whether `q` lies in the closed triangle `t`.
fn in_triangle(q: Uv, t: [Uv; 3]) -> bool {
    let s = orient(t[0], t[1], t[2]).signum();
    s != 0.0 && (0..3).all(|i| orient(t[i], t[(i + 1) % 3], q) * s >= 0.0)
}

/// Whether the hull of the quadratic net (corners `w`, bent edges'
/// controls `controls`, straight edges' midpoints) lies in `w`'s
/// triangle and those across its bent edges (`across[e]`): each control
/// off the triangle must lie in the one across its edge, and the
/// polygon of the corners and those controls must be convex, so the
/// hull is the triangle with an ear on each such edge.
fn net_inside(w: [Uv; 3], controls: [Option<Uv>; 3], across: [Option<[Uv; 3]>; 3]) -> bool {
    let s = orient(w[0], w[1], w[2]).signum();
    if s == 0.0 {
        return false;
    }
    let mut polygon: Vec<Uv> = Vec::with_capacity(6);
    for e in 0..3 {
        polygon.push(w[e]);
        let Some(c) = controls[e].filter(|&c| !in_triangle(c, w)) else {
            continue;
        };
        match across[e] {
            Some(t) if in_triangle(c, t) => polygon.push(c),
            _ => return false,
        }
    }
    let n = polygon.len();
    (0..n).all(|i| orient(polygon[i], polygon[(i + 1) % n], polygon[(i + 2) % n]) * s >= 0.0)
}

#[cfg(test)]
#[path = "chordal_bend_tests.rs"]
mod tests;
