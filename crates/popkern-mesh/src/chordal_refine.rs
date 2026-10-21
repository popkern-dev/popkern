//! Refined pairings for the triangles every single-piece pairing leaves
//! over the tolerance (PK-910).
//!
//! A two-triangle cap tiling a whole Gregory patch beside an activated
//! corner is within 2.2e-3 of its surface, but no affine, bilinear or
//! single quadratic correspondence over a triangle says so: the nearest
//! point map it would have to follow bends too far (its feet leave the
//! patch). A refined pairing follows it piecewise. Over a region `R` of
//! triangles, each triangle is cut into `n²` sub-triangles; a
//! sub-vertex `y` inside `R` is mapped to the surface parameter
//! nearest `L(y)`, and one on `R`'s boundary to itself. On each
//! sub-triangle `s` the map `Φ` is affine, so its pairing
//! `S(Φ(y)) ↔ L(y)` is the affine bound of the parameter triangle
//! `Φ(s)` against the planar triangle `L(s)`.
//!
//! The pairings cover `R` both ways when every `Φ(s)` keeps `s`'s
//! orientation: `Φ` is then continuous, the identity on `∂R`, and
//! (its linear homotopy to the identity fixing `∂R`) of degree one,
//! so every parameter of `R` has exactly one preimage and none
//! outside it does. Mesh to surface, `S(Φ(y))` lies on the face;
//! surface to mesh, every `S(x)`, `x ∈ R`, is some `S(Φ(y))`.
//!
//! Neighbours outside `R` keep their own pairings, so `R`'s boundary
//! must be an edge they treat as straight: a triangle the bend pass
//! bent toward a neighbour outside `R` leaves it. A region's bounds
//! are replaced together, and only when its worst bound improves.
//! Periodic charts are left to the caller (their triangles' parameters
//! are unwrapped per triangle, so shared edges need not agree).

use std::collections::{HashMap, HashSet};

use popkern_core::geometry::{Point3, SurfaceEval, SurfaceGeometry};

use super::chordal_bend::{BendTriangle, EdgeKey, key, orient};

type Uv = (f64, f64);

/// Subdivision levels tried, coarsest first (powers of two, so a
/// sub-vertex at a corner is the corner exactly). pentagon_tower's
/// caps certify at 4 (4.7e-3..5.0e-3 against a 5e-3 chord height;
/// 5.3e-3..6.4e-3 at 2).
const LEVELS: [usize; 2] = [2, 4];

/// Lower `bounds[k]` of the triangles over `tolerance` by refined
/// pairings. `bent` holds the edges the bend pass left bent;
/// `bound(w, p)` is the affine bound of a parameter triangle against a
/// planar one, `None` where it declines.
pub(super) fn refine_pass(
    surface: &SurfaceGeometry,
    tris: &[Option<BendTriangle>],
    bounds: &mut [f64],
    bent: &HashSet<EdgeKey>,
    tolerance: f64,
    mut bound: impl FnMut([Uv; 3], [Point3; 3]) -> Option<f64>,
) {
    let mut users: HashMap<EdgeKey, Vec<usize>> = HashMap::new();
    for (k, t) in tris.iter().enumerate() {
        let Some(t) = t else { continue };
        for e in 0..3 {
            users.entry(key(t, e)).or_default().push(k);
        }
    }
    let mut member: Vec<bool> = tris
        .iter()
        .zip(bounds.iter())
        .map(|(t, &b)| t.is_some() && b > tolerance)
        .collect();
    // The partner across edge `e` of `k` when it shares the edge's
    // parameters (one chart), else `None`.
    let across = |k: usize, e: usize| -> Option<usize> {
        let t = tris[k]?;
        let j = match users.get(&key(&t, e)).map(Vec::as_slice) {
            Some([a, b]) => {
                if *a == k {
                    *b
                } else {
                    *a
                }
            }
            _ => return None,
        };
        let o = tris[j]?;
        let same = |i: usize| {
            (0..3).any(|m| o.weld[m] == t.weld[i] && o.w[m].0 == t.w[i].0 && o.w[m].1 == t.w[i].1)
        };
        (same(e) && same((e + 1) % 3)).then_some(j)
    };
    // A bent edge to a triangle outside the region would have to stay
    // bent; drop such triangles until none is left.
    loop {
        let mut changed = false;
        for k in 0..tris.len() {
            let Some(t) = tris[k].filter(|_| member[k]) else {
                continue;
            };
            let leaks = (0..3)
                .any(|e| bent.contains(&key(&t, e)) && !across(k, e).is_some_and(|j| member[j]));
            if leaks {
                member[k] = false;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    // Regions: members joined across shared edges.
    let mut seen = vec![false; tris.len()];
    for start in 0..tris.len() {
        if !member[start] || seen[start] {
            continue;
        }
        let mut region = vec![start];
        seen[start] = true;
        let mut i = 0;
        while i < region.len() {
            let k = region[i];
            i += 1;
            for e in 0..3 {
                if let Some(j) = across(k, e).filter(|&j| member[j] && !seen[j]) {
                    seen[j] = true;
                    region.push(j);
                }
            }
        }
        let interior = |k: usize, e: usize| across(k, e).is_some_and(|j| member[j]);
        let current = region.iter().map(|&k| bounds[k]).fold(0.0, f64::max);
        let mut best: Option<(f64, Vec<f64>)> = None;
        for n in LEVELS {
            let Some(values) = refine_region(surface, tris, &region, &interior, n, &mut bound)
            else {
                continue;
            };
            let worst = values.iter().copied().fold(0.0, f64::max);
            if best.as_ref().is_none_or(|(b, _)| worst < *b) {
                best = Some((worst, values));
            }
            if worst <= tolerance {
                break;
            }
        }
        if let Some((_, values)) = best.filter(|(w, _)| *w < current) {
            for (&k, v) in region.iter().zip(values) {
                bounds[k] = v;
            }
        }
    }
}

/// The refined bounds of `region`'s triangles at level `n`, `None`
/// when a sub-triangle folds (loses its triangle's orientation) or a
/// piece declines. `interior(k, e)`: edge `e` of `tris[k]` is shared
/// with another region triangle in one chart, so its sub-vertices move.
fn refine_region(
    surface: &SurfaceGeometry,
    tris: &[Option<BendTriangle>],
    region: &[usize],
    interior: &dyn Fn(usize, usize) -> bool,
    n: usize,
    bound: &mut impl FnMut([Uv; 3], [Point3; 3]) -> Option<f64>,
) -> Option<Vec<f64>> {
    let shapes: Vec<[Uv; 3]> = region
        .iter()
        .filter_map(|&k| tris[k].map(|t| t.w))
        .collect();
    let mut lo = (f64::INFINITY, f64::INFINITY);
    let mut hi = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    for q in shapes.iter().flatten() {
        lo = (lo.0.min(q.0), lo.1.min(q.1));
        hi = (hi.0.max(q.0), hi.1.max(q.1));
    }
    let mut out = Vec::with_capacity(region.len());
    for &k in region {
        let t = tris[k]?;
        let sense = orient(t.w[0], t.w[1], t.w[2]).signum();
        if sense == 0.0 {
            return None;
        }
        // Sub-vertex (i, j) at barycentric (n − i − j, i, j) / n: its
        // surface parameter and mesh point.
        let mut at: HashMap<(usize, usize), (Uv, Point3)> = HashMap::new();
        for i in 0..=n {
            for j in 0..=n - i {
                let l = [n - i - j, i, j];
                let (y, p) = sub_vertex(&t, l, n);
                let zeros: Vec<usize> = (0..3).filter(|&m| l[m] == 0).collect();
                let moves = match zeros[..] {
                    // On edge e (vertices e, e + 1) the coordinate of
                    // vertex e + 2 is zero.
                    [z] => interior(k, (z + 1) % 3),
                    [] => true,
                    _ => false,
                };
                let x = if moves {
                    inside(y, nearest(surface, y, p, (lo, hi)), &shapes)
                } else {
                    y
                };
                at.insert((i, j), (x, p));
            }
        }
        let mut worst = 0.0_f64;
        for i in 0..n {
            for j in 0..n - i {
                let mut subs = vec![[(i, j), (i + 1, j), (i, j + 1)]];
                if i + j + 1 < n {
                    subs.push([(i + 1, j), (i + 1, j + 1), (i, j + 1)]);
                }
                for s in subs {
                    let mut v = [((0.0, 0.0), Point3::new(0.0, 0.0, 0.0)); 3];
                    for (slot, q) in v.iter_mut().zip(s) {
                        *slot = *at.get(&q)?;
                    }
                    let x = v.map(|q| q.0);
                    if orient(x[0], x[1], x[2]).signum() != sense {
                        return None;
                    }
                    worst = worst.max(bound(x, v.map(|q| q.1))?);
                }
            }
        }
        out.push(worst);
    }
    Some(out)
}

/// The parameter and mesh point at barycentric `l / n` of `t`. On an
/// edge they are interpolated from its lower-welded end, so the two
/// triangles sharing it compute the same floats.
fn sub_vertex(t: &BendTriangle, l: [usize; 3], n: usize) -> (Uv, Point3) {
    let nf = n as f64;
    if let Some(z) = l.iter().position(|&x| x == 0) {
        let (mut a, mut b) = ((z + 1) % 3, (z + 2) % 3);
        if t.weld[b] < t.weld[a] {
            std::mem::swap(&mut a, &mut b);
        }
        let s = l[b] as f64;
        let r = nf - s;
        let mix = |x: f64, y: f64| (x * r + y * s) / nf;
        let (pa, pb) = (t.p[a], t.p[b]);
        return (
            (mix(t.w[a].0, t.w[b].0), mix(t.w[a].1, t.w[b].1)),
            Point3::new(mix(pa.x, pb.x), mix(pa.y, pb.y), mix(pa.z, pb.z)),
        );
    }
    let f = l.map(|x| x as f64);
    let mix = |x: [f64; 3]| (x[0] * f[0] + x[1] * f[1] + x[2] * f[2]) / nf;
    (
        (mix(t.w.map(|q| q.0)), mix(t.w.map(|q| q.1))),
        Point3::new(
            mix(t.p.map(|q| q.x)),
            mix(t.p.map(|q| q.y)),
            mix(t.p.map(|q| q.z)),
        ),
    )
}

/// `x`, or the point toward it from `y` (halving the step) strictly
/// inside one of the region's triangles `shapes`; `y` when none is (a
/// foot on the region's boundary would flatten its piece).
fn inside(y: Uv, x: Uv, shapes: &[[Uv; 3]]) -> Uv {
    let strictly = |q: Uv, t: [Uv; 3]| {
        let s = orient(t[0], t[1], t[2]).signum();
        s != 0.0 && (0..3).all(|i| orient(t[i], t[(i + 1) % 3], q) * s > 0.0)
    };
    let mut t = 1.0;
    for _ in 0..24 {
        let q = (y.0 + t * (x.0 - y.0), y.1 + t * (x.1 - y.1));
        if shapes.iter().any(|&s| strictly(q, s)) {
            return q;
        }
        t *= 0.5;
    }
    y
}

/// The parameter in the box `(lo, hi)` nearest `c` on the surface, by
/// projected Gauss–Newton from `q`; `q` itself when it does not settle
/// to a finite point.
fn nearest(surface: &SurfaceGeometry, q: Uv, c: Point3, (lo, hi): (Uv, Uv)) -> Uv {
    let mut x = q;
    for _ in 0..12 {
        let r = surface.point_at(x.0, x.1) - c;
        let (su, sv) = (surface.partial_u(x.0, x.1), surface.partial_v(x.0, x.1));
        let (g, h, k) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
        let det = g * k - h * h;
        // NaN fails too.
        if det.abs().partial_cmp(&f64::MIN_POSITIVE) != Some(std::cmp::Ordering::Greater) {
            break;
        }
        x = (
            (x.0 - (k * su.dot(&r) - h * sv.dot(&r)) / det).clamp(lo.0, hi.0),
            (x.1 - (g * sv.dot(&r) - h * su.dot(&r)) / det).clamp(lo.1, hi.1),
        );
    }
    if x.0.is_finite() && x.1.is_finite() {
        x
    } else {
        q
    }
}

#[cfg(test)]
#[path = "chordal_refine_tests.rs"]
mod tests;
