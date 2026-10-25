//! Triangulation of planar regions bounded by loops: counter-clockwise
//! outer loops, clockwise holes. Each hole is bridged to the loop
//! around it (its rightmost vertex to a visible vertex of the outer
//! loop, after Eberly's "Triangulation by Ear Clipping"), then the
//! merged polygon is ear-clipped. Only strictly convex ears are cut, and
//! a boundary vertex on an ear's diagonal blocks it, so every boundary
//! vertex (including points on straight runs, which later cuts must
//! share) is a triangle corner: no T-junctions.

/// Signed doubled area of triangle (a, b, c).
fn cross(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
}

fn signed_area(l: &[u32], at: &dyn Fn(u32) -> [f64; 2]) -> f64 {
    (0..l.len())
        .map(|k| {
            let (a, b) = (at(l[k]), at(l[(k + 1) % l.len()]));
            a[0] * b[1] - b[0] * a[1]
        })
        .sum::<f64>()
        * 0.5
}

/// Whether `p` is inside loop `l` (even-odd).
fn contains(l: &[u32], at: &dyn Fn(u32) -> [f64; 2], p: [f64; 2]) -> bool {
    let mut inside = false;
    for k in 0..l.len() {
        let (a, b) = (at(l[k]), at(l[(k + 1) % l.len()]));
        if (a[1] > p[1]) != (b[1] > p[1]) {
            let x = a[0] + (p[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]);
            if x > p[0] {
                inside = !inside;
            }
        }
    }
    inside
}

/// Triangles (counter-clockwise) covering the region bounded by
/// `loops`, with vertex positions from `at`; `None` when it fails.
pub(super) fn triangulate(
    loops: &[Vec<u32>],
    at: &dyn Fn(u32) -> [f64; 2],
) -> Option<Vec<[u32; 3]>> {
    let areas: Vec<f64> = loops.iter().map(|l| signed_area(l, at)).collect();
    let outer: Vec<usize> = (0..loops.len()).filter(|&k| areas[k] > 0.0).collect();
    let mut holes: Vec<Vec<usize>> = vec![Vec::new(); outer.len()];
    for k in (0..loops.len()).filter(|&k| areas[k] <= 0.0) {
        let p = at(*loops[k].first()?);
        // The smallest outer loop around the hole.
        let host = outer
            .iter()
            .enumerate()
            .filter(|(_, o)| contains(&loops[**o], at, p))
            .min_by(|a, b| areas[*a.1].total_cmp(&areas[*b.1]))?
            .0;
        holes[host].push(k);
    }
    let mut out = Vec::new();
    for (i, &o) in outer.iter().enumerate() {
        let mut poly = loops[o].clone();
        // Rightmost holes first.
        let mut hs = holes[i].clone();
        let right = |k: usize| {
            loops[k]
                .iter()
                .map(|&v| at(v)[0])
                .fold(f64::NEG_INFINITY, f64::max)
        };
        hs.sort_by(|a, b| right(*b).total_cmp(&right(*a)));
        for k in hs {
            poly = bridge(&poly, &loops[k], at)?;
        }
        out.extend(ear_clip(poly, at)?);
    }
    Some(out)
}

/// `poly` with `hole` joined in by a bridge from the hole's rightmost
/// vertex to a vertex of `poly` it can see.
fn bridge(poly: &[u32], hole: &[u32], at: &dyn Fn(u32) -> [f64; 2]) -> Option<Vec<u32>> {
    let m = (0..hole.len()).max_by(|&a, &b| {
        let (pa, pb) = (at(hole[a]), at(hole[b]));
        pa[0].total_cmp(&pb[0]).then(pb[1].total_cmp(&pa[1]))
    })?;
    let pm = at(hole[m]);
    // Nearest crossing of the ray +x from M with an edge of `poly`.
    let mut best: Option<(f64, usize)> = None;
    for k in 0..poly.len() {
        let (a, b) = (at(poly[k]), at(poly[(k + 1) % poly.len()]));
        if (a[1] > pm[1]) == (b[1] > pm[1]) {
            continue;
        }
        let x = a[0] + (pm[1] - a[1]) / (b[1] - a[1]) * (b[0] - a[0]);
        if x >= pm[0] && best.is_none_or(|(bx, _)| x < bx) {
            best = Some((x, k));
        }
    }
    let (x, k) = best?;
    let i = [x, pm[1]];
    // The edge's end with the larger x is visible unless a reflex
    // vertex lies inside triangle (M, I, P); then the one of those
    // making the smallest angle with the ray.
    let (ka, kb) = (k, (k + 1) % poly.len());
    let mut p = if at(poly[ka])[0] > at(poly[kb])[0] {
        ka
    } else {
        kb
    };
    let pp = at(poly[p]);
    let n = poly.len();
    let mut best_angle = f64::INFINITY;
    let tri = [pm, i, pp];
    let positive = cross(tri[0], tri[1], tri[2]) >= 0.0;
    for j in 0..n {
        let q = at(poly[j]);
        let (prev, next) = (at(poly[(j + n - 1) % n]), at(poly[(j + 1) % n]));
        let reflex = cross(prev, q, next) < 0.0;
        if !reflex || j == p {
            continue;
        }
        let s = [
            cross(tri[0], tri[1], q),
            cross(tri[1], tri[2], q),
            cross(tri[2], tri[0], q),
        ];
        let inside = if positive {
            s.iter().all(|v| *v >= 0.0)
        } else {
            s.iter().all(|v| *v <= 0.0)
        };
        if inside {
            let d = [q[0] - pm[0], q[1] - pm[1]];
            let angle = (d[1] / (d[0] * d[0] + d[1] * d[1]).sqrt()).abs();
            if angle < best_angle {
                best_angle = angle;
                p = j;
            }
        }
    }
    let mut merged = Vec::with_capacity(poly.len() + hole.len() + 2);
    merged.extend_from_slice(&poly[..=p]);
    for t in 0..=hole.len() {
        merged.push(hole[(m + t) % hole.len()]);
    }
    merged.push(poly[p]);
    merged.extend_from_slice(&poly[p + 1..]);
    Some(merged)
}

/// Ear clipping of the simple (or bridged) counter-clockwise polygon.
fn ear_clip(mut poly: Vec<u32>, at: &dyn Fn(u32) -> [f64; 2]) -> Option<Vec<[u32; 3]>> {
    let mut out = Vec::new();
    let scale = poly
        .iter()
        .map(|&v| {
            let p = at(v);
            p[0].abs().max(p[1].abs())
        })
        .fold(0.0, f64::max)
        .max(1.0);
    let eps = 1.0e-12 * scale * scale;
    while poly.len() > 3 {
        let n = poly.len();
        let mut clipped = false;
        for i in 0..n {
            let (a, b, c) = (poly[(i + n - 1) % n], poly[i], poly[(i + 1) % n]);
            let (pa, pb, pc) = (at(a), at(b), at(c));
            let turn = cross(pa, pb, pc);
            if turn <= eps {
                continue;
            }
            let blocked = poly.iter().any(|&v| {
                let q = at(v);
                if [pa, pb, pc].iter().any(|r| r[0] == q[0] && r[1] == q[1]) {
                    return false;
                }
                cross(pa, pb, q) >= 0.0 && cross(pb, pc, q) >= 0.0 && cross(pc, pa, q) >= 0.0
            });
            if !blocked {
                out.push([a, b, c]);
                poly.remove(i);
                clipped = true;
                break;
            }
        }
        if !clipped {
            // What is left has no area (a straight run, or a bridge
            // folded on itself): nothing to cover.
            let area: f64 = (0..poly.len())
                .map(|k| {
                    let (p, q) = (at(poly[k]), at(poly[(k + 1) % poly.len()]));
                    p[0] * q[1] - q[0] * p[1]
                })
                .sum();
            if area.abs() <= eps * poly.len() as f64 {
                return Some(out);
            }
            return None;
        }
    }
    if poly.len() == 3 {
        let (pa, pb, pc) = (at(poly[0]), at(poly[1]), at(poly[2]));
        if cross(pa, pb, pc) > eps {
            out.push([poly[0], poly[1], poly[2]]);
        }
    }
    Some(out)
}
