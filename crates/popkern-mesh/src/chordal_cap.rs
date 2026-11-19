//! Cap-chart pairings for the chordal certificate (chart-service.md §6a,
//! PK-1267).
//!
//! A triangle the atlas meshed in a cap (the pole projection about a
//! chart singularity) is flat in that chart, not in the face chart: two
//! vertices a few degrees from a sphere's pole sit far apart in longitude,
//! and the face-parameter segment between them runs round the latitude
//! circle, far from their chord. Paired affinely in face parameters, such
//! a triangle's bound is that arc's sagitta (2.7e-2 on a full sphere whose
//! sampled deviation is 2.3e-3), however good the mesh.
//!
//! A cap pairing follows the cap chart instead. Each triangle `T` of a
//! region `R` (the cap triangles every other pairing leaves over the
//! tolerance) is cut into `n²` sub-triangles in the chart's coordinates
//! `y` (the projection of its corners onto the plane at the singular
//! point). A sub-vertex maps to the face parameters of the surface point
//! projecting to it (`RegularChart::to_surface_near`), and to the point
//! `L(y)` of the emitted triangle; on each sub-triangle the pairing is the
//! affine bound of its parameter triangle against `L` of its corners, with
//! the corners' discrepancies (a sub-triangle at the singular point is the
//! row line's rectangle, `corner_pieces`). A sub-vertex on an edge `R`
//! shares with a triangle outside it is interpolated in face parameters,
//! so the neighbour's straight pairing meets it; shared edges are cut from
//! their lower-keyed corner, so both sides read the same parameters.
//!
//! The pairings tile `R` when every sub-triangle keeps the chart's
//! orientation (checked, as `chordal_refine` checks its maps), and as `n`
//! doubles the sub-triangles' arcs shrink by four: the bound tends to the
//! mesh's own deviation.

use popkern_math::hash::FastHashMap as HashMap;
use popkern_math::hash::FastHashSet as HashSet;

use popkern_core::geometry::Vec3;
use popkern_core::geometry::chart::atlas::{CAP_LIMIT, cap_chart, cap_row};
use popkern_core::geometry::chart::ops::singular_point_near;
use popkern_core::geometry::chart::{ChartDescriptor, SingularKind};
use popkern_core::geometry::chart::{ChartedSurface, RegularChart};
use popkern_core::geometry::{Point3, SurfaceEval, SurfaceGeometry};
use popkern_math::PortableMath;

use super::chordal_pieces::{Piece, unwrap};

type Uv = (f64, f64);

/// Subdivision levels, coarsest first.
const LEVELS: [usize; 3] = [2, 4, 8];

/// A triangle for the cap pass: weld ids, face parameters, positions.
#[derive(Debug, Clone, Copy)]
pub(super) struct CapTriangle {
    pub(super) weld: [usize; 3],
    pub(super) w: [Uv; 3],
    pub(super) p: [Point3; 3],
}

/// One sub-vertex: face parameters and the emitted point paired with it.
type Sub = (Uv, Point3);

/// Lower `bounds[k]` of the triangles over `tolerance` lying in a cap by
/// cap pairings; `bound(piece)` is a piece's certified bound, `None` where
/// it declines (the triangle then keeps its bound).
pub(super) fn cap_pass(
    surface: &SurfaceGeometry,
    descriptor: &ChartDescriptor,
    tris: &[CapTriangle],
    bounds: &mut [f64],
    (tolerance, bent): (f64, &HashSet<(usize, usize)>),
    mut bound: impl FnMut(&Piece) -> Option<f64>,
) {
    let rows: Vec<RegularChart> = descriptor
        .singular_rows
        .iter()
        .filter(|r| r.kind == SingularKind::Chart)
        .filter_map(|r| cap_chart(descriptor, surface, &cap_row(descriptor, surface, r).ok()?))
        .collect();
    if rows.is_empty() {
        return;
    }
    let periods = ChartDescriptor::evaluation_periods_of(surface);
    // A triangle faces the nearest chart when its corners' normals lie
    // within the cap limit of the chart's (the projection a graph there,
    // its Jacobian at least cos 60°).
    let chart_of = |t: &CapTriangle| -> Option<usize> {
        let c = Point3::new(
            (t.p[0].x + t.p[1].x + t.p[2].x) / 3.0,
            (t.p[0].y + t.p[1].y + t.p[2].y) / 3.0,
            (t.p[0].z + t.p[1].z + t.p[2].z) / 3.0,
        );
        let k = (0..rows.len()).min_by(|&a, &b| {
            let d = |i: usize| {
                rows[i]
                    .origin()
                    .map_or(f64::INFINITY, |o| o.distance_to(&c))
            };
            d(a).total_cmp(&d(b))
        })?;
        let n = rows[k].plane_normal()?;
        // The service's normal (the row's at the singular point, where the
        // raw cross product vanishes), read within rounding of the point.
        let facing = t.w.iter().zip(&t.p).all(|(&uv, p)| {
            let exact = 1e-12 * (1.0 + p.x.abs().max(p.y.abs()).max(p.z.abs()));
            popkern_core::geometry::chart::ops::chart_normal(descriptor, surface, uv, *p, exact)
                .is_some_and(|m| m.dot(&n) >= CAP_LIMIT.pcos())
        });
        facing.then_some(k)
    };
    let facing: Vec<Option<usize>> = tris.iter().map(chart_of).collect();
    for (c, chart) in rows.iter().enumerate() {
        // A cap's region: every triangle facing its chart, when one of
        // them is over the tolerance (its pairings replace the region's
        // together: neighbours on one curved edge must read it alike).
        if !(0..tris.len()).any(|k| facing[k] == Some(c) && bounds[k] > tolerance) {
            continue;
        }
        let mut member: Vec<bool> = facing.iter().map(|f| *f == Some(c)).collect();
        // A triangle on the region's boundary whose edge there the bend
        // pass bent leaves it (the neighbour's pairing reads that edge
        // bent; the region's would read it straight).
        loop {
            let users = edge_users(tris, &member);
            let mut changed = false;
            for (k, t) in tris.iter().enumerate() {
                if member[k]
                    && (0..3).any(|e| {
                        let key = edge_key(t, e);
                        bent.contains(&key) && users.get(&key).is_none_or(|&n| n < 2)
                    })
                {
                    member[k] = false;
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
        let users = edge_users(tris, &member);
        // Each member's unit facet normal, facing the chart's, and each
        // edge's sum of its members' (the direction both triangles on it
        // pair its points along, so they read one point).
        let Some(axis) = chart.plane_normal() else {
            continue;
        };
        let facet = |t: &CapTriangle| -> Vec3 {
            let m = (t.p[1] - t.p[0]).cross(&(t.p[2] - t.p[0]));
            let m = m.normalized().unwrap_or(axis);
            if m.dot(&axis) < 0.0 { -m } else { m }
        };
        let mut edge_normal: HashMap<(usize, usize), Vec3> = HashMap::default();
        for (_, t) in tris.iter().enumerate().filter(|(k, _)| member[*k]) {
            let m = facet(t);
            for e in 0..3 {
                let sum = edge_normal
                    .entry(edge_key(t, e))
                    .or_insert(Vec3::new(0.0, 0.0, 0.0));
                *sum = *sum + m;
            }
        }
        let old = (0..tris.len())
            .filter(|&k| member[k])
            .map(|k| bounds[k])
            .fold(0.0, f64::max);
        for &n in &LEVELS {
            let mut level: Vec<(usize, f64)> = Vec::new();
            let mut declined = false;
            for (k, t) in tris.iter().enumerate().filter(|(k, _)| member[*k]) {
                // An edge follows the chart only between two members;
                // else it stays straight in face parameters.
                let curved = [0, 1, 2].map(|e| users.get(&edge_key(t, e)).is_some_and(|&u| u == 2));
                let normals = (
                    [0, 1, 2].map(|e| edge_normal.get(&edge_key(t, e)).copied().unwrap_or(axis)),
                    facet(t),
                );
                match triangle_at_level(
                    (surface, descriptor),
                    (chart, periods),
                    (t, normals),
                    curved,
                    n,
                    &mut bound,
                ) {
                    Some(b) => level.push((k, b)),
                    None => {
                        declined = true;
                        break;
                    }
                }
            }
            if declined {
                break;
            }
            let worst = level.iter().map(|x| x.1).fold(0.0, f64::max);
            if worst < old {
                for (k, b) in level {
                    bounds[k] = b;
                }
            }
            if worst <= tolerance {
                break;
            }
        }
    }
}

fn edge_key(t: &CapTriangle, e: usize) -> (usize, usize) {
    let (a, b) = (t.weld[e], t.weld[(e + 1) % 3]);
    (a.min(b), a.max(b))
}

fn edge_users(tris: &[CapTriangle], member: &[bool]) -> HashMap<(usize, usize), usize> {
    let mut users: HashMap<(usize, usize), usize> = HashMap::default();
    for (_, t) in tris.iter().enumerate().filter(|(k, _)| member[*k]) {
        for e in 0..3 {
            *users.entry(edge_key(t, e)).or_default() += 1;
        }
    }
    users
}

/// The cap pairing of `t` cut `n × n`, or `None` when a sub-vertex has no
/// chart preimage, a sub-triangle folds, or a piece declines.
///
/// A triangle with a corner at the singular point is cut on a polar grid
/// (rows parallel to the far edge, `n + 1` points on each): its
/// sub-triangles at the point then span `1/n` of its angle there, where a
/// barycentric grid's keep the whole angle at every level and their face
/// parameters' arcs never shrink. Along each edge both grids place the
/// same `n + 1` points, so neighbours of either kind meet.
fn triangle_at_level(
    (surface, descriptor): (&SurfaceGeometry, &ChartDescriptor),
    (chart, periods): (&RegularChart, (Option<f64>, Option<f64>)),
    (t, (edge_normals, facet_normal)): (&CapTriangle, ([Vec3; 3], Vec3)),
    curved: [bool; 3],
    n: usize,
    bound: &mut impl FnMut(&Piece) -> Option<f64>,
) -> Option<f64> {
    let charted = ChartedSurface::new(
        surface,
        chart.clone(),
        ((f64::MIN, f64::MAX), (f64::MIN, f64::MAX)),
    );
    let y: [Uv; 3] = [0, 1, 2].map(|i| chart.local_at(t.w[i], &t.p[i]));
    // Corner parameters on one branch, about corner 0.
    let w = unwrap(t.w, periods);
    let origin = chart.origin();
    let at_origin = |q: &Point3| {
        origin.is_some_and(|o| {
            let scale = 1.0 + o.x.abs().max(o.y.abs()).max(o.z.abs());
            q.distance_to(&o) <= 1e-12 * scale
        })
    };
    let lerp2 = |a: Uv, b: Uv, s: f64| (a.0 + (b.0 - a.0) * s, a.1 + (b.1 - a.1) * s);
    let lerp3 = |a: Point3, b: Point3, s: f64| a + (b - a) * s;
    // The surface point paired with the facet point `l`: along `d` from
    // it, found from the chart point `yy` (`along_normal`).
    let preimage = |l: Point3, d: Vec3, yy: Uv, seed: Uv| {
        along_normal((surface, chart, &charted), l, d, (yy, seed))
    };
    // A point of edge `e` (corner e → e + 1) at fraction `s` from e, cut
    // from its lower-keyed corner so both triangles on it read one point.
    let edge_point = |e: usize, s: f64| -> Option<Sub> {
        let (a, b) = (e, (e + 1) % 3);
        let (lo, hi, f) = if t.weld[a] <= t.weld[b] {
            (a, b, s)
        } else {
            (b, a, 1.0 - s)
        };
        let l = lerp3(t.p[lo], t.p[hi], f);
        // A corner on the singular point carries an arbitrary collapsed
        // coordinate: an edge into it is read from the row at its
        // partner's (a corner ball's boundary meridian lerped from the
        // pole's stored longitude spiralled off the meridian, and the
        // region declined at every level).
        let on_row = |at: usize, partner: usize| -> Uv {
            let exact = 1e-12 * (1.0 + t.p[at].x.abs().max(t.p[at].y.abs()).max(t.p[at].z.abs()));
            match singular_point_near(descriptor, &t.p[at], exact) {
                Some((row, true)) => (w[partner].0, row.at),
                Some((row, false)) => (row.at, w[partner].1),
                None => w[at],
            }
        };
        let wl = lerp2(on_row(lo, hi), on_row(hi, lo), f);
        if curved[e] {
            Some((preimage(l, edge_normals[e], lerp2(y[lo], y[hi], f), wl)?, l))
        } else {
            Some((wl, l))
        }
    };
    // An interior point with weights `k` on the corners.
    let interior = |k: [f64; 3]| -> Option<Sub> {
        let mix2 = |q: [Uv; 3]| {
            (
                k[0] * q[0].0 + k[1] * q[1].0 + k[2] * q[2].0,
                k[0] * q[0].1 + k[1] * q[1].1 + k[2] * q[2].1,
            )
        };
        let l = Point3::new(
            k[0] * t.p[0].x + k[1] * t.p[1].x + k[2] * t.p[2].x,
            k[0] * t.p[0].y + k[1] * t.p[1].y + k[2] * t.p[2].y,
            k[0] * t.p[0].z + k[1] * t.p[1].z + k[2] * t.p[2].z,
        );
        Some((preimage(l, facet_normal, mix2(y), mix2(w))?, l))
    };
    let nf = n as f64;
    let mut subs: HashMap<(usize, usize), Sub> = HashMap::default();
    let mut shapes: Vec<[(usize, usize); 3]> = Vec::new();
    if let Some(c) = (0..3).find(|&c| at_origin(&t.p[c])) {
        // Polar grid: (i, j) is row i from the point (i = 0 the point
        // itself), column j from corner A toward corner B.
        let (a, b) = ((c + 1) % 3, (c + 2) % 3);
        for i in 0..=n {
            for j in 0..=n {
                let sub = match (i, j) {
                    (0, _) => (w[c], t.p[c]),
                    (_, 0) if i == n => (w[a], t.p[a]),
                    (_, _) if i == n && j == n => (w[b], t.p[b]),
                    (_, 0) => edge_point(c, i as f64 / nf)?,
                    (_, _) if j == n => edge_point(b, 1.0 - i as f64 / nf)?,
                    (_, _) if i == n => edge_point(a, j as f64 / nf)?,
                    _ => {
                        let (r, s) = (i as f64 / nf, j as f64 / nf);
                        let mut k = [0.0; 3];
                        k[c] = 1.0 - r;
                        k[a] = r * (1.0 - s);
                        k[b] = r * s;
                        interior(k)?
                    }
                };
                subs.insert((i, j), sub);
            }
        }
        for j in 0..n {
            shapes.push([(0, 0), (1, j), (1, j + 1)]);
        }
        for i in 1..n {
            for j in 0..n {
                shapes.push([(i, j), (i + 1, j), (i + 1, j + 1)]);
                shapes.push([(i, j), (i + 1, j + 1), (i, j + 1)]);
            }
        }
        // Row 0's points are all the singular point.
        for j in 1..=n {
            subs.insert((0, j), (w[c], t.p[c]));
        }
    } else {
        // Barycentric grid: (i, j) with weights (n − i − j, i, j).
        for i in 0..=n {
            for j in 0..=(n - i) {
                let k0 = n - i - j;
                let sub = match (k0, i, j) {
                    (_, 0, 0) => (w[0], t.p[0]),
                    (0, _, 0) => (w[1], t.p[1]),
                    (0, 0, _) => (w[2], t.p[2]),
                    (_, _, 0) => edge_point(0, i as f64 / nf)?,
                    (0, _, _) => edge_point(1, j as f64 / nf)?,
                    (_, 0, _) => edge_point(2, k0 as f64 / nf)?,
                    _ => interior([k0 as f64 / nf, i as f64 / nf, j as f64 / nf])?,
                };
                subs.insert((i, j), sub);
            }
        }
        for i in 0..n {
            for j in 0..(n - i) {
                shapes.push([(i, j), (i + 1, j), (i, j + 1)]);
                if i + j + 1 < n {
                    shapes.push([(i + 1, j), (i + 1, j + 1), (i, j + 1)]);
                }
            }
        }
    }
    // Every sub-triangle off the singular point keeps one orientation
    // against the chart's (the pairings then tile the triangle).
    let cap_area = cross(y[0], y[1], y[2]);
    let mut sense: Option<f64> = None;
    let mut worst = 0.0_f64;
    for s in shapes {
        let corners = s.map(|ij| subs[&ij]);
        let ws = unwrap(corners.map(|c| c.0), periods);
        let ls = corners.map(|c| c.1);
        if !corners
            .iter()
            .any(|c| at_origin(&surface.point_at(c.0.0, c.0.1)))
        {
            let a = (cross(ws[0], ws[1], ws[2]) * cap_area).signum();
            if a == 0.0 || sense.is_some_and(|sg| sg != a) {
                return None;
            }
            sense = Some(a);
        }
        for piece in corner_pieces(descriptor, ws, ls, 0.0) {
            worst = worst.max(bound(&piece)?);
        }
    }
    Some(worst)
}

/// The face parameters of the surface point on the line through `l`
/// along `d`, by Newton in the chart's coordinates from `yy` (face
/// parameters `seed` near it). Paired along the chart's axis instead (the
/// chart point's preimage), a facet tilted from the axis read its
/// deviation over the cosine of the tilt: 5.5e-3 for a 4.1e-3 witness by a
/// 45° row, at every level. Unconverged, the last iterate stands: any
/// surface point is a sound pairing, only a looser one.
fn along_normal(
    (surface, chart, charted): (&SurfaceGeometry, &RegularChart, &ChartedSurface<'_>),
    l: Point3,
    d: Vec3,
    (yy, seed): (Uv, Uv),
) -> Option<Uv> {
    let d = d.normalized()?;
    let scale = 1.0 + l.x.abs().max(l.y.abs()).max(l.z.abs());
    let mut y = yy;
    let mut w = chart.to_surface_near(surface, y, Some(seed)).ok()?;
    for _ in 0..24 {
        let r = surface.point_at(w.0, w.1) - l;
        let rp = r - d * r.dot(&d);
        if rp.length() <= 1e-14 * scale {
            break;
        }
        let (px, py) = charted.first_partials(y.0, y.1);
        let (ax, ay) = (px - d * px.dot(&d), py - d * py.dot(&d));
        let (a11, a12, a22) = (ax.dot(&ax), ax.dot(&ay), ay.dot(&ay));
        let (b1, b2) = (-ax.dot(&rp), -ay.dot(&rp));
        let det = a11 * a22 - a12 * a12;
        if !det.is_finite() || det.abs() <= 1e-300 {
            break;
        }
        let next = (
            y.0 + (b1 * a22 - b2 * a12) / det,
            y.1 + (a11 * b2 - a12 * b1) / det,
        );
        let Ok(wn) = chart.to_surface_near(surface, next, Some(w)) else {
            break;
        };
        (y, w) = (next, wn);
    }
    Some(w)
}

/// A sub-triangle's parameter pieces. A corner on a singular point (to
/// rounding, read by the service in space) carries an arbitrary collapsed
/// coordinate: the triangle `P A B` stands for the region from the row
/// through `P` to `A`'s and `B`'s lines, mapped with both row corners at
/// `P` (a rectangle when `A` and `B` share the across coordinate, else two
/// triangles), so the pieces round the point tile it.
/// A triangle's pairing pieces when a corner lies at a singular row (read
/// from the descriptor in space within `reach`, or to rounding): the
/// corner stands for the whole row line, so the triangle is the rectangle
/// (or the quad split in two) from the row to its other corners, both row
/// corners at the singular point. The one authority for a corner at a
/// chart singularity in the certificate, for its first pass and its cap
/// pass alike.
pub(super) fn corner_pieces(
    d: &ChartDescriptor,
    w: [Uv; 3],
    p: [Point3; 3],
    reach: f64,
) -> Vec<Piece> {
    for k in 0..3 {
        let exact = 1e-12 * (1.0 + p[k].x.abs().max(p[k].y.abs()).max(p[k].z.abs()));
        let Some((row, collapses_u)) = singular_point_near(d, &p[k], exact.max(reach)) else {
            continue;
        };
        let (a, b) = ((k + 1) % 3, (k + 2) % 3);
        let on_row = |q: Uv| {
            if collapses_u {
                (q.0, row.at)
            } else {
                (row.at, q.1)
            }
        };
        let across = |q: Uv| if collapses_u { q.1 } else { q.0 };
        let along = |q: Uv| if collapses_u { q.0 } else { q.1 };
        let (pa, pb) = (on_row(w[a]), on_row(w[b]));
        if across(w[a]) == across(w[b]) {
            let span = (along(w[a]), along(w[b]));
            let height = (row.at, across(w[a]));
            return vec![if collapses_u {
                Piece::Rectangle(span, height, [p[k], p[k], p[a], p[b]])
            } else {
                Piece::Rectangle(height, span, [p[k], p[a], p[k], p[b]])
            }];
        }
        return vec![
            Piece::Triangle([pa, w[a], w[b]], [p[k], p[a], p[b]]),
            Piece::Triangle([pa, w[b], pb], [p[k], p[b], p[k]]),
        ];
    }
    vec![Piece::Triangle(w, p)]
}

fn cross(a: Uv, b: Uv, c: Uv) -> f64 {
    (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
}
