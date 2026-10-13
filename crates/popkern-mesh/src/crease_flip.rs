//! Crease-reducing edge flips — the second arm of
//! [`crate::TessellationOptions::max_facet_angle`] after
//! [`crate::angular_refine`].
//!
//! Angular refinement splits edges whose vertex normals spread past
//! the bound. It cannot see a crease whose vertex normals agree: two
//! facets that tilt off their own vertex normals in opposite senses
//! because their shared diagonal runs the wrong way across an
//! anisotropic cell (PK-690: thread runout corners crease 14–18° with
//! every vertex normal within 2°). Splitting there degenerates, since
//! longest-edge bisection keeps the bad diagonal and fans slivers onto
//! it. Swapping the diagonal is what removes the crease.
//!
//! The pass flips an interior edge `ab` of the quad `a d b c` to `cd`
//! when the largest crease over the quad's five edges strictly drops,
//! and only when the flip is admissible:
//! - both new triangles keep the chart orientation of the old pair,
//!   with a nonzero chart area (no fold in the parameter domain);
//! - each new facet faces the same way as its vertex normals;
//! - the new diagonal sags off the surface no more than the chord
//!   height (or the old diagonal, whichever is larger), so the chord
//!   contract holds;
//! - no vertex of the quad is chart-singular (its normal is undefined);
//! - the worse new triangle keeps at least half the shape quality of
//!   the worse old one, so a crease is never traded for a needle;
//! - neither new triangle's vertex normals spread past the angular
//!   bound (or past the old pair's spread, when that was larger), so
//!   the refinement's bound survives.
//!
//! Boundary edges have one triangle in the face mesh and are never
//! flipped, so the canonical edge samples neighbouring faces share are
//! untouched and watertightness holds by construction. Every accepted
//! flip lowers a local maximum; a flip budget bounds the work.

use std::collections::HashMap;

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

use crate::mesh::FaceMesh;
use crate::vertex_normal::singular_without_unique_normal;

/// Flips allowed per input triangle.
const FLIP_BUDGET_PER_TRIANGLE: usize = 4;

/// Strict-improvement margin (radians) a flip must clear.
const IMPROVEMENT: f64 = 1.0e-9;

pub(crate) type EdgeKey = (usize, usize);

pub(crate) fn key(a: usize, b: usize) -> EdgeKey {
    (a.min(b), a.max(b))
}

pub(crate) fn angle_between(a: Vec3, b: Vec3) -> f64 {
    let (la, lb) = (a.length(), b.length());
    if la < 1e-300 || lb < 1e-300 {
        return 0.0;
    }
    (a.dot(&b) / (la * lb)).clamp(-1.0, 1.0).pacos()
}

pub(crate) fn facet_normal(mesh: &FaceMesh, t: [usize; 3]) -> Vec3 {
    let p = |k: usize| mesh.vertices[t[k]].position;
    (p(1) - p(0)).cross(&(p(2) - p(0)))
}

pub(crate) fn chart_area(mesh: &FaceMesh, t: [usize; 3]) -> f64 {
    let uv = |k: usize| mesh.vertices[t[k]].uv;
    let (a, b, c) = (uv(0), uv(1), uv(2));
    0.5 * ((b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0))
}

/// The vertex of `t` that is neither `a` nor `b`.
pub(crate) fn apex(t: [usize; 3], a: usize, b: usize) -> Option<usize> {
    t.iter().copied().find(|&v| v != a && v != b)
}

/// Whether `t` traverses `a → b`.
pub(crate) fn runs(t: [usize; 3], a: usize, b: usize) -> bool {
    (0..3).any(|k| t[k] == a && t[(k + 1) % 3] == b)
}

/// Shape quality of a triangle: `area / longest edge²` (an equilateral
/// triangle scores √3/4, a needle or cap scores near 0).
fn quality(mesh: &FaceMesh, t: [usize; 3]) -> f64 {
    let p = |k: usize| mesh.vertices[t[k]].position;
    let longest = (0..3)
        .map(|k| p(k).distance_to(&p((k + 1) % 3)))
        .fold(0.0, f64::max);
    if longest < 1e-300 {
        return 0.0;
    }
    0.5 * facet_normal(mesh, t).length() / (longest * longest)
}

/// Largest angle between a triangle's vertex normals.
pub(crate) fn spread(mesh: &FaceMesh, t: [usize; 3]) -> f64 {
    (0..3)
        .map(|k| {
            angle_between(
                mesh.vertices[t[k]].normal,
                mesh.vertices[t[(k + 1) % 3]].normal,
            )
        })
        .fold(0.0, f64::max)
}

/// Crease across interior edge `e`; zero on a boundary edge.
fn crease(mesh: &FaceMesh, edge_tris: &HashMap<EdgeKey, Vec<usize>>, e: EdgeKey) -> f64 {
    match edge_tris.get(&e).map(Vec::as_slice) {
        Some(&[t0, t1]) => angle_between(
            facet_normal(mesh, mesh.triangles[t0]),
            facet_normal(mesh, mesh.triangles[t1]),
        ),
        _ => 0.0,
    }
}

/// Distance between a chord's midpoint and the surface point over its
/// chart midpoint; `None` when that point leaves the chord's
/// neighbourhood (the chart midpoint crossed a seam or singularity).
pub(crate) fn chord_sag(
    mesh: &FaceMesh,
    surface: &SurfaceGeometry,
    a: usize,
    b: usize,
) -> Option<f64> {
    let (va, vb) = (&mesh.vertices[a], &mesh.vertices[b]);
    let uv = (0.5 * (va.uv.0 + vb.uv.0), 0.5 * (va.uv.1 + vb.uv.1));
    let mid = Point3::new(
        0.5 * (va.position.x + vb.position.x),
        0.5 * (va.position.y + vb.position.y),
        0.5 * (va.position.z + vb.position.z),
    );
    let sag = surface.point_at(uv.0, uv.1).distance_to(&mid);
    (sag <= 0.5 * va.position.distance_to(&vb.position)).then_some(sag)
}

/// Flip interior edges of `mesh` while a flip strictly lowers the
/// largest crease over its quad (see the module docs). Only quads
/// creased past `crease_bound` are candidates (a crease within the
/// bound is what the facet bound promises already); `spread_bound`
/// caps new triangles' vertex-normal spread. Returns the
/// triangles the flips rewrote, sorted and deduplicated.
pub(crate) fn flip_creases(
    mesh: &mut FaceMesh,
    surface: &SurfaceGeometry,
    chord_height: f64,
    crease_bound: f64,
    spread_bound: f64,
) -> Vec<usize> {
    if mesh.triangles.len() < 2 || matches!(surface, SurfaceGeometry::Plane { .. }) {
        return Vec::new();
    }
    let singular: Vec<bool> = mesh
        .vertices
        .iter()
        .map(|v| singular_without_unique_normal(surface, v.uv.0, v.uv.1))
        .collect();
    let mut edge_tris: HashMap<EdgeKey, Vec<usize>> = HashMap::new();
    for (ti, t) in mesh.triangles.iter().enumerate() {
        for k in 0..3 {
            edge_tris
                .entry(key(t[k], t[(k + 1) % 3]))
                .or_default()
                .push(ti);
        }
    }
    let mut queue: Vec<EdgeKey> = edge_tris
        .iter()
        .filter(|(_, ts)| ts.len() == 2)
        .map(|(&e, _)| e)
        .collect();
    queue.sort_unstable();
    let mut budget = mesh.triangles.len() * FLIP_BUDGET_PER_TRIANGLE;
    let mut flipped: Vec<usize> = Vec::new();
    while let Some(e) = queue.pop() {
        if budget == 0 {
            break;
        }
        let Some(&[t0, t1]) = edge_tris.get(&e).map(Vec::as_slice) else {
            continue;
        };
        // Orient so t0 runs a → b and t1 runs b → a.
        let (a, b) = if runs(mesh.triangles[t0], e.0, e.1) {
            (e.0, e.1)
        } else {
            (e.1, e.0)
        };
        let (tri0, tri1) = (mesh.triangles[t0], mesh.triangles[t1]);
        if !runs(tri0, a, b) || !runs(tri1, b, a) {
            continue; // inconsistent winding: leave it alone
        }
        let (Some(c), Some(d)) = (apex(tri0, a, b), apex(tri1, a, b)) else {
            continue;
        };
        if c == d || edge_tris.contains_key(&key(c, d)) || [a, b, c, d].iter().any(|&v| singular[v])
        {
            continue;
        }
        let quad = [key(a, b), key(a, c), key(b, c), key(a, d), key(b, d)];
        let before = quad
            .iter()
            .map(|&q| crease(mesh, &edge_tris, q))
            .fold(0.0, f64::max);
        if before <= crease_bound {
            continue;
        }
        // New pair: (a, d, c) and (d, b, c).
        let (new0, new1) = ([a, d, c], [d, b, c]);
        let sense = chart_area(mesh, tri0).signum();
        let quad_area = chart_area(mesh, tri0).abs() + chart_area(mesh, tri1).abs();
        let chart_ok = [new0, new1].iter().all(|&t| {
            let area = chart_area(mesh, t);
            area.signum() == sense && area.abs() > 1e-9 * quad_area
        });
        let facing_ok = [new0, new1].iter().all(|&t| {
            let n = facet_normal(mesh, t);
            t.iter().all(|&v| n.dot(&mesh.vertices[v].normal) > 0.0)
        });
        // Shape guard: trading a crease for a needle only moves the
        // defect (measured: unguarded flips left 120 near-zero-area
        // triangles on mixed_radius_corner's caps).
        let worst_old = quality(mesh, tri0).min(quality(mesh, tri1));
        let worst_new = quality(mesh, new0).min(quality(mesh, new1));
        // Spread guard: the refinement already bounded every
        // triangle's vertex-normal spread; a flip must not undo it
        // (measured: an unguarded flip on mixed_radius_corner joined
        // two vertices whose normals spread 19°).
        let spread_old = spread(mesh, tri0).max(spread(mesh, tri1));
        let spread_new = spread(mesh, new0).max(spread(mesh, new1));
        if !chart_ok
            || !facing_ok
            || worst_new < 0.5 * worst_old
            || spread_new > spread_bound.max(spread_old)
        {
            continue;
        }
        let old_sag = chord_sag(mesh, surface, a, b).unwrap_or(0.0);
        if chord_sag(mesh, surface, c, d).is_none_or(|s| s > chord_height.max(old_sag)) {
            continue;
        }
        // Trial: apply, measure, keep or revert.
        mesh.triangles[t0] = new0;
        mesh.triangles[t1] = new1;
        let moved = |edge_tris: &mut HashMap<EdgeKey, Vec<usize>>, e: EdgeKey, from, to| {
            if let Some(list) = edge_tris.get_mut(&e) {
                for x in list.iter_mut() {
                    if *x == from {
                        *x = to;
                    }
                }
            }
        };
        edge_tris.remove(&key(a, b));
        edge_tris.insert(key(c, d), vec![t0, t1]);
        moved(&mut edge_tris, key(a, d), t1, t0);
        moved(&mut edge_tris, key(b, c), t0, t1);
        let after_quad = [key(c, d), key(a, c), key(b, c), key(a, d), key(b, d)];
        let after = after_quad
            .iter()
            .map(|&q| crease(mesh, &edge_tris, q))
            .fold(0.0, f64::max);
        if after < before - IMPROVEMENT {
            budget -= 1;
            flipped.extend([t0, t1]);
            queue.extend([key(a, c), key(b, c), key(a, d), key(b, d)]);
            continue;
        }
        mesh.triangles[t0] = tri0;
        mesh.triangles[t1] = tri1;
        edge_tris.remove(&key(c, d));
        edge_tris.insert(key(a, b), vec![t0, t1]);
        moved(&mut edge_tris, key(a, d), t0, t1);
        moved(&mut edge_tris, key(b, c), t1, t0);
    }
    flipped.sort_unstable();
    flipped.dedup();
    flipped
}

#[cfg(test)]
#[path = "crease_flip_tests.rs"]
mod tests;
