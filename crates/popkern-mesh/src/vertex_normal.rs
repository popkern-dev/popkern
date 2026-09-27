//! Mesh vertex normals at chart singularities.
//!
//! `SurfaceEval::normal_at` is the raw `∂u × ∂v`: it vanishes where
//! the chart degenerates (sphere and cone poles, collapsed NURBS or
//! Gregory edges) and FLIPS just past a pole, where the degenerate
//! partial changes sign. Boundary vertices take their uv from point
//! inversion, which can land a hair outside the natural range —
//! measured: the polyhedral dice's sphere-corner poles sat at
//! v = π + 4e-8, so their stored normals pointed INTO the ball (174°
//! from their neighbours), a dark spot in shaded views and a false
//! stripe knot in zebra views, at every tessellation setting.
//!
//! [`repair_singular_normals`] re-derives the normal of every vertex
//! that is chart-singular or outside the natural parameter range as
//! the one-sided LIMIT normal ([`limit_normal`]): clamp into the
//! range, then step a hair off the degenerate side. Healthy vertices
//! keep their normals bit-for-bit.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use crate::mesh::FaceMesh;

/// Relative degeneracy threshold: `|∂u × ∂v|` below this fraction of
/// the larger partial's squared length reads as a singular chart
/// point (a pole approached from either side, a collapsed edge).
const SINGULAR_RATIO: f64 = 1.0e-6;

/// Fraction of a finite parameter span stepped off a singular point
/// to read its limit normal.
const LIMIT_STEP: f64 = 1.0e-5;

/// Whether the chart is singular at `(u, v)`: the normal there is a
/// limit, not a direction the raw cross product can supply.
pub(crate) fn chart_singular(surface: &SurfaceGeometry, u: f64, v: f64) -> bool {
    if matches!(surface, SurfaceGeometry::Plane { .. }) {
        return false;
    }
    let (pu, pv) = (surface.partial_u(u, v), surface.partial_v(u, v));
    let big = pu.length().max(pv.length());
    big < 1e-300 || pu.cross(&pv).length() < SINGULAR_RATIO * big * big
}

fn outside_range(surface: &SurfaceGeometry, u: f64, v: f64) -> bool {
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    let (pu, pv) = surface.parameter_periodicity();
    (pu.is_none() && (u < u0 || u > u1)) || (pv.is_none() && (v < v0 || v > v1))
}

fn unit(n: Vec3) -> Option<Vec3> {
    let l = n.length();
    (l > 1e-300 && l.is_finite()).then(|| Vec3::new(n.x / l, n.y / l, n.z / l))
}

/// The surface's unit normal at `(u, v)` in the chart's own sense,
/// taken as the one-sided limit at singular or out-of-range points:
/// clamp non-periodic coordinates into the natural range, then, while
/// the chart stays singular, step toward the range interior along the
/// coordinate whose partial collapsed. `None` when no nearby regular
/// point exists.
pub(crate) fn limit_normal(surface: &SurfaceGeometry, u: f64, v: f64) -> Option<Vec3> {
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    let (per_u, per_v) = surface.parameter_periodicity();
    let clamp = |x: f64, lo: f64, hi: f64, periodic: bool| {
        if periodic || !lo.is_finite() || !hi.is_finite() {
            x
        } else {
            x.clamp(lo, hi)
        }
    };
    let (mut u, mut v) = (
        clamp(u, u0, u1, per_u.is_some()),
        clamp(v, v0, v1, per_v.is_some()),
    );
    let step = |lo: f64, hi: f64, x: f64| -> f64 {
        let span = if lo.is_finite() && hi.is_finite() {
            hi - lo
        } else {
            1.0
        };
        let toward = if lo.is_finite() && hi.is_finite() {
            if x > 0.5 * (lo + hi) { -1.0 } else { 1.0 }
        } else if lo.is_finite() {
            1.0
        } else {
            -1.0
        };
        toward * LIMIT_STEP * span
    };
    for _ in 0..4 {
        if !chart_singular(surface, u, v) {
            return unit(surface.normal_at(u, v));
        }
        // Step along the coordinate that moves AWAY from the collapse:
        // a vanishing ∂u (sphere pole row) is cured by moving in v.
        let (pu, pv) = (surface.partial_u(u, v), surface.partial_v(u, v));
        if pu.length() <= pv.length() {
            v += step(v0, v1, v);
        } else {
            u += step(u0, u1, u);
        }
    }
    unit(surface.normal_at(u, v))
}

/// Whether the chart is singular at `(u, v)` AND its limit normal
/// depends on the approach direction (a cone apex, a collapsed NURBS
/// edge): such a vertex has no single normal a facet could track. A
/// sphere pole is singular but direction-independent — its limit
/// normal is the radial direction from every meridian.
pub(crate) fn singular_without_unique_normal(surface: &SurfaceGeometry, u: f64, v: f64) -> bool {
    if !chart_singular(surface, u, v) {
        return false;
    }
    let (pu, pv) = (surface.partial_u(u, v), surface.partial_v(u, v));
    // Sweep the collapsed coordinate: every point of the collapsed
    // row is the same 3D point, so its limit normals must agree.
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    let span = |lo: f64, hi: f64| {
        if lo.is_finite() && hi.is_finite() {
            hi - lo
        } else {
            1.0
        }
    };
    let probes: Vec<(f64, f64)> = if pu.length() <= pv.length() {
        let d = span(u0, u1);
        [0.0, 0.25, 0.5, 0.75]
            .iter()
            .map(|f| (u + f * d, v))
            .collect()
    } else {
        let d = span(v0, v1);
        [0.0, 0.25, 0.5, 0.75]
            .iter()
            .map(|f| (u, v + f * d))
            .collect()
    };
    let normals: Vec<Vec3> = probes
        .iter()
        .filter_map(|&(a, b)| limit_normal(surface, a, b))
        .collect();
    normals.len() < probes.len() || normals.iter().any(|n| n.dot(&normals[0]) < 0.9999)
}

/// Re-derive the normals of `mesh`'s chart-singular and out-of-range
/// vertices as limit normals, in the face's facing (the majority sign
/// of its healthy vertices' stored normals against the chart normal).
pub(crate) fn repair_singular_normals(mesh: &mut FaceMesh, surface: &SurfaceGeometry) {
    if matches!(surface, SurfaceGeometry::Plane { .. }) {
        return;
    }
    let suspect: Vec<bool> = mesh
        .vertices
        .iter()
        .map(|v| chart_singular(surface, v.uv.0, v.uv.1) || outside_range(surface, v.uv.0, v.uv.1))
        .collect();
    if !suspect.iter().any(|&s| s) {
        return;
    }
    let mut vote = 0i64;
    for (v, &s) in mesh.vertices.iter().zip(&suspect) {
        if s {
            continue;
        }
        let d = surface.normal_at(v.uv.0, v.uv.1).dot(&v.normal);
        if d > 0.0 {
            vote += 1;
        } else if d < 0.0 {
            vote -= 1;
        }
    }
    let sign = if vote < 0 { -1.0 } else { 1.0 };
    for (v, &s) in mesh.vertices.iter_mut().zip(&suspect) {
        if !s {
            continue;
        }
        if let Some(n) = limit_normal(surface, v.uv.0, v.uv.1) {
            v.normal = Vec3::new(sign * n.x, sign * n.y, sign * n.z);
        }
    }
}

#[cfg(test)]
#[path = "vertex_normal_tests.rs"]
mod tests;
