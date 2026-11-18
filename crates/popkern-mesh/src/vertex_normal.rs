//! Mesh vertex normals at chart singularities, from the chart service.
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
//! The normal is the chart service's (`chart::ops::chart_normal`): at a
//! chart singularity the singular row's normal, which every approach
//! shares; at a geometric singularity (a cone's apex) none, the limit
//! depending on the approach; elsewhere `∂u × ∂v` with a coordinate
//! that has no period clamped into the domain. A point is at a singular
//! row within rounding of its point in space (`SINGULAR_EXACT` of its
//! coordinate scale), never by a ratio of partials.
//! [`repair_singular_normals`] re-derives the normal of every vertex at
//! a singular row or outside the natural parameter range; healthy
//! vertices keep their normals bit-for-bit.

use popkern_core::geometry::chart::ops::{SINGULAR_EXACT, chart_normal, evaluate};
use popkern_core::geometry::chart::{ChartDescriptor, SingularKind};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};

use crate::mesh::FaceMesh;

/// Rounding at `p`: a point within this of a singular row's point is on it.
fn exact(p: &Point3) -> f64 {
    SINGULAR_EXACT * p.x.abs().max(p.y.abs()).max(p.z.abs()).max(1.0)
}

/// The singular row whose normal `p` takes, if any: on its point to
/// rounding, or (a chart singularity) within 1e-5 of its collapse radius,
/// where the raw `S_u × S_v` is a cancellation and the service answers the
/// row's normal (as `chart_normal` reads it; a vanishing fillet strip's
/// near-tip vertices kept raw normals and display refinement chased them
/// to its growth budget).
fn row_at<'a>(
    chart: &'a ChartDescriptor,
    p: &Point3,
) -> Option<&'a popkern_core::geometry::chart::SingularRow> {
    let e = exact(p);
    chart.singular_rows.iter().find(|row| {
        let reach = match row.kind {
            SingularKind::Chart => e.max(1e-5 * row.distortion),
            _ => e,
        };
        row.point.distance_to(p) <= reach
    })
}

fn outside_range(chart: &ChartDescriptor, u: f64, v: f64) -> bool {
    let ((u0, u1), (v0, v1)) = chart.domain;
    let free = |axis| chart.period(axis).is_none();
    (free(popkern_core::geometry::chart::Axis::U) && (u < u0 || u > u1))
        || (free(popkern_core::geometry::chart::Axis::V) && (v < v0 || v > v1))
}

/// The surface's unit normal at `uv` (its point `p`) in the chart's own
/// sense, from the chart service: the singular row's normal at a chart
/// singularity, `None` at an apex or a degenerate chart point.
pub fn limit_normal(
    chart: &ChartDescriptor,
    surface: &SurfaceGeometry,
    uv: (f64, f64),
    p: Point3,
) -> Option<Vec3> {
    chart_normal(chart, surface, uv, p, exact(&p))
}

/// [`limit_normal`] at `uv`, its point evaluated there.
pub fn limit_normal_at(
    chart: &ChartDescriptor,
    surface: &SurfaceGeometry,
    uv: (f64, f64),
) -> Option<Vec3> {
    let p = evaluate(chart, surface, uv);
    if singular_without_unique_normal(chart, &p) {
        return None;
    }
    limit_normal(chart, surface, uv, p)
}

/// Whether `p` lies on a geometric singularity (a cone's apex): a vertex
/// there has no single normal a facet could track. A sphere's pole is a
/// chart singularity with one normal (the radial direction from every
/// meridian).
///
/// Within an apex's neighbourhood too, not only on its point: there
/// refinement can never meet a facet bound (every facet touching the tip
/// spans its sweep, however small), so a vertex within 1e-5 of the apex's
/// collapse radius (where its speed along the row is a millionth of its
/// speed away from it, as the old partial-ratio test read it) is exempt.
/// On the point alone, refinement chased a vanishing fillet strip's tip to
/// the growth budget (the tee saddle 3.4x slower).
pub(crate) fn singular_without_unique_normal(chart: &ChartDescriptor, p: &Point3) -> bool {
    let e = exact(p);
    chart.singular_rows.iter().any(|row| {
        row.kind != SingularKind::Chart && row.point.distance_to(p) <= e.max(1e-5 * row.distortion)
    })
}

/// Re-derive the normals of `mesh`'s vertices at a singular row or
/// outside the natural range as the service's normal, in the face's
/// facing (the majority sign of its healthy vertices' stored normals
/// against the chart normal).
pub(crate) fn repair_singular_normals(
    mesh: &mut FaceMesh,
    surface: &SurfaceGeometry,
    chart: &ChartDescriptor,
) {
    if matches!(surface, SurfaceGeometry::Plane { .. }) {
        return;
    }
    let suspect: Vec<bool> = mesh
        .vertices
        .iter()
        .map(|v| row_at(chart, &v.position).is_some() || outside_range(chart, v.uv.0, v.uv.1))
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
        if let Some(n) = limit_normal(chart, surface, v.uv, v.position) {
            v.normal = Vec3::new(sign * n.x, sign * n.y, sign * n.z);
        }
    }
}

#[cfg(test)]
#[path = "vertex_normal_tests.rs"]
mod tests;
