//! Arc-length sampling and surface-normal lookup for blend boundaries.
//!
//! Extracted from `boundary.rs` to keep the parent file under the project's
//! 500-line per-file structural cap. The types and functions here remain
//! re-exported from `super::boundary` so existing call sites are unaffected.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::BRep;

use super::boundary::{BlendBoundary, BoundarySegment};
use crate::fillet_surface::FilletSurfaceResult;

/// A sampled point on the blend boundary, with the adjacent surface normal.
#[derive(Debug, Clone)]
pub struct BoundarySample {
    /// 3D position on the boundary curve.
    pub position: Point3,
    /// Surface normal of the adjacent fillet or face surface at this position.
    pub normal: Vec3,
}

/// Sample the blend boundary at uniform arc-length intervals.
///
/// For each sample, projects the position onto the adjacent surface (torus for
/// `FilletEnd` segments, plane for `FaceEdge` segments) and evaluates the
/// surface normal. Falls back to `(0, 0, 1)` if UV projection fails.
///
/// Returns an empty `Vec` if `n_samples` is zero or the boundary has no segments.
pub fn sample_boundary_by_arc_length(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
    n_samples: usize,
) -> Vec<BoundarySample> {
    if n_samples == 0 || boundary.segments.is_empty() {
        return Vec::new();
    }

    const ARC_SUBDIVISIONS: usize = 32;

    // --- Step 1: compute per-segment arc lengths ---
    let mut seg_lengths: Vec<f64> = Vec::with_capacity(boundary.segments.len());

    for seg in &boundary.segments {
        let curve = match seg {
            BoundarySegment::FilletEnd { curve, .. } => curve,
            BoundarySegment::FaceEdge { curve, .. } => curve,
        };
        let (t_min, t_max) = curve.parameter_range();
        let mut length = 0.0;
        let mut prev = curve.point_at(t_min);
        for i in 1..=ARC_SUBDIVISIONS {
            let t = t_min + (t_max - t_min) * (i as f64) / ARC_SUBDIVISIONS as f64;
            let p = curve.point_at(t);
            let dx = p.x - prev.x;
            let dy = p.y - prev.y;
            let dz = p.z - prev.z;
            length += (dx * dx + dy * dy + dz * dz).sqrt();
            prev = p;
        }
        seg_lengths.push(length);
    }

    let total_length: f64 = seg_lengths.iter().sum();
    if total_length < f64::EPSILON {
        return Vec::new();
    }

    // --- Step 2: build cumulative breakpoints ---
    let mut breakpoints: Vec<f64> = Vec::with_capacity(boundary.segments.len() + 1);
    breakpoints.push(0.0);
    let mut acc = 0.0;
    for &l in &seg_lengths {
        acc += l;
        breakpoints.push(acc);
    }

    // --- Step 3: sample at uniform arc-length intervals ---
    let mut samples = Vec::with_capacity(n_samples);

    for i in 0..n_samples {
        let target = total_length * (i as f64) / (n_samples as f64);

        // Find which segment contains this arc-length position.
        let seg_idx = breakpoints
            .windows(2)
            .position(|w| target >= w[0] && target < w[1])
            .unwrap_or(boundary.segments.len().saturating_sub(1));

        let seg = &boundary.segments[seg_idx];
        let curve = match seg {
            BoundarySegment::FilletEnd { curve, .. } => curve,
            BoundarySegment::FaceEdge { curve, .. } => curve,
        };

        // Interpolate t within the segment proportionally.
        let seg_start = breakpoints[seg_idx];
        let seg_len = seg_lengths[seg_idx];
        let local_frac = if seg_len > f64::EPSILON {
            (target - seg_start) / seg_len
        } else {
            0.0
        };

        let (t_min, t_max) = curve.parameter_range();
        let t = t_min + (t_max - t_min) * local_frac.clamp(0.0, 1.0);
        let position = curve.point_at(t);

        // --- Step 4: look up the adjacent surface normal ---
        let normal = surface_normal_at_point(seg, &position, fillet_surfaces, brep);

        samples.push(BoundarySample { position, normal });
    }

    samples
}

/// Evaluate the surface normal of the surface adjacent to a boundary segment
/// at the given 3D position.
///
/// Uses UV projection to find the parameter-space location, then queries the
/// surface normal. Returns `(0, 0, 1)` if the surface or projection is unavailable.
fn surface_normal_at_point(
    seg: &BoundarySegment,
    position: &Point3,
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
) -> Vec3 {
    let fallback = Vec3::new(0.0, 0.0, 1.0);

    match seg {
        BoundarySegment::FilletEnd { fillet_index, .. } => {
            let Some(fillet) = fillet_surfaces.get(*fillet_index) else {
                return fallback;
            };
            let Some((u, v)) = project_point_to_uv(&fillet.surface, position) else {
                return fallback;
            };
            let n = fillet.surface.normal_at(u, v);
            n.normalized().unwrap_or(fallback)
        }
        BoundarySegment::FaceEdge { face, .. } => {
            let Some(face_data) = brep.faces.get(*face) else {
                return fallback;
            };
            let Some(surf_on_face) = &face_data.surface else {
                return fallback;
            };
            let surface = &surf_on_face.geometry;
            let Some((u, v)) = project_point_to_uv(surface, position) else {
                return fallback;
            };
            let n = surface.normal_at(u, v);
            n.normalized().unwrap_or(fallback)
        }
    }
}
