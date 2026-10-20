//! Canonical-station trim-loop sampling (PK-299).
//!
//! A face that declines boundary-first tessellates on the uniform-grid
//! grid, whose trim polylines were sampled adaptively from each
//! half-edge's PCURVE — its boundary vertices then lie on chords of
//! that private polyline, up to the UV chord-height off the true
//! edge curve, which is beyond the crack zipper's 2e-4 matching band.
//! Every shared edge between a fallback face and a served
//! boundary-first neighbour therefore reads open (measured on the
//! PK-295 fixture mid-campaign: 57 open mesh edges along a spiric).
//!
//! This module samples a trim loop at the CANONICAL per-edge stations
//! (`sample_edges` — the same samples every boundary-first face
//! uses), pairing each canonical 3-D point with the half-edge
//! pcurve's UV at the co-parameterized station. The uniform-grid face's
//! trim polyline then shares every vertex — bit-identical — with the
//! neighbours' chains, clip-introduced T-vertices land exactly ON
//! the neighbours' segments, and the tolerance-scoped zipper closes
//! the weld.
//!
//! Verification gate: each station's UV must round-trip through the
//! surface to its canonical point within a band; any violation
//! declines the WHOLE loop back to the caller's adaptive path (a
//! pcurve authored off-contract must not silently warp the trim).

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, EdgeHandle, LoopHandle};

use crate::boundary_first::EdgeSamples;

/// A sampled loop: paired UV points and canonical 3-D points.
pub(crate) type SampledLoop = (Vec<(f64, f64)>, Vec<Point3>);

/// Sample `loop_handle` at canonical edge stations. Returns
/// `(uv_points, points_3d)` or `None` when any half-edge lacks a
/// pcurve, an edge lacks canonical samples, or the co-param
/// verification fails.
pub(crate) fn sample_loop_canonical(
    brep: &BRep,
    loop_handle: LoopHandle,
    surface: &SurfaceGeometry,
    edge_samples: &HashMap<EdgeHandle, EdgeSamples>,
    chord_height: f64,
) -> Option<SampledLoop> {
    let trace = popkern_core::trace::trace_enabled("trim");
    let loop_ = brep.loops.get(loop_handle)?;
    let start_he = loop_.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return None;
    }
    // The UV image may legitimately sit a chord-height off the exact
    // surface point (fitted pcurves), plus the canonical samples'
    // own 3-D chord tolerance.
    let band = (chord_height * 10.0).max(1.0e-4);

    let mut uv_points = Vec::new();
    let mut points_3d = Vec::new();
    let mut current = start_he;
    loop {
        let he = &brep.half_edges[current];
        let pc = he.pcurve.as_ref();
        let Some(es) = edge_samples.get(&he.edge) else {
            if trace {
                eprintln!(
                    "  [trim-canonical] loop {loop_handle:?}: edge {:?} has no canonical samples",
                    he.edge
                );
            }
            return None;
        };
        let n = es.points.len();
        if n < 2 {
            return None;
        }
        let origin = brep.vertices.get(he.origin)?.position;
        // Traversal direction along the canonical stations. For a
        // closed rim (coincident endpoints) direction is decided by
        // the verification probe below, not the endpoints.
        let d_first = origin.distance_to(&es.points[0]);
        let d_last = origin.distance_to(&es.points[n - 1]);
        let endpoints_distinct = es.points[0].distance_to(&es.points[n - 1]) > 1.0e-9;
        // UV source: the co-parameterized pcurve when present, else
        // the surface's closed-form projection of the canonical 3-D
        // point (the missing-pcurve decline class); both are held to
        // the same round-trip band.
        let uv_at = |k: usize, forward: bool| -> Option<(f64, f64)> {
            match pc {
                Some(pc) => {
                    let f = es.fractions[k];
                    let dir_f = if forward { f } else { 1.0 - f };
                    let uv = pc.point_at(pc.t_start + dir_f * (pc.t_end - pc.t_start));
                    Some((uv.x, uv.y))
                }
                None => popkern_core::geometry::uv_projection::project_point_to_uv(
                    surface,
                    &es.points[k],
                ),
            }
        };
        let deviation = |forward: bool| -> Option<f64> {
            let mut worst = 0.0_f64;
            for k in [0, n / 2, n - 1] {
                let (u, v) = uv_at(k, forward)?;
                let img = surface.point_at(u, v);
                worst = worst.max(img.distance_to(&es.points[k]));
            }
            Some(worst)
        };
        let forward = if endpoints_distinct {
            d_first <= d_last
        } else {
            let (df, dr) = (
                deviation(true).unwrap_or(f64::MAX),
                deviation(false).unwrap_or(f64::MAX),
            );
            df <= dr
        };
        if deviation(forward)? > band {
            if trace {
                eprintln!(
                    "  [trim-canonical] loop {loop_handle:?}: edge {:?} co-param deviation {:?} > band {band:.6}",
                    he.edge,
                    deviation(forward)
                );
            }
            return None;
        }
        // Emit stations in traversal order; corner duplicates are
        // cleaned by the caller's dedup, matching the adaptive path.
        let order: Box<dyn Iterator<Item = usize>> = if forward {
            Box::new(0..n)
        } else {
            Box::new((0..n).rev())
        };
        for k in order {
            let (u, v) = uv_at(k, forward)?;
            let img = surface.point_at(u, v);
            if img.distance_to(&es.points[k]) > band {
                if trace {
                    eprintln!(
                        "  [trim-canonical] loop {loop_handle:?}: edge {:?} station {k} round-trip {:.6} > band {band:.6}",
                        he.edge,
                        img.distance_to(&es.points[k])
                    );
                }
                return None;
            }
            uv_points.push((u, v));
            points_3d.push(es.points[k]);
        }
        current = he.next;
        if current == start_he {
            break;
        }
    }
    if uv_points.len() < 3 {
        return None;
    }
    if trace {
        eprintln!(
            "  [trim-canonical] loop {loop_handle:?}: canonical sampling SERVED {} stations",
            uv_points.len()
        );
    }
    Some((uv_points, points_3d))
}
