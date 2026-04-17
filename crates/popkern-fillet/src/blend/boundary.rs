//! Blend region boundary identification.
//!
//! At a vertex where multiple fillets meet, this module computes the closed
//! boundary loop around the vertex. The boundary alternates between fillet-end
//! arcs (the transverse edge of each fillet surface near the vertex) and
//! original face edge segments between adjacent fillets.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use super::boundary_hex::{build_face_edge_bridge_between, reorient_arcs_for_hex_closure};
use crate::fillet_surface::FilletSurfaceResult;

/// Boundary of a vertex blend region.
#[derive(Debug)]
pub struct BlendBoundary {
    /// The vertex being blended.
    pub vertex: VertexHandle,
    /// Position of the vertex (saved before vertex removal by blend insertion).
    pub vertex_position: Point3,
    /// Boundary curve segments (alternating fillet-end arcs and face edge segments).
    pub segments: Vec<BoundarySegment>,
}

/// A single segment of a blend boundary.
#[derive(Debug, Clone)]
pub enum BoundarySegment {
    /// End arc of a fillet surface at the vertex.
    FilletEnd {
        /// Index of the fillet in the global fillet list.
        fillet_index: usize,
        /// The boundary curve segment.
        curve: NurbsCurve,
    },
    /// Edge segment on an original face between two adjacent fillets.
    FaceEdge {
        /// The face this segment lies on.
        face: FaceHandle,
        /// The boundary curve segment.
        curve: NurbsCurve,
    },
}

/// Compute the blend boundary at a vertex where multiple fillets meet.
///
/// Collects the end arcs of each fillet surface near the vertex, plus the
/// original face edge segments between fillets. These form a closed loop
/// around the vertex that defines the blend region.
///
/// # Errors
///
/// Returns an error if any fillet index is out of range or if the vertex
/// position cannot be determined.
pub fn compute_blend_boundary(
    vertex: VertexHandle,
    fillet_indices: &[usize],
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
    tolerance: &Tolerance,
) -> Result<BlendBoundary, String> {
    if fillet_indices.is_empty() {
        return Err("no fillet indices provided for blend boundary".into());
    }

    // Resolve the vertex position.
    let vertex_data = brep
        .vertices
        .get(vertex)
        .ok_or_else(|| format!("vertex {vertex:?} not found in BRep"))?;
    let vertex_pos = vertex_data.position;

    // Build each fillet-end arc and remember its two 3D endpoints.
    let n = fillet_indices.len();
    let mut fillet_arcs: Vec<NurbsCurve> = Vec::with_capacity(n);
    let mut fillet_endpoints: Vec<[Point3; 2]> = Vec::with_capacity(n);
    for &fi in fillet_indices {
        let fillet = fillet_surfaces
            .get(fi)
            .ok_or_else(|| format!("fillet index {fi} out of range"))?;
        let end_curve = extract_fillet_end_near_vertex(fillet, &vertex_pos, tolerance)?;
        let (t0, t1) = end_curve.parameter_range();
        let ep = [end_curve.point_at(t0), end_curve.point_at(t1)];
        fillet_arcs.push(end_curve);
        fillet_endpoints.push(ep);
    }

    // Find the planar faces at this vertex. §5.3 Round 8 C1 face-
    // edge bridges are only emitted when every adjacent face at
    // the vertex is a `Plane` — the curved-support case is a
    // follow-up round (needs surface-isocurve bridge construction).
    let faces_at_vertex = faces_at_vertex_set(brep, vertex);
    let planar_faces: Vec<(FaceHandle, Point3, Vec3)> = faces_at_vertex
        .iter()
        .filter_map(|&fh| {
            let face = brep.faces.get(fh)?;
            let sof = face.surface.as_ref()?;
            if let popkern_core::geometry::surface::SurfaceGeometry::Plane {
                origin, normal, ..
            } = &sof.geometry
            {
                Some((fh, *origin, *normal))
            } else {
                None
            }
        })
        .collect();

    // Emit alternating FilletEnd / FaceEdge segments. For
    // hexagonal closure the bridge between arc[k] and arc[k+1]
    // MUST connect `arc[k].end` to `arc[k+1].start` (not
    // arbitrary endpoint pair) — anything else breaks the
    // closed-loop invariant. Each fillet's arc endpoints come
    // from the trim-curve closest-point construction in
    // `extract_fillet_end_near_vertex`: index 0 is the arc start
    // (UV-side from trim_curve_1's nearest point) and index 1 is
    // the arc end (trim_curve_2's nearest point).
    //
    // `try_reorient_arcs_for_closure` re-orders fillet indexing
    // so consecutive arcs actually share a nearby endpoint via a
    // cube face — the caller's `fillet_indices` order may not be
    // CCW around the vertex.
    let mut segments = Vec::with_capacity(n * 2);
    let all_planar = !planar_faces.is_empty() && planar_faces.len() == faces_at_vertex.len();

    if all_planar && n >= 2 {
        let reorder = reorient_arcs_for_hex_closure(
            fillet_indices,
            &fillet_arcs,
            &fillet_endpoints,
            &planar_faces,
            tolerance.linear().max(1e-6),
        );
        if let Some((order, flip, arc_by_k)) = reorder {
            // `arc_by_k[k]` is the arc to emit at position `k`,
            // already re-oriented; `flip[k]` records whether we
            // reversed arc `order[k]` (so the bridge builder can
            // use the correct endpoints).
            for k in 0..n {
                let fi = fillet_indices[order[k]];
                segments.push(BoundarySegment::FilletEnd {
                    fillet_index: fi,
                    curve: arc_by_k[k].clone(),
                });
                let next_k = (k + 1) % n;
                // Bridge connects current arc's end to next arc's
                // start (post-reorientation endpoints).
                let current_endpoints = if flip[k] {
                    [fillet_endpoints[order[k]][1], fillet_endpoints[order[k]][0]]
                } else {
                    fillet_endpoints[order[k]]
                };
                let next_endpoints = if flip[next_k] {
                    [
                        fillet_endpoints[order[next_k]][1],
                        fillet_endpoints[order[next_k]][0],
                    ]
                } else {
                    fillet_endpoints[order[next_k]]
                };
                let arc_end = current_endpoints[1];
                let next_start = next_endpoints[0];
                if let Some(bridge) = build_face_edge_bridge_between(
                    &arc_end,
                    &next_start,
                    &planar_faces,
                    tolerance.linear().max(1e-6),
                ) {
                    segments.push(bridge);
                }
            }
            return Ok(BlendBoundary {
                vertex,
                vertex_position: vertex_pos,
                segments,
            });
        }
        // Fall through to the FilletEnd-only path if reorientation
        // couldn't determine a valid hex order (e.g., some arcs
        // don't share faces with neighbours — non-planar geometry,
        // degenerate corner, etc.). This preserves pre-Round 8
        // behaviour.
    }

    for k in 0..n {
        let fi = fillet_indices[k];
        segments.push(BoundarySegment::FilletEnd {
            fillet_index: fi,
            curve: fillet_arcs[k].clone(),
        });
    }

    Ok(BlendBoundary {
        vertex,
        vertex_position: vertex_pos,
        segments,
    })
}

/// Collect the set of face handles incident to `vertex` by walking
/// all half-edges originating from it.
fn faces_at_vertex_set(brep: &BRep, vertex: VertexHandle) -> Vec<FaceHandle> {
    let mut seen: Vec<FaceHandle> = Vec::new();
    for he in brep.half_edges.values() {
        if he.origin != vertex {
            continue;
        }
        let Some(loop_data) = brep.loops.get(he.loop_) else {
            continue;
        };
        let face = loop_data.face;
        if !seen.contains(&face) {
            seen.push(face);
        }
    }
    seen
}

/// Try to build a face-edge bridge between two consecutive fillet-
/// end arcs. The bridge is a straight-line NURBS on whichever
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

/// Extract the end curve of a fillet surface closest to a vertex.
///
/// Takes the trim curves of the fillet and finds the endpoint nearest to the
/// vertex position, then builds a short transverse arc between the two trim
/// curve endpoints.
fn extract_fillet_end_near_vertex(
    fillet: &FilletSurfaceResult,
    vertex_pos: &Point3,
    _tolerance: &Tolerance,
) -> Result<NurbsCurve, String> {
    // §5.3 Round 7 C4b.1a — find the closest POINT on each trim
    // curve to the vertex (not just one of the two endpoints).
    // The fillet's trim curves typically span the entire filleted
    // edge including some internal extended-rolling-ball overshoot
    // (e.g., from `x = -1` to `x = 2` on a unit cube at `x ∈ [0, 1]`),
    // so the nearest trim-curve points to a vertex at the edge
    // endpoint are in the INTERIOR of the trim curves, not at
    // their parametric endpoints. The legacy `d_start` / `d_end`
    // selection picked whichever pair was "less far" and would
    // land on the far-overshoot plane, producing boundary arcs
    // geometrically disconnected from the vertex.
    let p1 = closest_point_on_nurbs_curve(&fillet.trim_curve_1, vertex_pos, 64);
    let p2 = closest_point_on_nurbs_curve(&fillet.trim_curve_2, vertex_pos, 64);

    // Project both points onto the fillet surface so we can trace
    // an ARC ALONG THE SURFACE between them (the geometrically
    // correct fillet-end arc) instead of a straight line through
    // space.
    let Some((u1, v1)) = project_point_to_uv(&fillet.surface, &p1) else {
        return Err(format!("failed to project p1 = {p1:?} onto fillet surface"));
    };
    let Some((u2, v2)) = project_point_to_uv(&fillet.surface, &p2) else {
        return Err(format!("failed to project p2 = {p2:?} onto fillet surface"));
    };

    // Periodicity handling: for Cylinder / Cone / Torus the `u`
    // parameter is angular (period `TAU`). Walking linearly from
    // `u1` to `u2` may take the long way around the cylinder; pick
    // the shorter arc direction.
    let (u_periodicity, _) = fillet.surface.periodicity();
    let mut du = u2 - u1;
    if let Some(period) = u_periodicity
        && du.abs() > period * 0.5
    {
        du -= period * du.signum();
    }

    super::boundary_hex::fit_hermite_cubic_arc_on_surface(&fillet.surface, u1, v1, u2, v2, du)
}

/// Find the closest point on `curve` to `target` by sampling at
/// `n_samples + 1` uniform parameter values and picking the
/// minimum-distance sample. Chord-length-accurate to about
/// `curve_length / n_samples`; for our N=3-corner fixtures with
/// `n_samples = 64`, accuracy is sub-millimetre on unit-scale
/// models — well under the rolling-ball radius.
///
/// Could be improved with Newton polishing on the best sample,
/// but the extra precision hasn't been needed in practice.
fn closest_point_on_nurbs_curve(curve: &NurbsCurve, target: &Point3, n_samples: usize) -> Point3 {
    let (t0, t1) = curve.parameter_range();
    let mut best_d = f64::INFINITY;
    let mut best_p = curve.point_at(t0);
    let samples = n_samples.max(1);
    for i in 0..=samples {
        let t = t0 + (t1 - t0) * (i as f64) / (samples as f64);
        let p = curve.point_at(t);
        let d = p.distance_to(target);
        if d < best_d {
            best_d = d;
            best_p = p;
        }
    }
    best_p
}
