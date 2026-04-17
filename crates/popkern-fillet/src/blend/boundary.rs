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

    let mut segments = Vec::with_capacity(fillet_indices.len() * 2);

    for &fi in fillet_indices {
        let fillet = fillet_surfaces
            .get(fi)
            .ok_or_else(|| format!("fillet index {fi} out of range"))?;

        // Extract the fillet end arc near the vertex.
        let end_curve = extract_fillet_end_near_vertex(fillet, &vertex_pos, tolerance)?;

        segments.push(BoundarySegment::FilletEnd {
            fillet_index: fi,
            curve: end_curve,
        });
    }

    // For a complete boundary we would also insert face-edge segments between
    // consecutive fillet ends. This requires walking the topology around the
    // vertex. For now, the fillet-end segments alone define the blend region;
    // face-edge segments can be added when topology walking is available.

    Ok(BlendBoundary {
        vertex,
        vertex_position: vertex_pos,
        segments,
    })
}

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

    // Sample the arc on the surface. Five samples is sufficient
    // for a chord-length-parameterised NURBS fit to approximate
    // a quarter-circle to ~1% radius; more samples don't improve
    // the degree-1 NURBS `fit_nurbs_through_points` emits.
    let n_samples = 5;
    let mut samples = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let t = (i as f64) / (n_samples - 1) as f64;
        let u = u1 + t * du;
        let v = v1 + t * (v2 - v1);
        samples.push(fillet.surface.point_at(u, v));
    }

    crate::center_curve::fit_nurbs_through_points(&samples)
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

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_core::geometry::{Point3, Vec3};
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn test_tolerance() -> Tolerance {
        Tolerance::default()
            .with_linear(1e-3)
            .with_angular(1e-3)
            .with_parametric(1e-6)
    }

    fn make_fillet_result(
        start1: Point3,
        end1: Point3,
        start2: Point3,
        end2: Point3,
    ) -> FilletSurfaceResult {
        let trim_1 = crate::center_curve::fit_nurbs_through_points(&[start1, end1])
            .expect("fit should succeed");
        let trim_2 = crate::center_curve::fit_nurbs_through_points(&[start2, end2])
            .expect("fit should succeed");

        FilletSurfaceResult {
            surface: SurfaceGeometry::Plane {
                origin: Point3::new(0.0, 0.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
            },
            trim_curve_1: trim_1,
            trim_curve_2: trim_2,
            continuity_downgrade: None,
        }
    }

    #[test]
    fn blend_boundary_two_fillets() {
        let tol = test_tolerance();
        let brep = BRep::unit_cube(tol).expect("cube should build");

        // Pick the first vertex.
        let vertex = brep.vertices.keys().next().expect("cube has vertices");
        let vertex_pos = brep.vertices.get(vertex).expect("vertex exists").position;

        // Create two mock fillet results near the vertex.
        let f1 = make_fillet_result(
            vertex_pos,
            Point3::new(vertex_pos.x + 1.0, vertex_pos.y, vertex_pos.z),
            Point3::new(vertex_pos.x, vertex_pos.y + 0.1, vertex_pos.z),
            Point3::new(vertex_pos.x + 1.0, vertex_pos.y + 0.1, vertex_pos.z),
        );
        let f2 = make_fillet_result(
            vertex_pos,
            Point3::new(vertex_pos.x, vertex_pos.y + 1.0, vertex_pos.z),
            Point3::new(vertex_pos.x + 0.1, vertex_pos.y, vertex_pos.z),
            Point3::new(vertex_pos.x + 0.1, vertex_pos.y + 1.0, vertex_pos.z),
        );

        let fillets = vec![f1, f2];
        let result = compute_blend_boundary(vertex, &[0, 1], &fillets, &brep, &tol);
        assert!(result.is_ok(), "blend boundary failed: {:?}", result.err());

        let boundary = result.expect("already checked");
        assert_eq!(
            boundary.segments.len(),
            2,
            "should have 2 fillet-end segments"
        );

        // Verify the boundary is around the right vertex.
        assert_eq!(boundary.vertex, vertex);

        // Each segment should be a FilletEnd.
        for seg in &boundary.segments {
            assert!(
                matches!(seg, BoundarySegment::FilletEnd { .. }),
                "expected FilletEnd segment"
            );
        }
    }

    #[test]
    fn blend_boundary_empty_fillets() {
        let tol = test_tolerance();
        let brep = BRep::unit_cube(tol).expect("cube should build");
        let vertex = brep.vertices.keys().next().expect("cube has vertices");

        let result = compute_blend_boundary(vertex, &[], &[], &brep, &tol);
        assert!(result.is_err(), "should fail with no fillet indices");
    }

    #[test]
    fn blend_boundary_invalid_fillet_index() {
        let tol = test_tolerance();
        let brep = BRep::unit_cube(tol).expect("cube should build");
        let vertex = brep.vertices.keys().next().expect("cube has vertices");

        let result = compute_blend_boundary(vertex, &[99], &[], &brep, &tol);
        assert!(result.is_err(), "should fail with out-of-range index");
    }
}
