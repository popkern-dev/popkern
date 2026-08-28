//! Vertex blend surface construction (mitred and smooth modes).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ContinuityLevel, CornerMethod, VertexMode};
use crate::vertex_classifier::VertexBlendClassification;

use super::boundary::{BlendBoundary, BoundarySegment};

/// Result of vertex blend computation.
#[derive(Debug, Clone)]
pub struct BlendSurfaceResult {
    /// The constructed blend surface.
    pub surface: SurfaceGeometry,
    /// Boundary curves of the blend region.
    pub boundary_curves: Vec<NurbsCurve>,
    /// §5.3 method used to build `surface`. Emitted on
    /// [`crate::types::CornerPatchInfo::method`] so UI consumers
    /// know which construction path fired.
    pub method: CornerMethod,
    /// §7 continuity the construction delivers with the incident
    /// blend strips. G1 is the baseline for Gregory; G2 is reached
    /// when every incident fillet shares one radius and the
    /// curvature-matching interior ring exactly matches 1/r at
    /// each boundary sample.
    pub continuity: ContinuityLevel,
}

/// Build a vertex blend surface.
///
/// For [`VertexMode::Mitred`], fits a least-squares plane through the boundary
/// points and returns a planar surface.
///
/// For [`VertexMode::Smooth`] with two fillets, attempts a spherical cap
/// tangent to both fillet surfaces. If G1 continuity cannot be achieved,
/// falls back to NURBS patch fitting. For three or more fillets, degrades
/// gracefully to the mitred approach.
///
/// # Errors
///
/// Returns an error if the boundary has no segments or if the geometry is
/// degenerate (e.g., all boundary points coincide).
pub fn build_blend_surface(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
    vertex_mode: VertexMode,
    classification: VertexBlendClassification,
    tolerance: &Tolerance,
) -> Result<BlendSurfaceResult, String> {
    if boundary.segments.is_empty() {
        return Err("blend boundary has no segments".into());
    }

    // Collect all boundary points.
    let boundary_curves = collect_boundary_curves(boundary);
    let points = sample_boundary_points(&boundary_curves, 10);

    if points.len() < 3 {
        return Err("need at least 3 boundary points for blend surface".into());
    }

    match vertex_mode {
        VertexMode::Mitred => build_mitred_blend(&points, boundary_curves),
        VertexMode::Smooth => build_smooth_blend(
            boundary,
            fillet_surfaces,
            &points,
            boundary_curves,
            classification,
            tolerance,
            brep,
        ),
    }
}

/// Build a mitred (planar) blend surface.
///
/// Fits a least-squares plane through the boundary sample points.
fn build_mitred_blend(
    points: &[Point3],
    boundary_curves: Vec<NurbsCurve>,
) -> Result<BlendSurfaceResult, String> {
    let (origin, normal) = fit_plane_least_squares(points)?;

    // Construct orthonormal frame.
    let u_axis = pick_perpendicular(&normal);
    let v_axis = normal.cross(&u_axis);
    let v_axis = v_axis
        .normalized()
        .ok_or_else(|| "degenerate v_axis in blend plane".to_string())?;

    let surface = SurfaceGeometry::Plane {
        origin,
        normal,
        u_axis,
        v_axis,
    };

    Ok(BlendSurfaceResult {
        surface,
        boundary_curves,
        method: {
            if std::env::var_os("PK_SPHEREFIT_TRACE").is_some() {
                eprintln!("[corner-method] mitred");
            }
            CornerMethod::Mitred
        },
        continuity: ContinuityLevel::G0,
    })
}

/// Build a smooth blend surface.
///
/// Fallback chain: sphere → polar NURBS → mitred plane.
///
/// For two-fillet corners, the torus axes define a valid sphere bisector,
/// so the sphere is tried first. If it fails, polar NURBS is tried as an
/// intermediate fallback before the mitred plane.
///
/// For three or more fillets (e.g., a cube corner), the geometry is more
/// complex. Polar NURBS is tried first; the mitred plane is the final
/// fallback.
fn build_smooth_blend(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    points: &[Point3],
    boundary_curves: Vec<NurbsCurve>,
    classification: VertexBlendClassification,
    tolerance: &Tolerance,
    brep: &BRep,
) -> Result<BlendSurfaceResult, String> {
    let n_fillets = boundary
        .segments
        .iter()
        .filter(|s| matches!(s, BoundarySegment::FilletEnd { .. }))
        .count();

    // Fallback chain: (0) tangent-plane cap for Tangent classification
    // [§5.3 item #7]; (1) rolling-ball sphere for SameSide* [§5.3 C6];
    // (2) Gregory G1 polar-NURBS; (3) mitred plane G0.

    if matches!(classification, VertexBlendClassification::Tangent)
        && let Some(result) = super::tangent_corner::try_tangent_plane_blend(
            boundary,
            fillet_surfaces,
            points,
            &boundary_curves,
        )
    {
        return Ok(result);
    }

    let sphere_allowed = matches!(
        classification,
        VertexBlendClassification::SameSideConvex | VertexBlendClassification::SameSideConcave,
    );
    if sphere_allowed
        && let Some(result) = try_spherical_blend(
            boundary,
            fillet_surfaces,
            points,
            &boundary_curves,
            tolerance,
        )
    {
        return Ok(result);
    }

    // §5.3 polar-NURBS fallback — works for any N ≥ 2.
    // Rational-Gregory dispatch happens upstream in
    // `blend_insert::insert_vertex_blend` before this function.
    let _ = n_fillets;
    if let Some(result) =
        build_gregory_g1_blend(boundary, fillet_surfaces, brep, boundary_curves.clone())
    {
        return Ok(result);
    }

    build_mitred_blend(points, boundary_curves)
}

/// Build a polar NURBS blend surface.
///
/// Uses a polar parameterization centered at the vertex position.
/// The boundary points form the outer ring, and the vertex position
/// is the degenerate pole. Interior control points are adjusted for
/// G1 continuity with adjacent surfaces.
/// §5.3 Gregory G1/G2 N-sided corner patch — polar NURBS with a
/// curvature-matching interior ring stepped inward along each
/// boundary sample's surface normal. The degenerate pole sits at
/// the vertex position; the outer ring traces the N boundary
/// curves at arc-length-uniform spacing.
///
/// **Interior step (§5.3).** For a degree-2 Bezier in
/// the radial direction with control points `[pole P_0, interior
/// P_1, boundary P_2]`, the surface curvature at the boundary is
/// `κ = |n̂ × (v_pos − P_2)| / (2 · step²)` where `step = |P_2 − P_1|`.
/// Solving for the step that makes `κ = 1/r` (matching the incident
/// fillet strip's arc curvature) gives
/// `step = sqrt(r · |n̂ × (v_pos − P_2)| / 2)`. Per-sample step
/// delivers G2 at each boundary point when the incident fillets
/// share one radius; for mixed-radii corners we use the first
/// fillet's radius as the target and report G1 (the deviation at
/// mismatched samples is bounded by the radius ratio).
///
/// For `N = 2` this is the classical two-fillet smooth corner.
/// For `N ≥ 3` it's a generalisation of the same construction —
/// not a true rational Gregory patch, but a commercial-grade
/// single-face approximation that delivers G2 in the equal-radius
/// common case and G1 otherwise. True degree-3×3 rational Gregory
/// for `N = 3` / `N = 4` is a follow-up round.
fn build_gregory_g1_blend(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
    boundary_curves: Vec<NurbsCurve>,
) -> Option<BlendSurfaceResult> {
    use super::boundary::sample_boundary_by_arc_length;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let n_samples = 24;
    let samples = sample_boundary_by_arc_length(boundary, fillet_surfaces, brep, n_samples);

    if samples.len() < 3 {
        return None;
    }

    let v_pos = boundary.vertex_position;
    let m = samples.len();

    // Collect the incident fillet minor radii to decide the
    // curvature target and whether this corner qualifies as G2.
    let radii: Vec<f64> = fillet_surfaces
        .iter()
        .filter_map(|fr| match &fr.surface {
            SurfaceGeometry::Torus { minor_radius, .. } => Some(*minor_radius),
            SurfaceGeometry::Cylinder { radius, .. } => Some(*radius),
            _ => None,
        })
        .collect();
    let target_radius = radii.first().copied().unwrap_or(0.1);
    // G2 upgrade requires a confirmed uniform curvature target
    // from the incident analytic fillet surfaces. If the strips
    // are NURBS (G2 edge fillets produce NURBS, not Torus /
    // Cylinder), we cannot verify the target radius matches the
    // strip's actual boundary curvature at each sample — honest
    // answer is G1 at these corners. Likewise if N > radii.len()
    // (some incident strips are analytic, some NURBS): we don't
    // have a full set of radii to verify, so stay at G1.
    let radii_uniform = if radii.len() >= fillet_surfaces.len() && !radii.is_empty() {
        let (lo, hi) = radii.iter().fold((f64::INFINITY, 0.0_f64), |(lo, hi), &r| {
            (lo.min(r), hi.max(r))
        });
        (hi - lo) < 1e-6
    } else {
        false
    };

    // Build the num_u × num_v polar control point grid.
    // num_u = m samples around the boundary (u direction).
    // num_v = 3 layers in the radial direction (v direction):
    //   v=0: degenerate pole (all at vertex position)
    //   v=0.5: G2-adjusted interior ring (step per sample)
    //   v=1: boundary sample positions
    let num_u = m;
    let num_v = 3usize;
    let degree_u = 3usize.min(m - 1);
    let degree_v = 2usize;

    let mut control_points = Vec::with_capacity(num_u * num_v);
    for sample in &samples {
        // v=0: pole
        control_points.push(v_pos);
        // v=1: G2 interior — step = sqrt(r · |n̂ × (v_pos − P_2)| / 2).
        // Falls back to the prior radial_dist/3 heuristic if the
        // cross product vanishes (sample.normal aligned with the
        // sample-to-vertex direction — a degenerate edge case).
        let n_cross_r = sample.normal.cross(&(v_pos - sample.position));
        let k_mag = n_cross_r.length();
        let step = if k_mag > 1e-12 && target_radius > 0.0 {
            (target_radius * k_mag / 2.0).sqrt()
        } else {
            sample.position.distance_to(&v_pos) / 3.0
        };
        control_points.push(Point3::new(
            sample.position.x - step * sample.normal.x,
            sample.position.y - step * sample.normal.y,
            sample.position.z - step * sample.normal.z,
        ));
        // v=2: boundary
        control_points.push(sample.position);
    }

    // Clamped uniform knot vector in u: num_u + degree_u + 1 elements,
    // with degree_u+1 repeated values at each end and uniform interior knots.
    let n_interior_u = num_u - degree_u - 1;
    let mut knots_u: Vec<f64> = Vec::with_capacity(num_u + degree_u + 1);
    knots_u.extend(std::iter::repeat_n(0.0, degree_u + 1));
    for i in 1..=n_interior_u {
        knots_u.push(i as f64 / (n_interior_u + 1) as f64);
    }
    knots_u.extend(std::iter::repeat_n(1.0, degree_u + 1));

    // Clamped knot vector in v (degree 2, 3 control points): [0,0,0,1,1,1]
    let knots_v = vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let weights = vec![1.0; num_u * num_v];

    let nurbs = popkern_core::geometry::NurbsSurface {
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        control_points,
        num_u,
        num_v,
        weights,
    };

    // G2 when all incident fillets share one radius — the
    // curvature-matching step exactly matches 1/r at each
    // sample. Otherwise the target radius is wrong at mismatched
    // samples and we report G1.
    let continuity = if radii_uniform {
        ContinuityLevel::G2
    } else {
        ContinuityLevel::G1
    };
    Some(BlendSurfaceResult {
        surface: SurfaceGeometry::Nurbs(nurbs),
        boundary_curves,
        method: {
            if std::env::var_os("PK_SPHEREFIT_TRACE").is_some() {
                eprintln!("[corner-method] gregory");
            }
            CornerMethod::Gregory
        },
        continuity,
    })
}

/// Try to build a NURBS patch blend surface.
///
/// Creates a simple bilinear or Coons-style patch through the boundary
/// points. Returns `None` if the construction fails.
#[allow(dead_code)]
fn try_nurbs_blend(
    points: &[Point3],
    boundary_curves: &[NurbsCurve],
) -> Option<BlendSurfaceResult> {
    if points.len() < 4 {
        return None;
    }

    // Build a simple bilinear patch from the first 4 boundary points.
    let n_u = 2;
    let n_v = 2;
    let control_points = vec![points[0], points[1], points[2], points[3]];
    let weights = vec![1.0; 4];

    let nurbs = NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points,
        num_u: n_u,
        num_v: n_v,
        weights,
    };

    Some(BlendSurfaceResult {
        surface: SurfaceGeometry::Nurbs(nurbs),
        boundary_curves: boundary_curves.to_vec(),
        method: {
            if std::env::var_os("PK_SPHEREFIT_TRACE").is_some() {
                eprintln!("[corner-method] polar");
            }
            CornerMethod::PolarNurbs
        },
        continuity: ContinuityLevel::G0,
    })
}

/// Collect boundary curves from blend boundary segments.
fn collect_boundary_curves(boundary: &BlendBoundary) -> Vec<NurbsCurve> {
    boundary
        .segments
        .iter()
        .map(|seg| match seg {
            BoundarySegment::FilletEnd { curve, .. } => curve.clone(),
            BoundarySegment::FaceEdge { curve, .. } => curve.clone(),
        })
        .collect()
}

/// Sample points from a set of boundary curves.
fn sample_boundary_points(curves: &[NurbsCurve], samples_per_curve: usize) -> Vec<Point3> {
    let mut points = Vec::with_capacity(curves.len() * samples_per_curve);

    for curve in curves {
        let (t_min, t_max) = curve.parameter_range();
        for i in 0..samples_per_curve {
            let t = t_min + (t_max - t_min) * (i as f64) / ((samples_per_curve - 1).max(1) as f64);
            points.push(curve.point_at(t));
        }
    }

    points
}

use super::surface_math::{fit_plane_least_squares, pick_perpendicular};
use super::surface_sphere::try_spherical_blend;

#[cfg(test)]
#[path = "surface_tests.rs"]
mod tests;
