//! §5.3 item #11 — reflection-line regression diagnostics for
//! blend faces.
//!
//! Reflection lines are the Class-A surface-quality tool: a
//! parallel-stripe light source is reflected off the surface, and
//! the stripe boundaries trace lines whose smoothness reveals G2
//! continuity defects that curvature-deviation metrics miss.
//!
//! * [`reflection_lines_on_face`] — compute all reflection lines
//!   on a single face (wraps
//!   [`popkern_core::healing::compute_reflection_lines`]).
//! * [`reflection_line_boundary_kink`] — measure the max
//!   tangent-angle discontinuity between reflection lines on two
//!   adjacent faces at their shared boundary.
//! * [`blend_reflection_line_kink_max`] — the one-call Class-A
//!   continuity metric for a blend face against all its adjacent
//!   support faces.
//!
//! Companion to [`crate::highlight_line`], which uses specular
//! highlight peaks (N·H ≈ 1 contours) rather than reflected
//! stripes. Highlight lines concentrate near a point light and are
//! cheaper to compute; reflection lines cover the surface more
//! uniformly and are more sensitive to subtle G2 breaks.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::healing::{ReflectionLine, compute_reflection_lines};
use popkern_core::topology::{BRep, FaceHandle};

use crate::geometry_compute::face_param_range;

/// Result of running reflection-line analysis on a single face.
#[derive(Debug, Clone)]
pub struct ReflectionLineReport {
    /// The face that was analysed.
    pub face: FaceHandle,
    /// Grid density used for the sampling.
    pub grid_density: usize,
    /// Raw reflection-line polylines.
    pub lines: Vec<ReflectionLine>,
}

/// Compute reflection lines on a single face of the BRep.
///
/// Resolves the face's [`SurfaceGeometry`] and a bounded
/// `(u_range, v_range)` parameter box from the face's outer-loop
/// vertices, then delegates to
/// [`compute_reflection_lines`].
///
/// Returns an empty report on faces without surface geometry or
/// stale handles — not an error (faces may intentionally have no
/// surface during construction).
pub fn reflection_lines_on_face(
    brep: &BRep,
    face: FaceHandle,
    view_dir: Vec3,
    stripe_dir: Vec3,
    stripe_spacing: f64,
    grid_density: usize,
) -> ReflectionLineReport {
    let face_ref = match brep.faces.get(face) {
        Some(f) => f,
        None => {
            return ReflectionLineReport {
                face,
                grid_density,
                lines: Vec::new(),
            };
        }
    };
    let Some(surface_on_face) = face_ref.surface.as_ref() else {
        return ReflectionLineReport {
            face,
            grid_density,
            lines: Vec::new(),
        };
    };
    let (u_range, v_range) = face_param_range(brep, face, &surface_on_face.geometry);
    let lines = compute_reflection_lines_on_surface(
        &surface_on_face.geometry,
        u_range,
        v_range,
        view_dir,
        stripe_dir,
        stripe_spacing,
        grid_density,
    );
    ReflectionLineReport {
        face,
        grid_density,
        lines,
    }
}

/// Compute reflection lines on a standalone [`SurfaceGeometry`] —
/// delegating wrapper for callers that already have a surface +
/// UV range and don't want the BRep-face plumbing.
pub fn compute_reflection_lines_on_surface(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    view_dir: Vec3,
    stripe_dir: Vec3,
    stripe_spacing: f64,
    grid_density: usize,
) -> Vec<ReflectionLine> {
    compute_reflection_lines(
        surface,
        u_range,
        v_range,
        view_dir,
        stripe_dir,
        stripe_spacing,
        grid_density,
    )
}

/// §5.3 item #11 — reflection-line Class-A continuity diagnostic.
///
/// Sample the blend face's outer-loop vertices (boundary sample
/// points) and, at each one, compute the **reflection-line
/// tangent direction** on face_a and face_b via finite
/// differences of the reflected-stripe scalar field. The angle
/// between the two tangents IS the reflection-line kink at that
/// boundary point. Return the maximum angle across all sampled
/// boundary points.
///
/// The reflected-stripe scalar field at `(u, v)` is
/// `f(u, v) = (2·(N·view_dir)·N − view_dir) · stripe_dir`. Its
/// level sets are the reflection lines. The tangent to a level
/// set at `(u, v)` is perpendicular to `∇f` in the surface
/// tangent plane: `T = −(∂f/∂v)·S_u + (∂f/∂u)·S_v` (3D vector).
///
/// **Interpretation.** G2-continuous blend returns `≈ 0`; G1
/// typically O(radians·Δκ·h); G0 up to π/2.
///
/// Returns `None` if no valid sample pairs could be computed.
#[allow(clippy::too_many_arguments)]
pub fn reflection_line_boundary_kink(
    brep: &BRep,
    face_a: FaceHandle,
    face_b: FaceHandle,
    view_dir: Vec3,
    stripe_dir: Vec3,
    _stripe_spacing: f64,
    _grid_density: usize,
    boundary_tol: f64,
) -> Option<f64> {
    use popkern_core::geometry::surface::SurfaceEval;
    let geom_a = brep.faces.get(face_a)?.surface.as_ref()?.geometry.clone();
    let geom_b = brep.faces.get(face_b)?.surface.as_ref()?.geometry.clone();

    // Use the outer-loop 3D polyline AND a dense resampling at
    // the edge midpoints so the boundary is covered even when
    // outer-loop vertices land on different supports.
    let loop_pts = brep.face_outer_loop_positions(face_a);
    if loop_pts.len() < 2 {
        return None;
    }
    let mut boundary_pts: Vec<popkern_core::geometry::Point3> = loop_pts.clone();
    for i in 0..loop_pts.len() {
        let a = loop_pts[i];
        let b = loop_pts[(i + 1) % loop_pts.len()];
        for t in [0.25_f64, 0.5, 0.75] {
            boundary_pts.push(popkern_core::geometry::Point3::new(
                a.x + t * (b.x - a.x),
                a.y + t * (b.y - a.y),
                a.z + t * (b.z - a.z),
            ));
        }
    }

    let tol_sq = boundary_tol * boundary_tol;
    let mut max_angle: Option<f64> = None;
    for pt in &boundary_pts {
        // Locate the point on both surfaces.
        let Some((ua, va)) = geom_a.closest_point(pt) else {
            continue;
        };
        let Some((ub, vb)) = geom_b.closest_point(pt) else {
            continue;
        };
        // Validate the closest-point projections actually land
        // near the boundary vertex (face_b is the support — its
        // closest point to a blend-boundary vertex should be
        // that vertex itself within tolerance).
        let pa = geom_a.point_at(ua, va);
        let pb = geom_b.point_at(ub, vb);
        if (pa.x - pt.x).hypot(pa.y - pt.y).hypot(pa.z - pt.z).powi(2) > tol_sq {
            continue;
        }
        if (pb.x - pt.x).hypot(pb.y - pt.y).hypot(pb.z - pt.z).powi(2) > tol_sq {
            continue;
        }
        let Some(tan_a) = reflection_tangent_at(&geom_a, ua, va, view_dir, stripe_dir) else {
            continue;
        };
        let Some(tan_b) = reflection_tangent_at(&geom_b, ub, vb, view_dir, stripe_dir) else {
            continue;
        };
        let cos = tan_a.dot(&tan_b).clamp(-1.0, 1.0);
        // Tangents may be parallel OR anti-parallel across the
        // boundary (orientation is per-face-arbitrary).
        let angle = cos.abs().acos();
        max_angle = Some(match max_angle {
            Some(m) => m.max(angle),
            None => angle,
        });
    }
    max_angle
}

/// Compute the reflection-line tangent vector on `surface` at
/// parameter `(u, v)` for the given `view_dir` / `stripe_dir`.
/// Returns a unit 3D vector along the level set of the reflected-
/// stripe scalar field through `(u, v)`. Finite-difference
/// gradient on the scalar; tangent = `−f_v·S_u + f_u·S_v` (in
/// 3D, via the surface partials S_u, S_v).
fn reflection_tangent_at(
    surface: &SurfaceGeometry,
    u: f64,
    v: f64,
    view_dir: Vec3,
    stripe_dir: Vec3,
) -> Option<Vec3> {
    use popkern_core::geometry::surface::SurfaceEval;
    let h = 1e-4_f64;
    let f = |u: f64, v: f64| -> f64 {
        let n = surface
            .normal_at(u, v)
            .normalized()
            .unwrap_or(Vec3::new(0.0, 0.0, 1.0));
        let r = n * (2.0 * n.dot(&view_dir)) - view_dir;
        r.dot(&stripe_dir)
    };
    let f_u = (f(u + h, v) - f(u - h, v)) * (0.5 / h);
    let f_v = (f(u, v + h) - f(u, v - h)) * (0.5 / h);
    let s_u = surface.partial_u(u, v);
    let s_v = surface.partial_v(u, v);
    // Tangent = -f_v * S_u + f_u * S_v.
    let t = Vec3::new(
        -f_v * s_u.x + f_u * s_v.x,
        -f_v * s_u.y + f_u * s_v.y,
        -f_v * s_u.z + f_u * s_v.z,
    );
    t.normalized()
}

/// One-call Class-A continuity metric for a blend face: max
/// reflection-line kink at any boundary with any adjacent
/// **non-Gregory** support face.
///
/// Walks the outer loop of `blend_face`, collects each unique
/// adjacent face whose surface is not itself a Gregory sub-patch,
/// and calls [`reflection_line_boundary_kink`] pairwise. Gregory
/// sub-patch → Gregory sub-patch boundaries (shared spokes at a
/// corner) are **excluded** — those are measured under item #3's
/// spoke-G2 follow-up, not this blend-vs-support metric.
/// Returns the max kink across all remaining pairs.
#[allow(clippy::too_many_arguments)]
pub fn blend_reflection_line_kink_max(
    brep: &BRep,
    blend_face: FaceHandle,
    view_dir: Vec3,
    stripe_dir: Vec3,
    stripe_spacing: f64,
    grid_density: usize,
    boundary_tol: f64,
) -> Option<f64> {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let supports = adjacent_support_faces(brep, blend_face)?;
    let mut worst: Option<f64> = None;
    for support in supports {
        // Skip Gregory-vs-Gregory boundaries (spoke continuity is
        // item #3's scope, not the blend-vs-support metric).
        if brep
            .faces
            .get(support)
            .and_then(|f| f.surface.as_ref())
            .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Gregory(_)))
        {
            continue;
        }
        if let Some(angle) = reflection_line_boundary_kink(
            brep,
            blend_face,
            support,
            view_dir,
            stripe_dir,
            stripe_spacing,
            grid_density,
            boundary_tol,
        ) {
            worst = Some(match worst {
                Some(m) => m.max(angle),
                None => angle,
            });
        }
    }
    worst
}

/// Collect unique adjacent faces (with surface geometry) of a
/// blend face by walking the outer loop's half-edge twins.
pub(crate) fn adjacent_support_faces(
    brep: &BRep,
    blend_face: FaceHandle,
) -> Option<Vec<FaceHandle>> {
    let face = brep.faces.get(blend_face)?;
    let outer_loop = brep.loops.get(face.outer_loop)?;
    let start_he = outer_loop.half_edge;
    let mut supports: Vec<FaceHandle> = Vec::new();
    let mut cursor = start_he;
    let bound = brep.half_edges.len().max(1);
    for _ in 0..=bound {
        let he = brep.half_edges.get(cursor)?;
        let twin_h = he.twin;
        if twin_h != popkern_core::topology::HalfEdgeHandle::default()
            && let Some(twin) = brep.half_edges.get(twin_h)
            && let Some(loop_) = brep.loops.get(twin.loop_)
        {
            let f = loop_.face;
            if f != blend_face
                && !supports.contains(&f)
                && brep
                    .faces
                    .get(f)
                    .and_then(|ff| ff.surface.as_ref())
                    .is_some()
            {
                supports.push(f);
            }
        }
        cursor = he.next;
        if cursor == start_he {
            break;
        }
    }
    Some(supports)
}
