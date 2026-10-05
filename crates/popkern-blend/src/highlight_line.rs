//! §14.8 highlight-line preview.
//!
//! Thin wrapper over `popkern_core::healing::compute_highlight_lines`
//! that maps it onto a single BRep face — resolving the surface
//! geometry and finite parameter range automatically from the face
//! handle. This turns the core primitive into a one-call fillet
//! preview: "give me all the highlight lines on this blend face for
//! the given light and view."
//!
//! For a blend-vs-support boundary continuity diagnostic (§15 Class-A
//! "highlight lines cross boundaries without kink"), see
//! [`highlight_boundary_kink`], which reports the max tangent-angle
//! discontinuity between isophotes on two adjacent faces.
//!
//! # Example
//!
//! ```no_run
//! use popkern_core::geometry::Vec3;
//! use popkern_blend::highlight_line::{highlight_lines_on_face, HighlightLineReport};
//! # use popkern_core::topology::{BRep, FaceHandle};
//! # let brep: BRep = unimplemented!();
//! # let fillet_face: FaceHandle = unimplemented!();
//! let report = highlight_lines_on_face(
//!     &brep,
//!     fillet_face,
//!     Vec3::new(0.0, 0.0, 1.0),
//!     Vec3::new(0.0, 0.0, 1.0),
//!     50,
//! );
//! // `report.lines` is ready to render as a set of polylines.
//! ```

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::healing::{Isophote, compute_highlight_lines};
use popkern_core::topology::{BRep, FaceHandle};
use popkern_math::PortableMath;

use crate::geometry_compute::face_param_range;

/// Result of running highlight-line analysis on a single face.
#[derive(Debug, Clone)]
pub struct HighlightLineReport {
    /// The face that was analysed.
    pub face: FaceHandle,
    /// Grid density used for the sampling (passed through for
    /// logging or re-rendering at the same resolution).
    pub grid_density: usize,
    /// Raw isophote polylines from `compute_highlight_lines`.
    /// Each [`Isophote`] is a single contour at a specific intensity
    /// level; rendering a highlight-line preview means drawing all
    /// of them as thin white lines.
    pub lines: Vec<Isophote>,
}

impl HighlightLineReport {
    /// Total number of highlight line contours across all intensity
    /// levels.
    pub fn line_count(&self) -> usize {
        self.lines.len()
    }

    /// Total number of sample points across every highlight line.
    /// A good proxy for visualisation density.
    pub fn total_samples(&self) -> usize {
        self.lines.iter().map(|l| l.samples.len()).sum()
    }
}

/// Compute highlight lines on a single face of the BRep.
///
/// Resolves the face's [`SurfaceGeometry`] and a bounded `(u_range,
/// v_range)` parameter box from the face's outer-loop vertices, then
/// delegates to `popkern_core::healing::compute_highlight_lines`.
///
/// Returns a report with an empty `lines` vector when the face has no
/// attached surface — not an error (chamfer faces, freshly-split
/// faces without surface re-assignment). Callers that want to
/// guarantee a valid surface should `.unwrap_or_default()` the result
/// or check `face.surface.is_some()` up front.
///
/// Arguments:
/// - `brep` — the post-fillet BRep
/// - `face` — the face to analyse (typically a blend face from
///   `identify_blend` or a freshly-created one)
/// - `light_dir` — directional light vector (need not be unit)
/// - `view_dir` — view/camera direction (need not be unit). For an
///   orthographic preview pass `Vec3::new(0, 0, 1)`.
/// - `grid_density` — number of sample points in each parameter
///   direction. `50` is a good default; below `20` risks missing
///   contours.
pub fn highlight_lines_on_face(
    brep: &BRep,
    face: FaceHandle,
    light_dir: Vec3,
    view_dir: Vec3,
    grid_density: usize,
) -> HighlightLineReport {
    let face_ref = match brep.faces.get(face) {
        Some(f) => f,
        None => {
            return HighlightLineReport {
                face,
                grid_density,
                lines: Vec::new(),
            };
        }
    };
    let Some(surface_on_face) = face_ref.surface.as_ref() else {
        return HighlightLineReport {
            face,
            grid_density,
            lines: Vec::new(),
        };
    };
    let (u_range, v_range) = face_param_range(brep, face, &surface_on_face.geometry);
    let lines = compute_highlight_lines_on_surface(
        &surface_on_face.geometry,
        u_range,
        v_range,
        light_dir,
        view_dir,
        grid_density,
    );
    HighlightLineReport {
        face,
        grid_density,
        lines,
    }
}

/// Compute highlight lines on a standalone [`SurfaceGeometry`] — the
/// function the face-aware wrapper delegates to. Exposed so callers
/// that already have a `SurfaceGeometry` (e.g. from `preview_fillet`
/// before BRep application) can run the analysis without plumbing a
/// full BRep face handle through.
pub fn compute_highlight_lines_on_surface(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    light_dir: Vec3,
    view_dir: Vec3,
    grid_density: usize,
) -> Vec<Isophote> {
    compute_highlight_lines(surface, u_range, v_range, light_dir, view_dir, grid_density)
}

/// §15 Class-A continuity diagnostic: measure the max tangent-angle
/// discontinuity between isophotes on two adjacent faces.
///
/// Walks every pair of (sample on face_a, sample on face_b) whose 3D
/// points are within `boundary_tol` of each other — i.e. both land
/// near the shared edge — and computes the angle between the two
/// isophotes' local tangent directions. Returns the maximum angle
/// found (in radians), or `None` if no pairs were identified.
///
/// **Interpretation.** A G2-continuous blend produces `≈ 0` (the
/// tangent direction is continuous across the boundary). A G1 blend
/// produces values up to a few degrees. A G0 (positional-only) blend
/// produces large jumps, potentially π/2 or more.
///
/// **Cost.** O(density² × n_samples_per_isophote² × n_levels²) —
/// brute-force pair matching. Call with small grid densities
/// (`20–30`) for interactive feedback and larger values for offline
/// reports.
pub fn highlight_boundary_kink(
    brep: &BRep,
    face_a: FaceHandle,
    face_b: FaceHandle,
    light_dir: Vec3,
    view_dir: Vec3,
    grid_density: usize,
    boundary_tol: f64,
) -> Option<f64> {
    let report_a = highlight_lines_on_face(brep, face_a, light_dir, view_dir, grid_density);
    let report_b = highlight_lines_on_face(brep, face_b, light_dir, view_dir, grid_density);

    let mut max_angle: Option<f64> = None;
    let tol_sq = boundary_tol * boundary_tol;

    for line_a in &report_a.lines {
        let tangents_a = isophote_tangents(line_a);
        for (i_a, (p_a, tan_a)) in line_a.samples.iter().zip(tangents_a.iter()).enumerate() {
            // Skip endpoint tangents — they're one-sided.
            if i_a == 0 || i_a == line_a.samples.len() - 1 {
                continue;
            }
            let (_, _, pt_a) = p_a;
            for line_b in &report_b.lines {
                let tangents_b = isophote_tangents(line_b);
                for (i_b, (p_b, tan_b)) in line_b.samples.iter().zip(tangents_b.iter()).enumerate()
                {
                    if i_b == 0 || i_b == line_b.samples.len() - 1 {
                        continue;
                    }
                    let (_, _, pt_b) = p_b;
                    let dx = pt_a.x - pt_b.x;
                    let dy = pt_a.y - pt_b.y;
                    let dz = pt_a.z - pt_b.z;
                    let d_sq = dx * dx + dy * dy + dz * dz;
                    if d_sq > tol_sq {
                        continue;
                    }
                    let cos = tan_a.dot(tan_b).clamp(-1.0, 1.0);
                    // Take the absolute value — isophote tangents
                    // may be pointed in opposite directions on each
                    // face because the sampling orientation is per-
                    // face-arbitrary. A G2 cross should either
                    // align the tangents or anti-align them; we
                    // treat both as "zero angle."
                    let angle = (1.0 - cos.abs()).max(0.0).sqrt(); // chord approximation for small θ
                    let angle_rad = 2.0 * angle.pasin(); // exact for |cos| ≤ 1
                    let _ = angle;
                    let current = match max_angle {
                        Some(m) => m.max(angle_rad),
                        None => angle_rad,
                    };
                    max_angle = Some(current);
                }
            }
        }
    }

    max_angle
}

/// Compute per-sample tangent unit vectors for an [`Isophote`].
///
/// Uses central differences over the 3D points of the polyline. The
/// first and last samples use forward/backward differences.
///
/// Returns a vector the same length as `isophote.samples`; each entry
/// is a unit tangent at that sample. Degenerate (zero-length) tangent
/// vectors fall back to `(1, 0, 0)`.
fn isophote_tangents(isophote: &Isophote) -> Vec<Vec3> {
    let n = isophote.samples.len();
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![Vec3::new(1.0, 0.0, 0.0)];
    }
    let mut tangents = Vec::with_capacity(n);
    for i in 0..n {
        let (lo, hi) = if i == 0 {
            (0, 1)
        } else if i == n - 1 {
            (n - 2, n - 1)
        } else {
            (i - 1, i + 1)
        };
        let (_, _, p_lo) = isophote.samples[lo];
        let (_, _, p_hi) = isophote.samples[hi];
        let v = Vec3::new(p_hi.x - p_lo.x, p_hi.y - p_lo.y, p_hi.z - p_lo.z);
        tangents.push(v.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0)));
    }
    tangents
}

#[cfg(test)]
#[path = "highlight_line_tests.rs"]
mod tests;
