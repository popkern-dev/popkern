//! Post-operation blend quality metrics (§10.5).
//!
//! [`fillet_quality`] walks the new-face delta produced by a
//! fillet or chamfer operation and returns a [`BlendQuality`]
//! record with `R_min`, patch count, surface area, and a composite
//! quality score in `[0, 1]`.
//!
//! Broken out from `precheck` so each module stays under the
//! 500-line structural limit.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::BRep;

use crate::offset::self_intersect::check_offset_self_intersection;
use crate::types::ContinuityLevel;

/// Post-operation quality metrics for a successful blend (§10.5).
///
/// Produced by [`fillet_quality`] after a successful fillet. A UI
/// can display `composite_score` as a per-blend "health" indicator
/// and drill into the raw metrics for a detailed analysis panel.
///
/// All distance values are in model units; `surface_area` is the
/// sum of area of the new blend faces; `composite_score` is in
/// `[0, 1]` with 1.0 meaning "analytically exact on simple geometry".
#[derive(Debug, Clone, Copy)]
pub struct BlendQuality {
    /// Minimum radius of curvature anywhere on the blend surface.
    /// Manufacturing-critical: the CNC tool radius must be smaller
    /// than this value to machine the blend without gouging.
    pub r_min: f64,
    /// Number of new faces created by the operation.
    pub patch_count: usize,
    /// Total area of the new blend faces.
    pub surface_area: f64,
    /// Composite quality score in `[0, 1]`. 1.0 for analytical
    /// fast-path output on simple geometry; decreases monotonically
    /// with geometric difficulty. See [`ComputeQualityInput`] for
    /// the scoring formula.
    pub composite_score: f64,
    /// The highest [`ContinuityLevel`] delivered by any new blend
    /// face, inferred from the output surface geometry.
    ///
    /// A degree-5 v-direction NURBS sweep is the fingerprint of
    /// the G2 path (`g2_sweep_from_contact_solutions` produces
    /// exactly that shape); any other surface — analytic Torus /
    /// Cylinder / Sphere or the degree-1 arc sweep — is
    /// classified as G1 (the rolling ball gives G1 by construction).
    ///
    /// **Why max, not min.** A fillet on a cube edge routed through
    /// G2 produces a degree-5 NURBS strip (G2) flanked by small
    /// corner cleanup planes (G1) from the multi-chord pipeline.
    /// The user asked "did G2 work?" — the answer is yes, on the
    /// primary surface. The small corner planes are expected and
    /// unrelated to the spine's continuity. Max across faces
    /// captures that intent; inspect individual `face.surface`
    /// entries if a stricter per-face report is needed.
    ///
    /// UI surfaces should display this alongside a "requested"
    /// badge from [`crate::types::FilletOptions::continuity`] so
    /// callers can see at a glance when the engine fell back. The
    /// fallback also emits a
    /// [`crate::error::BlendWarning::ContinuityDowngraded`]
    /// warning on the [`crate::BlendResult::edge_warnings`] list.
    pub achieved_continuity: ContinuityLevel,
}

/// Inputs to [`fillet_quality`]. The function walks the BRep delta
/// (new faces vs old face set) so callers can either pass the raw
/// pre/post face handle sets (the typical case) or precomputed
/// entity lists.
#[derive(Debug, Clone)]
pub struct ComputeQualityInput<'a> {
    /// The BRep after the fillet operation completed.
    pub post_brep: &'a BRep,
    /// The set of face handles that existed *before* the operation,
    /// used to distinguish new blend faces from surviving original
    /// faces. A face handle that is in `post_brep.faces` but not in
    /// this set is counted as a new blend face.
    pub pre_faces: &'a std::collections::HashSet<popkern_core::topology::FaceHandle>,
    /// The radius that was requested (for relative quality scoring).
    pub requested_radius: f64,
}

/// Compute blend quality metrics for a completed fillet operation.
pub fn fillet_quality(input: &ComputeQualityInput) -> BlendQuality {
    let mut patch_count = 0usize;
    let mut surface_area = 0.0f64;
    let mut r_min = f64::INFINITY;
    // Achieved continuity is the MAX across all new faces — see the
    // BlendQuality::achieved_continuity docs for the rationale.
    let mut achieved_continuity = ContinuityLevel::G1;

    for (fh, face) in input.post_brep.faces.iter() {
        if input.pre_faces.contains(&fh) {
            continue;
        }
        patch_count += 1;
        let face_continuity = face
            .surface
            .as_ref()
            .map(|s| classify_face_continuity(&s.geometry))
            .unwrap_or(ContinuityLevel::G1);
        if face_continuity > achieved_continuity {
            achieved_continuity = face_continuity;
        }

        let area = face_area_estimate(input.post_brep, fh);
        surface_area += area;

        // Estimate the per-face minimum curvature radius from the
        // surface geometry. Analytic surfaces return exact values.
        // `face.surface` is `Option<SurfaceOnFace>` — skip faces
        // with no attached surface geometry.
        let Some(surface_on_face) = face.surface.as_ref() else {
            continue;
        };
        let (u_range, v_range) = crate::geometry_compute::face_param_range(
            input.post_brep,
            fh,
            &surface_on_face.geometry,
        );
        let face_r_min = surface_min_curvature_radius(&surface_on_face.geometry, u_range, v_range);
        if face_r_min < r_min {
            r_min = face_r_min;
        }
    }

    if !r_min.is_finite() {
        // No new faces or all flat — fall back to the requested radius.
        r_min = input.requested_radius;
    }

    // Composite scoring heuristic:
    //
    // * 1.0 if r_min ≥ requested_radius (analytical — the blend
    //   surface matches or exceeds the curvature we asked for),
    // * linear ramp to 0.0 as r_min drops toward tol_linear.
    let requested = input.requested_radius.max(f64::EPSILON);
    let composite_score = if r_min >= requested {
        1.0
    } else {
        (r_min / requested).clamp(0.0, 1.0)
    };

    BlendQuality {
        r_min,
        patch_count,
        surface_area,
        composite_score,
        achieved_continuity,
    }
}

/// Classify a blend face's surface geometry into the continuity
/// level it delivers at the blend boundaries.
///
/// - Degree-5 v-direction NURBS (the fingerprint of
///   [`crate::g2_cross_section::g2_sweep_from_contact_solutions`]) → G2
/// - Any other surface (Torus / Cylinder / Sphere / Plane / lower-
///   degree NURBS sweep) → G1 (rolling ball gives G1 by
///   construction; analytic fast paths represent rolling-ball
///   geometry exactly).
fn classify_face_continuity(surface: &SurfaceGeometry) -> ContinuityLevel {
    if let SurfaceGeometry::Nurbs(ns) = surface
        && ns.degree_v == 5
    {
        return ContinuityLevel::G2;
    }
    ContinuityLevel::G1
}

pub(crate) fn face_area_estimate(brep: &BRep, face: popkern_core::topology::FaceHandle) -> f64 {
    // Approximate area from the outer loop polygon (Newell's method
    // for 3D polygons).
    let positions = brep.face_outer_loop_positions(face);
    if positions.len() < 3 {
        return 0.0;
    }

    let mut nx = 0.0_f64;
    let mut ny = 0.0_f64;
    let mut nz = 0.0_f64;
    let n_pts = positions.len();
    for i in 0..n_pts {
        let a = positions[i];
        let b = positions[(i + 1) % n_pts];
        nx += (a.y - b.y) * (a.z + b.z);
        ny += (a.z - b.z) * (a.x + b.x);
        nz += (a.x - b.x) * (a.y + b.y);
    }
    0.5 * (nx * nx + ny * ny + nz * nz).sqrt()
}

/// Minimum radius of curvature for a surface over the given UV range.
///
/// Analytic surfaces return exact values:
/// - **Plane** — infinite (zero curvature everywhere).
/// - **Cylinder** — `radius` (circumferential curvature `1/r`, axial
///   curvature `0`). This covers straight-edge plane-plane fillets,
///   which now produce Cylinder surfaces instead of Torus.
/// - **Sphere** — `radius` (uniform curvature).
/// - **Torus** — `minor_radius` (the cross-section curvature `1/r`
///   of the rolling ball). This covers curved-edge fillets.
///
/// Cone and NURBS surfaces fall back to the sampling-based estimate
/// via the offset self-intersection probe.
fn surface_min_curvature_radius(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> f64 {
    match surface {
        SurfaceGeometry::Plane { .. } => f64::INFINITY,
        SurfaceGeometry::Cylinder { radius, .. } => *radius,
        SurfaceGeometry::Sphere { radius, .. } => *radius,
        SurfaceGeometry::Torus { minor_radius, .. } => *minor_radius,
        // Cone and NURBS: fall back to the sampling-based estimate.
        _ => {
            let probe = 1e9;
            check_offset_self_intersection(surface, probe, true, u_range, v_range)
                .or_else(|| check_offset_self_intersection(surface, probe, false, u_range, v_range))
                .unwrap_or(f64::INFINITY)
        }
    }
}
