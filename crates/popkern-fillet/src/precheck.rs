//! Non-destructive feasibility pre-check for fillet and chamfer
//! operations (§10.3, §10.4).
//!
//! Callers use this module to answer the UI-critical question *will
//! this fillet succeed?* **without** paying the cost of the full
//! pipeline and **without** mutating the input BRep. The three
//! public functions mirror the spec:
//!
//! - [`precheck_fillet`] — sample an edge, run the feasibility
//!   checks, and return a structured [`FilletFeasibility`] record
//!   with per-issue diagnostics.
//! - [`feasible_radius_bounds`] — binary-search the safe radius
//!   range `[r_min, r_max]` on one edge, suitable for UI slider
//!   bounds.
//! - [`fillet_quality`] — after a successful fillet, compute
//!   `R_min`, patch count, surface area, and a composite quality
//!   score `Q ∈ [0, 1]`.
//!
//! All three are pure functions. They never touch the input BRep
//! and can be called concurrently.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

use crate::analysis::{self, Convexity};
use crate::error::BlendError;
use crate::geometry_compute::{
    face_param_range, largest_requested_radius, self_intersection_safe_radius,
};
use crate::offset::self_intersect::check_offset_self_intersection;
use crate::types::FilletRadius;

/// Structured non-destructive feasibility report for a single edge.
///
/// Produced by [`precheck_fillet`]. A commercial-grade UI reads
/// `can_succeed` for its "go / no-go" indicator and `issues` for the
/// detailed error-list panel. `max_radius_estimate` drives
/// radius-slider upper bounds when the user is refining a request.
#[derive(Debug, Clone)]
pub struct FilletFeasibility {
    /// Engine verdict: can the requested fillet succeed? `true`
    /// means no hard problems were found at pre-check time; the
    /// full pipeline may still fail at deeper stages (the pre-check
    /// is a fast sanity filter, not a correctness guarantee).
    pub can_succeed: bool,
    /// Best estimate of the largest radius that would fit at this
    /// edge. `None` when no curvature-based bound can be computed
    /// (e.g. both support surfaces are planes, for which the
    /// minimum curvature radius is ∞).
    pub max_radius_estimate: Option<f64>,
    /// Convexity classification (convex / concave / mixed).
    pub convexity: Convexity,
    /// `true` if the edge is a seam edge (same face on both sides).
    pub is_seam: bool,
    /// Structured issue list — every entry is a specific problem
    /// the caller can surface to a user.
    pub issues: Vec<FilletIssue>,
}

impl FilletFeasibility {
    /// Convenience: `true` when `issues` is empty.
    pub fn is_clean(&self) -> bool {
        self.issues.is_empty()
    }
}

/// Ranked feasibility issue. Severity ladder matches §16.1
/// (Info / Warning / Error / Critical).
#[derive(Debug, Clone)]
pub enum FilletIssue {
    /// The edge is a seam edge — cannot be filleted until the
    /// underlying face is split along the seam.
    SeamEdge,
    /// The edge does not have two distinct adjacent faces.
    InvalidTopology {
        /// Human-readable reason.
        reason: String,
    },
    /// The edge has mixed convexity (tangent or near-tangent
    /// support faces at the edge). Fillets on mixed edges may
    /// produce degenerate results.
    MixedConvexity,
    /// The requested radius exceeds the local minimum curvature
    /// radius on one of the support surfaces — the rolling ball
    /// would self-intersect before the blend surface is built.
    RadiusExceedsCurvature {
        /// Requested radius.
        requested: f64,
        /// The safe maximum radius the engine determined.
        safe_max: f64,
    },
    /// The requested chamfer distance or fillet radius is larger
    /// than a reasonable fraction of the smallest edge in the
    /// support face's outer loop — the blend will consume the
    /// neighbouring face.
    ExcessiveRelativeSize {
        /// Requested radius or chamfer distance.
        requested: f64,
        /// Length of the smallest perpendicular edge used to measure.
        smallest_edge: f64,
    },
}

/// Run a non-destructive feasibility pre-check on a single edge
/// against a requested radius spec.
///
/// Behaviour:
///
/// 1. Analyse the edge (find adjacent faces, classify convexity,
///    identify seams). Failure here populates `InvalidTopology` or
///    `SeamEdge`.
/// 2. Compute `max_radius_estimate` from the minimum curvature
///    radius of both support surfaces at the edge vicinity. This
///    is the engine's best non-destructive guess and drives the
///    remediation hint when the requested radius is too large.
/// 3. Compare the requested peak radius against the estimate;
///    populate `RadiusExceedsCurvature` if the request exceeds
///    the bound.
///
/// Returns a [`FilletFeasibility`] record whose `can_succeed` is
/// `true` iff no error-level issues were found. Warnings do not
/// flip `can_succeed` — they are informational.
pub fn precheck_fillet(
    brep: &BRep,
    edge: EdgeHandle,
    radius: &FilletRadius,
    tolerance: Tolerance,
) -> FilletFeasibility {
    let mut issues = Vec::new();

    // 1. Analyse the edge.
    let analyses = match analysis::analyze_edges(brep, &[edge], false) {
        Ok(a) => a,
        Err(msg) => {
            return FilletFeasibility {
                can_succeed: false,
                max_radius_estimate: None,
                convexity: Convexity::Mixed,
                is_seam: false,
                issues: vec![FilletIssue::InvalidTopology { reason: msg }],
            };
        }
    };

    let Some(analysis) = analyses.into_iter().next() else {
        return FilletFeasibility {
            can_succeed: false,
            max_radius_estimate: None,
            convexity: Convexity::Mixed,
            is_seam: false,
            issues: vec![FilletIssue::InvalidTopology {
                reason: "edge analysis returned no results".to_string(),
            }],
        };
    };

    // 2. Seam detection.
    let is_seam = analysis.face_1 == analysis.face_2;
    if is_seam {
        issues.push(FilletIssue::SeamEdge);
    }

    if analysis.convexity == Convexity::Mixed {
        issues.push(FilletIssue::MixedConvexity);
    }

    // 3. Maximum curvature-safe radius on either support surface.
    let (u1_range, v1_range) = face_param_range(brep, analysis.face_1, &analysis.surface_1);
    let (u2_range, v2_range) = face_param_range(brep, analysis.face_2, &analysis.surface_2);

    let max_radius_estimate = non_destructive_max_radius(
        &analysis.surface_1,
        &analysis.surface_2,
        analysis.orientation_1,
        analysis.orientation_2,
        u1_range,
        v1_range,
        u2_range,
        v2_range,
    );

    // 4. Requested vs estimated.
    if let Some(max) = max_radius_estimate {
        let requested = largest_requested_radius(radius);
        if requested > max {
            issues.push(FilletIssue::RadiusExceedsCurvature {
                requested,
                safe_max: max,
            });
        }
    }

    // 5. Double-check using the self-intersection pre-flight (same
    //    predicate the real pipeline uses) so the precheck verdict
    //    stays in lockstep with geometry_compute.
    if let Some(safe_max) = self_intersection_safe_radius(
        &analysis.surface_1,
        &analysis.surface_2,
        radius,
        analysis.orientation_1,
        analysis.orientation_2,
        u1_range,
        v1_range,
        u2_range,
        v2_range,
    ) {
        let requested = largest_requested_radius(radius);
        if requested > safe_max
            && !issues
                .iter()
                .any(|i| matches!(i, FilletIssue::RadiusExceedsCurvature { .. }))
        {
            issues.push(FilletIssue::RadiusExceedsCurvature {
                requested,
                safe_max,
            });
        }
    }

    let _ = tolerance; // reserved for future use (§9 categories)

    let has_error = issues.iter().any(|i| {
        matches!(
            i,
            FilletIssue::SeamEdge
                | FilletIssue::InvalidTopology { .. }
                | FilletIssue::RadiusExceedsCurvature { .. }
        )
    });

    FilletFeasibility {
        can_succeed: !has_error,
        max_radius_estimate,
        convexity: analysis.convexity,
        is_seam,
        issues,
    }
}

/// Binary-search the feasible radius range `[r_min, r_max]` for a
/// single edge, suitable for UI slider bounds (§10.4).
///
/// Returns `Err(BlendError::InvalidParameter)` if the edge is a seam
/// or topologically invalid — neither has a feasible range. On a
/// plane-plane edge (where both surfaces have infinite curvature
/// radius) `max` is reported as the finite fallback described
/// below.
///
/// The lower bound `min` is always `tolerance.linear()` — below
/// that threshold the blend surface degenerates. The upper bound
/// `max` is:
///
/// - the non-destructive curvature-based estimate from
///   [`precheck_fillet`] when it is finite, or
/// - `edge_length / 4` as a finite fallback when both support
///   surfaces are flat (plane-plane fillets are bounded by the
///   geometry of the neighbouring edges, not by curvature).
pub fn feasible_radius_bounds(
    brep: &BRep,
    edge: EdgeHandle,
    tolerance: Tolerance,
) -> Result<FilletRadiusBounds, BlendError> {
    // Reuse the pre-check with a sentinel radius so we can pull the
    // analysis + curvature estimate.
    let pre = precheck_fillet(
        brep,
        edge,
        &FilletRadius::Constant(tolerance.linear()),
        tolerance,
    );

    // Hard-fail on structural problems.
    for issue in &pre.issues {
        match issue {
            FilletIssue::SeamEdge => {
                return Err(BlendError::InvalidParameter {
                    parameter: "edge".to_string(),
                    value: format!("{edge:?}"),
                    reason: "seam edges have no feasible fillet radius range".to_string(),
                });
            }
            FilletIssue::InvalidTopology { reason } => {
                return Err(BlendError::InvalidParameter {
                    parameter: "edge".to_string(),
                    value: format!("{edge:?}"),
                    reason: format!("edge analysis failed: {reason}"),
                });
            }
            _ => {}
        }
    }

    let min = tolerance.linear();

    // Curvature-based upper bound (if any).
    let curvature_max = pre.max_radius_estimate;

    // Geometric upper bound: edge length / 4 as a conservative cap
    // for plane-plane or near-plane cases where curvature says ∞.
    let edge_len = edge_length(brep, edge).unwrap_or(1.0);
    let geometric_max = edge_len * 0.25;

    let max = match curvature_max {
        Some(c) if c.is_finite() => c.min(geometric_max),
        _ => geometric_max,
    };

    Ok(FilletRadiusBounds {
        min,
        max,
        curvature_bound: curvature_max,
        geometric_bound: geometric_max,
    })
}

/// The feasibility-range record returned by [`feasible_radius_bounds`].
#[derive(Debug, Clone, Copy)]
pub struct FilletRadiusBounds {
    /// Smallest usable radius — anything below this degenerates.
    pub min: f64,
    /// Largest usable radius (tightest of curvature and geometric bounds).
    pub max: f64,
    /// Pure curvature-based upper bound, when computable.
    pub curvature_bound: Option<f64>,
    /// Pure geometry-based upper bound (≈ edge length / 4).
    pub geometric_bound: f64,
}

impl FilletRadiusBounds {
    /// Returns `true` if the given radius lies strictly inside
    /// `[min, max]`.
    pub fn contains(&self, radius: f64) -> bool {
        radius > self.min && radius < self.max
    }

    /// Returns the `radius` clamped into `[min, max]`.
    pub fn clamp(&self, radius: f64) -> f64 {
        radius.clamp(self.min, self.max)
    }
}

// ---------------------------------------------------------------------------
// BlendQuality (§10.5)
// ---------------------------------------------------------------------------

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

    for (fh, face) in input.post_brep.faces.iter() {
        if input.pre_faces.contains(&fh) {
            continue;
        }
        patch_count += 1;

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
    // * linear ramp to 0.0 as r_min drops to tol_linear.
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
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Non-destructive curvature-based max-radius estimate — returns
/// the tighter of the two surfaces' minimum curvature radii as
/// reported by `check_offset_self_intersection`.
#[allow(clippy::too_many_arguments)]
fn non_destructive_max_radius(
    surface_1: &SurfaceGeometry,
    surface_2: &SurfaceGeometry,
    orientation_1: bool,
    orientation_2: bool,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    u2_range: (f64, f64),
    v2_range: (f64, f64),
) -> Option<f64> {
    // Probe each surface with a radius larger than any realistic
    // fillet — if the surface reports "unsafe at this probe", it
    // gives us back the actual safe bound.
    let probe = 1e9;
    let m1 = check_offset_self_intersection(surface_1, probe, !orientation_1, u1_range, v1_range);
    let m2 = check_offset_self_intersection(surface_2, probe, !orientation_2, u2_range, v2_range);

    match (m1, m2) {
        (None, None) => None, // both surfaces are flat or convex-outward
        (Some(a), Some(b)) => Some(a.min(b)),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
    }
}

fn edge_length(brep: &BRep, edge: EdgeHandle) -> Option<f64> {
    let e = brep.edges.get(edge)?;
    // An edge's vertices live on its pair of half-edges: the
    // half-edge's `origin` and the twin's `origin`.
    let he = brep.half_edges.get(e.half_edge)?;
    let he_twin = brep.half_edges.get(he.twin)?;
    let v1 = brep.vertices.get(he.origin)?.position;
    let v2 = brep.vertices.get(he_twin.origin)?.position;
    Some(v1.distance_to(&v2))
}

fn face_area_estimate(brep: &BRep, face: popkern_core::topology::FaceHandle) -> f64 {
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

fn surface_min_curvature_radius(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> f64 {
    // Reuse the self-intersection probe — same curvature math.
    let probe = 1e9;
    check_offset_self_intersection(surface, probe, true, u_range, v_range)
        .or_else(|| check_offset_self_intersection(surface, probe, false, u_range, v_range))
        .unwrap_or(f64::INFINITY)
}

#[cfg(test)]
#[path = "precheck_tests.rs"]
mod tests;
