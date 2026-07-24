//! Fillet surface construction.
//!
//! Given a center curve and footprint points, builds the rolling-ball fillet
//! surface. For two-plane cases an analytic torus/cylinder section is returned;
//! otherwise a NURBS swept surface is constructed by skinning cross-section arcs.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};

use crate::center_curve::CenterCurveResult;
use crate::fillet_surface_analytic::{
    try_circular_center_torus, try_parallel_cylinder_cylinder, try_plane_plane_analytic,
};
use crate::types::ContinuityLevel;

/// Result of a fillet surface construction.
#[derive(Debug, Clone)]
pub struct FilletSurfaceResult {
    /// The constructed fillet surface.
    pub surface: SurfaceGeometry,
    /// Trim curve on the first surface (3D curve along footprint_1).
    pub trim_curve_1: NurbsCurve,
    /// Trim curve on the second surface (3D curve along footprint_2).
    pub trim_curve_2: NurbsCurve,
    /// §7 continuity downgrade: `Some(reason)` when the caller
    /// requested G2 but the G2 sweep errored and the builder fell
    /// back to the legacy G1 arc sweep. `None` when the delivered
    /// surface matches the requested continuity.
    pub continuity_downgrade: Option<String>,
}

/// Build the fillet surface from center-curve data and the two adjacent surfaces.
///
/// Attempts an analytic construction for the plane-plane case, falling back
/// to a NURBS swept surface for general surface pairs. The `is_constant_radius`
/// flag controls whether the analytic torus path is eligible: a torus has a
/// single minor radius, so it can only represent a constant-radius fillet —
/// passing `false` for a variable-radius spec forces the NURBS sweep path,
/// which honours the per-sample radii baked into `center_result.footprint_*`.
///
/// When `continuity` is `G2` AND `center_result.contact_solutions`
/// is `Some` (i.e., the spine walker produced per-sample contact solutions),
/// the sweep is built via [`crate::g2_cross_section::g2_sweep_from_contact_solutions`]
/// — a degree-5 Bezier cross-section sweep that matches the surface
/// normal curvature at the boundary (G2). Otherwise the legacy
/// degree-1 arc sweep is used. The G2 path needs `orientation_*`
/// to translate the popkern face-normal convention into the
/// `BallSide` convention the contact solver uses.
#[allow(clippy::too_many_arguments)]
pub fn build_fillet_surface(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    center_result: &CenterCurveResult,
    radius: f64,
    arc_samples: usize,
    is_constant_radius: bool,
    continuity: ContinuityLevel,
    _orientation_1: bool,
    _orientation_2: bool,
) -> Result<FilletSurfaceResult, String> {
    // G2 path: route through the degree-5 Bezier sweep when the
    // caller asked for it AND the center-curve producer (the spine
    // walker) gave us per-sample contact solutions. If the G2
    // sweep errors on a specific cross-section (ZeroChord,
    // DegenerateNormal, TangentDegenerate), fall back to the
    // legacy G1 arc sweep and record the reason so the caller
    // emits a ContinuityDowngraded warning.
    let mut g2_fallback_reason: Option<String> = None;
    if continuity == ContinuityLevel::G2
        && let Some(solutions) = center_result.contact_solutions.as_deref()
        && solutions.len() >= 2
    {
        use crate::g2_cross_section::g2_sweep_from_contact_solutions;

        // The cross-section builder derives its contact normals from
        // each converged solution's ball center, so no orientation →
        // BallSide translation is needed here (the previous mapping
        // hardcoded the concave-edge case and mirrored every convex
        // blend — PK-091).
        match g2_sweep_from_contact_solutions(solutions, surface1, surface2) {
            Ok(nurbs) => {
                let trim_1 =
                    crate::center_curve::fit_nurbs_through_points(&center_result.footprint_1)?;
                let trim_2 =
                    crate::center_curve::fit_nurbs_through_points(&center_result.footprint_2)?;
                return Ok(FilletSurfaceResult {
                    surface: SurfaceGeometry::Nurbs(nurbs),
                    trim_curve_1: trim_1,
                    trim_curve_2: trim_2,
                    continuity_downgrade: None,
                });
            }
            Err(e) => {
                g2_fallback_reason = Some(format!("G2 sweep failed: {e}"));
                // Fall through to the legacy arc-sweep path.
            }
        }
    }

    // Try analytic torus fast paths for constant-radius fillets.
    // Variable-radius specs cannot use an analytic torus — its
    // single minor radius would silently average the per-sample
    // radii — so we must build a true NURBS sweep that follows
    // the varying footprints.
    if is_constant_radius {
        // §11.1 plane-plane → torus (dihedral fillet).
        if let (
            SurfaceGeometry::Plane { normal: n1, .. },
            SurfaceGeometry::Plane { normal: n2, .. },
        ) = (surface1, surface2)
            && let Some(result) = try_plane_plane_analytic(n1, n2, center_result, radius)
        {
            return Ok(result);
        }
        // §11.1 analytic torus for fillets with a circular center
        // curve: plane-sphere (always) and plane-cylinder (when
        // the plane is perpendicular to the cylinder axis).
        if let Some(result) = try_circular_center_torus(surface1, surface2, center_result, radius) {
            return Ok(result);
        }
        // §11.1 analytic cylinder for parallel-axes cylinder-cylinder.
        if let Some(result) =
            try_parallel_cylinder_cylinder(surface1, surface2, center_result, radius)
        {
            return Ok(result);
        }
    }

    // General case: NURBS swept surface (G1 arc). If we fell here
    // after a G2 sweep failure, attach the downgrade reason.
    let mut result = build_nurbs_sweep(
        center_result,
        radius,
        arc_samples,
        is_constant_radius.then_some((surface1, surface2)),
    )?;
    if g2_fallback_reason.is_some() {
        result.continuity_downgrade = g2_fallback_reason;
    }
    Ok(result)
}

/// Build a NURBS swept surface by skinning cross-section arcs along the center curve.
///
/// Each cross-section is a polyline arc from footprint_1[i] to footprint_2[i],
/// centered at the corresponding center-curve sample. The arcs are sampled
/// at `arc_samples` points and skinned into a tensor-product surface.
fn build_nurbs_sweep(
    center_result: &CenterCurveResult,
    radius: f64,
    arc_samples: usize,
    supports: Option<(&SurfaceGeometry, &SurfaceGeometry)>,
) -> Result<FilletSurfaceResult, String> {
    // Prefer the edge-span sweep channel when the producer trimmed a
    // grossly overhanging carrier (see CenterCurveResult docs): the
    // legacy over-extended artifacts stay for the face-trim
    // splitter; sections here must stay inside the edge span.
    let (center_curve, footprint_1, footprint_2) = match (
        center_result.sweep_center_curve.as_ref(),
        center_result.sweep_footprint_1.as_ref(),
        center_result.sweep_footprint_2.as_ref(),
    ) {
        (Some(c), Some(f1), Some(f2)) => (c, f1, f2),
        _ => (
            &center_result.center_curve,
            &center_result.footprint_1,
            &center_result.footprint_2,
        ),
    };
    let center_result_curve = center_curve;
    // Footprints can differ by one sample on closed spines (the
    // closure sample lands in one list only); indexing beyond the
    // shorter list panicked. Use the common prefix; a larger
    // mismatch is a walker defect surfaced as an error.
    let n_sections = footprint_1.len().min(footprint_2.len());
    if footprint_1.len().abs_diff(footprint_2.len()) > 1 {
        return Err(format!(
            "footprint length mismatch: {} vs {}",
            footprint_1.len(),
            footprint_2.len()
        ));
    }
    if n_sections < 2 {
        return Err("need at least 2 cross-sections for sweep".into());
    }
    let n_arc = arc_samples.max(3); // at least 3 points per arc

    // Sample center curve at n_sections points.
    let (t_min, t_max) =
        popkern_core::geometry::curve::CurveEval::parameter_range(center_result_curve);

    // Section-exact skinning: each cross-section is a circular arc,
    // exactly representable as ONE rational quadratic row (weight
    // cos(θ/2)) for sweeps up to ~115°. A degree-1 loft of exact
    // rows makes every STATION iso the true circle — the polyline
    // skin's chord sagitta (measured 4.3e-4 at r=0.10 with the
    // default sampling) was the floor for terminal settlement and
    // seam conformity on variable strips. Sections beyond the single
    // -segment range fall back to the legacy sampled skin, recorded
    // as a continuity downgrade.
    let mut rational_rows: Option<Vec<([Point3; 3], f64)>> = Some(Vec::with_capacity(n_sections));
    let mut sampled_points = Vec::with_capacity(n_sections * n_arc);
    for i in 0..n_sections {
        let t = t_min + (t_max - t_min) * (i as f64) / ((n_sections - 1).max(1) as f64);
        let center = popkern_core::geometry::curve::CurveEval::point_at(center_result_curve, t);
        let fp1 = &footprint_1[i];
        let fp2 = &footprint_2[i];
        if let Some(rows) = rational_rows.as_mut() {
            match rational_arc_row(&center, fp1, fp2) {
                Some(row) => rows.push(row),
                None => rational_rows = None,
            }
        }
        let arc_pts = sample_arc(&center, fp1, fp2, radius, n_arc);
        if std::env::var_os("PK_SWEEP_TRACE").is_some() {
            eprintln!(
                "[sweep] sec {i}/{n_sections} t {t:.3} c ({:.3},{:.3},{:.3}) fp1 ({:.3},{:.3},{:.3}) fp2 ({:.3},{:.3},{:.3}) mid ({:.3},{:.3},{:.3})",
                center.x,
                center.y,
                center.z,
                fp1.x,
                fp1.y,
                fp1.z,
                fp2.x,
                fp2.y,
                fp2.z,
                arc_pts[arc_pts.len() / 2].x,
                arc_pts[arc_pts.len() / 2].y,
                arc_pts[arc_pts.len() / 2].z
            );
        }
        sampled_points.extend_from_slice(&arc_pts);
    }

    // U = along edge (n_sections), V = around arc.
    let degree_u = 1.min(n_sections - 1);
    let (nurbs, continuity_downgrade) = if let Some(rows) = rational_rows {
        let mut control_points = Vec::with_capacity(n_sections * 3);
        let mut weights = Vec::with_capacity(n_sections * 3);
        for (cps, w_mid) in &rows {
            control_points.extend_from_slice(cps);
            weights.extend_from_slice(&[1.0, *w_mid, 1.0]);
        }
        let loft = NurbsSurface {
            degree_u,
            degree_v: 2,
            knots_u: build_clamped_uniform_knots(n_sections, degree_u),
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points,
            num_u: n_sections,
            num_v: 3,
            weights,
        };
        (loft, None)
    } else {
        let degree_v = 1.min(n_arc - 1);
        (
            NurbsSurface {
                degree_u,
                degree_v,
                knots_u: build_clamped_uniform_knots(n_sections, degree_u),
                knots_v: build_clamped_uniform_knots(n_arc, degree_v),
                control_points: sampled_points,
                num_u: n_sections,
                num_v: n_arc,
                weights: vec![1.0; n_sections * n_arc],
            },
            Some(
                "section sweep exceeds the single-segment rational range; sampled polyline skin"
                    .to_string(),
            ),
        )
    };

    // Chart overlap (PK-141/PK-177): extend the loft by exact
    // extrapolated station rows at each end on curved-support
    // strips, until the chart reaches past the corner-ball station
    // (clamped charts left the corner chord endpoint OUTSIDE both
    // charts — measured 2.4e-2/3.7e-2 clamp distances at W). Both
    // the rational 3-row loft and the sampled-skin fallback extend;
    // rows are re-projected onto the supports and re-arced, so the
    // extension is exact station geometry; trims and anchors are
    // untouched (walk-based window extension was falsified three
    // ways — see the ledger).
    let (mut nurbs, continuity_downgrade) = (nurbs, continuity_downgrade);
    if let Some((s1, s2)) = supports {
        let curved = !matches!(s1, SurfaceGeometry::Plane { .. })
            || !matches!(s2, SurfaceGeometry::Plane { .. });
        if curved && nurbs.num_u >= 3 {
            crate::loft_chart_extension::extend_loft_rows(&mut nurbs, s1, s2, radius);
        }
    }
    let trim_1 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_1)?;
    let trim_2 = crate::center_curve::fit_nurbs_through_points(&center_result.footprint_2)?;

    Ok(FilletSurfaceResult {
        surface: SurfaceGeometry::Nurbs(nurbs),
        trim_curve_1: trim_1,
        trim_curve_2: trim_2,
        continuity_downgrade,
    })
}

/// One EXACT rational-quadratic arc row from `start` to `end` about
/// `center`: control points [start, tangent-intersection, end] with
/// mid weight cos(θ/2). Returns `None` when the section is outside
/// the well-conditioned single-segment range (θ > ~115°, radii
/// mismatched, or degenerate).
fn rational_arc_row(center: &Point3, start: &Point3, end: &Point3) -> Option<([Point3; 3], f64)> {
    let v1 = *start - *center;
    let v2 = *end - *center;
    let r1 = v1.length();
    let r2 = v2.length();
    if r1 < 1e-12 || r2 < 1e-12 {
        return None;
    }
    if (r1 - r2).abs() > 1e-6 * r1.max(r2) {
        return None;
    }
    let u1 = v1 * (1.0 / r1);
    let u2 = v2 * (1.0 / r2);
    let cos_theta = u1.dot(&u2).clamp(-1.0, 1.0);
    let theta = cos_theta.acos();
    if !(1e-9..=2.0).contains(&theta) {
        return None;
    }
    let half = 0.5 * theta;
    let w_mid = half.cos();
    let mid_dir = (u1 + u2).normalized()?;
    let r = 0.5 * (r1 + r2);
    let mid = *center + mid_dir * (r / w_mid);
    Some(([*start, mid, *end], w_mid))
}

/// Sample a circular arc from `start` to `end` centered at `center`.
///
/// The arc lies in the plane defined by the three points. If `start` and
/// `end` are equidistant from `center` at approximately `radius`, the arc
/// sweeps the shorter angle between them.
pub(crate) fn sample_arc(
    center: &Point3,
    start: &Point3,
    end: &Point3,
    _radius: f64,
    n: usize,
) -> Vec<Point3> {
    let v1 = *start - *center;
    let v2 = *end - *center;
    let r1 = v1.length();
    let r2 = v2.length();

    // Handle degenerate cases.
    if r1 < 1e-15 || r2 < 1e-15 || n < 2 {
        return vec![*start; n];
    }

    // Normalize to unit vectors.
    let u = Vec3::new(v1.x / r1, v1.y / r1, v1.z / r1);
    let w_raw = v2 - u * v2.dot(&u);
    let w_len = w_raw.length();

    if w_len < 1e-15 {
        // start and end are nearly coincident or antipodal — just linearly interpolate.
        let mut pts = Vec::with_capacity(n);
        for i in 0..n {
            let t = i as f64 / (n - 1).max(1) as f64;
            pts.push(Point3::new(
                start.x + (end.x - start.x) * t,
                start.y + (end.y - start.y) * t,
                start.z + (end.z - start.z) * t,
            ));
        }
        return pts;
    }

    let v = Vec3::new(w_raw.x / w_len, w_raw.y / w_len, w_raw.z / w_len);

    // Compute sweep angle.
    let cos_theta = (v1.dot(&v2) / (r1 * r2)).clamp(-1.0, 1.0);
    let theta = cos_theta.acos();

    // Use average radius for the arc.
    let r = (r1 + r2) * 0.5;

    let mut pts = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f64 / (n - 1).max(1) as f64;
        let angle = theta * t;
        let x = r * angle.cos();
        let y = r * angle.sin();
        pts.push(Point3::new(
            center.x + u.x * x + v.x * y,
            center.y + u.y * x + v.y * y,
            center.z + u.z * x + v.z * y,
        ));
    }

    pts
}

/// Build clamped uniform knot vector for `n` control points at given degree.
pub(crate) fn build_clamped_uniform_knots(n: usize, degree: usize) -> Vec<f64> {
    let n_knots = n + degree + 1;
    let n_interior = n.saturating_sub(degree);
    let mut knots = Vec::with_capacity(n_knots);

    knots.extend(std::iter::repeat_n(0.0, degree + 1));
    for i in 1..n_interior {
        knots.push(i as f64 / n_interior as f64);
    }
    knots.extend(std::iter::repeat_n(1.0, degree + 1));

    knots
}

#[cfg(test)]
#[path = "fillet_surface_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "fillet_surface_analytic_tests.rs"]
mod analytic_tests;

#[cfg(test)]
#[path = "fillet_surface_analytic_cylinder_tests.rs"]
mod analytic_cylinder_tests;
