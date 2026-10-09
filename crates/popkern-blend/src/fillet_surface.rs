//! Fillet surface construction.
//!
//! Given a center curve and footprint points, builds the rolling-ball fillet
//! surface. For two-plane cases an analytic torus/cylinder section is returned;
//! otherwise a NURBS swept surface is constructed by skinning cross-section arcs.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3};

use crate::center_curve::CenterCurveResult;
use crate::fillet_surface_analytic::{
    try_circular_center_torus, try_cylindrical_analytic, try_plane_plane_analytic,
};
use crate::types::ContinuityLevel;

/// Result of a fillet surface construction.
#[derive(Debug, Clone, PartialEq)]
pub struct FilletSurfaceResult {
    /// The constructed fillet surface.
    pub surface: SurfaceGeometry,
    /// Trim curve on the first surface (3D curve along footprint_1).
    pub trim_curve_1: NurbsCurve,
    /// Trim curve on the second surface (3D curve along footprint_2).
    pub trim_curve_2: NurbsCurve,
    /// §7 continuity downgrade: `Some(reason)` when the caller
    /// requested G2 but the G2 sweep errored and the builder fell
    /// back to the G1 arc sweep. `None` when the delivered
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
/// normal curvature at the boundary (G2). Otherwise the
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
    profile: crate::types::BlendProfile,
) -> Result<FilletSurfaceResult, String> {
    // PK-189 conic-by-rho profile: same contacts and end tangents as
    // the rolling-ball arc (the tangent triangle is profile-
    // independent), but the section interior follows the requested
    // shoulder factor. The analytic torus/cylinder fast paths and
    // the G2 sweep represent CIRCULAR/G2 sections only, so a conic
    // request routes straight to the exact rational section sweep
    // with the conic middle weight; a station outside the single-
    // segment rational range is a structured error, never a silent
    // circular downgrade (rho at the circular value must reproduce
    // the arc to roundoff, so the conic path subsumes the arc).
    let conic_rho = match profile {
        crate::types::BlendProfile::Conic { rho } => Some(rho),
        _ => None,
    };
    if let Some(rho) = conic_rho {
        return build_nurbs_sweep(
            center_result,
            radius,
            arc_samples,
            is_constant_radius.then_some((surface1, surface2)),
            Some(rho),
        );
    }
    // G2 path: route through the degree-5 Bezier sweep when the
    // caller asked for it AND the center-curve producer (the spine
    // walker) gave us per-sample contact solutions. If the G2
    // sweep errors on a specific cross-section (ZeroChord,
    // DegenerateNormal, TangentDegenerate), fall back to the
    // G1 arc sweep and record the reason so the caller
    // emits a ContinuityDowngraded warning.
    let mut g2_fallback_reason: Option<String> = None;
    if continuity == ContinuityLevel::G2
        && let Some(solutions) = center_result.contact_solutions.as_deref()
        && solutions.len() >= 2
    {
        use crate::g2_cross_section::g2_sweep_from_contact_solutions;

        // Circular spine (coaxial plane × cylinder/cone/sphere): the
        // exact surface of revolution of one G2 section — its rows
        // are the rails, where the cubic loft below left them by
        // 2.6e-4 (PK-389). Certified per station; `None` keeps the
        // loft.
        // The strip's boundary rows are the contact solutions
        // themselves (each G2 section starts and ends at them), so the
        // trims are fitted through the same points. The clipped
        // centre-curve footprints are a different construction and
        // left the base rail of a G2 pentagonal-pyramid strip 1.25e-2
        // off the strip (3.9e-2 check 9, PK-531).
        let contacts_0: Vec<_> = solutions.iter().map(|s| s.contact0).collect();
        let contacts_1: Vec<_> = solutions.iter().map(|s| s.contact1).collect();
        if let Some(nurbs) = crate::g2_revolution::g2_revolution_from_contact_solutions(
            solutions, surface1, surface2,
        ) {
            let trim_1 = crate::center_curve::fit_nurbs_through_points(&contacts_0)?;
            let trim_2 = crate::center_curve::fit_nurbs_through_points(&contacts_1)?;
            return Ok(FilletSurfaceResult {
                surface: SurfaceGeometry::Nurbs(nurbs),
                trim_curve_1: trim_1,
                trim_curve_2: trim_2,
                continuity_downgrade: None,
            });
        }

        // The cross-section builder derives its contact normals from
        // each converged solution's ball center, so no orientation →
        // BallSide translation is needed here (the previous mapping
        // hardcoded the concave-edge case and mirrored every convex
        // blend — PK-091).
        match g2_sweep_from_contact_solutions(solutions, surface1, surface2) {
            Ok(nurbs) => {
                let trim_1 = crate::center_curve::fit_nurbs_through_points(&contacts_0)?;
                let trim_2 = crate::center_curve::fit_nurbs_through_points(&contacts_1)?;
                return Ok(FilletSurfaceResult {
                    surface: SurfaceGeometry::Nurbs(nurbs),
                    trim_curve_1: trim_1,
                    trim_curve_2: trim_2,
                    continuity_downgrade: None,
                });
            }
            Err(e) => {
                g2_fallback_reason = Some(format!("G2 sweep failed: {e}"));
                // Fall through to the arc-sweep path.
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
        // §11.1 analytic cylinder: parallel-axes cylinder pairs, and a
        // plane parallel to a cylinder axis with the edge along it.
        if let Some(result) = try_cylindrical_analytic(surface1, surface2, center_result, radius) {
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
        None,
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
    conic_rho: Option<f64>,
) -> Result<FilletSurfaceResult, String> {
    // Prefer the edge-span sweep channel when the producer trimmed a
    // grossly overhanging carrier (see CenterCurveResult docs): the
    // historical over-extended artifacts stay for the face-trim
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

    // Station refinement (PK-177): re-solve stations against the
    // supports (~1e-10; walker ~1e-5); declines keep walker values.
    // Refined stations keep their solved centres (exact rows, PK-734).
    let mut station_centres: Vec<Option<Point3>> = vec![None; n_sections];
    let refined: Option<Vec<(Point3, Point3)>> = supports.map(|(s1, s2)| {
        (0..n_sections)
            .map(|i| {
                match crate::station_geometry::refine_station(
                    s1,
                    s2,
                    &footprint_1[i],
                    &footprint_2[i],
                    radius,
                ) {
                    Some((a, b, c)) => {
                        station_centres[i] = Some(c);
                        (a, b)
                    }
                    None => (footprint_1[i], footprint_2[i]),
                }
            })
            .collect()
    });
    let footprint_1: Vec<Point3> = match &refined {
        Some(r) => r.iter().map(|(a, _)| *a).collect(),
        None => footprint_1[..n_sections].to_vec(),
    };
    let footprint_2: Vec<Point3> = match &refined {
        Some(r) => r.iter().map(|(_, b)| *b).collect(),
        None => footprint_2[..n_sections].to_vec(),
    };
    let footprint_1 = &footprint_1;
    let footprint_2 = &footprint_2;

    // Sample center curve at n_sections points.
    let (t_min, t_max) =
        popkern_core::geometry::curve::CurveEval::parameter_range(center_result_curve);

    // Section-exact skinning: one rational quadratic row per
    // station (<=~115°); two-segment bi-rational beyond; sampled
    // skin as the last resort (continuity downgrade recorded).
    let mut rational_rows: Option<Vec<([Point3; 3], f64)>> = Some(Vec::with_capacity(n_sections));
    let mut bi_rows: Option<Vec<([Point3; 5], f64)>> = Some(Vec::with_capacity(n_sections));
    let mut sampled_points = Vec::with_capacity(n_sections * n_arc);
    for i in 0..n_sections {
        let t = t_min + (t_max - t_min) * (i as f64) / ((n_sections - 1).max(1) as f64);
        let center = station_centres[i].unwrap_or_else(|| {
            popkern_core::geometry::curve::CurveEval::point_at(center_result_curve, t)
        });
        let fp1 = &footprint_1[i];
        let fp2 = &footprint_2[i];
        // A collapsed section is the strip's vanishing end (PK-764).
        if fp1.distance_to(fp2) <= radius * 1e-9 {
            push_collapsed_row(
                fp1,
                &mut rational_rows,
                &mut bi_rows,
                &mut sampled_points,
                n_arc,
            );
            continue;
        }
        if let Some(rows) = rational_rows.as_mut() {
            match rational_arc_row(&center, fp1, fp2) {
                Some((cps, w_circular)) => {
                    // Conic profile: identical tangent-triangle
                    // control points, shoulder weight rho/(1-rho)
                    // instead of cos(theta/2).
                    let w = match conic_rho {
                        Some(rho) => crate::types::BlendProfile::conic_weight(rho),
                        None => w_circular,
                    };
                    rows.push((cps, w));
                }
                None if conic_rho.is_some() => {
                    return Err(format!(
                        "conic profile station {i}/{n_sections} is outside the                          single-segment rational section range (opening angle too                          wide or degenerate contacts); conic sections do not                          downgrade silently"
                    ));
                }
                None => rational_rows = None,
            }
        }
        if let Some(rows) = bi_rows.as_mut() {
            match crate::station_geometry::bi_rational_arc_row(&center, fp1, fp2) {
                Some(row) => rows.push(row),
                None => bi_rows = None,
            }
        }
        let arc_pts = sample_arc(&center, fp1, fp2, radius, n_arc);
        if popkern_core::trace::trace_enabled("sweep") {
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
    } else if let Some(rows) = bi_rows {
        // Two-segment piecewise-rational rows: EXACT arcs for
        // sweeps beyond the single-segment range (up to ~170° per
        // station). The sampled polyline skin's chordal deviation
        // (measured 8.7e-5 at r = 0.05 with 8-point arcs) was the
        // wall x strip mutual-disagreement floor on the boss
        // corner-adjacent rail.
        let mut control_points = Vec::with_capacity(n_sections * 5);
        let mut weights = Vec::with_capacity(n_sections * 5);
        for (cps, w) in &rows {
            control_points.extend_from_slice(cps);
            weights.extend_from_slice(&[1.0, *w, 1.0, *w, 1.0]);
        }
        (
            NurbsSurface {
                degree_u,
                degree_v: 2,
                knots_u: build_clamped_uniform_knots(n_sections, degree_u),
                knots_v: vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0],
                control_points,
                num_u: n_sections,
                num_v: 5,
                weights,
            },
            None,
        )
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

    // Chart overlap: extend curved-support lofts by exact rows past
    // the corner-ball station (all layouts; trims untouched).
    let (mut nurbs, continuity_downgrade) = (nurbs, continuity_downgrade);
    if let Some((s1, s2)) = supports {
        let curved = !matches!(s1, SurfaceGeometry::Plane { .. })
            || !matches!(s2, SurfaceGeometry::Plane { .. });
        // Conic sections skip the overhang extension: its rows are
        // derived from the circular contact construction and would
        // kink a conic loft at the edge-span boundary.
        if curved && nurbs.num_u >= 3 && conic_rho.is_none() {
            crate::loft_chart_extension::extend_loft_rows(&mut nurbs, s1, s2, radius);
        }
    }
    // PK-177 smooth loft for constant radius (boss 1.17e-5 -> 1.01e-5);
    // variable radius keeps the degree-1 loft, whose trims ARE its rails
    // (apply anchors need that parameterization: 52 opens otherwise).
    if supports.is_some()
        && nurbs.degree_u == 1
        && let Some((smooth, t1, t2)) = crate::station_geometry::smooth_reloft_with_rails(&nurbs)
    {
        return Ok(FilletSurfaceResult {
            surface: SurfaceGeometry::Nurbs(smooth),
            trim_curve_1: t1,
            trim_curve_2: t2,
            continuity_downgrade,
        });
    }
    let (trim_1, trim_2) = if degree_u == 1 {
        (
            crate::station_geometry::station_polyline(footprint_1)?,
            crate::station_geometry::station_polyline(footprint_2)?,
        )
    } else {
        (
            crate::center_curve::fit_nurbs_through_points(footprint_1)?,
            crate::center_curve::fit_nurbs_through_points(footprint_2)?,
        )
    };

    Ok(FilletSurfaceResult {
        surface: SurfaceGeometry::Nurbs(nurbs),
        trim_curve_1: trim_1,
        trim_curve_2: trim_2,
        continuity_downgrade,
    })
}

#[path = "fillet_surface_rows.rs"]
mod fillet_surface_rows;
pub(crate) use fillet_surface_rows::{
    build_clamped_uniform_knots, push_collapsed_row, rational_arc_row, sample_arc, strip_row_points,
};

#[cfg(test)]
#[path = "fillet_surface_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "fillet_surface_analytic_tests.rs"]
mod analytic_tests;

#[cfg(test)]
#[path = "fillet_surface_analytic_cylinder_tests.rs"]
mod analytic_cylinder_tests;
