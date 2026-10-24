//! The pre-marching coincidence fast path — split from
//! `surface_surface.rs` (500-line cap rule).

use super::{SurfaceSurfaceIntersection, SurfaceSurfaceIntersectionResult};
use crate::error::CoreError;
use crate::geometry::intersection::coincident::{
    closest_point_distance_pub, detect_coincident_regions,
};
use crate::geometry::intersection::coincident_analytic::are_surfaces_coincident_analytic;
use crate::geometry::intersection::plane_plane::try_plane_plane;
use crate::geometry::intersection::{
    cylinder_cone, cylinder_cylinder, cylinder_sphere, plane_cone, plane_cylinder, plane_sphere,
    plane_torus, radial_level_curve, torus_coaxial,
};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

/// Attempt to detect coincident surfaces before the expensive marching solver.
///
/// **Tier A** — analytic O(1) comparison for same-type analytic pairs.
/// **Tier B** — sample-based: 4-corner quick rejection followed by the
/// existing 8×8 grid check for NURBS / mixed-type pairs.
#[allow(clippy::too_many_arguments)]
pub(super) fn try_coincident_fast(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Option<SurfaceSurfaceIntersectionResult> {
    // Tier A: analytic comparison for same-type analytic surfaces.
    if are_surfaces_coincident_analytic(surface1, surface2, tolerance) {
        return Some(SurfaceSurfaceIntersectionResult {
            intersections: vec![SurfaceSurfaceIntersection::CoincidentRegion {
                u1_range,
                v1_range,
                u2_range,
                v2_range,
            }],
        });
    }

    // Tier B: sample-based detection for NURBS / mixed types.
    // Quick rejection: sample the 4 corners of surface1's parameter domain and
    // project each onto surface2. If any corner is far away, bail out early.
    let corners = [
        (u1_range.0, v1_range.0),
        (u1_range.0, v1_range.1),
        (u1_range.1, v1_range.0),
        (u1_range.1, v1_range.1),
    ];
    let rejection_threshold = tolerance.linear() * 10.0;
    for (u, v) in &corners {
        let p = surface1.point_at(*u, *v);
        let dist = closest_point_distance_pub(p, surface2, u2_range, v2_range);
        if dist > rejection_threshold {
            return None;
        }
    }

    // All corners are close — run the full 8×8 grid check.
    let regions = detect_coincident_regions(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    );
    if regions.is_empty() {
        None
    } else {
        Some(SurfaceSurfaceIntersectionResult {
            intersections: regions,
        })
    }
}

/// Analytic pair dispatch: every exact SSI fast path, tried in order,
/// ending with the coincidence detector. `Ok(Some(..))` is the exact
/// answer; `Ok(None)` falls through to the marching solver. Split
/// from `intersect_surfaces_inner` (500-line cap rule) — the chain's
/// ordering and per-arm rationale comments moved verbatim.
#[allow(clippy::too_many_arguments)]
pub(super) fn try_analytic_pair(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Result<Option<SurfaceSurfaceIntersectionResult>, CoreError> {
    // Fast path: both surfaces are planes → analytic solver.
    if let (
        SurfaceGeometry::Plane {
            origin: o1,
            normal: n1,
            u_axis: ua1,
            v_axis: va1,
        },
        SurfaceGeometry::Plane {
            origin: o2,
            normal: n2,
            u_axis: ua2,
            v_axis: va2,
        },
    ) = (surface1, surface2)
    {
        let intersections = try_plane_plane(
            *o1, *n1, *ua1, *va1, u1_range, v1_range, *o2, *n2, *ua2, *va2, u2_range, v2_range,
            tolerance,
        )?;
        return Ok(Some(SurfaceSurfaceIntersectionResult { intersections }));
    }

    // Fast path: plane-sphere pair → analytic solver.
    if let Some(intersections) = plane_sphere::try_plane_sphere(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    ) {
        return Ok(Some(SurfaceSurfaceIntersectionResult { intersections }));
    }

    // Fast path: plane-cylinder pair → analytic solver.
    if let Some(intersections) = plane_cylinder::try_plane_cylinder(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    ) {
        return Ok(Some(SurfaceSurfaceIntersectionResult { intersections }));
    }

    // Fast path: plane-cone pair → analytic solver.
    if let Some(intersections) = plane_cone::try_plane_cone(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    ) {
        return Ok(Some(SurfaceSurfaceIntersectionResult { intersections }));
    }

    // Fast path: COAXIAL cylinder-cone pair → exact circle(s); the
    // marching path fragmented this everyday class (PK-147).
    if let Some(intersections) = cylinder_cone::try_cylinder_cone(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    ) {
        return Ok(Some(SurfaceSurfaceIntersectionResult { intersections }));
    }

    // Fast path: COAXIAL cylinder-sphere pair → exact circles, or an
    // exact empty answer for the tangent / enclosed regimes. The
    // marching path re-marched a fillet strip × corner sphere tangent
    // circle from every one of ~2100 tangent seeds (PK-390).
    if let Some(intersections) = cylinder_sphere::try_cylinder_sphere(
        surface1,
        (u1_range, v1_range),
        surface2,
        (u2_range, v2_range),
        tolerance,
    ) {
        return Ok(Some(SurfaceSurfaceIntersectionResult { intersections }));
    }

    // Fast path: cylinder-cylinder with non-parallel axes → exact
    // ruling-swept branches. The marching path fragments the
    // cross-drill class into overlapping partial traces and loses
    // whole branches (PK-170).
    if let Some(intersections) =
        cylinder_cylinder::try_cylinder_cylinder(surface1, v1_range, surface2, v2_range, tolerance)
    {
        return Ok(Some(SurfaceSurfaceIntersectionResult { intersections }));
    }

    // PK-279 torus analytic arms: plane×torus in any orientation
    // (exact circles / sampled-analytic spiric sections) and coaxial
    // {cylinder, cone, sphere, torus}×torus (exact circles). A `None`
    // (non-torus pair, skew/offset axes, numeric-guard decline) keeps
    // the marching fallback; coincident coaxial tori fall through to
    // the coincident detector. POPKERN_ANALYTIC_TORUS_SSI=0 is the
    // diagnostic kill-switch (PK-239 precedent) pinning the marching
    // path for A/B measurement and containment.
    // Deliberately a per-call env read (NOT OnceLock'd like the trace
    // gates): the torus-SSI acceptance ratio test toggles this
    // kill-switch mid-process to A/B the marching path, and a cached
    // read froze its second measurement on the analytic arm. One read
    // per intersect call is negligible next to the solve itself.
    let torus_arms_enabled = std::env::var_os("POPKERN_ANALYTIC_TORUS_SSI")
        .map(|v| v != "0")
        .unwrap_or(true);
    if torus_arms_enabled {
        if let Some(intersections) = plane_torus::try_plane_torus(
            surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
        ) {
            return Ok(Some(SurfaceSurfaceIntersectionResult { intersections }));
        }
        if let Some(intersections) = torus_coaxial::try_torus_coaxial(
            surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
        ) {
            return Ok(Some(SurfaceSurfaceIntersectionResult { intersections }));
        }
    }

    // Radial level-curve arm: NURBS band × cylinder pairs whose
    // distance-to-axis field crosses the cylinder radius as a
    // transversal graph, scoped to multi-wrap traces (thread-groove
    // flank × shank — PK-461). The marched path priced an M6 cut at
    // 10–30 s and its 200k-step budget bounded one cut to ~35 turns.
    // Every gate failure falls through to marching;
    // POPKERN_ANALYTIC_HELICAL_SSI=0 pins marching for A/B.
    if let Some(intersections) = radial_level_curve::try_radial_level_curve(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    ) {
        return Ok(Some(SurfaceSurfaceIntersectionResult { intersections }));
    }

    // Fast path: coincident surfaces → return CoincidentRegion immediately.
    if let Some(result) = try_coincident_fast(
        surface1, u1_range, v1_range, surface2, u2_range, v2_range, tolerance,
    ) {
        return Ok(Some(result));
    }

    Ok(None)
}
