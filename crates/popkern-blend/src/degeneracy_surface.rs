//! Surface-level §13.1 degeneracy detectors, extracted from
//! `degeneracy.rs` to stay under the 500-line structural cap.
//!
//! Contains:
//!
//! - **§13.1(3) nearly-tangent supports** — dihedral angle at the
//!   edge midpoint below `2 * tolerance.angular()`.
//! - **§13.1(4) parameterization singularity** — NaN/∞ in
//!   `point_at` / `partial_u` / `partial_v` at any sample on the
//!   surface's UV grid.
//! - **§13.1(6) degenerate support normal + limit-normal repair**
//!   — vanishing normal at the edge midpoint, with a
//!   limit-normal-from-neighbours repair that averages compass-
//!   point neighbour normals.
//!
//! The structural detectors for collapsed edges, endpoint cusps,
//! and zero-area patches live in the parent `degeneracy.rs`.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_math::PortableMath;

use crate::analysis::EdgeAnalysis;
use crate::degeneracy::{DegeneracyKind, SupportSide};

/// §13.1(6) detection + **repair**: vanishing surface normal at
/// the edge midpoint, with a limit-normal-from-neighbours repair
/// applied transparently when the exact normal is degenerate but
/// nearby samples produce a coherent limit.
///
/// **Spec repair** (§13.1(6)): "limit-normal from neighbors". When
/// the exact `normal_at(u, v)` is degenerate on a support surface
/// (the edge midpoint projects to a sphere pole, cone apex, or
/// other isolated parameterisation singularity), [`limit_normal_at`]
/// averages normals from a small UV neighbourhood and returns the
/// result if coherent. When the limit exists, the detector treats
/// the support as usable and does NOT fire — the contact solver
/// will step away from the singular point during Newton iteration
/// and use the well-defined normals of the neighbourhood.
///
/// Only when **even the limit is degenerate** (a genuinely broken
/// surface or a sharp corner where neighbouring normals flip
/// direction and cancel) does the detector reject.
pub(crate) fn check_degenerate_support_normal(
    brep: &BRep,
    analysis: &EdgeAnalysis,
    tolerance: &Tolerance,
) -> Result<(), DegeneracyKind> {
    let midpoint = match crate::geometry_compute::get_edge_midpoint(brep, analysis) {
        Ok(p) => p,
        Err(_) => return Ok(()),
    };

    let check_side = |surface: &SurfaceGeometry, side: SupportSide| -> Result<(), DegeneracyKind> {
        let Some((u, v)) = surface.closest_point(&midpoint) else {
            return Ok(());
        };
        let normal = surface.normal_at(u, v);
        let mag = normal.length();
        if mag >= tolerance.linear() {
            return Ok(());
        }
        // §13.1(6) repair: try the limit-normal-from-neighbours.
        let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
        let u_span = if u_min.is_finite() && u_max.is_finite() {
            u_max - u_min
        } else {
            1.0
        };
        let v_span = if v_min.is_finite() && v_max.is_finite() {
            v_max - v_min
        } else {
            1.0
        };
        let step_u = u_span * 0.01;
        let step_v = v_span * 0.01;
        if limit_normal_at(surface, u, v, step_u, step_v, tolerance.linear()).is_some() {
            return Ok(());
        }
        Err(DegeneracyKind::DegenerateSupportNormal {
            side,
            magnitude: mag,
        })
    };

    check_side(&analysis.surface_1, SupportSide::Face1)?;
    check_side(&analysis.surface_2, SupportSide::Face2)?;
    Ok(())
}

/// Compute a limit normal at `(u, v)` by averaging unit normals from
/// compass-point neighbours in a small UV neighbourhood, clamped to
/// the surface's valid parameter range.
///
/// Used when the exact surface normal at `(u, v)` is degenerate
/// (||N|| < `tol`) but nearby points have well-defined normals.
/// The canonical case is an isolated parameterisation singularity
/// like a sphere pole or cone apex, where `partial_u × partial_v`
/// vanishes at the singular point itself but the surrounding
/// surface is smooth and well-normed.
///
/// **Algorithm.** Samples at offsets `(±step_u, 0)`, `(0, ±step_v)`,
/// and `(±step_u, ±step_v)` — eight points on a rectangle centred
/// at `(u, v)`. Samples that fall outside the surface's own
/// `parameter_range()` are **skipped** rather than evaluated,
/// because outside-range evaluation on a non-periodic surface gives
/// a mirror position whose normal flips sign and cancels the
/// in-range neighbours when summed (canonical example: sphere
/// pole at `v = 0`, where sampling at `v = -step` crosses into a
/// different chart). For each in-range sample whose own normal
/// magnitude is above `tol`, add the unit normal to the
/// accumulator. Returns `Some(unit_normal)` when the accumulated
/// coherence `|∑unit_n| / count ≥ 0.8` on at least one contributing
/// sample. Otherwise returns `None`.
///
/// **Why a coherence ratio?** The number of in-range neighbours
/// varies: 8 for an interior point on a closed surface, as few as
/// 3 at a sphere pole (where half the compass falls outside the
/// `v ≥ 0` half). A fixed `|∑| ≥ 0.5` threshold would fail at
/// boundaries even when the surviving neighbours agree perfectly.
/// Normalising by the count gives a threshold that scales to the
/// number of contributing samples: `|∑unit| / N = 1` means all
/// unit normals point in one direction (perfect coherence);
/// `|∑unit| / N = 0` means they cancel out. The 0.8 threshold
/// accepts neighbourhoods where neighbours agree within ≈ 37° on
/// average, which covers the typical pole/apex case while
/// rejecting sharp-corner cancellations.
pub(crate) fn limit_normal_at(
    surface: &SurfaceGeometry,
    u: f64,
    v: f64,
    step_u: f64,
    step_v: f64,
    tol: f64,
) -> Option<Vec3> {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    let offsets: [(f64, f64); 8] = [
        (step_u, 0.0),
        (-step_u, 0.0),
        (0.0, step_v),
        (0.0, -step_v),
        (step_u, step_v),
        (step_u, -step_v),
        (-step_u, step_v),
        (-step_u, -step_v),
    ];
    let mut accum = Vec3::new(0.0, 0.0, 0.0);
    let mut count = 0usize;
    for (du, dv) in offsets {
        let uu = u + du;
        let vv = v + dv;
        if u_min.is_finite() && uu < u_min {
            continue;
        }
        if u_max.is_finite() && uu > u_max {
            continue;
        }
        if v_min.is_finite() && vv < v_min {
            continue;
        }
        if v_max.is_finite() && vv > v_max {
            continue;
        }
        let n = surface.normal_at(uu, vv);
        let mag = n.length();
        if mag >= tol {
            accum = accum + n * (1.0 / mag);
            count += 1;
        }
    }
    if count == 0 {
        return None;
    }
    let sum_mag = accum.length();
    let coherence = sum_mag / (count as f64);
    if coherence < 0.8 {
        return None;
    }
    Some(accum * (1.0 / sum_mag))
}

/// §13.1(3) detection: nearly-tangent supports — the two face
/// normals are anti-parallel within `2 * tolerance.angular()`.
///
/// popkern face normals point INTO the solid. A dihedral angle of
/// π (flat, no edge) corresponds to `n1 · n2 = 1` (parallel); a
/// dihedral angle of 0 (knife-edge) corresponds to `n1 · n2 = -1`
/// (anti-parallel). We flag the knife-edge case because it is what
/// makes a rolling-ball fillet impossible to place: as the dihedral
/// shrinks, the rolling-ball centre escapes to infinity.
///
/// The flat case (dihedral ≈ π, `n1 · n2 ≈ 1`) is handled by the
/// existing `Convexity::Mixed` classification in
/// `analysis::analyze_edges` and is not re-checked here to avoid
/// double-reporting.
///
/// **Spec repair** (§13.1(3): "increase damping, regularize,
/// warn/reject"). The *reject* arm is this detector. The
/// *increase damping, regularize* arm is implemented by the
/// adaptive Levenberg-Marquardt regularization in
/// [`crate::contact::solve_rolling_ball_contact`]: when the
/// Jacobian condition number climbs above 1e8 (which it does for
/// tight-but-nonzero dihedrals well before this detector's hard
/// threshold fires), the solver automatically blends toward LM
/// gradient descent. The walker's adaptive step-size shrinkage
/// in [`crate::walker_adaptive`] complements it by taking
/// smaller dt steps on near-tangent spine samples. Together
/// these handle every tight geometry that *can* be filleted;
/// the hard reject here covers only the cases where the
/// dihedral literally vanishes.
pub(crate) fn check_nearly_tangent_supports(
    brep: &BRep,
    analysis: &EdgeAnalysis,
    tolerance: &Tolerance,
) -> Result<(), DegeneracyKind> {
    let midpoint = match crate::geometry_compute::get_edge_midpoint(brep, analysis) {
        Ok(p) => p,
        Err(_) => return Ok(()),
    };
    let Some((u1, v1)) = analysis.surface_1.closest_point(&midpoint) else {
        return Ok(());
    };
    let Some((u2, v2)) = analysis.surface_2.closest_point(&midpoint) else {
        return Ok(());
    };
    let n1 = analysis.surface_1.normal_at(u1, v1);
    let n2 = analysis.surface_2.normal_at(u2, v2);
    let m1 = n1.length();
    let m2 = n2.length();
    if m1 < tolerance.linear() || m2 < tolerance.linear() {
        return Ok(()); // already caught upstream by check_degenerate_support_normal
    }
    let n1u = n1 * (1.0 / m1);
    let n2u = n2 * (1.0 / m2);
    let dot = n1u.dot(&n2u);

    let sign_1 = if analysis.orientation_1 { 1.0 } else { -1.0 };
    let sign_2 = if analysis.orientation_2 { 1.0 } else { -1.0 };
    let face_dot = sign_1 * sign_2 * dot;

    let threshold = -(1.0 - 2.0 * tolerance.angular() * tolerance.angular());
    if face_dot < threshold {
        let angle_between_face_normals = face_dot.clamp(-1.0, 1.0).pacos();
        let dihedral = std::f64::consts::PI - angle_between_face_normals;
        return Err(DegeneracyKind::NearlyTangentSupports {
            angle_rad: dihedral,
        });
    }
    Ok(())
}

/// §13.1(4): NaN/∞ in surface evaluation at sampled (u, v) points.
///
/// Sweeps each support surface on an 8×8 UV grid over the natural
/// parameter range (infinite dimensions clamp to `[-10, 10]` so
/// the sweep terminates). At each sample checks `point_at`,
/// `partial_u`, and `partial_v` for non-finite coordinates. Fires
/// on the first non-finite sample found, so the cost is bounded
/// at 8×8×9 = 576 float checks in the clean case.
///
/// This detector catches genuinely broken surfaces — pathological
/// NURBS control polygons, bad knot vectors, arithmetic singularities
/// in derived formulas. popkern's analytic primitives return finite
/// values across their natural parameter range, so they will only
/// fire on the NURBS case.
pub(crate) fn check_parameterization_singularity(
    analysis: &EdgeAnalysis,
) -> Result<(), DegeneracyKind> {
    check_surface_singularity(&analysis.surface_1, SupportSide::Face1)?;
    check_surface_singularity(&analysis.surface_2, SupportSide::Face2)?;
    Ok(())
}

fn check_surface_singularity(
    surface: &SurfaceGeometry,
    side: SupportSide,
) -> Result<(), DegeneracyKind> {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    let u_lo = if u_min.is_finite() { u_min } else { -10.0 };
    let u_hi = if u_max.is_finite() { u_max } else { 10.0 };
    let v_lo = if v_min.is_finite() { v_min } else { -10.0 };
    let v_hi = if v_max.is_finite() { v_max } else { 10.0 };

    const N: usize = 8;
    let denom = (N - 1) as f64;
    for i in 0..N {
        let u = u_lo + (u_hi - u_lo) * (i as f64) / denom;
        for j in 0..N {
            let v = v_lo + (v_hi - v_lo) * (j as f64) / denom;
            let p = surface.point_at(u, v);
            let pu = surface.partial_u(u, v);
            let pv = surface.partial_v(u, v);
            if !p.x.is_finite()
                || !p.y.is_finite()
                || !p.z.is_finite()
                || !pu.x.is_finite()
                || !pu.y.is_finite()
                || !pu.z.is_finite()
                || !pv.x.is_finite()
                || !pv.y.is_finite()
                || !pv.z.is_finite()
            {
                return Err(DegeneracyKind::ParameterizationSingularity { side, uv: (u, v) });
            }
        }
    }
    Ok(())
}
