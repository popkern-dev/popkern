//! Co-parameterized chart pcurve fitting for analytic surfaces.
//!
//! Given a 3D curve and a surface with a closed-form chart inverse,
//! builds a pcurve `fit` with `surface(fit(t)) ≈ curve(t)` for the
//! SAME `t` — the co-parameterization contract downstream chart
//! embeddings, UV clipping, and the STEP exporter rely on. This is
//! the `plane_sphere` fitting pattern promoted to a shared API
//! (PK-279 introduced it for the torus intersection arms; PK-292
//! replaced the closed-fillet rim's 17-sample polyline author with it
//! after the polyline's chart sag put a boolean trim endpoint 1.8e-5
//! off the rim arc and the split stage rejected the junction).

use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::uv_projection::project_point_to_uv;
use crate::tolerance::Tolerance;

/// Fit the chart image of `curve` over `range` on `surface` as a
/// co-parameterized degree-1 NURBS.
///
/// Samples the curve at parameters spanning `range`, charts each
/// sample with the surface's closed-form inverse, unwraps periodic
/// chart axes by continuity (period branch within half a period of
/// the previous sample), and interpolates with knots AT the sample
/// parameters — so the pcurve keeps the curve's own
/// parameterization. Sample density doubles until every mid-knot
/// deviation passes the linear tolerance or the 4097-point cap is
/// reached (chord sag shrinks quadratically per doubling; the cap
/// matches the `plane_sphere` precedent and its recorded
/// piecewise-quadratic upgrade candidate).
///
/// Chart-singularity guard: a sample whose inverse fails (pole, apex
/// generator) keeps its neighbour's chart value; if the FIRST sample
/// fails, or any chart value is non-finite, the fit returns `None`
/// and the caller keeps its fallback.
pub fn fit_chart_pcurve(
    surface: &SurfaceGeometry,
    curve: &dyn CurveEval,
    range: (f64, f64),
    tolerance: &Tolerance,
) -> Option<NurbsCurve> {
    let (t0, t1) = range;
    if !(t1 - t0).is_finite() || (t1 - t0).abs() < tolerance.numerical_floor() {
        return None;
    }
    let gate = tolerance.linear().max(1.0e-9);
    let (period_u, period_v) = surface.parameter_periodicity();

    let unwrap = |x: f64, prev: Option<f64>, period: Option<f64>| -> f64 {
        match (prev, period) {
            (Some(p), Some(period)) => x - period * ((x - p) / period).round(),
            _ => x,
        }
    };

    let mut n = 65usize;
    loop {
        let mut params = Vec::with_capacity(n);
        let mut uvs: Vec<(f64, f64)> = Vec::with_capacity(n);
        for i in 0..n {
            let t = t0 + (t1 - t0) * i as f64 / (n - 1) as f64;
            let prev = uvs.last().copied();
            let (u, v) = match project_point_to_uv(surface, &curve.point_at(t)) {
                Some((u, v)) => (
                    unwrap(u, prev.map(|(u, _)| u), period_u),
                    unwrap(v, prev.map(|(_, v)| v), period_v),
                ),
                None => prev?,
            };
            if !u.is_finite() || !v.is_finite() {
                return None;
            }
            params.push(t);
            uvs.push((u, v));
        }
        // Mid-knot verification: the polyline evaluated halfway
        // between knots must land back on the 3D curve.
        let mut worst = 0.0f64;
        for i in 0..n - 1 {
            let t_mid = 0.5 * (params[i] + params[i + 1]);
            let (u, v) = (
                0.5 * (uvs[i].0 + uvs[i + 1].0),
                0.5 * (uvs[i].1 + uvs[i + 1].1),
            );
            let on_surface = surface.point_at(u, v);
            worst = worst.max(on_surface.distance_to(&curve.point_at(t_mid)));
        }
        if worst <= gate || n >= 4097 {
            let mut knots = Vec::with_capacity(n + 2);
            knots.push(params[0]);
            knots.extend(params.iter().copied());
            knots.push(params[n - 1]);
            return Some(NurbsCurve {
                degree: 1,
                knots,
                control_points: uvs
                    .into_iter()
                    .map(|(u, v)| Point3::new(u, v, 0.0))
                    .collect(),
                weights: vec![1.0; n],
            });
        }
        n = (n - 1) * 2 + 1;
    }
}

/// Exact pcurve candidate for a curve lying in a Plane face. Dispatches
/// per curve form; returns `None` when no exact construction applies
/// (the caller falls back to the sampled fit, and the survival rule
/// measures whichever candidate is produced — exactness is verified,
/// not assumed).
pub fn exact_plane_pcurve(
    surface: &SurfaceGeometry,
    curve: &CurveGeometry,
    t_lo: f64,
    t_hi: f64,
) -> Option<Pcurve> {
    exact_plane_arc_pcurve(surface, curve, t_lo, t_hi)
        .or_else(|| exact_plane_nurbs_pcurve(surface, curve, t_lo, t_hi))
}

/// Exact angle-parameterized elliptical pcurve for an Arc or full
/// Circle lying in a Plane face — the boss-ring rail class.
///
/// `Arc::point_at` IS angle-linear, so an angle-parameterized
/// elliptical pcurve is rate-consistent with the curve; the generic
/// 23-sample cubic interpolation leaves ~`r·(Δθ)⁴/384` of sagitta
/// instead (measured 2.5e-4..9.2e-4 on the corpus boss-ring rails —
/// above the 1e-5 commercial trim gate, which read three watertight
/// parts as achieved-G0 with silent requested-G1 downgrades).
/// Projection uses the same 2×2 SVD decomposition as the plane×torus
/// analytic edge triple.
fn exact_plane_arc_pcurve(
    surface: &SurfaceGeometry,
    curve: &CurveGeometry,
    t_lo: f64,
    t_hi: f64,
) -> Option<Pcurve> {
    let SurfaceGeometry::Plane {
        origin,
        u_axis,
        v_axis,
        ..
    } = surface
    else {
        return None;
    };
    let (center, normal, radius) = match curve {
        CurveGeometry::Arc {
            center,
            normal,
            radius,
            ..
        }
        | CurveGeometry::Circle {
            center,
            normal,
            radius,
        } => (center, normal, radius),
        _ => return None,
    };
    let (bu, bv) = crate::geometry::analytic_curves::circle_basis(*normal);
    let ul2 = u_axis.dot(u_axis);
    let vl2 = v_axis.dot(v_axis);
    if ul2 < 1e-30 || vl2 < 1e-30 {
        return None;
    }
    let cen: crate::geometry::Vec3 = *center - Point3::new(origin.x, origin.y, origin.z);
    let (a_uu, a_vu) = (radius * bu.dot(u_axis) / ul2, radius * bu.dot(v_axis) / vl2);
    let (a_uv, a_vv) = (radius * bv.dot(u_axis) / ul2, radius * bv.dot(v_axis) / vl2);
    let (u_radius, v_radius, rotation, phase) =
        crate::geometry::edge_triple::decompose_2x2(a_uu, a_uv, a_vu, a_vv);
    Some(Pcurve::elliptical(
        (cen.dot(u_axis) / ul2, cen.dot(v_axis) / vl2),
        u_radius,
        v_radius,
        rotation,
        t_lo - phase,
        t_hi - phase,
        t_lo,
        t_hi,
    ))
}

/// Exact NURBS pcurve for a NURBS curve lying in a Plane face — the
/// boolean hole-rim class (PK-287).
///
/// A plane's UV map is affine, and affine maps commute with rational
/// combinations: mapping every control point into UV while keeping
/// degree, knots, and weights yields a pcurve whose trace and
/// parameterization BOTH match the 3-D curve exactly. The generic
/// 23-sample projection walk + cubic interpolation leaves ~1.4e-4 of
/// interpolation sagitta on a r=0.3 rim circle instead — above the
/// stale-pcurve gate, so every repair pass re-flagged the rim forever
/// (`fix_stale_pcurves` could never measure its own recompute better).
/// Same orthogonal-axes projection convention as the arc path; the
/// caller's survival rule verifies the candidate by measurement, so a
/// non-orthogonal plane frame degrades to the sampled fit rather than
/// shipping a wrong exact-looking pcurve.
fn exact_plane_nurbs_pcurve(
    surface: &SurfaceGeometry,
    curve: &CurveGeometry,
    t_lo: f64,
    t_hi: f64,
) -> Option<Pcurve> {
    let (
        SurfaceGeometry::Plane {
            origin,
            u_axis,
            v_axis,
            ..
        },
        CurveGeometry::Nurbs(nurbs),
    ) = (surface, curve)
    else {
        return None;
    };
    let ul2 = u_axis.dot(u_axis);
    let vl2 = v_axis.dot(v_axis);
    if ul2 < 1e-30 || vl2 < 1e-30 {
        return None;
    }
    let mapped: Vec<Point3> = nurbs
        .control_points
        .iter()
        .map(|p| {
            let rel: crate::geometry::Vec3 = *p - Point3::new(origin.x, origin.y, origin.z);
            Point3::new(rel.dot(u_axis) / ul2, rel.dot(v_axis) / vl2, 0.0)
        })
        .collect();
    let full_domain = Pcurve {
        geometry: crate::geometry::pcurve::PcurveGeometry::Nurbs(
            crate::geometry::nurbs_curve::NurbsCurve {
                degree: nurbs.degree,
                knots: nurbs.knots.clone(),
                control_points: mapped,
                weights: nurbs.weights.clone(),
            },
        ),
        t_start: t_lo,
        t_end: t_hi,
    };
    // Normalize to the authored-pcurve structural contract: knot
    // domain == parameter range, forward along the half-edge. The
    // mapped NURBS spans the WHOLE curve domain and may run backward
    // (twin side); downstream consumers (UV walker, split matching)
    // are built against sliced-forward pcurves, and shipping the
    // full-circle domain broke the flush-union chain's face splits
    // (measured: flush_arm/flush_gusset SplitFailed). `to_nurbs`
    // slices and reverses by exact knot operations.
    Some(Pcurve::from_nurbs(full_domain.to_nurbs()))
}

/// Analytic-first chart pcurve author: exact constructions where the
/// chart admits them, `None` otherwise (callers fall back to
/// [`fit_chart_pcurve`]).
///
/// Covered exactly: plane charts — angle-parameterized elliptical
/// arcs and affine NURBS images ([`exact_plane_pcurve`], the PK-210
/// machinery). Iso-line images on periodic charts need NO analytic
/// arm: the adaptive fit interpolates collinear chart samples with
/// zero error at its first density, so a fitted iso-line is already
/// exact. (An explicit `Pcurve::Linear` arm was tried and REMOVED:
/// full-period angle-ranged Linear pcurves hung a boolean UV-clip
/// consumer — recorded under PK-292 — while buying no accuracy.)
///
/// Exactness matters beyond absolute accuracy: two seam sides
/// authored exactly agree to machine precision, while independently
/// FITTED sides carry uncorrelated chord sag and the pcurve-true
/// two-sided continuity evidence reads the difference (measured:
/// replacing the closed-fillet author's correlated 17-sample
/// polylines with per-side adaptive fits moved the corpus
/// bossed_housing G0 seam evidence from ~1e-10 to 1.977e-6 — the
/// exact plane arm restores machine-level agreement AND absolute
/// accuracy, PK-292).
pub fn analytic_chart_pcurve(
    surface: &SurfaceGeometry,
    curve_geometry: &CurveGeometry,
    range: (f64, f64),
) -> Option<Pcurve> {
    exact_plane_pcurve(surface, curve_geometry, range.0, range.1)
}
