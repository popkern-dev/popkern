//! Saddle station pushback for the blend boundary (PK-234).
//!
//! At an `OppositeSide` saddle the strips' tangency rails CROSS on
//! their shared support faces: each neighbor strip consumes the
//! support band up to the crossing, so the closest-to-vertex arc
//! stations used by `extract_fillet_end_near_vertex` land past the
//! carved rail ends and the emitted bridges run through consumed
//! territory (measured on the notched-cube saddle: default station
//! x=0.9875 vs rail end x=0.9; bridge point 3.3e-2 off the neighbor
//! strip that owns the region). Convex corners never hit this — the
//! bridge region there lies beyond the rails, which is why the
//! default stations serve every cube/pyramid fixture.
//!
//! The pushback moves each fillet's arc station back past every
//! IN-STRIP rail crossing with its neighbors, plus a setback margin,
//! then rebuilds the arc through the pushed rail points (exact chart
//! iso where available). Bridges downstream rebuild automatically
//! from the new arc endpoints and cross real face regions.
//!
//! Crossings on the OUT-OF-STRIP side of a trim (the reflex pair's
//! tangent-line crossing — in-domain only because the trims carry
//! rolling-ball overshoot padding) impose nothing: the strips never
//! reach them, and their junction is a genuine bridge.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_math::PortableMath;

use crate::fillet_surface::FilletSurfaceResult;

/// A crossing between two trims can only push a station when its gap
/// is a genuine contact, not a near miss. Rail crossings measure
/// ≤ 1e-8 on the saddle; the wobble class on variable-radius corners
/// reaches ~3e-3.
const CROSSING_GAP_BAND: f64 = 1.0e-3;

/// Setback margin past the farthest crossing, as a fraction of the
/// crossing's arclength from the vertex-side station (the crossing
/// distance IS the neighbor's local setback scale).
const SETBACK_MARGIN_FRACTION: f64 = 0.25;

/// Push each fillet's end-arc station back past every in-strip rail
/// crossing with its neighbors and rebuild the arcs. No-op for
/// fillets with no in-strip crossings.
pub(super) fn apply_ring_station_pushback(
    fillet_indices: &[usize],
    fillet_surfaces: &[FilletSurfaceResult],
    vertex_pos: &Point3,
    fillet_arcs: &mut [NurbsCurve],
    fillet_endpoints: &mut [[Point3; 2]],
) -> Result<(), String> {
    let n = fillet_indices.len();
    // Required pushback arclength per (local fillet, trim slot).
    let mut required: Vec<[f64; 2]> = vec![[0.0; 2]; n];
    for a in 0..n {
        for b in (a + 1)..n {
            let fa = &fillet_surfaces[fillet_indices[a]];
            let fb = &fillet_surfaces[fillet_indices[b]];
            for (sa, ca) in [&fa.trim_curve_1, &fa.trim_curve_2].into_iter().enumerate() {
                for (sb, cb) in [&fb.trim_curve_1, &fb.trim_curve_2].into_iter().enumerate() {
                    let Some((gap, ta, tb)) = closest_pair_params(ca, cb) else {
                        continue;
                    };
                    if gap > CROSSING_GAP_BAND {
                        continue;
                    }
                    if let Some(s) = in_strip_arclength(ca, ta, vertex_pos) {
                        required[a][sa] = required[a][sa].max(s);
                    }
                    if let Some(s) = in_strip_arclength(cb, tb, vertex_pos) {
                        required[b][sb] = required[b][sb].max(s);
                    }
                }
            }
        }
    }

    // Two-phase: compute every pushed arc first, commit only when
    // all succeed — a partial pushback would ship a boundary whose
    // stations disagree with its bridges.
    let mut pushed: Vec<Option<(NurbsCurve, [Point3; 2])>> = vec![None; n];
    for k in 0..n {
        // The station must clear BOTH rails' worst crossings — the
        // arc spans rail to rail at one station.
        let s_needed = required[k][0].max(required[k][1]);
        if s_needed <= 0.0 {
            continue;
        }
        let s_pushed = s_needed * (1.0 + SETBACK_MARGIN_FRACTION);
        let fillet = &fillet_surfaces[fillet_indices[k]];
        let p1 = point_at_arclength_from_vertex(&fillet.trim_curve_1, vertex_pos, s_pushed)
            .ok_or_else(|| {
                format!("saddle pushback: fillet {k} trim 1 too short for {s_pushed:.4}")
            })?;
        let p2 = point_at_arclength_from_vertex(&fillet.trim_curve_2, vertex_pos, s_pushed)
            .ok_or_else(|| {
                format!("saddle pushback: fillet {k} trim 2 too short for {s_pushed:.4}")
            })?;
        // Exact station cross-section first (analytic cylinder strips
        // — the straight-edge fillet class), then the chart-edge iso
        // (NURBS sweeps whose station lands on a chart boundary),
        // then the interior iso (PK-239: NURBS sweep strips at a
        // pushed interior station — the exact section, since a
        // surface iso lies on the surface by identity), then the
        // Hermite on-surface fit (general strips, ~1e-3 κ class —
        // same fallback the default extraction uses; fails the ring
        // carve's 1e-5 on-surface postcondition, so reaching it
        // means the ring declines with the constructor named).
        let arc = exact_cylinder_station_arc(&fillet.surface, &p1, &p2)
            .or_else(|| super::boundary_hex::exact_boundary_iso_between(&fillet.surface, &p1, &p2))
            .or_else(|| exact_interior_iso_station_arc(&fillet.surface, &p1, &p2))
            .or_else(|| hermite_station_arc(&fillet.surface, &p1, &p2))
            .ok_or_else(|| {
                format!("saddle pushback: no station arc for fillet {k} at {s_pushed:.4}")
            })?;
        let (t0, t1) = arc.parameter_range();
        let endpoints = [arc.point_at(t0), arc.point_at(t1)];
        if popkern_core::trace::trace_enabled("corner-carve") {
            eprintln!(
                "  [ring-stations] fillet {k} pushed to arclength {s_pushed:.4} \
                 (crossings demanded {s_needed:.4}); arc ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4})",
                endpoints[0].x,
                endpoints[0].y,
                endpoints[0].z,
                endpoints[1].x,
                endpoints[1].y,
                endpoints[1].z,
            );
        }
        pushed[k] = Some((arc, endpoints));
    }
    for (k, entry) in pushed.into_iter().enumerate() {
        if let Some((arc, endpoints)) = entry {
            fillet_arcs[k] = arc;
            fillet_endpoints[k] = endpoints;
        }
    }
    Ok(())
}

/// Grid + window-refined closest approach between two trims,
/// unrestricted. Returns `(gap, t_on_a, t_on_b)`.
fn closest_pair_params(ca: &NurbsCurve, cb: &NurbsCurve) -> Option<(f64, f64, f64)> {
    const S: usize = 96;
    let (a0, a1) = ca.parameter_range();
    let (b0, b1) = cb.parameter_range();
    let mut best = (f64::INFINITY, a0, b0);
    for i in 0..=S {
        let ta = a0 + (a1 - a0) * i as f64 / S as f64;
        let pa = ca.point_at(ta);
        for j in 0..=S {
            let tb = b0 + (b1 - b0) * j as f64 / S as f64;
            let d = pa.distance_to(&cb.point_at(tb));
            if d < best.0 {
                best = (d, ta, tb);
            }
        }
    }
    let mut wa = (a1 - a0).abs() / S as f64;
    let mut wb = (b1 - b0).abs() / S as f64;
    for _ in 0..12 {
        let (ta_lo, ta_hi) = ((best.1 - wa).max(a0.min(a1)), (best.1 + wa).min(a0.max(a1)));
        for i in 0..=16 {
            let ta = ta_lo + (ta_hi - ta_lo) * i as f64 / 16.0;
            let d = ca.point_at(ta).distance_to(&cb.point_at(best.2));
            if d < best.0 {
                best = (d, ta, best.2);
            }
        }
        let (tb_lo, tb_hi) = ((best.2 - wb).max(b0.min(b1)), (best.2 + wb).min(b0.max(b1)));
        for j in 0..=16 {
            let tb = tb_lo + (tb_hi - tb_lo) * j as f64 / 16.0;
            let d = ca.point_at(best.1).distance_to(&cb.point_at(tb));
            if d < best.0 {
                best = (d, best.1, tb);
            }
        }
        wa *= 0.35;
        wb *= 0.35;
    }
    best.0.is_finite().then_some(best)
}

/// Arclength from the trim's closest-to-vertex station (polished,
/// shared with every station caller, PK-548) to `t_x`,
/// but only when `t_x` lies on the STRIP-INTERIOR side of that
/// station (toward the trim's domain midpoint — trims run the swept
/// edge's length plus overshoot, so the midpoint side is the strip).
/// Out-of-strip crossings return `None`.
fn in_strip_arclength(curve: &NurbsCurve, t_x: f64, vertex_pos: &Point3) -> Option<f64> {
    let (t0, t1) = curve.parameter_range();
    let t_v = super::boundary_hex::closest_param_on_nurbs_curve(curve, vertex_pos, 128);
    let t_mid = 0.5 * (t0 + t1);
    let interior_sign = (t_mid - t_v).signum();
    if ((t_x - t_v) * interior_sign) <= 0.0 {
        return None;
    }
    Some(polyline_arclength(curve, t_v, t_x))
}

/// Point on the trim at `arclength` from its closest-to-vertex
/// station, walking toward the domain midpoint. `None` when the trim
/// ends first.
fn point_at_arclength_from_vertex(
    curve: &NurbsCurve,
    vertex_pos: &Point3,
    arclength: f64,
) -> Option<Point3> {
    let (t0, t1) = curve.parameter_range();
    const S: usize = 512;
    let t_v = super::boundary_hex::closest_param_on_nurbs_curve(curve, vertex_pos, 128);
    let t_mid = 0.5 * (t0 + t1);
    let t_end = if t_mid >= t_v { t1.max(t0) } else { t1.min(t0) };
    let mut acc = 0.0;
    let mut prev = curve.point_at(t_v);
    for i in 1..=S {
        let t = t_v + (t_end - t_v) * i as f64 / S as f64;
        let p = curve.point_at(t);
        let step = prev.distance_to(&p);
        if acc + step >= arclength {
            let f = if step > 0.0 {
                (arclength - acc) / step
            } else {
                0.0
            };
            return Some(Point3::new(
                prev.x + (p.x - prev.x) * f,
                prev.y + (p.y - prev.y) * f,
                prev.z + (p.z - prev.z) * f,
            ));
        }
        acc += step;
        prev = p;
    }
    None
}

/// Chordal arclength between two params (256 samples).
fn polyline_arclength(curve: &NurbsCurve, t_a: f64, t_b: f64) -> f64 {
    const S: usize = 256;
    let mut acc = 0.0;
    let mut prev = curve.point_at(t_a);
    for i in 1..=S {
        let t = t_a + (t_b - t_a) * i as f64 / S as f64;
        let p = curve.point_at(t);
        acc += prev.distance_to(&p);
        prev = p;
    }
    acc
}

/// Exact rational-quadratic cross-section arc of an analytic
/// cylinder strip at one station, from `p1` to `p2` (both on the
/// cylinder at the same axial station, minor turn < π). Standard
/// conic construction: apex control point on the bisector at
/// `r / cos(θ/2)`, weights `[1, cos(θ/2), 1]`.
fn exact_cylinder_station_arc(
    surface: &popkern_core::geometry::SurfaceGeometry,
    p1: &Point3,
    p2: &Point3,
) -> Option<NurbsCurve> {
    use popkern_core::geometry::SurfaceGeometry;
    let SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    } = surface
    else {
        return None;
    };
    let axis = axis.normalized()?;
    let station = (*p1 - *origin).dot(&axis);
    let station_2 = (*p2 - *origin).dot(&axis);
    let band = 1.0e-6 * radius.max(1.0);
    if (station - station_2).abs() > band {
        return None; // not a single cross-section
    }
    let center = *origin + axis * (0.5 * (station + station_2));
    let a = *p1 - center;
    let b = *p2 - center;
    if (a.length() - radius).abs() > band || (b.length() - radius).abs() > band {
        return None;
    }
    let a_hat = a.normalized()?;
    let b_hat = b.normalized()?;
    let cos_theta = a_hat.dot(&b_hat).clamp(-1.0, 1.0);
    let theta = cos_theta.pacos();
    if !(1.0e-9..std::f64::consts::PI - 1.0e-6).contains(&theta) {
        return None; // degenerate or major arc
    }
    let half_cos = (0.5 * theta).pcos();
    let bisector = (a_hat + b_hat).normalized()?;
    let apex = center + bisector * (*radius / half_cos);
    Some(NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![*p1, apex, *p2],
        weights: vec![1.0, half_cos, 1.0],
    })
}

/// Exact interior-iso station arc for NURBS sweep strips (PK-239):
/// when the two pushed rail points project onto a common iso line of
/// the strip chart — same u at opposite v borders, or same v at
/// opposite u borders — the surface's exact iso curve at that
/// station IS the cross-section, on-surface by identity. The conic
/// fillet class needs this: its strips are rational sweeps whose
/// stations land between chart rows, where the Hermite fit's ~1e-3
/// curvature error fails the ring carve's 1e-5 postcondition.
fn exact_interior_iso_station_arc(
    surface: &popkern_core::geometry::SurfaceGeometry,
    p1: &Point3,
    p2: &Point3,
) -> Option<NurbsCurve> {
    use popkern_core::geometry::CurveEval as _;
    use popkern_core::geometry::SurfaceEval as _;
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    let popkern_core::geometry::SurfaceGeometry::Nurbs(n) = surface else {
        return None;
    };
    let (u1, v1) = project_point_to_uv(surface, p1)?;
    let (u2, v2) = project_point_to_uv(surface, p2)?;
    let ((u_lo, u_hi), (v_lo, v_hi)) = n.parameter_range();
    let u_span = (u_hi - u_lo).abs().max(1.0e-12);
    let v_span = (v_hi - v_lo).abs().max(1.0e-12);
    // Same-iso decision band: parametric, relative to the chart span.
    // Loose on purpose — a wrong candidate is rejected by the
    // endpoint identity gate below, never returned.
    const ISO_PARAM_BAND: f64 = 1.0e-3;
    let mut candidates: Vec<NurbsCurve> = Vec::new();
    if (u1 - u2).abs() <= ISO_PARAM_BAND * u_span {
        candidates.push(n.iso_curve_at_u(0.5 * (u1 + u2)));
    }
    if (v1 - v2).abs() <= ISO_PARAM_BAND * v_span {
        candidates.push(n.iso_curve_at_v(0.5 * (v1 + v2)));
    }
    let scale = p1.distance_to(p2).max(1e-9);
    let gate = scale * 1.0e-6;
    for candidate in candidates {
        let (t0, t1) = candidate.parameter_range();
        let (c0, c1) = (candidate.point_at(t0), candidate.point_at(t1));
        if c0.distance_to(p1) <= gate && c1.distance_to(p2) <= gate {
            return Some(candidate);
        }
        if c0.distance_to(p2) <= gate && c1.distance_to(p1) <= gate {
            return Some(NurbsCurve {
                degree: candidate.degree,
                knots: super::boundary_hex::reversed_clamped_knots(&candidate.knots),
                control_points: candidate.control_points.iter().rev().copied().collect(),
                weights: candidate.weights.iter().rev().copied().collect(),
            });
        }
    }
    None
}

/// Hermite on-surface station arc for general strips: project the
/// pushed rail points into UV and reuse the boundary extraction's
/// on-surface cubic fit (periodicity-unwrapped `du`, same as the
/// default path).
fn hermite_station_arc(
    surface: &popkern_core::geometry::SurfaceGeometry,
    p1: &Point3,
    p2: &Point3,
) -> Option<NurbsCurve> {
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    let (u1, v1) = project_point_to_uv(surface, p1)?;
    let (u2, v2) = project_point_to_uv(surface, p2)?;
    let (u_periodicity, _) = surface.periodicity();
    let mut du = u2 - u1;
    if let Some(period) = u_periodicity
        && du.abs() > period * 0.5
    {
        du -= period * du.signum();
    }
    super::boundary_hex::fit_hermite_cubic_arc_on_surface(surface, u1, v1, u2, v2, du).ok()
}
