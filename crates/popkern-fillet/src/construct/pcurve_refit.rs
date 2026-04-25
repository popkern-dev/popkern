//! Global pcurve re-fit pass.
//!
//! After topology construction (face splits + merges + surface
//! reassignments + edge reconciliation), some half-edge pcurves are
//! stale: they were authored against an earlier version of the face's
//! surface (or the parent face before a poke / split), and the same
//! 3-D curve traced through them no longer hits the right surface
//! points. Symptom: OCCT-strict slicers tessellate the pcurve into
//! a polyline that diverges from the actual edge geometry, producing
//! non-manifold edges in the imported mesh.
//!
//! This pass walks every edge with a 3-D curve, samples it densely,
//! and re-fits each adjacent half-edge's pcurve as a degree-3 NURBS
//! through the projected (u, v) samples on the *current* face surface.
//! Only replaces the pcurve if the new fit is geometrically closer
//! to the 3-D curve than the old one — a conservative heal that
//! never makes things worse.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::BRep;

/// Walk every edge in `brep`, refit pcurves on each adjacent half-
/// edge from the 3-D curve via closest-point projection. Only
/// installs the new pcurve if its surface-roundtrip error is smaller
/// than the existing pcurve's.
pub(crate) fn refit_all_pcurves(brep: &mut BRep) {
    let edge_handles: Vec<popkern_core::topology::EdgeHandle> = brep.edges.keys().collect();
    let n_samples = 12;
    for edge in edge_handles {
        let Some(curve_on_edge) = brep.edges.get(edge).and_then(|e| e.curve.clone()) else {
            continue;
        };
        let primary = brep.edges[edge].half_edge;
        let Some(he) = brep.half_edges.get(primary) else {
            continue;
        };
        let twin = he.twin;

        for &half in &[primary, twin] {
            let Some(he_data) = brep.half_edges.get(half) else {
                continue;
            };
            let face = brep.loops[he_data.loop_].face;
            let Some(surface) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            let surface = surface.geometry.clone();
            // Direction: does the half-edge's origin match the
            // curve's t_start or t_end? If it matches t_start, the
            // half-edge runs forward; otherwise reversed.
            let v_origin = he_data.origin;
            let p_origin = brep.vertices[v_origin].position;
            let p_t_start = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
            let p_t_end = curve_on_edge.geometry.point_at(curve_on_edge.t_end);
            let forward = p_origin.distance_to(&p_t_start) <= p_origin.distance_to(&p_t_end);

            let new_pcurve = match fit_pcurve_for_half_edge(
                &curve_on_edge.geometry,
                curve_on_edge.t_start,
                curve_on_edge.t_end,
                &surface,
                forward,
                n_samples,
            ) {
                Some(p) => p,
                None => continue,
            };
            // Compare new vs old pcurve roundtrip error against the
            // 3-D curve. Keep whichever is closer.
            let old_pcurve = brep.half_edges[half].pcurve.clone();
            let new_pcurve_pc: popkern_core::geometry::pcurve::Pcurve = new_pcurve.clone().into();
            let new_err =
                pcurve_roundtrip_error(&curve_on_edge.geometry, &surface, &new_pcurve_pc, forward);
            let old_err = match &old_pcurve {
                Some(pc) => pcurve_roundtrip_error(&curve_on_edge.geometry, &surface, pc, forward),
                None => f64::INFINITY,
            };
            // Only replace if the new pcurve is *meaningfully*
            // better (≥10% reduction) — avoids churn from random
            // numerical jitter.
            if new_err < old_err * 0.9 {
                brep.half_edges[half].pcurve = Some(new_pcurve.into());
            }
        }
        let _ = twin;
    }
}

fn fit_pcurve_for_half_edge(
    curve_3d: &popkern_core::geometry::curve::CurveGeometry,
    t_start: f64,
    t_end: f64,
    surface: &SurfaceGeometry,
    forward: bool,
    n_samples: usize,
) -> Option<NurbsCurve> {
    let mut points: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut params: Vec<f64> = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let s = i as f64 / (n_samples - 1) as f64;
        let t_s = if forward { s } else { 1.0 - s };
        let t = t_start + t_s * (t_end - t_start);
        let p_3d = curve_3d.point_at(t);
        let (u, v) = surface.closest_point(&p_3d)?;
        // Clamp to surface parameter range (Gregory etc. extrapolate
        // unreliably outside).
        let ((u_lo, u_hi), (v_lo, v_hi)) = surface.parameter_range();
        let u = if u_lo.is_finite() && u_hi.is_finite() {
            u.clamp(u_lo, u_hi)
        } else {
            u
        };
        let v = if v_lo.is_finite() && v_hi.is_finite() {
            v.clamp(v_lo, v_hi)
        } else {
            v
        };
        points.push(Point3::new(u, v, 0.0));
        params.push(s);
    }
    let degree = 3;
    let n_cps = (degree + 2).min(points.len());
    popkern_core::geometry::fitting::fit_curve(&points, &params, degree, n_cps).ok()
}

/// Sample the pcurve at 9 evenly-spaced parameters, evaluate the
/// surface at each (u, v), and return the max distance to the
/// corresponding 3-D curve point. `forward = false` means the
/// pcurve runs backwards relative to the 3-D curve.
fn pcurve_roundtrip_error(
    curve_3d: &popkern_core::geometry::curve::CurveGeometry,
    surface: &SurfaceGeometry,
    pcurve: &popkern_core::geometry::pcurve::Pcurve,
    forward: bool,
) -> f64 {
    let (pa_lo, pa_hi) = pcurve.parameter_range();
    let mut max_err = 0.0_f64;
    for i in 0..9 {
        let s = i as f64 / 8.0;
        let t_s = if forward { s } else { 1.0 - s };
        let t_p = pa_lo + s * (pa_hi - pa_lo);
        let uv = pcurve.point_at(t_p);
        let p_surf = surface.point_at(uv.x, uv.y);
        // The 3-D curve uses its own parameter range — sample by `t_s`.
        let (t0, t1) = curve_3d_param_range(curve_3d);
        let t = t0 + t_s * (t1 - t0);
        let p_3d = curve_3d.point_at(t);
        let d = p_surf.distance_to(&p_3d);
        if d > max_err {
            max_err = d;
        }
    }
    max_err
}

fn curve_3d_param_range(curve: &popkern_core::geometry::curve::CurveGeometry) -> (f64, f64) {
    use popkern_core::geometry::curve::CurveGeometry as CG;
    match curve {
        CG::Line { .. } => (0.0, 1.0),
        CG::Circle { .. } => (0.0, std::f64::consts::TAU),
        CG::Arc {
            start_angle,
            end_angle,
            ..
        } => (*start_angle, *end_angle),
        CG::Nurbs(n) => n.parameter_range(),
    }
}
