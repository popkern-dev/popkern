//! Face-level parameter boxes with period-aware unwrapping.
//!
//! Vertex-projected boxes ([`super::uv_projection::bounded_parameter_range`])
//! fail two ways on periodic surfaces: vertices can project into
//! different period branches (inflating the box to the full period
//! for a narrow band), and a full-wrap face's few vertices can all
//! sit AT the seam (collapsing the box to a sliver while the face
//! spans 2π — measured: a rod's full cylinder face collapsed to the
//! seam angle and a face blend seeded on the rod's FAR side).
//! Sampling the boundary curves and unwrapping consecutive angles
//! yields the true covered arc in both cases.

use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::uv_projection::{
    bounded_parameter_range, exact_bounded_parameter_range, project_point_to_uv,
};
use crate::geometry::{Point3, SurfaceEval};

type UvBox = ((f64, f64), (f64, f64));

/// Padded face box (50 %, min 1.0 — the intersection-solver headroom
/// convention of `bounded_parameter_range`; periodic u is never
/// padded past one period). Non-periodic surfaces with finite native
/// ranges return those unchanged.
pub fn face_bounded_parameter_range(
    brep: &crate::topology::BRep,
    face: crate::topology::FaceHandle,
    surface: &SurfaceGeometry,
) -> UvBox {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    let (u_period, _) = surface.parameter_periodicity();
    if u_period.is_none()
        && u_min.is_finite()
        && u_max.is_finite()
        && v_min.is_finite()
        && v_max.is_finite()
    {
        return ((u_min, u_max), (v_min, v_max));
    }
    let Some(((u_lo, u_hi), (v_lo, v_hi))) = sampled_face_box(brep, face, surface) else {
        return bounded_parameter_range(&brep.face_outer_loop_positions(face), surface);
    };
    let pad = |lo: f64, hi: f64| -> (f64, f64) {
        let w = (hi - lo).abs();
        let p = (0.5 * w).max(1.0);
        (lo - p, hi + p)
    };
    let (v0p, v1p) = pad(v_lo, v_hi);
    let (u0p, u1p) = if let Some(period) = u_period {
        if u_hi - u_lo >= period - 1e-9 {
            (u_lo, u_lo + period)
        } else {
            let w = u_hi - u_lo;
            let p = (0.5 * w).max(0.2).min(0.5 * (period - w));
            (u_lo - p, u_hi + p)
        }
    } else {
        pad(u_lo, u_hi)
    };
    ((u0p, u1p), (v0p, v1p))
}

/// Exact (unpadded) face box from sampled boundary curves — the true
/// UV extent of the face, period-unwrapped. Use for seeding and trim
/// gating where padding-band points must not count as on-face. Falls
/// back to the vertex-projected exact box when boundary sampling
/// yields nothing.
pub fn face_exact_parameter_range(
    brep: &crate::topology::BRep,
    face: crate::topology::FaceHandle,
    surface: &SurfaceGeometry,
) -> UvBox {
    sampled_face_box(brep, face, surface).unwrap_or_else(|| {
        exact_bounded_parameter_range(&brep.face_outer_loop_positions(face), surface)
    })
}

/// Raw sampled box: 8 samples per boundary edge curve (vertex
/// fallback), per-sample UV (cylinders analytically,
/// branch-canonical), consecutive-u unwrapping across the period.
fn sampled_face_box(
    brep: &crate::topology::BRep,
    face: crate::topology::FaceHandle,
    surface: &SurfaceGeometry,
) -> Option<UvBox> {
    let (u_period, _) = surface.parameter_periodicity();
    let mut samples: Vec<Point3> = Vec::new();
    let face_data = brep.faces.get(face)?;
    let start = brep.loops[face_data.outer_loop].half_edge;
    let mut he = start;
    for _ in 0..10_000 {
        let data = &brep.half_edges[he];
        if let Some(curve) = brep.edges[data.edge].curve.as_ref() {
            let ceval: &dyn CurveEval = &curve.geometry;
            let origin = brep.vertices[data.origin].position;
            let (t0, t1) = (curve.t_start, curve.t_end);
            let forward =
                ceval.point_at(t0).distance_to(&origin) <= ceval.point_at(t1).distance_to(&origin);
            for k in 0..8 {
                let f = k as f64 / 8.0;
                let t = if forward {
                    t0 + (t1 - t0) * f
                } else {
                    t1 - (t1 - t0) * f
                };
                samples.push(ceval.point_at(t));
            }
        } else {
            samples.push(brep.vertices[data.origin].position);
        }
        he = data.next;
        if he == start {
            break;
        }
    }
    if samples.is_empty() {
        return None;
    }

    let uv_of = |p: &Point3| -> Option<(f64, f64)> {
        if let SurfaceGeometry::Cylinder { origin, axis, .. } = surface {
            let (perp1, perp2) = crate::geometry::analytic_surfaces::orthonormal_basis(*axis);
            let d = *p - *origin;
            Some((d.dot(&perp2).atan2(d.dot(&perp1)), d.dot(axis)))
        } else {
            project_point_to_uv(surface, p)
        }
    };
    let mut us: Vec<f64> = Vec::with_capacity(samples.len());
    let mut v_lo = f64::INFINITY;
    let mut v_hi = f64::NEG_INFINITY;
    for p in &samples {
        if let Some((u, v)) = uv_of(p) {
            us.push(u);
            v_lo = v_lo.min(v);
            v_hi = v_hi.max(v);
        }
    }
    if us.is_empty() || !v_lo.is_finite() {
        return None;
    }
    let (mut u_lo, mut u_hi) = (us[0], us[0]);
    if let Some(period) = u_period {
        let mut walked = us[0];
        for &u in us.iter().skip(1) {
            let mut w = u;
            while w - walked > 0.5 * period {
                w -= period;
            }
            while w - walked < -0.5 * period {
                w += period;
            }
            walked = w;
            u_lo = u_lo.min(w);
            u_hi = u_hi.max(w);
        }
        // Full-wrap faces: clamp to exactly one period.
        if u_hi - u_lo > period {
            u_hi = u_lo + period;
        }
    } else {
        for &u in &us {
            u_lo = u_lo.min(u);
            u_hi = u_hi.max(u);
        }
    }
    Some(((u_lo, u_hi), (v_lo, v_hi)))
}
