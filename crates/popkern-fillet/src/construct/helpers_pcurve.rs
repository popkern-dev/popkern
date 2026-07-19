//! Pcurve projection helpers shared by the construction phases:
//! straight and trim-following degree-1 pcurves, plus the UV
//! projection fallback for points off the analytic inverse.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{BRep, FaceHandle};

/// Build a simple linear NURBS curve between two points.
pub(crate) fn line_nurbs(a: &Point3, b: &Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![*a, *b],
        weights: vec![1.0, 1.0],
    }
}

/// Build a pcurve (UV-space line) between two points projected onto a surface.
pub(crate) fn compute_pcurve_for_trim(
    surface: &SurfaceGeometry,
    pt_a: &Point3,
    pt_b: &Point3,
) -> NurbsCurve {
    let (u_a, v_a) = project_to_uv(surface, pt_a);
    let (u_b, v_b) = project_to_uv(surface, pt_b);
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(u_a, v_a, 0.0), Point3::new(u_b, v_b, 0.0)],
        weights: vec![1.0, 1.0],
    }
}

/// Project a 3D point onto a surface's UV parameter space.
pub(crate) fn project_to_uv(surface: &SurfaceGeometry, point: &Point3) -> (f64, f64) {
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    match project_point_to_uv(surface, point) {
        Some(uv) => uv,
        None => {
            // Fallback: grid search.
            use popkern_core::geometry::surface::SurfaceEval;
            let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
            let u_lo = if u_min.is_finite() { u_min } else { -10.0 };
            let u_hi = if u_max.is_finite() { u_max } else { 10.0 };
            let v_lo = if v_min.is_finite() { v_min } else { -10.0 };
            let v_hi = if v_max.is_finite() { v_max } else { 10.0 };
            let n = 20;
            let mut best = ((u_lo + u_hi) * 0.5, (v_lo + v_hi) * 0.5);
            let mut best_d = f64::INFINITY;
            for iu in 0..n {
                let u = u_lo + (u_hi - u_lo) * (iu as f64 / (n - 1) as f64);
                for iv in 0..n {
                    let v = v_lo + (v_hi - v_lo) * (iv as f64 / (n - 1) as f64);
                    let d = surface.point_at(u, v).distance_to(point);
                    if d < best_d {
                        best_d = d;
                        best = (u, v);
                    }
                }
            }
            best
        }
    }
}

/// Compute a pcurve for a trim curve on the given face surface.
/// Endpoint-projected straight pcurve (the trim-aware variant is
/// [`pcurve_for_face_along`]; this one remains for exact 2-point
/// trims and its unit tests).
#[cfg_attr(not(test), allow(dead_code))]
pub fn pcurve_for_face(
    brep: &BRep,
    face: FaceHandle,
    trim_start: &Point3,
    trim_end: &Point3,
) -> Result<NurbsCurve, String> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("face {face:?} not found"))?;
    let surface = face_data
        .surface
        .as_ref()
        .ok_or_else(|| format!("face {face:?} has no surface"))?;
    Ok(compute_pcurve_for_trim(
        &surface.geometry,
        trim_start,
        trim_end,
    ))
}

/// Pcurve for a face-splitting trim that follows a CURVED 3-D trim:
/// samples the curve and projects each sample into the face's chart
/// (degree-1 through the projections, uniform knots — matching the
/// polyline sub-trim authored by `multi_chord`). Straight trims (2
/// control points) delegate to the endpoint projection so analytic
/// splits keep their exact linear pcurves.
pub fn pcurve_for_face_along(
    brep: &BRep,
    face: FaceHandle,
    trim_3d: &NurbsCurve,
) -> Result<NurbsCurve, String> {
    use popkern_core::geometry::CurveEval;
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("face {face:?} not found"))?;
    let surface = face_data
        .surface
        .as_ref()
        .ok_or_else(|| format!("face {face:?} has no surface"))?;
    let (t0, t1) = trim_3d.parameter_range();
    if trim_3d.control_points.len() <= 2 {
        let a = trim_3d.point_at(t0);
        let b = trim_3d.point_at(t1);
        return Ok(compute_pcurve_for_trim(&surface.geometry, &a, &b));
    }
    let n = trim_3d.control_points.len().max(3);
    let mut uv_points = Vec::with_capacity(n);
    for k in 0..n {
        let t = t0 + (t1 - t0) * k as f64 / (n - 1) as f64;
        let p = trim_3d.point_at(t);
        let (u, v) =
            popkern_core::geometry::uv_projection::project_point_to_uv(&surface.geometry, &p)
                .ok_or_else(|| format!("trim sample does not project onto face {face:?}"))?;
        uv_points.push(Point3::new(u, v, 0.0));
    }
    let degree = 1;
    let n_pts = uv_points.len();
    let mut knots = Vec::with_capacity(n_pts + 2);
    knots.push(0.0);
    for k in 0..n_pts {
        knots.push(k as f64 / (n_pts - 1) as f64);
    }
    knots.push(1.0);
    Ok(NurbsCurve {
        degree,
        knots,
        control_points: uv_points,
        weights: vec![1.0; n_pts],
    })
}
