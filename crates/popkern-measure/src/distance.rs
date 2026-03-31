//! Minimum distance between entity pairs.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::tight_bounded_parameter_range;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::MeasureError;

#[path = "distance_newton.rs"]
mod newton;
use newton::{refine_curve_param, refine_surface_params};

/// A tagged union for any topology handle type.
#[derive(Debug, Clone, Copy)]
pub enum TopoHandle {
    /// A vertex handle.
    Vertex(VertexHandle),
    /// An edge handle.
    Edge(EdgeHandle),
    /// A face handle.
    Face(FaceHandle),
}

/// Compute the minimum distance between two topology entities.
///
/// Dispatches to specialized sampling-based sub-functions depending on the
/// entity-pair type.  The two entities may belong to different `BRep` solids.
///
/// # Errors
///
/// Returns [`MeasureError::InvalidHandle`] if either handle is stale, or
/// [`MeasureError::MissingCurve`] / [`MeasureError::MissingSurface`] if the
/// edge or face lacks geometric data.
pub fn min_distance(
    brep_a: &BRep,
    handle_a: TopoHandle,
    brep_b: &BRep,
    handle_b: TopoHandle,
    tolerance: &Tolerance,
) -> Result<f64, MeasureError> {
    match (handle_a, handle_b) {
        (TopoHandle::Vertex(va), TopoHandle::Vertex(vb)) => vertex_vertex(brep_a, va, brep_b, vb),
        (TopoHandle::Vertex(v), TopoHandle::Edge(e)) => vertex_edge(brep_a, v, brep_b, e),
        (TopoHandle::Edge(e), TopoHandle::Vertex(v)) => vertex_edge(brep_b, v, brep_a, e),
        (TopoHandle::Vertex(v), TopoHandle::Face(f)) => {
            vertex_face(brep_a, v, brep_b, f, tolerance)
        }
        (TopoHandle::Face(f), TopoHandle::Vertex(v)) => {
            vertex_face(brep_b, v, brep_a, f, tolerance)
        }
        (TopoHandle::Edge(ea), TopoHandle::Edge(eb)) => edge_edge(brep_a, ea, brep_b, eb),
        (TopoHandle::Edge(e), TopoHandle::Face(f)) => edge_face(brep_a, e, brep_b, f, tolerance),
        (TopoHandle::Face(f), TopoHandle::Edge(e)) => edge_face(brep_b, e, brep_a, f, tolerance),
        (TopoHandle::Face(fa), TopoHandle::Face(fb)) => {
            face_face(brep_a, fa, brep_b, fb, tolerance)
        }
    }
}

// ---------------------------------------------------------------------------
// Dispatch helpers
// ---------------------------------------------------------------------------

/// Minimum distance between two vertices.
fn vertex_vertex(
    brep_a: &BRep,
    va: VertexHandle,
    brep_b: &BRep,
    vb: VertexHandle,
) -> Result<f64, MeasureError> {
    let pa = brep_a
        .vertices
        .get(va)
        .ok_or_else(|| MeasureError::InvalidHandle("vertex".into()))?
        .position;
    let pb = brep_b
        .vertices
        .get(vb)
        .ok_or_else(|| MeasureError::InvalidHandle("vertex".into()))?
        .position;
    Ok(pa.distance_to(&pb))
}

/// Minimum distance from a vertex to an edge.
///
/// Samples the edge curve at 64 uniformly-spaced parameter values, then
/// refines the closest sample with Newton iteration for a precise result.
fn vertex_edge(
    brep_v: &BRep,
    v: VertexHandle,
    brep_e: &BRep,
    e: EdgeHandle,
) -> Result<f64, MeasureError> {
    let p = brep_v
        .vertices
        .get(v)
        .ok_or_else(|| MeasureError::InvalidHandle("vertex".into()))?
        .position;
    let coe = brep_e
        .edges
        .get(e)
        .ok_or_else(|| MeasureError::InvalidHandle("edge".into()))?
        .curve
        .as_ref()
        .ok_or(MeasureError::MissingCurve)?;

    let t0 = coe.t_start;
    let t1 = coe.t_end;
    const N: usize = 64;
    let mut min_d = f64::MAX;
    let mut t_best = t0;
    for i in 0..=N {
        let t = t0 + (t1 - t0) * (i as f64 / N as f64);
        let q = coe.geometry.point_at(t);
        let d = p.distance_to(&q);
        if d < min_d {
            min_d = d;
            t_best = t;
        }
    }

    // Newton refinement: find t where (C(t) - P) · C'(t) = 0.
    t_best = refine_curve_param(&coe.geometry, t_best, t0, t1, [p.x, p.y, p.z]);
    let final_pt = coe.geometry.point_at(t_best);
    Ok(p.distance_to(&final_pt))
}

/// Minimum distance from a vertex to a face.
///
/// Samples the outer boundary vertices plus a 16×16 UV grid over the face's
/// bounded parameter domain, then refines with Newton iteration.
fn vertex_face(
    brep_v: &BRep,
    v: VertexHandle,
    brep_f: &BRep,
    f: FaceHandle,
    _tolerance: &Tolerance,
) -> Result<f64, MeasureError> {
    let p = brep_v
        .vertices
        .get(v)
        .ok_or_else(|| MeasureError::InvalidHandle("vertex".into()))?
        .position;

    let face_data = brep_f
        .faces
        .get(f)
        .ok_or_else(|| MeasureError::InvalidHandle("face".into()))?;
    let sof = face_data
        .surface
        .as_ref()
        .ok_or(MeasureError::MissingSurface)?;

    let positions = brep_f.face_outer_loop_positions(f);
    let (u_range, v_range) = tight_bounded_parameter_range(&positions, &sof.geometry);
    let (u0, u1) = u_range;
    let (v0, v1) = v_range;

    let mut min_d = f64::MAX;
    let mut u_best = u0;
    let mut v_best = v0;

    // Sample boundary positions directly.
    for q in &positions {
        let d = p.distance_to(q);
        if d < min_d {
            min_d = d;
            // Boundary positions don't map back to (u,v) easily; keep grid seed.
        }
    }

    // Sample 16×16 UV grid.
    const N: usize = 16;
    for i in 0..=N {
        for j in 0..=N {
            let u = u0 + (u1 - u0) * (i as f64 / N as f64);
            let vv = v0 + (v1 - v0) * (j as f64 / N as f64);
            let q = sof.geometry.point_at(u, vv);
            let d = p.distance_to(&q);
            if d < min_d {
                min_d = d;
                u_best = u;
                v_best = vv;
            }
        }
    }

    // Newton refinement: find (u, v) where the foot-point condition holds.
    (u_best, v_best) = refine_surface_params(
        &sof.geometry,
        u_best,
        v_best,
        [u0, u1, v0, v1],
        [p.x, p.y, p.z],
    );
    let final_sp = sof.geometry.point_at(u_best, v_best);
    Ok(p.distance_to(&final_sp))
}

/// Minimum distance between two edges.
///
/// Samples each curve at 64 points to find the closest pair, then refines
/// both parameter values with Newton iteration.
fn edge_edge(
    brep_a: &BRep,
    ea: EdgeHandle,
    brep_b: &BRep,
    eb: EdgeHandle,
) -> Result<f64, MeasureError> {
    let coe_a = brep_a
        .edges
        .get(ea)
        .ok_or_else(|| MeasureError::InvalidHandle("edge".into()))?
        .curve
        .as_ref()
        .ok_or(MeasureError::MissingCurve)?;
    let coe_b = brep_b
        .edges
        .get(eb)
        .ok_or_else(|| MeasureError::InvalidHandle("edge".into()))?
        .curve
        .as_ref()
        .ok_or(MeasureError::MissingCurve)?;

    const N: usize = 64;
    let mut min_d = f64::MAX;
    let mut ta_best = coe_a.t_start;
    let mut tb_best = coe_b.t_start;

    for i in 0..=N {
        let ta = coe_a.t_start + (coe_a.t_end - coe_a.t_start) * (i as f64 / N as f64);
        let pa = coe_a.geometry.point_at(ta);
        for j in 0..=N {
            let tb = coe_b.t_start + (coe_b.t_end - coe_b.t_start) * (j as f64 / N as f64);
            let pb = coe_b.geometry.point_at(tb);
            let d = pa.distance_to(&pb);
            if d < min_d {
                min_d = d;
                ta_best = ta;
                tb_best = tb;
            }
        }
    }

    // Newton refinement: minimise |C_a(ta) - C_b(tb)|^2 jointly.
    // Conditions: (C_a - C_b)·C_a' = 0 and (C_a - C_b)·C_b' = 0.
    for _ in 0..10 {
        let pa = coe_a.geometry.point_at(ta_best);
        let pb = coe_b.geometry.point_at(tb_best);
        let diff = pa - pb;
        let da = coe_a.geometry.tangent_at(ta_best);
        let db = coe_b.geometry.tangent_at(tb_best);
        let f0 = diff.dot(&da);
        let f1 = -diff.dot(&db);
        let j00 = da.dot(&da);
        let j01 = -da.dot(&db);
        let j11 = db.dot(&db);
        let det = j00 * j11 - j01 * j01;
        if det.abs() < f64::EPSILON {
            break;
        }
        let dta = (j11 * f0 - j01 * f1) / det;
        let dtb = (j00 * f1 - j01 * f0) / det;
        let ta_new = (ta_best - dta).clamp(coe_a.t_start, coe_a.t_end);
        let tb_new = (tb_best - dtb).clamp(coe_b.t_start, coe_b.t_end);
        if (ta_new - ta_best).abs() < 1e-12 && (tb_new - tb_best).abs() < 1e-12 {
            ta_best = ta_new;
            tb_best = tb_new;
            break;
        }
        ta_best = ta_new;
        tb_best = tb_new;
    }
    let pa = coe_a.geometry.point_at(ta_best);
    let pb = coe_b.geometry.point_at(tb_best);
    Ok(pa.distance_to(&pb))
}

/// Minimum distance between an edge and a face.
///
/// Samples the edge at 32 points and the face UV grid at 16×16 to seed the
/// best parameters, then refines with alternating Newton iteration between the
/// curve parameter and the surface (u, v) parameters.
fn edge_face(
    brep_e: &BRep,
    e: EdgeHandle,
    brep_f: &BRep,
    f: FaceHandle,
    _tolerance: &Tolerance,
) -> Result<f64, MeasureError> {
    let coe = brep_e
        .edges
        .get(e)
        .ok_or_else(|| MeasureError::InvalidHandle("edge".into()))?
        .curve
        .as_ref()
        .ok_or(MeasureError::MissingCurve)?;
    let face_data = brep_f
        .faces
        .get(f)
        .ok_or_else(|| MeasureError::InvalidHandle("face".into()))?;
    let surface_on_face = face_data
        .surface
        .as_ref()
        .ok_or(MeasureError::MissingSurface)?;

    let positions = brep_f.face_outer_loop_positions(f);
    let ((u0, u1), (v0, v1)) = tight_bounded_parameter_range(&positions, &surface_on_face.geometry);

    const NE: usize = 32;
    const NF: usize = 16;
    let dt = (coe.t_end - coe.t_start) / NE as f64;
    let du = (u1 - u0) / NF as f64;
    let dv = (v1 - v0) / NF as f64;

    let mut min_d = f64::MAX;
    let mut t_best = coe.t_start;
    let mut u_best = u0;
    let mut v_best = v0;

    for i in 0..=NE {
        let t = coe.t_start + dt * i as f64;
        let pe = coe.geometry.point_at(t);
        for ui in 0..=NF {
            for vi in 0..=NF {
                let u = u0 + du * ui as f64;
                let v = v0 + dv * vi as f64;
                let pf = surface_on_face.geometry.point_at(u, v);
                let d = pe.distance_to(&pf);
                if d < min_d {
                    min_d = d;
                    t_best = t;
                    u_best = u;
                    v_best = v;
                }
            }
        }
    }

    // Alternating Newton refinement.
    for _ in 0..5 {
        // Step 1: fix (u, v), refine t — solve (C(t) - S(u,v))·C'(t) = 0.
        let sp = surface_on_face.geometry.point_at(u_best, v_best);
        t_best = refine_curve_param(
            &coe.geometry,
            t_best,
            coe.t_start,
            coe.t_end,
            [sp.x, sp.y, sp.z],
        );

        // Step 2: fix t, refine (u, v) — solve Gram 2×2 system.
        let cp = coe.geometry.point_at(t_best);
        (u_best, v_best) = refine_surface_params(
            &surface_on_face.geometry,
            u_best,
            v_best,
            [u0, u1, v0, v1],
            [cp.x, cp.y, cp.z],
        );
    }

    let final_cp = coe.geometry.point_at(t_best);
    let final_sp = surface_on_face.geometry.point_at(u_best, v_best);
    Ok(final_cp.distance_to(&final_sp))
}

/// Minimum distance between two faces.
///
/// Samples a 16×16 UV grid on each face to seed the best parameters, then
/// refines with alternating Newton iteration between the two surfaces.
fn face_face(
    brep_a: &BRep,
    fa: FaceHandle,
    brep_b: &BRep,
    fb: FaceHandle,
    _tolerance: &Tolerance,
) -> Result<f64, MeasureError> {
    let face_a = brep_a
        .faces
        .get(fa)
        .ok_or_else(|| MeasureError::InvalidHandle("face".into()))?;
    let surf_a = face_a
        .surface
        .as_ref()
        .ok_or(MeasureError::MissingSurface)?;

    let face_b = brep_b
        .faces
        .get(fb)
        .ok_or_else(|| MeasureError::InvalidHandle("face".into()))?;
    let surf_b = face_b
        .surface
        .as_ref()
        .ok_or(MeasureError::MissingSurface)?;

    let positions_a = brep_a.face_outer_loop_positions(fa);
    let ((au0, au1), (av0, av1)) = tight_bounded_parameter_range(&positions_a, &surf_a.geometry);

    let positions_b = brep_b.face_outer_loop_positions(fb);
    let ((bu0, bu1), (bv0, bv1)) = tight_bounded_parameter_range(&positions_b, &surf_b.geometry);

    const N: usize = 16;
    let dua = (au1 - au0) / N as f64;
    let dva = (av1 - av0) / N as f64;
    let dub = (bu1 - bu0) / N as f64;
    let dvb = (bv1 - bv0) / N as f64;

    let mut min_d = f64::MAX;
    let mut ua_best = au0;
    let mut va_best = av0;
    let mut ub_best = bu0;
    let mut vb_best = bv0;

    for i in 0..=N {
        for j in 0..=N {
            let ua = au0 + dua * i as f64;
            let va = av0 + dva * j as f64;
            let pa = surf_a.geometry.point_at(ua, va);
            for k in 0..=N {
                for l in 0..=N {
                    let ub = bu0 + dub * k as f64;
                    let vb = bv0 + dvb * l as f64;
                    let pb = surf_b.geometry.point_at(ub, vb);
                    let d = pa.distance_to(&pb);
                    if d < min_d {
                        min_d = d;
                        ua_best = ua;
                        va_best = va;
                        ub_best = ub;
                        vb_best = vb;
                    }
                }
            }
        }
    }

    // Alternating Newton refinement.
    for _ in 0..5 {
        // Step 1: fix (ub, vb), project onto surface A.
        let pb = surf_b.geometry.point_at(ub_best, vb_best);
        (ua_best, va_best) = refine_surface_params(
            &surf_a.geometry,
            ua_best,
            va_best,
            [au0, au1, av0, av1],
            [pb.x, pb.y, pb.z],
        );

        // Step 2: fix (ua, va), project onto surface B.
        let pa = surf_a.geometry.point_at(ua_best, va_best);
        (ub_best, vb_best) = refine_surface_params(
            &surf_b.geometry,
            ub_best,
            vb_best,
            [bu0, bu1, bv0, bv1],
            [pa.x, pa.y, pa.z],
        );
    }

    let final_pa = surf_a.geometry.point_at(ua_best, va_best);
    let final_pb = surf_b.geometry.point_at(ub_best, vb_best);
    Ok(final_pa.distance_to(&final_pb))
}

#[cfg(test)]
#[path = "distance_tests.rs"]
mod tests;
