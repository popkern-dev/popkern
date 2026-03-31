//! Minimum distance between entity pairs.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::tight_bounded_parameter_range;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use crate::MeasureError;

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
/// Samples the edge curve at 64 uniformly-spaced parameter values and returns
/// the smallest distance found.
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
    let curve_on_edge = brep_e
        .edges
        .get(e)
        .ok_or_else(|| MeasureError::InvalidHandle("edge".into()))?
        .curve
        .as_ref()
        .ok_or(MeasureError::MissingCurve)?;

    let t0 = curve_on_edge.t_start;
    let t1 = curve_on_edge.t_end;
    const N: usize = 64;
    let mut min_d = f64::MAX;
    for i in 0..=N {
        let t = t0 + (t1 - t0) * (i as f64 / N as f64);
        let q = curve_on_edge.geometry.point_at(t);
        let d = p.distance_to(&q);
        if d < min_d {
            min_d = d;
        }
    }
    Ok(min_d)
}

/// Minimum distance from a vertex to a face.
///
/// Samples the outer boundary vertices plus a 16×16 UV grid over the face's
/// bounded parameter domain.
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
    let surface_on_face = face_data
        .surface
        .as_ref()
        .ok_or(MeasureError::MissingSurface)?;

    let positions = brep_f.face_outer_loop_positions(f);
    let ((u0, u1), (v0, v1)) = tight_bounded_parameter_range(&positions, &surface_on_face.geometry);

    let mut min_d = f64::MAX;

    // Sample boundary positions directly.
    for q in &positions {
        let d = p.distance_to(q);
        if d < min_d {
            min_d = d;
        }
    }

    // Sample 16×16 UV grid.
    const N: usize = 16;
    for i in 0..=N {
        for j in 0..=N {
            let u = u0 + (u1 - u0) * (i as f64 / N as f64);
            let vv = v0 + (v1 - v0) * (j as f64 / N as f64);
            let q = surface_on_face.geometry.point_at(u, vv);
            let d = p.distance_to(&q);
            if d < min_d {
                min_d = d;
            }
        }
    }
    Ok(min_d)
}

/// Minimum distance between two edges.
///
/// Samples each curve at 64 points and returns the closest pair distance.
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
    let pts_a: Vec<_> = (0..=N)
        .map(|i| {
            let t = coe_a.t_start + (coe_a.t_end - coe_a.t_start) * (i as f64 / N as f64);
            coe_a.geometry.point_at(t)
        })
        .collect();
    let pts_b: Vec<_> = (0..=N)
        .map(|i| {
            let t = coe_b.t_start + (coe_b.t_end - coe_b.t_start) * (i as f64 / N as f64);
            coe_b.geometry.point_at(t)
        })
        .collect();

    let mut min_d = f64::MAX;
    for pa in &pts_a {
        for pb in &pts_b {
            let d = pa.distance_to(pb);
            if d < min_d {
                min_d = d;
            }
        }
    }
    Ok(min_d)
}

/// Minimum distance between an edge and a face.
///
/// Samples the edge at 32 points and the face UV grid at 16×16, then returns
/// the closest pair distance.
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

    const NE: usize = 32;
    let edge_pts: Vec<_> = (0..=NE)
        .map(|i| {
            let t = coe.t_start + (coe.t_end - coe.t_start) * (i as f64 / NE as f64);
            coe.geometry.point_at(t)
        })
        .collect();

    let positions = brep_f.face_outer_loop_positions(f);
    let ((u0, u1), (v0, v1)) = tight_bounded_parameter_range(&positions, &surface_on_face.geometry);

    const NF: usize = 16;
    let face_pts: Vec<_> = (0..=NF)
        .flat_map(|i| {
            let u = u0 + (u1 - u0) * (i as f64 / NF as f64);
            (0..=NF).map(move |j| {
                let vv = v0 + (v1 - v0) * (j as f64 / NF as f64);
                surface_on_face.geometry.point_at(u, vv)
            })
        })
        .collect();

    let mut min_d = f64::MAX;
    for pe in &edge_pts {
        for pf in &face_pts {
            let d = pe.distance_to(pf);
            if d < min_d {
                min_d = d;
            }
        }
    }
    Ok(min_d)
}

/// Minimum distance between two faces.
///
/// Samples a 16×16 UV grid on each face and returns the closest pair distance.
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
    let pts_a: Vec<_> = (0..=N)
        .flat_map(|i| {
            let u = au0 + (au1 - au0) * (i as f64 / N as f64);
            (0..=N).map(move |j| {
                let v = av0 + (av1 - av0) * (j as f64 / N as f64);
                surf_a.geometry.point_at(u, v)
            })
        })
        .collect();
    let pts_b: Vec<_> = (0..=N)
        .flat_map(|i| {
            let u = bu0 + (bu1 - bu0) * (i as f64 / N as f64);
            (0..=N).map(move |j| {
                let v = bv0 + (bv1 - bv0) * (j as f64 / N as f64);
                surf_b.geometry.point_at(u, v)
            })
        })
        .collect();

    let mut min_d = f64::MAX;
    for pa in &pts_a {
        for pb in &pts_b {
            let d = pa.distance_to(pb);
            if d < min_d {
                min_d = d;
            }
        }
    }
    Ok(min_d)
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn make_cube() -> BRep {
        BRep::unit_cube(Tolerance::default()).expect("cube")
    }

    #[test]
    fn same_vertex_distance_is_zero() {
        let cube = make_cube();
        let v = cube.vertices.keys().next().unwrap();
        let d = min_distance(
            &cube,
            TopoHandle::Vertex(v),
            &cube,
            TopoHandle::Vertex(v),
            &Tolerance::default(),
        )
        .unwrap();
        assert!(d < 1e-10, "self-distance should be 0, got {d}");
    }

    #[test]
    fn vertex_to_vertex_distance() {
        let cube = make_cube();
        let verts: Vec<_> = cube.vertices.keys().collect();
        let d = min_distance(
            &cube,
            TopoHandle::Vertex(verts[0]),
            &cube,
            TopoHandle::Vertex(verts[1]),
            &Tolerance::default(),
        )
        .unwrap();
        assert!(d > 0.0 && d <= 1.74, "distance = {d}");
    }

    #[test]
    fn vertex_to_edge_distance() {
        let cube = make_cube();
        let v = cube.vertices.keys().next().unwrap();
        let e = cube.edges.keys().next().unwrap();
        let d = min_distance(
            &cube,
            TopoHandle::Vertex(v),
            &cube,
            TopoHandle::Edge(e),
            &Tolerance::default(),
        )
        .unwrap();
        assert!(d >= 0.0 && d <= 1.74, "distance = {d}");
    }

    #[test]
    fn edge_to_vertex_symmetry() {
        let cube = make_cube();
        let v = cube.vertices.keys().next().unwrap();
        let e = cube.edges.keys().next().unwrap();
        let tol = Tolerance::default();
        let d_ve = min_distance(
            &cube,
            TopoHandle::Vertex(v),
            &cube,
            TopoHandle::Edge(e),
            &tol,
        )
        .unwrap();
        let d_ev = min_distance(
            &cube,
            TopoHandle::Edge(e),
            &cube,
            TopoHandle::Vertex(v),
            &tol,
        )
        .unwrap();
        assert!(
            (d_ve - d_ev).abs() < 1e-10,
            "symmetry violated: {d_ve} vs {d_ev}"
        );
    }

    #[test]
    fn vertex_to_face_distance() {
        let cube = make_cube();
        let v = cube.vertices.keys().next().unwrap();
        let f = cube.faces.keys().next().unwrap();
        let tol = Tolerance::default();
        let d = min_distance(
            &cube,
            TopoHandle::Vertex(v),
            &cube,
            TopoHandle::Face(f),
            &tol,
        )
        .unwrap();
        assert!(d >= 0.0 && d <= 1.74, "distance = {d}");
    }

    #[test]
    fn edge_to_edge_distance() {
        let cube = make_cube();
        let edges: Vec<_> = cube.edges.keys().collect();
        let tol = Tolerance::default();
        let d = min_distance(
            &cube,
            TopoHandle::Edge(edges[0]),
            &cube,
            TopoHandle::Edge(edges[0]),
            &tol,
        )
        .unwrap();
        assert!(d < 1e-10, "self edge distance should be ~0, got {d}");
    }

    #[test]
    fn face_to_face_distance() {
        let cube = make_cube();
        let faces: Vec<_> = cube.faces.keys().collect();
        let tol = Tolerance::default();
        let d = min_distance(
            &cube,
            TopoHandle::Face(faces[0]),
            &cube,
            TopoHandle::Face(faces[0]),
            &tol,
        )
        .unwrap();
        assert!(d < 1e-10, "self face distance should be ~0, got {d}");
    }

    #[test]
    fn invalid_vertex_handle_returns_error() {
        let cube = make_cube();
        let bad_v = VertexHandle::default();
        let v = cube.vertices.keys().next().unwrap();
        let result = min_distance(
            &cube,
            TopoHandle::Vertex(bad_v),
            &cube,
            TopoHandle::Vertex(v),
            &Tolerance::default(),
        );
        assert!(result.is_err());
    }
}
