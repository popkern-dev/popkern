use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, FaceHandle, HalfEdge};

pub(super) fn dump_pcurve_trajectory(
    brep: &BRep,
    sof: &popkern_core::geometry::surface::SurfaceOnFace,
    he: &HalfEdge,
    twin_face: FaceHandle,
    step: usize,
    pa: &popkern_core::geometry::pcurve::Pcurve,
) {
    let edge = brep.edges.get(he.edge);
    let twin_kind = brep
        .faces
        .get(twin_face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| match &s.geometry {
            SurfaceGeometry::Plane { .. } => "Plane",
            SurfaceGeometry::Cylinder { .. } => "Cyl",
            _ => "Other",
        })
        .unwrap_or("none");
    let (pa_lo, pa_hi) = pa.parameter_range();
    eprintln!(
        "    step{step} edge={:?} twin={twin_kind} pa_range=[{pa_lo:.3},{pa_hi:.3}]",
        he.edge
    );
    if let Some(coe) = edge.and_then(|e| e.curve.as_ref()) {
        eprintln!("      curve3d_range=[{:.3},{:.3}]", coe.t_start, coe.t_end);
    }
    for s in 0..=4 {
        let f = s as f64 / 4.0;
        let t = pa_lo + f * (pa_hi - pa_lo);
        let uv = pa.point_at(t);
        let p_surf = sof.geometry.point_at(uv.x, uv.y);
        let p_curve = edge.and_then(|e| e.curve.as_ref()).map(|c| {
            let tt = c.t_start + f * (c.t_end - c.t_start);
            c.geometry.point_at(tt)
        });
        let dev = p_curve.map(|p| p.distance_to(&p_surf));
        eprintln!(
            "      s={f:.2}: pa=({:.3},{:.3}) surf=({:.3},{:.3},{:.3}) curve3d={:?} dev={dev:?}",
            uv.x,
            uv.y,
            p_surf.x,
            p_surf.y,
            p_surf.z,
            p_curve.map(|p| (p.x, p.y, p.z)),
        );
    }
}
