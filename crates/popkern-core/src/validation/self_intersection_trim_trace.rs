use crate::geometry::Point3;
use crate::geometry::surface::SurfaceGeometry;
use crate::topology::FaceHandle;

use super::{TrimContainment, UvPoint};

pub(super) struct TrimOverlapTraceSample<'a> {
    pub(super) face: FaceHandle,
    pub(super) surface: &'a SurfaceGeometry,
    pub(super) trim_contains: Option<TrimContainment>,
    pub(super) uv: UvPoint,
}

pub(super) fn trace_trim_overlap_sample(
    a: TrimOverlapTraceSample<'_>,
    b: TrimOverlapTraceSample<'_>,
) {
    if std::env::var_os("POPKERN_SELF_INTERSECTION_TRACE").is_none() {
        return;
    }
    eprintln!(
        "[self-intersection trim] {:?}/{} uv=({:.6},{:.6}) inside={:?} \
         {:?}/{} uv=({:.6},{:.6}) inside={:?}",
        a.face,
        surface_kind(a.surface),
        a.uv.u,
        a.uv.v,
        a.trim_contains,
        b.face,
        surface_kind(b.surface),
        b.uv.u,
        b.uv.v,
        b.trim_contains,
    );
}

pub(super) fn trace_trim_overlap_point(
    face_a: FaceHandle,
    surface_a: &SurfaceGeometry,
    trim_a_contains: Option<TrimContainment>,
    face_b: FaceHandle,
    surface_b: &SurfaceGeometry,
    trim_b_contains: Option<TrimContainment>,
    point: &Point3,
) {
    if std::env::var_os("POPKERN_SELF_INTERSECTION_TRACE").is_none() {
        return;
    }
    eprintln!(
        "[self-intersection trim] {face_a:?}/{} inside={trim_a_contains:?} \
         {face_b:?}/{} inside={trim_b_contains:?} point=({:.6},{:.6},{:.6})",
        surface_kind(surface_a),
        surface_kind(surface_b),
        point.x,
        point.y,
        point.z,
    );
}

fn surface_kind(surface: &SurfaceGeometry) -> &'static str {
    match surface {
        SurfaceGeometry::Plane { .. } => "Plane",
        SurfaceGeometry::Cylinder { .. } => "Cylinder",
        SurfaceGeometry::Cone { .. } => "Cone",
        SurfaceGeometry::Sphere { .. } => "Sphere",
        SurfaceGeometry::Torus { .. } => "Torus",
        SurfaceGeometry::Nurbs(_) => "Nurbs",
        SurfaceGeometry::Gregory(_) => "Gregory",
        SurfaceGeometry::PolyhedralSpline(_) => "PolyhedralSpline",
    }
}
