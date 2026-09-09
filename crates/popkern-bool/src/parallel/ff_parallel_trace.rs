//! Trace dumps for the parallel face-face dispatch (split from
//! `ff_parallel.rs` under the 500-line file cap).

use popkern_core::geometry::CurveEval;
use popkern_core::topology::FaceHandle;

pub(super) fn trace_pair_surfaces(
    face_a: FaceHandle,
    face_b: FaceHandle,
    surf_a: &popkern_core::geometry::surface::SurfaceGeometry,
    surf_b: &popkern_core::geometry::surface::SurfaceGeometry,
    ((u1_min, u1_max), (v1_min, v1_max)): ((f64, f64), (f64, f64)),
    ((u2_min, u2_max), (v2_min, v2_max)): ((f64, f64), (f64, f64)),
) {
    let kind = |s: &popkern_core::geometry::surface::SurfaceGeometry| match s {
        popkern_core::geometry::surface::SurfaceGeometry::Plane { .. } => "Plane",
        popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. } => "Cylinder",
        popkern_core::geometry::surface::SurfaceGeometry::Sphere { .. } => "Sphere",
        popkern_core::geometry::surface::SurfaceGeometry::Cone { .. } => "Cone",
        popkern_core::geometry::surface::SurfaceGeometry::Torus { .. } => "Torus",
        _ => "Nurbs/other",
    };
    eprintln!(
        "  [pair-surf] {face_a:?}({}) u[{u1_min:.1},{u1_max:.1}] v[{v1_min:.1},{v1_max:.1}] x {face_b:?}({}) u[{u2_min:.1},{u2_max:.1}] v[{v2_min:.1},{v2_max:.1}]",
        kind(surf_a),
        kind(surf_b)
    );
    if let popkern_core::geometry::surface::SurfaceGeometry::Plane { u_axis, v_axis, .. } = surf_a {
        eprintln!(
            "    [pair-axes] A ua=({:.2},{:.2},{:.2}) va=({:.2},{:.2},{:.2})",
            u_axis.x, u_axis.y, u_axis.z, v_axis.x, v_axis.y, v_axis.z
        );
    }
}

pub(super) fn trace_clip_input(
    face_a: FaceHandle,
    face_b: FaceHandle,
    surf_a: &popkern_core::geometry::surface::SurfaceGeometry,
    surf_b: &popkern_core::geometry::surface::SurfaceGeometry,
    curve: &popkern_core::geometry::curve_on_surface::CurveOnSurface,
    kept: usize,
) {
    let (t0, t1) = curve.curve_3d.parameter_range();
    let p0 = curve.curve_3d.point_at(t0);
    let p1 = curve.curve_3d.point_at(t1);
    let plane_dist = |s: &popkern_core::geometry::surface::SurfaceGeometry,
                      p: &popkern_core::geometry::Point3|
     -> f64 {
        match s {
            popkern_core::geometry::surface::SurfaceGeometry::Plane { origin, normal, .. } => {
                ((p.x - origin.x) * normal.x
                    + (p.y - origin.y) * normal.y
                    + (p.z - origin.z) * normal.z)
                    .abs()
            }
            _ => f64::NAN,
        }
    };
    eprintln!(
        "  [clip] A{face_a:?} x B{face_b:?}: t[{t0:.3},{t1:.3}] ({:.9},{:.9},{:.9})->({:.9},{:.9},{:.9}) dA={:.1e} dB={:.1e} kept={}",
        p0.x,
        p0.y,
        p0.z,
        p1.x,
        p1.y,
        p1.z,
        plane_dist(surf_a, &p0),
        plane_dist(surf_b, &p0),
        kept
    );
}

pub(super) fn trace_clip_output(
    face_a: FaceHandle,
    face_b: FaceHandle,
    clipped: &popkern_core::geometry::curve_on_surface::CurveOnSurface,
) {
    let (t0, t1) = clipped.curve_3d.parameter_range();
    let p0 = clipped.curve_3d.point_at(t0);
    let p1 = clipped.curve_3d.point_at(t1);
    eprintln!(
        "  [clip-out] A{face_a:?} x B{face_b:?}: t[{t0:.3},{t1:.3}] ({:.9},{:.9},{:.9})->({:.9},{:.9},{:.9}) closed={}",
        p0.x, p0.y, p0.z, p1.x, p1.y, p1.z, clipped.is_closed
    );
}
