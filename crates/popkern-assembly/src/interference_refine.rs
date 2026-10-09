//! Refinement of closest mesh points onto the true faces for the
//! interference check (split from `interference.rs` under the
//! 500-line budget).

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::topology::{BRep, FaceHandle};

use crate::mesh_bvh::{LocalMesh, WorldMesh};

/// Alternating projection of the closest mesh points onto the true
/// faces; `Some` only when both end points lie on their trimmed faces.
#[allow(clippy::too_many_arguments)]
pub(crate) fn refine(
    brep_a: &BRep,
    fa: FaceHandle,
    mesh_a: &WorldMesh,
    la: &LocalMesh,
    brep_b: &BRep,
    fb: FaceHandle,
    mesh_b: &WorldMesh,
    lb: &LocalMesh,
    here: Point3,
    there: Point3,
) -> Option<(Point3, Point3)> {
    let sa = &brep_a.faces.get(fa)?.surface.as_ref()?.geometry;
    let sb = &brep_b.faces.get(fb)?.surface.as_ref()?.geometry;
    let size = mesh_a
        .bounds
        .diagonal()
        .max(mesh_b.bounds.diagonal())
        .max(1e-12);
    let project = |s: &popkern_core::geometry::surface::SurfaceGeometry, p: &Point3| {
        let (u, v) = project_point_to_uv(s, p)?;
        Some(s.point_at(u, v))
    };
    let mut p = project(sa, &here)?;
    let mut q = project(sb, &there)?;
    for _ in 0..64 {
        let q_next = project(sb, &p)?;
        let p_next = project(sa, &q_next)?;
        let moved = (p_next - p).length() + (q_next - q).length();
        p = p_next;
        q = q_next;
        if moved <= 1e-14 * size {
            break;
        }
    }
    // On the face means inside its EXACT trim. A flat face's triangles
    // carry no surface deviation yet overshoot a concave curved
    // boundary by the chord sagitta — 0.08 mm past a cycloid disc's
    // flank, onto a ring pin 0.02 mm away, which then read as contact
    // at distance 0 (PK-747). The triangles answer only when the trim
    // test cannot.
    let on_trim = |brep: &BRep,
                   s: &popkern_core::geometry::surface::SurfaceGeometry,
                   face: FaceHandle,
                   x: &Point3| {
        let (u, v) = project_point_to_uv(s, x)?;
        let polygon = brep.face_uv_trim_containment(face, s, u, v).map(|c| {
            matches!(
                c,
                popkern_core::validation::FaceUvContainment::Interior
                    | popkern_core::validation::FaceUvContainment::Boundary
            )
        })?;
        // The polygon samples curved boundaries coarsely: near one, the
        // boundary curve itself decides the side.
        Some(
            polygon
                && crate::trim_side::material_side_near_boundary(brep, face, s, x).unwrap_or(true),
        )
    };
    let on_face = |mesh: &WorldMesh, local: &LocalMesh, face: FaceHandle, x: &Point3| {
        let mut nearest = f64::INFINITY;
        let mut eps: f64 = 0.0;
        for (k, tri) in mesh.tris.iter().enumerate() {
            let i = mesh.local[k];
            if local.faces[i] != face {
                continue;
            }
            eps = eps.max(local.eps[i]);
            nearest = nearest.min((crate::tri_geometry::closest_on_triangle(x, tri) - *x).length());
        }
        nearest <= eps + 1e-9 * size
    };
    let inside_a = on_trim(brep_a, sa, fa, &p).unwrap_or_else(|| on_face(mesh_a, la, fa, &p));
    let inside_b = on_trim(brep_b, sb, fb, &q).unwrap_or_else(|| on_face(mesh_b, lb, fb, &q));
    (inside_a && inside_b).then_some((p, q))
}

/// The faces whose triangles hold mesh point `x` (within `1e-9` of
/// `size`), `own` first: more than one when `x` lies on a mesh edge or
/// vertex between faces.
pub(crate) fn faces_at(
    mesh: &WorldMesh,
    local: &LocalMesh,
    x: &Point3,
    own: FaceHandle,
    size: f64,
) -> Vec<FaceHandle> {
    let mut faces = vec![own];
    for (k, tri) in mesh.tris.iter().enumerate() {
        let face = local.faces[mesh.local[k]];
        if !faces.contains(&face)
            && (crate::tri_geometry::closest_on_triangle(x, tri) - *x).length() <= 1e-9 * size
        {
            faces.push(face);
        }
    }
    faces
}

/// The diagonal of `brep`'s vertex bounds (the part size the
/// tessellation chord is a fraction of).
pub(crate) fn bounds_diagonal(brep: &BRep) -> f64 {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for v in brep.vertices.values() {
        for (k, x) in [v.position.x, v.position.y, v.position.z]
            .into_iter()
            .enumerate()
        {
            lo[k] = lo[k].min(x);
            hi[k] = hi[k].max(x);
        }
    }
    (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt()
}
