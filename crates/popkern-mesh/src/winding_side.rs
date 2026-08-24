//! Component-side decision for global mesh winding coherence.
//!
//! Parity votes orient the faces of one connected component
//! consistently, but the component's GLOBAL side needs a geometric
//! fact: an outer boundary encloses positive volume, while a VOID
//! boundary (the inner skin of a hollow body) must enclose negative
//! volume — its facing points into the void, away from the material.
//! The former rule ("every component positive") meshed every hollow
//! body's void skin inside-out: a thickened cube measured
//! outer+inner (2.24) instead of the wall volume (1.216), and any
//! STL export faced the void skin into the wall material.
//!
//! Enclosure is decided by ray-crossing parity against the OTHER
//! components' triangles — a pure position test, independent of the
//! (not-yet-fixed) triangle facings, so it can run before sides are
//! chosen. Odd crossings = the component sits inside another closed
//! component = void boundary.

use popkern_core::geometry::{Point3, Vec3};

use crate::mesh::FaceMesh;

/// Deterministic, axis-skew ray direction (no normalization needed for
/// parity counting; skew components dodge axis-aligned degeneracies).
const RAY_DIR: Vec3 = Vec3 {
    x: 0.577_350_1,
    y: 0.577_351_2,
    z: 0.577_352_3,
};

/// Whether `component` (indices into `face_meshes`) lies inside another
/// component of the same body, by ray-crossing parity from one of its
/// surface points against all triangles of faces OUTSIDE the component.
pub(crate) fn component_is_enclosed(component: &[usize], face_meshes: &[FaceMesh]) -> bool {
    let in_component = |fi: usize| component.contains(&fi);
    // Sample point: first triangle centroid of the component.
    let Some(origin) = component.iter().find_map(|&fi| {
        let fm = &face_meshes[fi];
        fm.triangles.first().map(|t| {
            let (a, b, c) = (
                fm.vertices[t[0]].position,
                fm.vertices[t[1]].position,
                fm.vertices[t[2]].position,
            );
            Point3::new(
                (a.x + b.x + c.x) / 3.0,
                (a.y + b.y + c.y) / 3.0,
                (a.z + b.z + c.z) / 3.0,
            )
        })
    }) else {
        return false;
    };

    let mut crossings = 0_usize;
    for (fi, fm) in face_meshes.iter().enumerate() {
        if in_component(fi) {
            continue;
        }
        for t in &fm.triangles {
            let (a, b, c) = (
                fm.vertices[t[0]].position,
                fm.vertices[t[1]].position,
                fm.vertices[t[2]].position,
            );
            if ray_hits_triangle(&origin, &RAY_DIR, &a, &b, &c) {
                crossings += 1;
            }
        }
    }
    crossings % 2 == 1
}

/// Möller–Trumbore ray/triangle test for `origin + t·dir`, `t > 0`.
fn ray_hits_triangle(origin: &Point3, dir: &Vec3, a: &Point3, b: &Point3, c: &Point3) -> bool {
    let e1 = *b - *a;
    let e2 = *c - *a;
    let p = dir.cross(&e2);
    let det = e1.dot(&p);
    if det.abs() < 1e-15 {
        return false; // ray parallel to triangle plane
    }
    let inv_det = 1.0 / det;
    let s = *origin - *a;
    let u = s.dot(&p) * inv_det;
    if !(0.0..=1.0).contains(&u) {
        return false;
    }
    let q = s.cross(&e1);
    let v = dir.dot(&q) * inv_det;
    if v < 0.0 || u + v > 1.0 {
        return false;
    }
    let t = e2.dot(&q) * inv_det;
    t > 1e-12
}
