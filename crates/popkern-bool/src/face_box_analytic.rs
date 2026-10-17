//! Trim-scale boxes for sphere and torus faces (PK-875).
//!
//! The whole-surface grid boxed a face by its untrimmed surface: the
//! six-axis arm's j1 housing dome, a 38.5 mm cap of a 157.7 mm sphere,
//! boxed a 315 mm cube and paired with the base's column wall and
//! connector boss 66-109 mm away, whose records on the untrimmed sphere
//! the clip then had to discard.
//!
//! A coordinate's extreme over a face lies on its boundary or at an
//! interior critical point of the coordinate on the surface. The
//! boundary curves' conservative extents are already in the box; the
//! critical points of a sphere (centre ± r along each axis) and of a
//! torus (four per axis: the tube circles facing the axis direction)
//! join it when their chart point lies in the face's padded trim
//! window. The window contains the face, so the box contains it.
//! A face with any boundary use lacking a pcurve keeps the grid: its
//! window would fall back to vertex positions and could under-cover.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

/// Pad of the sampled trim window, as a fraction of its span (as the
/// NURBS sub-patch hull pads it).
const WINDOW_PAD: f64 = 0.05;

/// The critical points of a sphere or torus face that lie in its trim
/// window, or `None` when the face keeps the whole-surface grid.
pub(crate) fn trim_scale_points(
    brep: &BRep,
    face: FaceHandle,
    surface: &SurfaceGeometry,
) -> Option<Vec<Point3>> {
    let candidates: Vec<(Point3, bool)> = match surface {
        SurfaceGeometry::Sphere {
            center,
            radius,
            axis,
        } => axes()
            .into_iter()
            .flat_map(|e| [e, -e])
            .map(|e| {
                // A pole's u is arbitrary: its window test is v alone.
                let pole = e.cross(axis).length() < 1e-9;
                (*center + e * *radius, pole)
            })
            .collect(),
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => {
            let a = axis.normalized()?;
            let mut out = Vec::new();
            for e in axes() {
                let perp = e - a * e.dot(&a);
                // Along the axis the extremes are whole circles: any u.
                let (radial, free) = match perp.normalized() {
                    Some(w) if perp.length() > 1e-9 => ([w, -w], false),
                    _ => {
                        let any = popkern_core::geometry::analytic_surfaces::orthonormal_basis(a).0;
                        ([any, -any], true)
                    }
                };
                for rho in radial {
                    let tube = *center + rho * *major_radius;
                    let n = (rho * e.dot(&rho) + a * e.dot(&a)).normalized()?;
                    for s in [1.0, -1.0] {
                        out.push((tube + n * (s * *minor_radius), free));
                    }
                }
            }
            out
        }
        _ => return None,
    };
    if !every_use_has_a_pcurve(brep, face) {
        return None;
    }
    let ((u0, u1), (v0, v1)) = brep.face_trim_uv_window(face)?;
    let pad = |lo: f64, hi: f64| {
        let p = WINDOW_PAD * (hi - lo);
        (lo - p, hi + p)
    };
    let ((u0, u1), (v0, v1)) = (pad(u0, u1), pad(v0, v1));
    let (pu, pv) = surface.periodicity();
    let within = |x: f64, lo: f64, hi: f64, period: Option<f64>| match period {
        Some(p) if p > 0.0 => (-3..=3).any(|k| {
            let y = x + f64::from(k) * p;
            y >= lo && y <= hi
        }),
        _ => x >= lo && x <= hi,
    };
    let mut points = Vec::new();
    for (p, free_u) in candidates {
        let Some((u, v)) = surface.closest_point(&p) else {
            // No chart point: keep it (over-covering is safe).
            points.push(p);
            continue;
        };
        if within(v, v0, v1, pv) && (free_u || within(u, u0, u1, pu)) {
            points.push(p);
        }
    }
    Some(points)
}

fn axes() -> [Vec3; 3] {
    [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ]
}

/// Whether every half-edge of the face's loops carries a pcurve (the
/// trim window reads them).
fn every_use_has_a_pcurve(brep: &BRep, face: FaceHandle) -> bool {
    let Some(f) = brep.faces.get(face) else {
        return false;
    };
    for lh in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let Some(start) = brep.loops.get(lh).map(|l| l.half_edge) else {
            return false;
        };
        let mut he = start;
        for _ in 0..=brep.half_edges.len() {
            let Some(h) = brep.half_edges.get(he) else {
                return false;
            };
            if h.pcurve.is_none() {
                return false;
            }
            he = h.next;
            if he == start {
                break;
            }
        }
    }
    true
}

#[cfg(test)]
#[path = "face_box_analytic_tests.rs"]
mod tests;
