//! Facet-deviation refinement for the structured tensor strips
//! (`boundary_first_ruled_strip`).
//!
//! A tensor cell whose iso-line chords each meet the chord height still
//! carries diagonal facets up to about twice it (the two directions'
//! sags add on a synclastic surface: slit torus fillet bands sat 0.0084
//! off at a 0.005 chord, a rounded washer 0.0067; PK-653). The strip's
//! stations ARE its boundary samples, shared with the neighbouring
//! faces, so the cure is the shared split mechanism: every boundary
//! segment spanning an over-deviation triangle's cell, both rails at
//! its row and both cross chains at its column, is split at its span
//! midpoint, and the next tessellation round meshes both faces of each
//! edge on the refined samples. Splitting every spanning segment keeps
//! the rails paired station for station.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::boundary_first_deviation::normal_gap;
use super::{BoundaryRing, SplitRequest};
use crate::mesh::FaceMesh;

/// Request splits of the ring segments that span any triangle of `fm`
/// whose centroid lies farther than `chord_height` from `surface`
/// (measured along the surface normal). Returns the request count.
pub(super) fn request_deviation_splits(
    fm: &FaceMesh,
    surface: &SurfaceGeometry,
    ring: &BoundaryRing,
    chord_height: f64,
    split_requests: &mut Vec<SplitRequest>,
) -> usize {
    let n = ring.uv.len();
    if n < 4 || chord_height <= 0.0 {
        return 0;
    }
    // Chart metric at the ring's centre, to classify each segment by
    // its dominant (running) direction.
    let (mut u0, mut u1, mut v0, mut v1) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for &(u, v) in &ring.uv {
        u0 = u0.min(u);
        u1 = u1.max(u);
        v0 = v0.min(v);
        v1 = v1.max(v);
    }
    let (um, vm) = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
    let p = surface.point_at(um, vm);
    let eps_u = ((u1 - u0) * 1e-4).max(1e-9);
    let eps_v = ((v1 - v0) * 1e-4).max(1e-9);
    let su = surface.point_at(um + eps_u, vm).distance_to(&p) / eps_u;
    let sv = surface.point_at(um, vm + eps_v).distance_to(&p) / eps_v;
    let mut spanning: Vec<bool> = vec![false; n];
    for t in &fm.triangles {
        let [a, b, c] = t.map(|i| &fm.vertices[i]);
        let centroid = Point3::new(
            (a.position.x + b.position.x + c.position.x) / 3.0,
            (a.position.y + b.position.y + c.position.y) / 3.0,
            (a.position.z + b.position.z + c.position.z) / 3.0,
        );
        let cuv = (
            (a.uv.0 + b.uv.0 + c.uv.0) / 3.0,
            (a.uv.1 + b.uv.1 + c.uv.1) / 3.0,
        );
        if normal_gap(surface, cuv, centroid) <= chord_height {
            continue;
        }
        for (k, flag) in spanning.iter_mut().enumerate() {
            let (pa, pb) = (ring.uv[k], ring.uv[(k + 1) % n]);
            let runs_u = (su * (pb.0 - pa.0)).abs() >= (sv * (pb.1 - pa.1)).abs();
            let (lo, hi, x) = if runs_u {
                (pa.0.min(pb.0), pa.0.max(pb.0), cuv.0)
            } else {
                (pa.1.min(pb.1), pa.1.max(pb.1), cuv.1)
            };
            if x >= lo && x <= hi {
                *flag = true;
            }
        }
    }
    let mut count = 0;
    for (k, &split) in spanning.iter().enumerate() {
        let (edge, f_lo, f_hi) = ring.sources[k];
        if split && f_hi != f_lo {
            split_requests.push((edge, 0.5 * (f_lo + f_hi)));
            count += 1;
        }
    }
    count
}
