//! Delaunay flips in 3D (E4 step 1 of
//! docs/architecture/display-tessellation-in-face-folds.md).
//!
//! Boundary-first faces are triangulated by a Delaunay of the chart
//! scaled by ONE metric taken at the chart centre. Where the surface's
//! metric changes across the chart (a coil's helix, a twisted blade) or
//! is sheared, triangles that are well-shaped in that scaled chart are
//! obtuse in 3D, and an obtuse facet's normal drifts from its vertex
//! normals: the leaning interior facets of the display contract sat at
//! the 12 deg spread limit with a ~6 deg smallest angle (coil) and
//! 3 deg (propeller), against an 18 deg median. This pass restores the
//! Delaunay property in 3D: an interior edge flips when the two angles
//! opposite it sum past pi (Lawson). Guards keep everything the earlier
//! passes established: both new triangles keep the old pair's chart
//! winding with area and face along their vertex normals, their
//! vertex-normal spread and facet tilt stay within the bound (or the old
//! pair's), and the new diagonal sags no more than the chord height (or
//! the old diagonal). Chart-singular vertices are left alone.

use std::collections::HashMap;

use popkern_core::geometry::surface::SurfaceGeometry;

use crate::crease_flip::{
    EdgeKey, angle_between, apex, chart_area, chord_sag, facet_normal, key, runs, spread,
};
use crate::mesh::FaceMesh;
use crate::vertex_normal::singular_without_unique_normal;

/// Flip budget per triangle: Lawson flipping terminates on its own; the
/// budget only bounds pathological guard cycles.
const BUDGET_PER_TRIANGLE: usize = 4;

/// Largest angle between a triangle's facet normal and its vertex
/// normals.
fn tilt(mesh: &FaceMesh, t: [usize; 3]) -> f64 {
    let n = facet_normal(mesh, t);
    t.iter()
        .map(|&v| angle_between(n, mesh.vertices[v].normal))
        .fold(0.0, f64::max)
}

/// The 3D angle at `c` in triangle `(a, b, c)`.
fn angle_at(mesh: &FaceMesh, a: usize, b: usize, c: usize) -> f64 {
    let p = |i: usize| mesh.vertices[i].position;
    angle_between(p(a) - p(c), p(b) - p(c))
}

/// Flip interior edges of `mesh` to the 3D Delaunay triangulation of
/// its vertices, under the guards in the module docs. Returns the
/// number of flips.
pub(crate) fn delaunay_flips_3d(
    mesh: &mut FaceMesh,
    surface: &SurfaceGeometry,
    chord_height: f64,
    spread_bound: f64,
) -> usize {
    if mesh.triangles.len() < 2 || matches!(surface, SurfaceGeometry::Plane { .. }) {
        return 0;
    }
    let singular: Vec<bool> = mesh
        .vertices
        .iter()
        .map(|v| singular_without_unique_normal(surface, v.uv.0, v.uv.1))
        .collect();
    let mut edge_tris: HashMap<EdgeKey, Vec<usize>> = HashMap::new();
    for (ti, t) in mesh.triangles.iter().enumerate() {
        for k in 0..3 {
            edge_tris
                .entry(key(t[k], t[(k + 1) % 3]))
                .or_default()
                .push(ti);
        }
    }
    let mut queue: Vec<EdgeKey> = edge_tris
        .iter()
        .filter(|(_, ts)| ts.len() == 2)
        .map(|(&e, _)| e)
        .collect();
    queue.sort_unstable();
    let mut budget = mesh.triangles.len() * BUDGET_PER_TRIANGLE;
    let mut flips = 0usize;
    while let Some(e) = queue.pop() {
        if budget == 0 {
            break;
        }
        let Some(&[t0, t1]) = edge_tris.get(&e).map(Vec::as_slice) else {
            continue;
        };
        let (a, b) = if runs(mesh.triangles[t0], e.0, e.1) {
            (e.0, e.1)
        } else {
            (e.1, e.0)
        };
        let (tri0, tri1) = (mesh.triangles[t0], mesh.triangles[t1]);
        if !runs(tri0, a, b) || !runs(tri1, b, a) {
            continue;
        }
        let (Some(c), Some(d)) = (apex(tri0, a, b), apex(tri1, a, b)) else {
            continue;
        };
        if c == d || edge_tris.contains_key(&key(c, d)) || [a, b, c, d].iter().any(|&v| singular[v])
        {
            continue;
        }
        // Lawson: flip when the opposite angles sum past pi (with a
        // margin so cocircular quads do not flip back and forth).
        if angle_at(mesh, a, b, c) + angle_at(mesh, b, a, d) <= std::f64::consts::PI + 1e-9 {
            continue;
        }
        let (new0, new1) = ([a, d, c], [d, b, c]);
        let sense = chart_area(mesh, tri0).signum();
        let quad_area = chart_area(mesh, tri0).abs() + chart_area(mesh, tri1).abs();
        let chart_ok = [new0, new1].iter().all(|&t| {
            let area = chart_area(mesh, t);
            area.signum() == sense && area.abs() > 1e-9 * quad_area
        });
        let facing_ok = [new0, new1].iter().all(|&t| {
            let n = facet_normal(mesh, t);
            t.iter().all(|&v| n.dot(&mesh.vertices[v].normal) > 0.0)
        });
        let spread_ok = spread(mesh, new0).max(spread(mesh, new1))
            <= spread_bound.max(spread(mesh, tri0).max(spread(mesh, tri1)));
        // Never trade shape for lean: the display contract is tilt.
        let tilt_ok = tilt(mesh, new0).max(tilt(mesh, new1))
            <= spread_bound.max(tilt(mesh, tri0).max(tilt(mesh, tri1)));
        if !chart_ok || !facing_ok || !spread_ok || !tilt_ok {
            continue;
        }
        let old_sag = chord_sag(mesh, surface, a, b).unwrap_or(0.0);
        if chord_sag(mesh, surface, c, d).is_none_or(|s| s > chord_height.max(old_sag)) {
            continue;
        }
        mesh.triangles[t0] = new0;
        mesh.triangles[t1] = new1;
        let moved = |edge_tris: &mut HashMap<EdgeKey, Vec<usize>>, e: EdgeKey, from, to| {
            if let Some(list) = edge_tris.get_mut(&e) {
                for x in list.iter_mut() {
                    if *x == from {
                        *x = to;
                    }
                }
            }
        };
        edge_tris.remove(&key(a, b));
        edge_tris.insert(key(c, d), vec![t0, t1]);
        moved(&mut edge_tris, key(a, d), t1, t0);
        moved(&mut edge_tris, key(b, c), t0, t1);
        queue.extend([key(a, c), key(b, c), key(a, d), key(b, d)]);
        budget -= 1;
        flips += 1;
    }
    flips
}
