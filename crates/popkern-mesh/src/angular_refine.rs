//! Conforming angular refinement of a boundary-first face mesh — the
//! enforcement arm of [`crate::TessellationOptions::max_facet_angle`].
//!
//! Every boundary-first builder (Delaunay lattice, ruled strip,
//! revolution, Gregory patches) sizes its interior by its own local
//! rule, and none of them can promise a per-triangle normal bound:
//! the Delaunay clearance band and protecting disks evacuate lattice
//! points next to the boundary, so the triangles that bridge the
//! boundary to the first surviving row span one to several cells
//! (measured before this pass: torus rounds 20–43°, setback fillet
//! strips 35° as boundary-only rail-to-rail chords, sphere caps 25°,
//! against a 12° bound).
//!
//! The pass is path-independent: while a triangle's vertex normals
//! spread past the bound, bisect its longest INTERIOR edge at the
//! chart midpoint (evaluated on the surface) and split both incident
//! triangles. Boundary edges are never split — they carry the
//! canonical edge samples both neighbouring faces share (already
//! angle-bounded by the edge sampler), so watertightness holds by
//! construction. Guards: a chart-singular vertex (pole, collapsed
//! edge — its normal direction is undefined) does not vote on spread;
//! a midpoint whose surface image leaves the edge's chord
//! neighbourhood (a fold across a chart singularity) marks the edge
//! unsplittable instead of inserting it; a triangle budget bounds the
//! work. Triangles left over the bound after the pass are boundary-
//! limited or guard-limited, never silently dropped.

use std::collections::HashMap;

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};

use crate::mesh::{FaceMesh, MeshVertex};
use crate::vertex_normal::{chart_singular, limit_normal, singular_without_unique_normal};

/// Growth budget: the refined mesh may hold at most this many times
/// the input triangle count (plus a floor for tiny inputs).
const GROWTH_FACTOR: usize = 16;
const GROWTH_FLOOR: usize = 4096;

/// A chart coordinate whose 3D travel across an edge's parameter span
/// is below this fraction of the edge length is collapsed at that
/// endpoint (a pole): the coordinate carries no position there.
const POLE_TRAVEL: f64 = 1.0e-3;

fn angle_between(a: Vec3, b: Vec3) -> f64 {
    let (la, lb) = (a.length(), b.length());
    if la < 1e-300 || lb < 1e-300 {
        return 0.0;
    }
    (a.dot(&b) / (la * lb)).clamp(-1.0, 1.0).acos()
}

/// Refine `mesh` in place until every triangle's vertex-normal spread
/// is within `max_angle` radians, or no admissible interior split
/// remains, or the growth budget is spent.
pub(crate) fn refine_face_mesh(mesh: &mut FaceMesh, surface: &SurfaceGeometry, max_angle: f64) {
    if mesh.triangles.is_empty() || matches!(surface, SurfaceGeometry::Plane { .. }) {
        return;
    }
    // The stored normals already carry the face's facing; recover the
    // sign once so inserted vertices match it.
    let sign = mesh
        .vertices
        .iter()
        .find_map(|v| {
            let n = surface.normal_at(v.uv.0, v.uv.1);
            let d = n.dot(&v.normal);
            (d.abs() > 1e-6).then_some(if d > 0.0 { 1.0 } else { -1.0 })
        })
        .unwrap_or(1.0);
    let mut singular: Vec<bool> = mesh
        .vertices
        .iter()
        .map(|v| singular_without_unique_normal(surface, v.uv.0, v.uv.1))
        .collect();
    let key = |a: usize, b: usize| (a.min(b), a.max(b));
    let mut edge_tris: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
    for (ti, t) in mesh.triangles.iter().enumerate() {
        for k in 0..3 {
            edge_tris
                .entry(key(t[k], t[(k + 1) % 3]))
                .or_default()
                .push(ti);
        }
    }
    let budget = (mesh.triangles.len() * GROWTH_FACTOR).max(GROWTH_FLOOR);
    let mut unsplittable: std::collections::HashSet<(usize, usize)> = Default::default();
    // Longest admissible interior edge of a triangle (3D length).
    // Edge measure for longest-edge selection: normal angle first (the
    // GAUSS-MAP metric the bound lives in), 3D length to break ties.
    // Plain 3D length refines the wrong direction on anisotropic
    // surfaces: a small-minor-radius torus round's violating edges run
    // across the short minor arc, and longest-edge bisection spent the
    // whole budget halving long azimuthal edges whose normals barely
    // differ (measured: revolved_knob rounds budget-capped at 16x with
    // 43° triangles left). Along rulings the angle is ~0, so those
    // edges split only when LEPP forces them — anisotropic by design.
    let measure = |mesh: &FaceMesh, singular: &[bool], a: usize, b: usize| -> (f64, f64) {
        let angle = if singular[a] || singular[b] {
            0.0
        } else {
            angle_between(mesh.vertices[a].normal, mesh.vertices[b].normal)
        };
        (
            angle,
            mesh.vertices[a]
                .position
                .distance_to(&mesh.vertices[b].position),
        )
    };
    let exceeds = |x: (f64, f64), y: (f64, f64)| -> bool {
        let tie = 1e-9 * x.0.max(y.0) + 1e-12;
        if (x.0 - y.0).abs() > tie {
            x.0 > y.0
        } else {
            x.1 > y.1 * (1.0 + 1e-12)
        }
    };
    let longest = |mesh: &FaceMesh,
                   singular: &[bool],
                   edge_tris: &HashMap<(usize, usize), Vec<usize>>,
                   unsplittable: &std::collections::HashSet<(usize, usize)>,
                   t: [usize; 3]|
     -> Option<((usize, usize), (f64, f64))> {
        let mut best: Option<((usize, usize), (f64, f64))> = None;
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            let e = key(a, b);
            if unsplittable.contains(&e) || edge_tris.get(&e).is_none_or(|v| v.len() != 2) {
                continue;
            }
            let m = measure(mesh, singular, a, b);
            if best.is_none_or(|(_, l)| exceeds(m, l)) {
                best = Some((e, m));
            }
        }
        best
    };
    // Rivara longest-edge propagation (LEPP): a split edge must be
    // the longest edge of BOTH its triangles, or the neighbour is
    // bisected across a short edge and repeated rounds degenerate it
    // (measured: collinear zero-area triangles after ~90 splits on an
    // anisotropic cylinder chart). A triangle whose chosen edge is not
    // its neighbour's longest first FORCES the neighbour's own
    // longest-edge split; lengths grow strictly along the path, so it
    // terminates.
    let mut queue: Vec<(usize, bool)> = (0..mesh.triangles.len()).map(|t| (t, false)).collect();
    while let Some((ti, forced)) = queue.pop() {
        if mesh.triangles.len() + 2 > budget {
            break;
        }
        let t = mesh.triangles[ti];
        if !forced {
            // Over-bound pairs of this triangle, classified by whether
            // their edge can be split. A triangle is refined only when
            // it has violations AND every one of them is splittable:
            // an over-bound BOUNDARY (or guarded) edge can never be
            // fixed from inside, and splitting the triangle's other
            // edges only fans ever-thinner slivers onto it (measured:
            // an interior vertex driven to v = 1 - 1e-16 against a 90°
            // boundary edge until the budget ran out).
            let mut violations = 0usize;
            let mut limited = false;
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                if singular[a]
                    || singular[b]
                    || angle_between(mesh.vertices[a].normal, mesh.vertices[b].normal) <= max_angle
                {
                    continue;
                }
                violations += 1;
                let e = key(a, b);
                if unsplittable.contains(&e) || edge_tris.get(&e).is_none_or(|v| v.len() != 2) {
                    limited = true;
                }
            }
            if violations == 0 || limited {
                continue;
            }
        }
        let Some(((a, b), edge_measure)) = longest(mesh, &singular, &edge_tris, &unsplittable, t)
        else {
            continue; // boundary- or guard-limited
        };
        let neighbour = edge_tris
            .get(&(a, b))
            .and_then(|v| v.iter().copied().find(|&x| x != ti));
        if let Some(tn) = neighbour
            && let Some((en, measure_n)) = longest(
                mesh,
                &singular,
                &edge_tris,
                &unsplittable,
                mesh.triangles[tn],
            )
            && en != (a, b)
            && exceeds(measure_n, edge_measure)
        {
            queue.push((ti, forced));
            queue.push((tn, true));
            continue;
        }
        let len = edge_measure.1;
        let (va, vb) = (&mesh.vertices[a], &mesh.vertices[b]);
        let mut uv = (0.5 * (va.uv.0 + vb.uv.0), 0.5 * (va.uv.1 + vb.uv.1));
        // An edge leaving a pole (or any point where one chart
        // coordinate stops moving the surface) splits along the OTHER
        // endpoint's iso-line in that coordinate: the pole's own value
        // there is arbitrary (every meridian meets it), so the plain uv
        // midpoint would land on a foreign meridian and the fold guard
        // would refuse it. "Collapsed" is judged at the edge's scale —
        // the coordinate's 3D travel over this edge's parameter span,
        // against the edge length — not by a fixed ratio: a pole vertex
        // inverted to v = π − 3e-6 is 1e-6-regular yet moves 1e-7 over
        // a whole meridian gap.
        let (du, dv) = ((va.uv.0 - vb.uv.0).abs(), (va.uv.1 - vb.uv.1).abs());
        for (s_end, r_end) in [(va, vb), (vb, va)] {
            let pu = surface.partial_u(s_end.uv.0, s_end.uv.1).length();
            let pv = surface.partial_v(s_end.uv.0, s_end.uv.1).length();
            if pu * du < POLE_TRAVEL * len && pv * dv >= POLE_TRAVEL * len {
                uv.0 = r_end.uv.0;
            } else if pv * dv < POLE_TRAVEL * len && pu * du >= POLE_TRAVEL * len {
                uv.1 = r_end.uv.1;
            }
        }
        let position = surface.point_at(uv.0, uv.1);
        let chord_mid = Point3::new(
            0.5 * (va.position.x + vb.position.x),
            0.5 * (va.position.y + vb.position.y),
            0.5 * (va.position.z + vb.position.z),
        );
        // Fold guard: the true midpoint of a sub-facet chord sits
        // within about half the chord of it; an image farther away
        // crossed a chart singularity or seam.
        let unit_normal = (!chart_singular(surface, uv.0, uv.1))
            .then(|| limit_normal(surface, uv.0, uv.1))
            .flatten();
        let Some(n) = unit_normal.filter(|_| position.distance_to(&chord_mid) <= 0.5 * len) else {
            unsplittable.insert((a, b));
            queue.push((ti, forced));
            continue;
        };
        let normal = Vec3::new(sign * n.x, sign * n.y, sign * n.z);
        let m = mesh.vertices.len();
        mesh.vertices.push(MeshVertex {
            position,
            normal,
            uv,
        });
        singular.push(false);
        let incident = edge_tris.remove(&(a, b)).unwrap_or_default();
        for &tj in &incident {
            let tri = mesh.triangles[tj];
            // Rotate so the split edge is (tri[0], tri[1]) in the
            // triangle's own winding, then halve it.
            let Some(r) = (0..3).find(|&k| key(tri[k], tri[(k + 1) % 3]) == (a, b)) else {
                continue; // unreachable: edge_tris lists only incident triangles
            };
            let (p, q, c) = (tri[r], tri[(r + 1) % 3], tri[(r + 2) % 3]);
            let new_index = mesh.triangles.len();
            mesh.triangles[tj] = [p, m, c];
            mesh.triangles.push([m, q, c]);
            // Edge bookkeeping: (p,m) and (m,c) belong to tj; (m,q)
            // and (m,c) to the new triangle; (q,c) moves from tj to
            // the new triangle.
            edge_tris.entry(key(p, m)).or_default().push(tj);
            edge_tris.entry(key(m, q)).or_default().push(new_index);
            let mc = edge_tris.entry(key(m, c)).or_default();
            mc.push(tj);
            mc.push(new_index);
            if let Some(list) = edge_tris.get_mut(&key(q, c)) {
                for x in list.iter_mut() {
                    if *x == tj {
                        *x = new_index;
                    }
                }
            }
            queue.push((tj, false));
            queue.push((new_index, false));
        }
    }
}

#[cfg(test)]
#[path = "angular_refine_tests.rs"]
mod tests;
