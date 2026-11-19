//! The lattice rung's band piece: a uniform lattice of interior points in
//! the face chart (the band holds no singular row) and the band's ring
//! samples, triangulated in the chart metric with the rings as
//! constraints. A ring segment the triangulation leaves uncovered is split
//! at its midpoint (on its edge's curve, or on its row, an iso-line in the
//! band chart) and the triangulation runs again: segment splitting
//! converges where flips alone cannot, so this rung has no coverage
//! decline of its own. Every ring sample is a vertex, so the caps, which
//! carry the band's vertices on their rows, share every vertex there.

use super::super::*;
use crate::mesh::MeshVertex;
use popkern_core::geometry::Vec3;

/// Split rounds before the band declines: each halves every uncovered
/// segment, so twelve resolve a segment 4096 times longer than its
/// neighbourhood's lattice spacing.
const SPLIT_ROUNDS: usize = 12;

/// Mesh the band piece (`rings` charted on one branch of the face chart,
/// ring 0 its outer boundary): `None` when a ring segment stays uncovered
/// after every split round or the triangulation fails.
pub(in crate::boundary_first) fn mesh_band(
    (face, surface, orientation): (FaceHandle, &SurfaceGeometry, bool),
    rings: &[BoundaryRing],
    (brep, descriptor): (&BRep, &ChartDescriptor),
    options: &TessellationOptions,
) -> Option<FaceMesh> {
    let mut rings = rings.to_vec();
    // Holes onto the outer ring's period branch, as the structured band
    // places them (a hole a period off reads as material).
    super::super::super::boundary_first_unwrap::place_rings_on_chart(&mut rings, surface);
    let (mut lo, mut hi) = (
        (f64::INFINITY, f64::INFINITY),
        (f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for &(u, v) in rings.iter().flat_map(|r| r.uv.iter()) {
        lo = (lo.0.min(u), lo.1.min(v));
        hi = (hi.0.max(u), hi.1.max(v));
    }
    if !(hi.0 > lo.0 && hi.1 > lo.1) {
        why(face, "the band's rings bound no box");
        return None;
    }
    let ((u0, u1), (v0, v1)) = ((lo.0, hi.0), (lo.1, hi.1));
    // The chart metric at the box's centre, its shear through the first
    // fundamental form's Cholesky factor (as the structured band reads it).
    let (su, sv, shear, rest) = {
        let (um, vm) = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
        let p = surface.point_at(um, vm);
        let (eu, ev) = (((u1 - u0) * 1e-4).max(1e-9), ((v1 - v0) * 1e-4).max(1e-9));
        let pu = (surface.point_at(um + eu, vm) - p) * (1.0 / eu);
        let pv = (surface.point_at(um, vm + ev) - p) * (1.0 / ev);
        let (du, dv) = (pu.length().max(1e-9), pv.length().max(1e-9));
        let shear = pu.dot(&pv) / du;
        let rest = (dv * dv - shear * shear).max((1e-3 * dv).powi(2)).sqrt();
        (du, dv, shear, rest)
    };
    // Uniform stations meeting the chord height across the box.
    let cap = options
        .max_subdivisions
        .max(options.min_subdivisions)
        .max(1);
    let n = crate::adaptive::chord_station_count(
        surface,
        (u0, u1),
        (v0, v1),
        options.chord_height,
        (options.min_subdivisions.max(1), cap),
    );
    let (wide, tall) = (su * (u1 - u0), sv * (v1 - v0));
    let cell = wide.max(tall) / n as f64;
    let (n_u, n_v) = (
        ((wide / cell).ceil() as usize).max(1),
        ((tall / cell).ceil() as usize).max(1),
    );
    let cells = (wide / n_u as f64, tall / n_v as f64);
    for _ in 0..SPLIT_ROUNDS {
        let index =
            super::super::super::boundary_first_ring_index::RingIndex::new(&rings, su, sv, cells);
        let mut uv_pts: Vec<(f64, f64)> = Vec::new();
        let mut pos_pts: Vec<Option<Point3>> = Vec::new();
        let mut ranges = Vec::new();
        for ring in &rings {
            let start = uv_pts.len();
            uv_pts.extend(ring.uv.iter().copied());
            pos_pts.extend(ring.pos.iter().copied().map(Some));
            ranges.push((start, uv_pts.len()));
        }
        // Stations half a cell in from the box's sides, which ring samples
        // reach: a station on a ring segment is a T-point.
        for i in 0..n_u {
            for j in 0..n_v {
                let p = (
                    u0 + (u1 - u0) * (i as f64 + 0.5) / n_u as f64,
                    v0 + (v1 - v0) * (j as f64 + 0.5) / n_v as f64,
                );
                if index.inside(p) && !index.within_clearance(p) {
                    uv_pts.push(p);
                    pos_pts.push(None);
                }
            }
        }
        let (mut triangles, alias) =
            super::super::super::boundary_first_triangulate::triangulate_region(
                face,
                &uv_pts,
                (su, shear, rest),
                &index,
                &ranges,
            )
            .or_else(|| {
                why(face, "the triangulation failed");
                None
            })?;
        let mut edges = popkern_math::hash::FastHashSet::default();
        for t in &triangles {
            for (x, y) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                edges.insert((x.min(y), x.max(y)));
            }
        }
        let uncovered: Vec<Vec<usize>> = ranges
            .iter()
            .map(|&(start, end)| {
                let m = end - start;
                (0..m)
                    .filter(|&k| {
                        let (a, b) = (alias[start + k], alias[start + (k + 1) % m]);
                        a != b && !edges.contains(&(a.min(b), a.max(b)))
                    })
                    .collect()
            })
            .collect();
        if uncovered.iter().all(Vec::is_empty) {
            // Facets held to the chord height, as every atlas piece is.
            let chart = super::super::super::boundary_first_deviation::Chart {
                face,
                surface,
                planar: false,
                scale: (su, sv),
                tri_scale: (su, shear, rest),
                ring_index: &index,
                ring_ranges: &ranges,
                alias: &alias,
            };
            super::super::super::boundary_first_deviation::refine_by_deviation(
                &chart,
                &mut uv_pts,
                &mut pos_pts,
                &mut triangles,
                options.chord_height,
            );
            return Some(materialize(
                face,
                (surface, orientation),
                uv_pts,
                pos_pts,
                triangles,
            ));
        }
        for (ring, ks) in rings.iter_mut().zip(&uncovered) {
            for &k in ks.iter().rev() {
                split_segment(brep, (surface, descriptor), ring, k);
            }
        }
    }
    why(face, "a ring segment stayed uncovered");
    None
}

fn why(face: FaceHandle, reason: &str) {
    if popkern_core::trace::trace_enabled("atlas") {
        eprintln!("    [atlas] face {face:?}: lattice band declined: {reason}");
    }
}

/// Split ring segment `k` at its midpoint: on its edge's curve at the
/// middle fraction (its parameters the surface's there, on the segment's
/// branch), or, for a row or slit segment (no edge), at the middle of its
/// parameters, an iso-line in the band chart.
fn split_segment(
    brep: &BRep,
    (surface, descriptor): (&SurfaceGeometry, &ChartDescriptor),
    ring: &mut BoundaryRing,
    k: usize,
) {
    let n = ring.uv.len();
    let (a, b) = (ring.uv[k], ring.uv[(k + 1) % n]);
    let mid = (0.5 * (a.0 + b.0), 0.5 * (a.1 + b.1));
    let (edge, fa, fb) = ring.sources[k];
    let curve = (edge != EdgeHandle::default())
        .then(|| brep.edges.get(edge).and_then(|e| e.curve.as_ref()))
        .flatten();
    let (uv, p, halves) = match curve {
        Some(curve) => {
            let f = 0.5 * (fa + fb);
            let p = curve
                .geometry
                .point_at(curve.t_start + (curve.t_end - curve.t_start) * f);
            let uv = surface
                .closest_point(&p)
                .map(|q| super::on_branch_of(q, mid, descriptor))
                .unwrap_or(mid);
            (uv, p, ((edge, fa, f), (edge, f, fb)))
        }
        None => (
            mid,
            surface.point_at(mid.0, mid.1),
            (ring.sources[k], ring.sources[k]),
        ),
    };
    ring.uv.insert(k + 1, uv);
    ring.pos.insert(k + 1, p);
    ring.sources[k] = halves.0;
    ring.sources.insert(k + 1, halves.1);
}

/// The band's face mesh: ring samples keep their positions, interior
/// points take theirs from the surface; normals face along `orientation`,
/// and triangles wind with it.
fn materialize(
    face: FaceHandle,
    (surface, orientation): (&SurfaceGeometry, bool),
    uv_pts: Vec<(f64, f64)>,
    pos_pts: Vec<Option<Point3>>,
    mut triangles: Vec<[usize; 3]>,
) -> FaceMesh {
    let vertices = uv_pts
        .iter()
        .zip(pos_pts)
        .map(|(&(u, v), pos)| {
            let mut normal = surface.normal_at(u, v);
            let len = normal.length();
            if len > 1e-300 {
                normal = Vec3::new(normal.x / len, normal.y / len, normal.z / len);
            }
            if !orientation {
                normal = -normal;
            }
            MeshVertex {
                position: pos.unwrap_or_else(|| surface.point_at(u, v)),
                normal,
                uv: (u, v),
            }
        })
        .collect();
    if !orientation {
        for t in &mut triangles {
            t.swap(1, 2);
        }
    }
    FaceMesh {
        face,
        vertices,
        triangles,
    }
}
