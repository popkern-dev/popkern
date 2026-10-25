//! A closed mesh cut to an axis-aligned box: the part of a solid inside
//! the box, closed by caps on the box's faces (the submodel's domain).
//!
//! The mesh is cut by one half-space at a time. Triangles crossing the
//! plane are split at their edges (each crossing edge's point shared by
//! both triangles, so the result stays watertight); the cut edges,
//! reversed, chain into the section's loops, which are triangulated in
//! the plane (holes bridged to their outer loop, then ear clipping) as
//! the cap. A plane within 1e-6 of the solid's extent of a mesh vertex
//! is moved off it by tiny, deterministic steps.

use super::{FaceId, SolidMesh};

/// Where a triangle of the cut mesh comes from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Origin {
    /// Part of a face of the original solid.
    Face(FaceId),
    /// A cap on the box face normal to `axis` (`upper`: the face at the
    /// box's upper bound).
    Cap {
        /// The box face's normal axis.
        axis: usize,
        /// Whether it is the upper face.
        upper: bool,
    },
}

/// A solid cut to a box.
#[derive(Debug, Clone, PartialEq)]
pub struct CutMesh {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Triangles, outward.
    pub triangles: Vec<[u32; 3]>,
    /// The origin of each triangle.
    pub origin: Vec<Origin>,
    /// The box actually cut to (each bound moved off mesh vertices).
    pub lo: [f64; 3],
    /// Upper corner of that box.
    pub hi: [f64; 3],
}

/// Why a cut failed.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum CutError {
    /// The solid has no material inside the box.
    #[error("the solid does not reach into the box")]
    Empty,
    /// A section's cut edges do not close into loops (the input is not
    /// closed, or the cut met round-off it cannot resolve).
    #[error("the section on box face {axis} does not close")]
    OpenSection {
        /// The plane's normal axis.
        axis: usize,
    },
    /// A section loop could not be triangulated.
    #[error("the section on box face {axis} cannot be triangulated")]
    Triangulation {
        /// The plane's normal axis.
        axis: usize,
    },
}

/// Cuts closed, outward `mesh` to the box [`lo`, `hi`].
pub fn cut_to_box(mesh: &SolidMesh, lo: [f64; 3], hi: [f64; 3]) -> Result<CutMesh, CutError> {
    let mut vertices = mesh.vertices.clone();
    let mut triangles = mesh.triangles.clone();
    let mut origin: Vec<Origin> = mesh
        .triangle_face
        .iter()
        .map(|f| Origin::Face(*f))
        .collect();
    let (mut blo, mut bhi) = (lo, hi);
    let extent = (0..3)
        .map(|a| {
            let (l, u) = bounds(&vertices, a);
            u - l
        })
        .fold(0.0, f64::max);
    for axis in 0..3 {
        for upper in [false, true] {
            let (vmin, vmax) = bounds(&vertices, axis);
            let mut c = if upper { bhi[axis] } else { blo[axis] };
            // A plane beyond the solid cuts nothing.
            if (upper && c >= vmax) || (!upper && c <= vmin) {
                continue;
            }
            c = off_vertices(&vertices, axis, c, extent);
            if upper {
                bhi[axis] = c;
            } else {
                blo[axis] = c;
            }
            let cut = cut_half(&vertices, &triangles, &origin, axis, c, upper, extent)?;
            vertices = cut.0;
            triangles = cut.1;
            origin = cut.2;
        }
    }
    if triangles.is_empty() {
        return Err(CutError::Empty);
    }
    // Drop unused vertices, keeping order.
    let mut used = vec![u32::MAX; vertices.len()];
    let mut kept = Vec::new();
    for t in &triangles {
        for &v in t {
            if used[v as usize] == u32::MAX {
                used[v as usize] = kept.len() as u32;
                kept.push(vertices[v as usize]);
            }
        }
    }
    let triangles = triangles
        .iter()
        .map(|t| t.map(|v| used[v as usize]))
        .collect();
    Ok(CutMesh {
        vertices: kept,
        triangles,
        origin,
        lo: blo,
        hi: bhi,
    })
}

fn bounds(vertices: &[[f64; 3]], axis: usize) -> (f64, f64) {
    vertices
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, u), v| {
            (l.min(v[axis]), u.max(v[axis]))
        })
}

/// `c`, stepped by 3e-6 of `extent` until no vertex lies within 1e-6
/// of `extent` of it (a crossing point nearly on a vertex makes a
/// near-degenerate section).
fn off_vertices(vertices: &[[f64; 3]], axis: usize, c: f64, extent: f64) -> f64 {
    let tol = 1.0e-6 * extent;
    let mut c = c;
    for _ in 0..64 {
        if vertices.iter().all(|v| (v[axis] - c).abs() > tol) {
            return c;
        }
        c += 3.0e-6 * extent;
    }
    c
}

type Cut = (Vec<[f64; 3]>, Vec<[u32; 3]>, Vec<Origin>);

/// Keeps the part of the mesh on the inner side of the plane x_axis = c
/// (below it when `upper`), capped.
fn cut_half(
    vertices: &[[f64; 3]],
    triangles: &[[u32; 3]],
    origin: &[Origin],
    axis: usize,
    c: f64,
    upper: bool,
    extent: f64,
) -> Result<Cut, CutError> {
    let inside = |v: [f64; 3]| if upper { v[axis] < c } else { v[axis] > c };
    let mut out_v = vertices.to_vec();
    let mut out_t = Vec::new();
    let mut out_o = Vec::new();
    // The crossing point of each cut edge, by its (lower, higher) ends,
    // welded by position: overlapping collinear edges (zero-area
    // triangles in a tessellation) cross at one point.
    let mut crossing: std::collections::BTreeMap<(u32, u32), u32> = Default::default();
    let mut at_position: std::collections::BTreeMap<[i64; 3], u32> = Default::default();
    let quantum = 1.0e-10 * extent;
    let mut point = |a: u32, b: u32, out_v: &mut Vec<[f64; 3]>| -> u32 {
        let key = (a.min(b), a.max(b));
        if let Some(&p) = crossing.get(&key) {
            return p;
        }
        let (pa, pb) = (vertices[key.0 as usize], vertices[key.1 as usize]);
        let t = (c - pa[axis]) / (pb[axis] - pa[axis]);
        let mut p = [0, 1, 2].map(|k| pa[k] + t * (pb[k] - pa[k]));
        p[axis] = c;
        let q = p.map(|x| (x / quantum).round() as i64);
        let index = *at_position.entry(q).or_insert_with(|| {
            out_v.push(p);
            (out_v.len() - 1) as u32
        });
        crossing.insert(key, index);
        index
    };
    // Cap boundary segments, directed for the cap (reverse of the kept
    // polygons' cut edges).
    let mut segments: Vec<(u32, u32)> = Vec::new();
    for (t, o) in triangles.iter().zip(origin) {
        let flags = t.map(|v| inside(vertices[v as usize]));
        match flags.iter().filter(|f| **f).count() {
            3 => {
                out_t.push(*t);
                out_o.push(*o);
            }
            0 => {}
            _ => {
                // Walk the triangle's edges, keeping inside corners and
                // inserting crossing points.
                let mut poly: Vec<u32> = Vec::new();
                let mut cut_edge: Vec<u32> = Vec::new();
                for k in 0..3 {
                    let (a, b) = (t[k], t[(k + 1) % 3]);
                    let (ia, ib) = (flags[k], flags[(k + 1) % 3]);
                    if ia {
                        poly.push(a);
                    }
                    if ia != ib {
                        let p = point(a, b, &mut out_v);
                        poly.push(p);
                        cut_edge.push(p);
                    }
                }
                // The polygon's edge along the plane runs from the point
                // where it leaves the inside to where it comes back.
                if let [p, q] = cut_edge[..] {
                    let (from, to) = if poly_edge(&poly, p, q) {
                        (p, q)
                    } else {
                        (q, p)
                    };
                    // A cut edge welded to a point bounds no cap.
                    if to != from {
                        segments.push((to, from));
                    }
                }
                for k in 1..poly.len().saturating_sub(1) {
                    out_t.push([poly[0], poly[k], poly[k + 1]]);
                    out_o.push(*o);
                }
            }
        }
    }
    let cap = cap_triangles(&out_v, &segments, axis, upper)?;
    for t in cap {
        out_t.push(t);
        out_o.push(Origin::Cap { axis, upper });
    }
    Ok((out_v, out_t, out_o))
}

/// Whether `p` is immediately followed by `q` in polygon `poly`.
fn poly_edge(poly: &[u32], p: u32, q: u32) -> bool {
    (0..poly.len()).any(|k| poly[k] == p && poly[(k + 1) % poly.len()] == q)
}

/// Triangles of the cap bounded by the directed `segments` in the plane
/// normal to `axis`, outward (along +axis when `upper`).
fn cap_triangles(
    vertices: &[[f64; 3]],
    segments: &[(u32, u32)],
    axis: usize,
    upper: bool,
) -> Result<Vec<[u32; 3]>, CutError> {
    if segments.is_empty() {
        return Ok(Vec::new());
    }
    // Chain into loops.
    let mut next: std::collections::BTreeMap<u32, u32> = Default::default();
    for &(a, b) in segments {
        if next.insert(a, b).is_some() {
            return Err(CutError::OpenSection { axis });
        }
    }
    let mut loops: Vec<Vec<u32>> = Vec::new();
    let mut seen: std::collections::BTreeSet<u32> = Default::default();
    for &(start, _) in segments {
        if seen.contains(&start) {
            continue;
        }
        let mut l = vec![start];
        seen.insert(start);
        let mut v = start;
        loop {
            let Some(&n) = next.get(&v) else {
                return Err(CutError::OpenSection { axis });
            };
            if n == start {
                break;
            }
            if !seen.insert(n) {
                return Err(CutError::OpenSection { axis });
            }
            l.push(n);
            v = n;
        }
        loops.push(l);
    }
    // 2-D frame (u, w) with the cap's outward normal toward the viewer:
    // counter-clockwise loops are outer boundaries.
    let (u, w) = ((axis + 1) % 3, (axis + 2) % 3);
    let flip = !upper;
    let at = |v: u32| {
        let p = vertices[v as usize];
        if flip { [p[w], p[u]] } else { [p[u], p[w]] }
    };
    super::cap_fill::triangulate(&loops, &at).ok_or(CutError::Triangulation { axis })
}

#[cfg(test)]
#[path = "box_cut_tests.rs"]
mod tests;
