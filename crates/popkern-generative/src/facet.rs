//! A design surface as a body: the closed, oriented triangle mesh of a
//! density design becomes a planar-faced BRep (PK-956).
//!
//! The mesh is welded within the weld tolerance, checked closed and
//! manifold at every edge and vertex, and its triangles merged into
//! planar faces: a patch grows from a seed triangle over neighbours of
//! the same owner (the input face they lie on, or the design's new
//! surface) whose corners lie on the seed's plane, and becomes one face
//! with its holes; a patch whose boundary pinches (a vertex it passes
//! twice) is grown again as disks. A vertex between two straight edges
//! of the same two faces is then dropped. Separate pieces become
//! separate solids. Nothing is approximated beyond the plane tolerance:
//! the body encloses the mesh's volume.

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_part::PartError;
use popkern_part::poly::{PolyFace, build_polyhedron_brep};

/// Why a design surface does not make a body. Every geometric decline
/// names where it happened, in the mesh's coordinates.
#[derive(Debug, thiserror::Error)]
pub enum FacetError {
    /// The surface has no triangles (no part of the design reaches the
    /// iso level).
    #[error("the design surface is empty")]
    Empty,
    /// An edge is not shared by exactly two triangles in opposite
    /// directions.
    #[error("the design surface is not closed at {at:?}: an edge used {uses} times")]
    NotClosed {
        /// The edge's midpoint.
        at: [f64; 3],
        /// Directed uses of the edge.
        uses: usize,
    },
    /// Several triangle fans meet at one vertex (two members touching
    /// at a point).
    #[error("the design surface pinches at {at:?}: {fans} triangle fans meet there")]
    PinchedVertex {
        /// The vertex.
        at: [f64; 3],
        /// Separate fans around it.
        fans: usize,
    },
    /// A triangle with no area whose corners are distinct (a sliver the
    /// weld cannot remove).
    #[error("degenerate triangle at {at:?}")]
    DegenerateTriangle {
        /// Its centroid.
        at: [f64; 3],
    },
    /// The polyhedron builder or its validation declined.
    #[error("the faceted body does not build: {0}")]
    Polyhedron(Box<PartError>),
}

/// Tolerances of the faceting, mm.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FacetOptions {
    /// Vertices closer than this are one.
    pub weld: f64,
    /// A triangle joins a face when its corners lie within this of the
    /// face's plane (and a vertex is dropped when it lies within this of
    /// the line through its neighbours).
    pub plane: f64,
}

impl FacetOptions {
    /// Tolerances for a body built with `tolerance`: welds at a tenth
    /// of its linear tolerance, planes at a hundredth (the design
    /// surface's crossings on a cell face are exact to rounding).
    pub fn for_tolerance(tolerance: &Tolerance) -> FacetOptions {
        FacetOptions {
            weld: 0.1 * tolerance.linear(),
            plane: 0.01 * tolerance.linear(),
        }
    }
}

/// What the faceting did, for the feature's report.
#[derive(Debug, Clone, Copy, PartialEq, Default, serde::Serialize, serde::Deserialize)]
pub struct FacetStats {
    /// Triangles in.
    pub triangles: usize,
    /// Vertices merged by the weld.
    pub welded: usize,
    /// Triangles dropped as collapsed by the weld.
    pub collapsed: usize,
    /// Faces out.
    pub faces: usize,
    /// Vertices dropped between collinear edges.
    pub collinear: usize,
    /// The volume the welded mesh encloses (the body's, to rounding:
    /// the faceting moves no vertex).
    pub mesh_volume: f64,
}

/// The faceted body: its BRep, each face's owner in face order, and the
/// stats.
#[derive(Debug)]
pub struct FacetedBody {
    /// The body (validated).
    pub brep: BRep,
    /// Per face, in the BRep's face order: the owner of its triangles.
    pub owners: Vec<u32>,
    /// What the faceting did.
    pub stats: FacetStats,
}

type V3 = [f64; 3];

fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn mid(a: V3, b: V3) -> V3 {
    [
        0.5 * (a[0] + b[0]),
        0.5 * (a[1] + b[1]),
        0.5 * (a[2] + b[2]),
    ]
}

/// Welds `vertices` within `tol`: each vertex maps to the first earlier
/// one within the tolerance (a hash grid of cell `tol`, so the search
/// is the 27 neighbouring cells). Returns the map and the merged count.
fn weld(vertices: &[V3], tol: f64) -> (Vec<u32>, usize) {
    let cell = tol.max(f64::MIN_POSITIVE);
    let key = |p: V3| p.map(|x| (x / cell).floor() as i64);
    let mut grid: HashMap<[i64; 3], Vec<u32>> = HashMap::new();
    let mut map = Vec::with_capacity(vertices.len());
    let mut merged = 0;
    for (i, &p) in vertices.iter().enumerate() {
        let k = key(p);
        let mut found = None;
        'search: for dz in -1..=1 {
            for dy in -1..=1 {
                for dx in -1..=1 {
                    if let Some(list) = grid.get(&[k[0] + dx, k[1] + dy, k[2] + dz]) {
                        for &j in list {
                            let d = sub(p, vertices[j as usize]);
                            if dot(d, d) <= tol * tol {
                                found = Some(j);
                                break 'search;
                            }
                        }
                    }
                }
            }
        }
        match found {
            Some(j) => {
                map.push(j);
                merged += 1;
            }
            None => {
                map.push(i as u32);
                grid.entry(k).or_default().push(i as u32);
            }
        }
    }
    (map, merged)
}

/// The mesh after the weld: triangles and owners with collapsed
/// triangles dropped.
struct Welded {
    triangles: Vec<[u32; 3]>,
    owners: Vec<u32>,
    collapsed: usize,
    welded: usize,
}

fn welded(vertices: &[V3], triangles: &[[u32; 3]], owners: &[u32], tol: f64) -> Welded {
    let (map, welded) = weld(vertices, tol);
    let mut out = Welded {
        triangles: Vec::with_capacity(triangles.len()),
        owners: Vec::with_capacity(triangles.len()),
        collapsed: 0,
        welded,
    };
    for (t, &o) in triangles.iter().zip(owners) {
        let t = t.map(|i| map[i as usize]);
        if t[0] == t[1] || t[1] == t[2] || t[2] == t[0] {
            out.collapsed += 1;
        } else {
            out.triangles.push(t);
            out.owners.push(o);
        }
    }
    out
}

/// Directed edge → triangle, after checking every edge is used once in
/// each direction.
fn edge_map(
    vertices: &[V3],
    triangles: &[[u32; 3]],
) -> Result<HashMap<(u32, u32), u32>, FacetError> {
    let mut edges: HashMap<(u32, u32), u32> = HashMap::with_capacity(3 * triangles.len());
    for (t, tri) in triangles.iter().enumerate() {
        for k in 0..3 {
            let e = (tri[k], tri[(k + 1) % 3]);
            if edges.insert(e, t as u32).is_some() {
                return Err(FacetError::NotClosed {
                    at: mid(vertices[e.0 as usize], vertices[e.1 as usize]),
                    uses: 2,
                });
            }
        }
    }
    for &(a, b) in edges.keys() {
        if !edges.contains_key(&(b, a)) {
            return Err(FacetError::NotClosed {
                at: mid(vertices[a as usize], vertices[b as usize]),
                uses: 1,
            });
        }
    }
    Ok(edges)
}

/// Declines a vertex where more than one fan of triangles meets: the
/// triangles around it, linked across the edges they share at it, must
/// form one cycle.
fn check_vertex_fans(
    vertices: &[V3],
    triangles: &[[u32; 3]],
    edges: &HashMap<(u32, u32), u32>,
) -> Result<(), FacetError> {
    let mut around: Vec<u32> = vec![0; vertices.len()];
    for tri in triangles {
        for &v in tri {
            around[v as usize] += 1;
        }
    }
    // Walking one fan from each vertex's first triangle: the fan's
    // length must be every triangle at the vertex.
    let mut first: Vec<Option<u32>> = vec![None; vertices.len()];
    for (t, tri) in triangles.iter().enumerate() {
        for &v in tri {
            first[v as usize].get_or_insert(t as u32);
        }
    }
    for (v, start) in first.iter().enumerate() {
        let Some(start) = *start else { continue };
        let v = v as u32;
        let mut t = start;
        let mut count = 0u32;
        loop {
            count += 1;
            let tri = triangles[t as usize];
            // The edge leaving v in t, (v, w); its twin (w, v) is in the
            // next triangle around v.
            let k = tri.iter().position(|&x| x == v).unwrap_or(0);
            let w = tri[(k + 1) % 3];
            match edges.get(&(w, v)) {
                Some(&next) if next != start => t = next,
                _ => break,
            }
            if count > around[v as usize] {
                break;
            }
        }
        if count != around[v as usize] {
            return Err(FacetError::PinchedVertex {
                at: vertices[v as usize],
                fans: (around[v as usize] as usize)
                    .div_ceil(count as usize)
                    .max(2),
            });
        }
    }
    Ok(())
}

/// Each triangle's unit normal and plane offset; a triangle whose
/// height is under the weld tolerance declines.
fn planes(vertices: &[V3], triangles: &[[u32; 3]], tol: f64) -> Result<Vec<(V3, f64)>, FacetError> {
    triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| vertices[i as usize]);
            let n = cross(sub(b, a), sub(c, a));
            let len = dot(n, n).sqrt();
            let longest = [sub(b, a), sub(c, b), sub(a, c)]
                .iter()
                .map(|e| dot(*e, *e).sqrt())
                .fold(0.0, f64::max);
            if len <= tol * longest {
                let at = [0, 1, 2].map(|d| (a[d] + b[d] + c[d]) / 3.0);
                return Err(FacetError::DegenerateTriangle { at });
            }
            let n = n.map(|x| x / len);
            Ok((n, dot(n, a)))
        })
        .collect()
}

#[path = "facet_regions.rs"]
mod regions;

/// The body enclosed by the closed mesh `vertices`/`triangles`
/// (counter-clockwise from outside), each triangle tagged with an
/// `owner`; faces never merge triangles of different owners.
///
/// # Errors
///
/// [`FacetError`]: an empty, open, pinched or degenerate mesh, or a
/// body the polyhedron builder declines.
pub fn facet(
    vertices: &[V3],
    triangles: &[[u32; 3]],
    owners: &[u32],
    options: &FacetOptions,
    tolerance: Tolerance,
) -> Result<FacetedBody, FacetError> {
    if triangles.is_empty() {
        return Err(FacetError::Empty);
    }
    let w = welded(vertices, triangles, owners, options.weld);
    if w.triangles.is_empty() {
        return Err(FacetError::Empty);
    }
    let edges = edge_map(vertices, &w.triangles)?;
    check_vertex_fans(vertices, &w.triangles, &edges)?;
    let planes = planes(vertices, &w.triangles, options.weld)?;
    let faces = regions::grow(
        vertices,
        &w.triangles,
        &w.owners,
        &planes,
        &edges,
        options.plane,
    );
    // Every loop, outer then holes per face, so the collinear pass sees
    // each edge from both sides.
    let mut loops: Vec<Vec<u32>> = Vec::new();
    let mut face_loops: Vec<std::ops::Range<usize>> = Vec::new();
    for f in &faces {
        let first = loops.len();
        loops.push(f.cycle.clone());
        loops.extend(f.holes.iter().cloned());
        face_loops.push(first..loops.len());
    }
    let collinear = regions::drop_collinear(vertices, &mut loops, options.plane);
    // Compact the vertices the loops use.
    let mut index: HashMap<u32, usize> = HashMap::new();
    let mut points: Vec<Point3> = Vec::new();
    let mut compact = |c: &[u32]| -> Vec<usize> {
        c.iter()
            .map(|&v| {
                *index.entry(v).or_insert_with(|| {
                    let p = vertices[v as usize];
                    points.push(Point3::new(p[0], p[1], p[2]));
                    points.len() - 1
                })
            })
            .collect()
    };
    let poly: Vec<PolyFace> = face_loops
        .iter()
        .map(|r| PolyFace {
            outer: compact(&loops[r.start]),
            inners: loops[r.start + 1..r.end]
                .iter()
                .map(|c| compact(c))
                .collect(),
        })
        .collect();
    let brep = build_polyhedron_brep(&points, &poly, tolerance)
        .map_err(|e| FacetError::Polyhedron(Box::new(e)))?;
    Ok(FacetedBody {
        brep,
        owners: faces.iter().map(|f| f.owner).collect(),
        stats: FacetStats {
            triangles: triangles.len(),
            welded: w.welded,
            collapsed: w.collapsed,
            faces: faces.len(),
            collinear,
            mesh_volume: w
                .triangles
                .iter()
                .map(|t| {
                    let [a, b, c] = t.map(|i| vertices[i as usize]);
                    dot(a, cross(b, c)) / 6.0
                })
                .sum(),
        },
    })
}

#[cfg(test)]
#[path = "facet_tests.rs"]
mod tests;
