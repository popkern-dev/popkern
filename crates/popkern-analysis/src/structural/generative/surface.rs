//! A design's surface: the iso-surface x̃ = ½ of its material as a
//! closed, consistently oriented triangle mesh (PH-117).
//!
//! The design body is the part ∩ {x̄ ≥ ½}. The field is sampled at the
//! cell centres (the dual grid), each cell's value the smaller of its
//! density (fragments weighted by volume) and its fill (an entry's
//! volume over its cell's), cells outside the part 0, so the surface
//! closes. Between a solid cell and an empty one the ½ crossing falls
//! on their shared face: a crisp design keeps its cells' outline,
//! rounded by half a cell at its edges, and a member meeting a face of
//! the part keeps its section there. (Averaging cells to the corner
//! nodes, as this first did, rounded a 2-cell member to a diamond of
//! half its section, and shrank it to a point where it met the part's
//! boundary: measured, the strength check's factor fell to 0.03 at such
//! a pinch.) Marching tetrahedra (every cube cut the same way into six
//! tetrahedra along its main diagonal) has no ambiguous case, and
//! vertices are shared by edge, so the mesh is watertight; each
//! triangle faces away from the material.

use super::cells::CellField;
use super::{Design, DesignEntry};

/// A triangle mesh, mm.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceMesh {
    /// Vertex positions.
    pub vertices: Vec<[f64; 3]>,
    /// Triangles as vertex indices, counter-clockwise seen from outside.
    pub triangles: Vec<[u32; 3]>,
}

impl SurfaceMesh {
    /// The enclosed volume, by the divergence theorem, mm³.
    pub fn volume(&self) -> f64 {
        self.triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| self.vertices[i as usize]);
                (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                    + a[2] * (b[0] * c[1] - b[1] * c[0]))
                    / 6.0
            })
            .sum()
    }

    /// Whether every edge is shared by exactly two triangles that use it
    /// in opposite directions (closed and consistently oriented).
    pub fn is_closed(&self) -> bool {
        let mut directed: Vec<(u32, u32)> = self
            .triangles
            .iter()
            .flat_map(|t| [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])])
            .collect();
        directed.sort_unstable();
        if directed.windows(2).any(|w| w[0] == w[1]) {
            return false;
        }
        directed
            .iter()
            .all(|&(a, b)| directed.binary_search(&(b, a)).is_ok())
    }
}

/// The six tetrahedra of a cube (corner bit 1: +x, 2: +y, 4: +z), all
/// sharing the diagonal 0–7, so neighbouring cubes cut their shared face
/// the same way.
const TETS: [[usize; 4]; 6] = [
    [0, 1, 3, 7],
    [0, 1, 5, 7],
    [0, 2, 3, 7],
    [0, 2, 6, 7],
    [0, 4, 5, 7],
    [0, 4, 6, 7],
];

/// A design's material on its grid: a dense block of cells over the
/// design's bounding box with a zero layer around it, each cell's value
/// the smaller of its density and its fill (0 outside the part). It is
/// what [`design_surface`] marches, kept apart so a host can march a
/// design it recorded mid-solve (a film of the design converging,
/// PK-959).
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialGrid {
    /// Cell size, mm.
    pub cell_size: f64,
    /// Cells along each axis, the zero layer included.
    pub dims: [usize; 3],
    /// The lower corner of cell (0, 0, 0), mm.
    pub corner: [f64; 3],
    /// Per cell (x fastest, then y, then z): its material, 0..1.
    pub material: Vec<f64>,
}

impl MaterialGrid {
    /// The grid of `entries` at `densities` on cells of size
    /// `cell_size`; `None` without entries or with a size not positive.
    pub fn new(entries: &[DesignEntry], densities: &[f64], cell_size: f64) -> Option<MaterialGrid> {
        let cells = CellField::new(entries, densities, cell_size)?;
        Some(MaterialGrid {
            cell_size,
            dims: cells.dims,
            corner: cells.corner,
            material: (0..cells.len()).map(|c| cells.material(c)).collect(),
        })
    }
}

/// The iso-surface of `design` at material fraction `iso` (½ for the
/// design's outline).
pub fn design_surface(design: &Design, iso: f64) -> SurfaceMesh {
    match MaterialGrid::new(&design.entries, &design.densities, design.cell_size) {
        Some(grid) => grid_surface(&grid, iso),
        None => SurfaceMesh {
            vertices: Vec::new(),
            triangles: Vec::new(),
        },
    }
}

/// The iso-surface of `grid` at material fraction `iso`: empty for a
/// grid under two cells along an axis or whose material does not match
/// its dimensions.
pub fn grid_surface(grid: &MaterialGrid, iso: f64) -> SurfaceMesh {
    let cdims = grid.dims;
    if cdims.iter().any(|&d| d < 2) || grid.material.len() != cdims[0] * cdims[1] * cdims[2] {
        return SurfaceMesh {
            vertices: Vec::new(),
            triangles: Vec::new(),
        };
    }
    // The field is sampled at the cell centres (the dual grid).
    let nodes = &grid.material;
    let node_index = |q: [usize; 3]| q[0] + cdims[0] * (q[1] + cdims[1] * q[2]);
    let h = grid.cell_size;
    let position = |q: [usize; 3]| [0, 1, 2].map(|a| grid.corner[a] + h * (q[a] as f64 + 0.5));
    let mut vertices: Vec<[f64; 3]> = Vec::new();
    // Edge (node a, node b), a < b, to its vertex: looked up only, never
    // iterated, so the vertex order is the cubes' order.
    let mut edges: std::collections::HashMap<(usize, usize), u32> =
        std::collections::HashMap::new();
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    // Each cube of eight neighbouring cell centres, by its lowest cell.
    for z in 0..cdims[2] - 1 {
        for y in 0..cdims[1] - 1 {
            for x in 0..cdims[0] - 1 {
                let corner = |c: usize| [x + (c & 1), y + ((c >> 1) & 1), z + ((c >> 2) & 1)];
                for tet in TETS {
                    let q = tet.map(corner);
                    let ids = q.map(node_index);
                    let v = ids.map(|i| nodes[i]);
                    let inside: Vec<usize> = (0..4).filter(|&i| v[i] >= iso).collect();
                    if inside.is_empty() || inside.len() == 4 {
                        continue;
                    }
                    let outside: Vec<usize> = (0..4).filter(|&i| v[i] < iso).collect();
                    let mut vertex = |i: usize, o: usize| -> u32 {
                        let k = (ids[i].min(ids[o]), ids[i].max(ids[o]));
                        match edges.get(&k) {
                            Some(&id) => id,
                            None => {
                                // From the outside node (below iso) to the
                                // inside one (at or above).
                                let (a, b) = (position(q[o]), position(q[i]));
                                let t = (iso - v[o]) / (v[i] - v[o]);
                                let p = [0, 1, 2].map(|d| a[d] + t * (b[d] - a[d]));
                                let id = vertices.len() as u32;
                                vertices.push(p);
                                edges.insert(k, id);
                                id
                            }
                        }
                    };
                    // The crossing edges, as (inside, outside) corners,
                    // in cyclic order around the polygon.
                    let crossing: Vec<(usize, usize)> = match (inside.len(), outside.len()) {
                        (1, 3) => outside.iter().map(|&o| (inside[0], o)).collect(),
                        (3, 1) => inside.iter().map(|&i| (i, outside[0])).collect(),
                        _ => {
                            let (a, b) = (inside[0], inside[1]);
                            let (c, d) = (outside[0], outside[1]);
                            vec![(a, c), (a, d), (b, d), (b, c)]
                        }
                    };
                    let polygon: Vec<u32> = crossing.iter().map(|&(i, o)| vertex(i, o)).collect();
                    // Orient on the edges' midpoints (never degenerate,
                    // where a crossing at a node makes the triangle
                    // itself flat): the normal points from the inside
                    // corners to the outside ones.
                    let mid = |(i, o): (usize, usize)| {
                        let (a, b) = (position(q[i]), position(q[o]));
                        [0, 1, 2].map(|d| 0.5 * (a[d] + b[d]))
                    };
                    let centroid = |set: &[usize]| {
                        let n = set.len() as f64;
                        [0, 1, 2].map(|d| set.iter().map(|&c| position(q[c])[d]).sum::<f64>() / n)
                    };
                    let (cin, cout) = (centroid(&inside), centroid(&outside));
                    let away = [0, 1, 2].map(|d| cout[d] - cin[d]);
                    let mut push = |t: [usize; 3]| {
                        let [a, b, c] = t.map(|k| mid(crossing[k]));
                        let u = [0, 1, 2].map(|d| b[d] - a[d]);
                        let w = [0, 1, 2].map(|d| c[d] - a[d]);
                        let n = [
                            u[1] * w[2] - u[2] * w[1],
                            u[2] * w[0] - u[0] * w[2],
                            u[0] * w[1] - u[1] * w[0],
                        ];
                        let facing = n[0] * away[0] + n[1] * away[1] + n[2] * away[2];
                        let ids = t.map(|k| polygon[k]);
                        triangles.push(if facing >= 0.0 {
                            ids
                        } else {
                            [ids[0], ids[2], ids[1]]
                        });
                    };
                    push([0, 1, 2]);
                    if polygon.len() == 4 {
                        push([0, 2, 3]);
                    }
                }
            }
        }
    }
    SurfaceMesh {
        vertices,
        triangles,
    }
}

#[cfg(test)]
#[path = "surface_tests.rs"]
mod tests;
