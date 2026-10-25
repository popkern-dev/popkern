//! Local wall thickness: from a surface point, a ray along the inward
//! normal is marched through the grid cells (Amanatides–Woo) and tested
//! against each cell's clipped facet pieces; the first hit is the
//! thickness. Used for the thin-wall `UnderResolved` policy (E4: walls
//! under 1.5 cells).
//!
//! A reading is a wall only when the surface the ray leaves through
//! faces back at the start (outward normal within 60° of the ray, as
//! the far side of a wall does). Beside an acute convex edge the ray
//! instead crosses the wedge and leaves through the other face, nearly
//! side-on (PH-055: 0.004 mm on tapered_fillet's end, at 2 cells).

use crate::clip::{ClippedSolid, FacetPiece};
use crate::grid::Grid;
use crate::model::SolidMesh;

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Least cosine between the ray and the outward normal of the surface
/// it leaves through for the reading to be a wall.
pub const WALL_FACING_COS: f64 = 0.5;

/// Distance along the ray `p + t d` (t > `t_min`) to `piece`, if hit.
fn hit(piece: &FacetPiece, p: [f64; 3], d: [f64; 3], t_min: f64) -> Option<f64> {
    let n = piece.normal;
    let denom = dot(d, n);
    // Grazing planes give round-off "hits" at any distance.
    if denom.abs() < 1e-6 {
        return None;
    }
    let v = piece.polygon.vertices();
    let t = dot(sub(v[0], p), n) / denom;
    if t <= t_min {
        return None;
    }
    let q = [p[0] + t * d[0], p[1] + t * d[1], p[2] + t * d[2]];
    let scale = piece.area.sqrt();
    for i in 0..v.len() {
        let (a, b) = (v[i], v[(i + 1) % v.len()]);
        if dot(cross(sub(b, a), sub(q, a)), n) < -1e-12 * scale * scale {
            return None;
        }
    }
    Some(t)
}

/// Distance from `p` along unit `d` to the solid's surface, if within
/// `max_dist`, the surface faces back at `p` and its triangle is not
/// `unreliable`.
pub fn ray_to_surface(
    grid: &Grid,
    clipped: &ClippedSolid,
    p: [f64; 3],
    d: [f64; 3],
    max_dist: f64,
    unreliable: &[bool],
) -> Option<f64> {
    let h = grid.cell_size();
    let o = grid.origin();
    let dims = grid.dims();
    // Hits closer than this are the starting surface itself; a real
    // wall that thin is far below any resolution anyway.
    let t_min = 1e-3 * h;
    let start = [0, 1, 2].map(|a| p[a] + t_min * d[a]);
    let mut ijk = [0usize; 3];
    let mut step = [0i64; 3];
    let mut t_next = [f64::INFINITY; 3];
    let mut t_delta = [f64::INFINITY; 3];
    for a in 0..3 {
        let c = ((start[a] - o[a]) / h).floor();
        if c < 0.0 || c >= dims[a] as f64 {
            return None;
        }
        ijk[a] = c as usize;
        if d[a] > 0.0 {
            step[a] = 1;
            t_next[a] = (o[a] + (c + 1.0) * h - p[a]) / d[a];
            t_delta[a] = h / d[a];
        } else if d[a] < 0.0 {
            step[a] = -1;
            t_next[a] = (o[a] + c * h - p[a]) / d[a];
            t_delta[a] = -h / d[a];
        }
    }
    let pieces = &clipped.pieces;
    loop {
        let cell = grid.cell_index(ijk)?;
        let lo = pieces.partition_point(|q| q.cell < cell);
        let hi = pieces.partition_point(|q| q.cell <= cell);
        let best = pieces[lo..hi]
            .iter()
            .filter_map(|q| {
                let bad = unreliable
                    .get(q.triangle as usize)
                    .copied()
                    .unwrap_or(false);
                hit(q, p, d, t_min).map(|t| (t, if bad { f64::NAN } else { dot(q.normal, d) }))
            })
            .fold((f64::INFINITY, 0.0), |a, b| if b.0 < a.0 { b } else { a });
        let exit = t_next.iter().copied().fold(f64::INFINITY, f64::min);
        if best.0 <= exit {
            // The first surface met ends the material along the ray: a
            // wall only when it faces back at the start (never when it
            // is unreliable: NaN fails the comparison).
            return (best.0 <= max_dist && best.1 >= WALL_FACING_COS).then_some(best.0);
        }
        if exit > max_dist {
            return None;
        }
        let axis = (0..3)
            .min_by(|&a, &b| t_next[a].total_cmp(&t_next[b]))
            .unwrap_or(0);
        let next = ijk[axis] as i64 + step[axis];
        if next < 0 || next >= dims[axis] as i64 {
            return None;
        }
        ijk[axis] = next as usize;
        t_next[axis] += t_delta[axis];
    }
}

/// Wall thickness at each triangle of `mesh` (from its centroid along
/// the inward normal), `None` where thicker than `max_dist`, where the
/// triangle is `unreliable`, or where the ray first meets an unreliable
/// one.
pub fn triangle_thickness(
    grid: &Grid,
    clipped: &ClippedSolid,
    mesh: &SolidMesh,
    unreliable: &[bool],
    max_dist: f64,
) -> Vec<Option<f64>> {
    mesh.triangles
        .iter()
        .zip(unreliable.iter().copied().chain(std::iter::repeat(false)))
        .map(|(t, bad)| {
            if bad {
                return None;
            }
            let [a, b, c] = t.map(|i| mesh.vertices[i as usize]);
            let n = cross(sub(b, a), sub(c, a));
            let len = dot(n, n).sqrt();
            if len == 0.0 {
                return None;
            }
            let inward = n.map(|x| -x / len);
            let centroid = [0, 1, 2].map(|k| (a[k] + b[k] + c[k]) / 3.0);
            ray_to_surface(grid, clipped, centroid, inward, max_dist, unreliable)
        })
        .collect()
}

/// Thin-wall verdicts of one solid's display mesh.
pub struct ThinWalls {
    /// Thinnest wall below the threshold, if any.
    pub thinnest: Option<f64>,
    /// Per triangle: its reading is below the threshold.
    pub triangle: Vec<bool>,
    /// Per vertex: on a thin triangle.
    pub vertex: Vec<bool>,
}

/// Triangles of `mesh` whose wall reading (searched to 4 cells) is under
/// `THIN_WALL_CELLS` cells, their vertices, and the thinnest reading.
pub fn thin_walls(
    grid: &Grid,
    clipped: &ClippedSolid,
    mesh: &SolidMesh,
    unreliable: &[bool],
) -> ThinWalls {
    let h = grid.cell_size();
    let threshold = super::THIN_WALL_CELLS * h;
    let thickness = triangle_thickness(grid, clipped, mesh, unreliable, 4.0 * h);
    let triangle: Vec<bool> = thickness
        .iter()
        .map(|t| t.is_some_and(|w| w < threshold))
        .collect();
    let mut vertex = vec![false; mesh.vertices.len()];
    let mut thinnest: Option<f64> = None;
    for ((t, tri), thin) in thickness.iter().zip(&mesh.triangles).zip(&triangle) {
        if let (Some(w), true) = (t, thin) {
            thinnest = Some(thinnest.map_or(*w, |m| m.min(*w)));
            for &v in tri {
                vertex[v as usize] = true;
            }
        }
    }
    ThinWalls {
        thinnest,
        triangle,
        vertex,
    }
}
