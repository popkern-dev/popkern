//! Local thickness of a voxel body (Hildebrand & Rüegsegger 1997): a
//! point's thickness is the diameter of the largest ball inside the
//! body that holds it. The body is the cells whose material is at least
//! ½; balls are centred on the half-cell lattice (cell centres, faces,
//! edges and corners), so a slab n cells thick reads n cells, and each
//! ball's radius is its centre's distance to the nearest empty cell.
//!
//! A convex corner is held only by small balls whatever the body's
//! size, so the thinnest single point says nothing, and the band along
//! a convex edge no member ball reaches is a cell wide on a coarse grid
//! (measured: 40% of the bored block at 1.25 mm cells under a 3 mm
//! member). What is reported is the volume beyond the corner reach of
//! every member ball and its typical (median) thickness.

use super::cells::CellField;

/// The body's volume thinner than a length.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThinRegion {
    /// Its share of the body's volume.
    pub fraction: f64,
    /// Its median local thickness (by volume), mm.
    pub thickness: f64,
}

/// The body of `cells` (material ≥ ½) thinner than `length`: the
/// volume farther than √3·length/2 from every centre of a ball of
/// diameter `length` inside the body (so a cube's corner beside such a
/// ball is not thin, and a plate thinner than `length` is). `None` for
/// an empty body.
pub(super) fn thin_region(cells: &CellField, length: f64) -> Option<ThinRegion> {
    let local = local_thickness(cells, length, false)?;
    let mut thin: Vec<(f64, f64)> = (0..cells.len())
        .filter(|&c| local.weight[c] > 0.0 && !local.thick[c])
        .map(|c| (local.thickness[c], local.weight[c]))
        .collect();
    let volume: f64 = thin.iter().map(|t| t.1).sum();
    Some(ThinRegion {
        fraction: volume / local.total,
        thickness: quantile(&mut thin, 0.5),
    })
}

/// The `share` quantile by volume of the body's member thickness (each
/// cell's, up to `cap`), mm: of a design, a low quantile is its
/// thinnest members. A cell's member thickness is the largest ball
/// that reaches it within √3 of its radius, the reach of a cube's
/// corner: a member's edges and corners read the member, not the small
/// balls that alone cover them (which read a 4-cell slab's edge cells
/// as 2 cells). `None` for an empty body.
pub(super) fn thickness_quantile(cells: &CellField, cap: f64, share: f64) -> Option<f64> {
    let local = local_thickness(cells, cap, true)?;
    let mut all: Vec<(f64, f64)> = (0..cells.len())
        .filter(|&c| local.weight[c] > 0.0)
        .map(|c| (local.thickness[c], local.weight[c]))
        .collect();
    Some(quantile(&mut all, share))
}

/// The `share` quantile of (value, weight) pairs by weight.
fn quantile(pairs: &mut [(f64, f64)], share: f64) -> f64 {
    let total: f64 = pairs.iter().map(|t| t.1).sum();
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut below = 0.0;
    let mut value = 0.0;
    for (t, w) in pairs.iter() {
        value = *t;
        below += w;
        if below >= share * total {
            break;
        }
    }
    value
}

/// Per cell of a body: its volume (0 outside), local thickness (capped
/// at `length`), and whether it lies within the corner reach of a ball
/// of diameter `length`.
struct Local {
    weight: Vec<f64>,
    total: f64,
    thickness: Vec<f64>,
    thick: Vec<bool>,
}

/// The local thickness of `cells`' body up to `length`, each ball
/// reaching √3 of its radius with `corner_reach`.
fn local_thickness(cells: &CellField, length: f64, corner_reach: bool) -> Option<Local> {
    let h = cells.h;
    let solid: Vec<bool> = (0..cells.len()).map(|c| cells.material(c) >= 0.5).collect();
    let weight: Vec<f64> = (0..cells.len())
        .map(|c| if solid[c] { cells.fill[c] } else { 0.0 })
        .collect();
    let total: f64 = weight.iter().sum();
    if total <= 0.0 {
        return None;
    }
    let radius = 0.5 * length;
    // Cells within `reach` of a point, along each axis.
    let span = |reach: f64| (reach / h).ceil() as i64 + 1;
    let dims = cells.dims;
    let lattice = dims.map(|d| 2 * d + 1);
    // Lattice point k sits at corner + k·h/2.
    let point = |k: [usize; 3]| [0, 1, 2].map(|a| cells.corner[a] + 0.5 * h * k[a] as f64);
    let cells_near = |p: [f64; 3], reach: f64| {
        let lo = [0, 1, 2].map(|a| {
            let c = ((p[a] - cells.corner[a]) / h).floor() as i64 - span(reach) + 1;
            c.clamp(0, dims[a] as i64) as usize
        });
        let hi = [0, 1, 2].map(|a| {
            let c = ((p[a] - cells.corner[a]) / h).floor() as i64 + span(reach);
            c.clamp(0, dims[a] as i64) as usize
        });
        (lo, hi)
    };
    // Each lattice point's ball: the distance to the nearest empty
    // cell, capped at the radius (a larger ball decides nothing here).
    let at = |k: [usize; 3]| k[0] + lattice[0] * (k[1] + lattice[1] * k[2]);
    let mut ball = vec![0.0f64; lattice[0] * lattice[1] * lattice[2]];
    for kz in 1..lattice[2] - 1 {
        for ky in 1..lattice[1] - 1 {
            for kx in 1..lattice[0] - 1 {
                let p = point([kx, ky, kz]);
                let (lo, hi) = cells_near(p, radius);
                let mut d = radius;
                for z in lo[2]..hi[2] {
                    for y in lo[1]..hi[1] {
                        for x in lo[0]..hi[0] {
                            let q = [x, y, z];
                            if !solid[cells.index(q)] {
                                d = d.min(cube_distance(p, cells.centre(q), h));
                            }
                        }
                    }
                }
                ball[at([kx, ky, kz])] = d;
            }
        }
    }
    // Each cell's thickness, and whether it lies within √3 r of a
    // centre whose ball is a full member (the reach of a cube's corner
    // from the member ball in it).
    let near = 3.0f64.sqrt() * radius;
    let mut thickness = vec![0.0f64; cells.len()];
    let mut thick = vec![false; cells.len()];
    for kz in 1..lattice[2] - 1 {
        for ky in 1..lattice[1] - 1 {
            for kx in 1..lattice[0] - 1 {
                let d = ball[at([kx, ky, kz])];
                if d <= 0.0 {
                    continue;
                }
                // A full ball whose six neighbours are full too reaches
                // nothing its neighbours do not.
                let member = d >= radius;
                let inner = member
                    && [[1, 0, 0], [0, 1, 0], [0, 0, 1]].iter().all(|o| {
                        let up = [kx + o[0], ky + o[1], kz + o[2]];
                        let down = [kx - o[0], ky - o[1], kz - o[2]];
                        ball[at(up)] >= radius && ball[at(down)] >= radius
                    });
                let p = point([kx, ky, kz]);
                let widen: f64 = if corner_reach { 3.0 } else { 1.0 };
                let own = widen.sqrt() * d;
                let reach = if member && !inner { near.max(own) } else { own };
                let (lo, hi) = cells_near(p, reach);
                for z in lo[2]..hi[2] {
                    for y in lo[1]..hi[1] {
                        for x in lo[0]..hi[0] {
                            let q = [x, y, z];
                            let c = cells.index(q);
                            if !solid[c] {
                                continue;
                            }
                            let o = cells.centre(q);
                            let r2: f64 = (0..3).map(|a| (o[a] - p[a]) * (o[a] - p[a])).sum();
                            if r2 <= widen * d * d * (1.0 + 1e-12) {
                                thickness[c] = thickness[c].max(2.0 * d);
                                thick[c] |= member;
                            }
                            if member && r2 <= near * near * (1.0 + 1e-12) {
                                thick[c] = true;
                            }
                        }
                    }
                }
            }
        }
    }
    Some(Local {
        weight,
        total,
        thickness,
        thick,
    })
}

/// The distance from `p` to the cube of side `h` centred at `c`.
fn cube_distance(p: [f64; 3], c: [f64; 3], h: f64) -> f64 {
    (0..3)
        .map(|a| {
            let d = ((p[a] - c[a]).abs() - 0.5 * h).max(0.0);
            d * d
        })
        .sum::<f64>()
        .sqrt()
}

#[cfg(test)]
#[path = "thickness_tests.rs"]
mod tests;
