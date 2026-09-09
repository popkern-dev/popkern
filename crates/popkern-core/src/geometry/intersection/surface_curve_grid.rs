//! Knot-span-aware surface AABB grid for the general surface-curve
//! fallback (split from `surface_curve_general.rs` under the
//! 500-line file budget).

use crate::geometry::Point3;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// Number of cells per axis in the surface bounding-box grid.
pub(super) const GRID_N: usize = 8;

// ---------------------------------------------------------------------------
// Surface AABB grid
// ---------------------------------------------------------------------------

pub(super) struct SurfaceCell {
    pub(super) u_range: (f64, f64),
    pub(super) v_range: (f64, f64),
    pub(super) min: Point3,
    pub(super) max: Point3,
}

/// Grid density per axis: knot-span-aware for NURBS surfaces. A fixed
/// `n` under-seeds a long wrapped chart — an 8-wide grid over a 6-turn
/// helical wall put Newton seeds up to 3/8 turn from the crossing, the
/// solves diverged, and every classification ray through the coil
/// counted ZERO wall crossings (interior points read Outside and the
/// double-grind's kept cap fragment was discarded — PK-441). Knot
/// spans measure where the surface actually bends; the cap keeps the
/// candidate set bounded.
pub(super) fn grid_density(surface: &SurfaceGeometry, base: usize) -> (usize, usize) {
    match surface {
        SurfaceGeometry::Nurbs(s) => {
            let spans = |knots: &[f64]| {
                knots
                    .windows(2)
                    .filter(|w| w[1] > w[0] + 1.0e-12)
                    .count()
                    .max(1)
            };
            let d = (
                spans(&s.knots_u).clamp(base, 48),
                spans(&s.knots_v).clamp(base, 48),
            );
            if crate::trace::trace_enabled("ray-classify") {
                eprintln!(
                    "  [scg-grid] nurbs spans u={} v={} -> grid {}x{}",
                    spans(&s.knots_u),
                    spans(&s.knots_v),
                    d.0,
                    d.1
                );
            }
            d
        }
        _ => (base, base),
    }
}

/// Build one cell whose AABB CONTAINS its surface patch: the box of
/// the four corners plus the centre evaluation, inflated by the
/// centre's deviation from the corner-only box. A corner-only box
/// under-covers a curved cell by the sagitta — a 48-cell grid over a
/// 6.5-turn helical thread wall left each cell's bulge 0.27 world
/// units outside its box, every classification ray through the coil
/// was AABB-culled before Newton ran, interior points read Outside,
/// and the fragment keep/discard verdicts scrambled per sample
/// (PK-460). The centre deviation measures the local sagitta; adding
/// it once more bounds the quarter-cell bulges the centre itself
/// cannot see (quadratic-safe), and a flat cell inflates by zero, so
/// planar/fine geometry keeps byte-identical boxes.
pub(super) fn padded_cell(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> SurfaceCell {
    let corners = [
        surface.point_at(u_range.0, v_range.0),
        surface.point_at(u_range.0, v_range.1),
        surface.point_at(u_range.1, v_range.0),
        surface.point_at(u_range.1, v_range.1),
    ];
    let mut bmin = corners[0];
    let mut bmax = corners[0];
    for c in &corners[1..] {
        bmin = Point3::new(bmin.x.min(c.x), bmin.y.min(c.y), bmin.z.min(c.z));
        bmax = Point3::new(bmax.x.max(c.x), bmax.y.max(c.y), bmax.z.max(c.z));
    }
    let centre = surface.point_at(0.5 * (u_range.0 + u_range.1), 0.5 * (v_range.0 + v_range.1));
    let dev = {
        let dx = (bmin.x - centre.x).max(centre.x - bmax.x).max(0.0);
        let dy = (bmin.y - centre.y).max(centre.y - bmax.y).max(0.0);
        let dz = (bmin.z - centre.z).max(centre.z - bmax.z).max(0.0);
        (dx * dx + dy * dy + dz * dz).sqrt()
    };
    bmin = Point3::new(
        bmin.x.min(centre.x) - dev,
        bmin.y.min(centre.y) - dev,
        bmin.z.min(centre.z) - dev,
    );
    bmax = Point3::new(
        bmax.x.max(centre.x) + dev,
        bmax.y.max(centre.y) + dev,
        bmax.z.max(centre.z) + dev,
    );
    SurfaceCell {
        u_range,
        v_range,
        min: bmin,
        max: bmax,
    }
}

/// Build a grid of surface cells with patch-containing AABBs (see
/// [`padded_cell`]). Density is per-axis (see [`grid_density`]).
pub(super) fn build_surface_grid(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    base_n: usize,
) -> Vec<SurfaceCell> {
    let (nu, nv) = grid_density(surface, base_n);
    let (u0, u1) = u_range;
    let (v0, v1) = v_range;
    let du = (u1 - u0) / nu as f64;
    let dv = (v1 - v0) / nv as f64;
    let mut cells = Vec::with_capacity(nu * nv);
    for iu in 0..nu {
        for iv in 0..nv {
            cells.push(padded_cell(
                surface,
                (u0 + iu as f64 * du, u0 + (iu + 1) as f64 * du),
                (v0 + iv as f64 * dv, v0 + (iv + 1) as f64 * dv),
            ));
        }
    }
    cells
}

/// Interior distinct knot values strictly inside `(a, b)` — the
/// per-span refinement stations for a hit cell.
pub(super) fn interior_knots_in(knots: &[f64], a: f64, b: f64) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    for &k in knots {
        if k > a + 1.0e-12
            && k < b - 1.0e-12
            && out.last().is_none_or(|&last| (k - last).abs() > 1.0e-12)
        {
            out.push(k);
        }
    }
    out
}
