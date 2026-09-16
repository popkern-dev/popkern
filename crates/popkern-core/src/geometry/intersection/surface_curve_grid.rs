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

/// Build one cell whose AABB CONTAINS its surface patch: the box of a
/// 3×3 station lattice, inflated by the lattice's own measured
/// sagitta. A corner-only box under-covers a curved cell by the
/// sagitta — a 48-cell grid over a 6.5-turn helical thread wall left
/// each cell's bulge 0.27 world units outside its box, every
/// classification ray through the coil was AABB-culled before Newton
/// ran, interior points read Outside, and the fragment keep/discard
/// verdicts scrambled per sample (PK-460). The earlier corner+centre
/// box inflated only by the CENTRE's deviation outside the corner box
/// — zero whenever the centre lies inside it — while a patch whose
/// extremum sits between centre and corner still bulges out by up to
/// the half-cell sag (measured, PK-460 tl=8: a 74°-azimuth helical
/// band cell boxed x ≤ 2.95 with its true bulge at 3.0625; the
/// classify ray's ±5% window pad had rephased the grid into exactly
/// that blind spot, the exit crossing was culled before Newton ran,
/// and the misclassified groove strip mis-wired a one-pitch seam
/// slit). The 3×3 lattice's row/column sags measure the local
/// curvature scale directly; for curvature-dominated (quadratic)
/// deviation the bulge beyond the 9-station box is bounded by the
/// geometric tail S·(1/4 + 1/16 + …) = S/3, and the box inflates by
/// exactly that bound — full-S inflation was measured too blunt (the
/// world-unit sag of coarse multi-turn cells fattened every box, the
/// cull lost its selectivity, and the LH two-start and M5 countersunk
/// serves regressed to new loop-gap declines). A flat cell measures
/// zero sag and keeps byte-identical boxes.
pub(super) fn padded_cell(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> SurfaceCell {
    let us = [u_range.0, 0.5 * (u_range.0 + u_range.1), u_range.1];
    let vs = [v_range.0, 0.5 * (v_range.0 + v_range.1), v_range.1];
    let mut stations = [[Point3::new(0.0, 0.0, 0.0); 3]; 3];
    for (i, &u) in us.iter().enumerate() {
        for (j, &v) in vs.iter().enumerate() {
            stations[i][j] = surface.point_at(u, v);
        }
    }
    let mut bmin = stations[0][0];
    let mut bmax = stations[0][0];
    for row in &stations {
        for p in row {
            bmin = Point3::new(bmin.x.min(p.x), bmin.y.min(p.y), bmin.z.min(p.z));
            bmax = Point3::new(bmax.x.max(p.x), bmax.y.max(p.y), bmax.z.max(p.z));
        }
    }
    // Measured sagitta: each 3-station row/column's middle vs its end
    // chord midpoint. No extra evaluations — the lattice is its own
    // curvature probe.
    let sag_of = |a: &Point3, m: &Point3, b: &Point3| -> f64 {
        let cx = 0.5 * (a.x + b.x);
        let cy = 0.5 * (a.y + b.y);
        let cz = 0.5 * (a.z + b.z);
        ((m.x - cx).powi(2) + (m.y - cy).powi(2) + (m.z - cz).powi(2)).sqrt()
    };
    let mut sag = 0.0_f64;
    for (i, row) in stations.iter().enumerate() {
        sag = sag.max(sag_of(&row[0], &row[1], &row[2]));
        sag = sag.max(sag_of(&stations[0][i], &stations[1][i], &stations[2][i]));
    }
    let pad = sag / 3.0;
    bmin = Point3::new(bmin.x - pad, bmin.y - pad, bmin.z - pad);
    bmax = Point3::new(bmax.x + pad, bmax.y + pad, bmax.z + pad);
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
