//! SAME-SURFACE patch exemption for the self-intersection validator
//! (PK-308): two faces on ONE analytic surface are subsets of a
//! single embedded surface — they can only meet where their trimmed
//! domains overlap modulo the chart lattice. When the coherent exact
//! boxes are disjoint (boundary touch allowed) the pair cannot
//! interpenetrate, and the SSI solver must not be consulted at all:
//! identical geometry reads as whole-surface coincidence regardless
//! of the boxes (measured: the raw torus primitive flagged its own
//! tube-opposite quarter-patches as SelfIntersection).
//! Domain-overlapping same-surface pairs (a genuinely folded shell)
//! still run the full check.

use crate::geometry::Point3;
use crate::geometry::intersection::coincident_analytic::are_surfaces_coincident_analytic;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

/// True when both faces lie on the same analytic surface AND their
/// coherent (corner-chained) uv boxes are disjoint modulo the chart
/// lattice — the pair is exempt from the self-intersection probe.
pub(super) fn same_surface_domains_disjoint(
    geom_i: &SurfaceGeometry,
    geom_j: &SurfaceGeometry,
    positions_i: &[Point3],
    positions_j: &[Point3],
    tolerance: &Tolerance,
) -> bool {
    if !are_surfaces_coincident_analytic(geom_i, geom_j, tolerance) {
        return false;
    }
    let (pu, pv) = geom_i.parameter_periodicity();
    // Coherent boxes by corner-chaining: the raw exact ranges wrap at
    // chart-branch corners (a quarter patch touching v = 2π projects
    // corners to v = 0 and its box balloons across the chart), so
    // each face's loop positions are unwrapped corner-to-corner
    // before boxing.
    let coherent_box = |positions: &[Point3]| -> Option<((f64, f64), (f64, f64))> {
        let mut prev: Option<(f64, f64)> = None;
        let (mut lo_u, mut hi_u, mut lo_v, mut hi_v) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
        for pnt in positions {
            let (mut u, mut v) = crate::geometry::uv_projection::project_point_to_uv(geom_i, pnt)?;
            if let Some((pu0, pv0)) = prev {
                if let Some(p) = pu {
                    u += ((pu0 - u) / p).round() * p;
                }
                if let Some(p) = pv {
                    v += ((pv0 - v) / p).round() * p;
                }
            }
            prev = Some((u, v));
            lo_u = lo_u.min(u);
            hi_u = hi_u.max(u);
            lo_v = lo_v.min(v);
            hi_v = hi_v.max(v);
        }
        Some(((lo_u, hi_u), (lo_v, hi_v)))
    };
    let (Some((bu_i, bv_i)), Some((bu_j, bv_j))) =
        (coherent_box(positions_i), coherent_box(positions_j))
    else {
        return false;
    };
    let overlap_1d = |a: (f64, f64), b: (f64, f64), period: Option<f64>| -> f64 {
        let (b0, b1) = match period {
            Some(p) if p.is_finite() && p > 0.0 => {
                let shift = (((a.0 + a.1) - (b.0 + b.1)) * 0.5 / p).round() * p;
                (b.0 + shift, b.1 + shift)
            }
            _ => b,
        };
        a.1.min(b1) - a.0.max(b0)
    };
    let margin = 1.0e-7 * (1.0 + (bu_i.1 - bu_i.0).abs().max((bv_i.1 - bv_i.0).abs()));
    let du = overlap_1d(bu_i, bu_j, pu);
    let dv = overlap_1d(bv_i, bv_j, pv);
    du <= margin || dv <= margin
}
