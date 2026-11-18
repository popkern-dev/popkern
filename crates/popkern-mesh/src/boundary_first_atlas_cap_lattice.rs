//! A cap's lattice in its pole projection (chart-service.md §6a, step 5):
//! the cell is the largest whose segments along the projection's lines,
//! inside the cap, keep their sagitta (and, with an angular bound, their
//! tangent turn) within the chord height, the criterion the band's
//! pipeline reads along its iso-lines, so the density does not jump
//! across the row. The cap's box is not the cap: a square round a 45°
//! sphere cap reaches the equator at its corners, and a deviation read
//! along its diagonal sized the cap for a curvature it does not hold
//! (twice the band's triangles, PK-1271).
//!
//! The lattice is uniform in arc length from the singular point, not in
//! the projection: a tangent-plane projection compresses radially toward
//! the row (by the cosine of the normal's turn, 0.71 at a 45° row), so a
//! lattice uniform in it, sized at the rim, was twice as dense as needed
//! at the pole. Stations are laid in warped coordinates `w` whose radius
//! is the arc length along the ray (measured on the surface, averaged over
//! azimuths) and mapped into the projection along their ray.

use popkern_math::PortableMath;

use popkern_core::geometry::Point3;
use popkern_core::geometry::chart::atlas::cap_chart;
use popkern_core::geometry::chart::{Axis, ChartDescriptor, ChartedSurface};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::super::Cut;
use crate::TessellationOptions;

/// Radial samples of the warp per cap.
const RADIAL: usize = 96;

/// Azimuths the radial arc length is averaged over.
const AZIMUTHS: usize = 8;

/// The radial warp of a cap's projection: projection radius against arc
/// length from the singular point, both increasing from zero.
#[derive(Debug, Clone)]
pub(super) struct RadialWarp {
    rho: Vec<f64>,
    arc: Vec<f64>,
}

impl RadialWarp {
    /// Measured on `charted` out to `reach` in the projection.
    fn of(charted: &ChartedSurface<'_>, reach: f64) -> Self {
        let rho: Vec<f64> = (0..=RADIAL)
            .map(|k| reach * k as f64 / RADIAL as f64)
            .collect();
        let mut arc = vec![0.0; RADIAL + 1];
        let rays: Vec<(f64, f64)> = (0..AZIMUTHS)
            .map(|a| {
                let phi = std::f64::consts::TAU * (a as f64 + 0.5) / AZIMUTHS as f64;
                (phi.pcos(), phi.psin())
            })
            .collect();
        let mut prev: Vec<Point3> = rays.iter().map(|_| charted.point_at(0.0, 0.0)).collect();
        for k in 1..=RADIAL {
            let mut step = 0.0;
            for (j, &(c, s)) in rays.iter().enumerate() {
                let p = charted.point_at(rho[k] * c, rho[k] * s);
                step += p.distance_to(&prev[j]);
                prev[j] = p;
            }
            arc[k] = arc[k - 1] + step / AZIMUTHS as f64;
        }
        // A degenerate profile (no length) is the identity.
        if !arc[RADIAL].is_finite() || arc[RADIAL] <= 0.0 || arc.windows(2).any(|w| w[1] <= w[0]) {
            arc.clone_from(&rho);
        }
        Self { rho, arc }
    }

    /// The arc length at the projection radius `r` (linear past the
    /// measured reach).
    fn arc_at(&self, r: f64) -> f64 {
        interpolate(&self.rho, &self.arc, r)
    }

    /// The projection point of the warped station `w`.
    pub(super) fn xy(&self, w: (f64, f64)) -> (f64, f64) {
        let s = w.0.phypot(w.1);
        if s == 0.0 {
            return (0.0, 0.0);
        }
        let r = interpolate(&self.arc, &self.rho, s);
        (w.0 * r / s, w.1 * r / s)
    }

    /// The largest ratio of arc length to projection radius: the scale
    /// that takes the projection's box over the warped stations' box.
    fn stretch(&self) -> f64 {
        self.rho
            .iter()
            .zip(&self.arc)
            .skip(1)
            .map(|(r, s)| s / r)
            .fold(1.0_f64, f64::max)
    }

    /// The projection's step per unit arc length at the reach, radially
    /// and round the ray, the smaller: the cell's size in the projection
    /// at the rim, per unit cell.
    fn rim_scale(&self) -> f64 {
        let k = self.rho.len() - 1;
        let radial = (self.rho[k] - self.rho[k - 1]) / (self.arc[k] - self.arc[k - 1]);
        let round = self.rho[k] / self.arc[k];
        radial.min(round).min(1.0)
    }
}

/// `ys` at `x` on the increasing table `xs`, linear between samples and
/// past the last.
fn interpolate(xs: &[f64], ys: &[f64], x: f64) -> f64 {
    let n = xs.len();
    let i = xs.partition_point(|&t| t <= x).clamp(1, n - 1);
    let (x0, x1, y0, y1) = (xs[i - 1], xs[i], ys[i - 1], ys[i]);
    y0 + (y1 - y0) * (x - x0) / (x1 - x0)
}

/// A cap's lattice: `n` stations a side of `cell` (arc length) from the
/// warped box's corner `lo`, mapped into the projection by `warp`.
#[derive(Debug, Clone)]
pub(super) struct CapLattice {
    pub(super) n: usize,
    pub(super) cell: f64,
    lo: (f64, f64),
    warp: RadialWarp,
    /// The cell's size in the projection at the rim, where the projection
    /// is most compressed: the ring index's clearance scale.
    pub(super) rim_cell: f64,
}

impl CapLattice {
    /// Station `(i, j)` in the projection, offset half a cell from the
    /// box's corners.
    pub(super) fn station(&self, i: usize, j: usize) -> (f64, f64) {
        self.warp.xy((
            self.lo.0 + self.cell * (i as f64 + 0.5),
            self.lo.1 + self.cell * (j as f64 + 0.5),
        ))
    }
}

/// The share of the chord height a lattice segment may sag: a facet whose
/// edges each sag `s` bulges `4s/3` at its centre on a sphere (an
/// equilateral one), so sized to the full chord a corner ball's cap
/// sagged 5.3e-3 under a 5e-3 chord, its centroids within it.
const FACET_SHARE: f64 = 0.75;

/// Cross stations per direction at which the lines are read, as
/// fractions of the box's extent off its centre.
const CROSS: [f64; 5] = [-0.4, -0.2, 0.0, 0.2, 0.4];

/// The worst segment error at `n` cells over the warped box (sagitta, or
/// the turn scaled onto the chord height), read within `reach` (arc
/// length) of the singular point along lines of both axes, a cell's
/// sides, and both diagonals, a cell's diagonal steps: the lattice's
/// triangles span a cell's diagonal, whose sagitta on a sphere is twice a
/// side's, and a triangle's centroid sits well inside its deviation (a
/// 45° cap sized on its sides alone bowed 5.3e-3 under a 5e-3 chord at
/// its diagonals, the centroid test blind).
fn worst_at(
    (charted, warp): (&ChartedSurface<'_>, &RadialWarp),
    (lo, hi): ((f64, f64), (f64, f64)),
    reach: f64,
    n: usize,
    options: &TessellationOptions,
) -> f64 {
    let inside = |p: (f64, f64)| p.0.phypot(p.1) <= reach;
    let extent = (hi.0 - lo.0).max(hi.1 - lo.1);
    let h = extent / n as f64;
    let centre = (0.5 * (lo.0 + hi.0), 0.5 * (lo.1 + hi.1));
    let max_turn = options.angular_bound();
    let point = |w: (f64, f64)| {
        let p = warp.xy(w);
        charted.point_at(p.0, p.1)
    };
    let mut worst = 0.0_f64;
    // Each direction's step (a cell's side or diagonal) and the unit
    // offset across it.
    let directions = [
        ((h, 0.0), (0.0, 1.0)),
        ((0.0, h), (1.0, 0.0)),
        ((h, h), (-1.0, 1.0)),
        ((h, -h), (1.0, 1.0)),
    ];
    let steps = n as i64;
    for (step, across) in directions {
        for c in CROSS {
            let base = (
                centre.0 + across.0 * c * extent,
                centre.1 + across.1 * c * extent,
            );
            let at = |t: f64| (base.0 + step.0 * t, base.1 + step.1 * t);
            for k in -steps..steps {
                let (a, m, b) = (at(k as f64), at(k as f64 + 0.5), at(k as f64 + 1.0));
                if !(inside(a) && inside(b)) {
                    continue;
                }
                let (pa, pm, pb) = (point(a), point(m), point(b));
                let chord_mid = Point3::new(
                    0.5 * (pa.x + pb.x),
                    0.5 * (pa.y + pb.y),
                    0.5 * (pa.z + pb.z),
                );
                worst = worst.max(pm.distance_to(&chord_mid));
                if let Some(max_turn) = max_turn {
                    let turn = crate::facet_angle::curve_turn(pa, pm, pb);
                    worst = worst.max(options.chord_height * turn / max_turn);
                }
            }
        }
    }
    worst
}

/// A cap's lattice over the projection box `(lo, hi)` of a cap reaching
/// `reach` from the singular point. The count grows by a quarter (not
/// doubling, which overshoots the bound by up to two per direction) until
/// the segments inside the cap meet the chord height, then while a cell's
/// normal spread at the centre exceeds the angular bound.
pub(super) fn cap_lattice(
    charted: &ChartedSurface<'_>,
    (lo, hi): ((f64, f64), (f64, f64)),
    reach: f64,
    options: &TessellationOptions,
) -> CapLattice {
    let warp = RadialWarp::of(charted, reach);
    let k = warp.stretch();
    let (lo, hi) = ((lo.0 * k, lo.1 * k), (hi.0 * k, hi.1 * k));
    let arc_reach = warp.arc_at(reach);
    let extent = (hi.0 - lo.0).max(hi.1 - lo.1);
    let max_n = options.max_subdivisions.max(1);
    let grow = |n: usize| (n + (n / 4).max(1)).min(max_n);
    let mut n = options.min_subdivisions.max(1);
    while n < max_n
        && worst_at((charted, &warp), (lo, hi), arc_reach, n, options)
            > FACET_SHARE * options.chord_height
    {
        n = grow(n);
    }
    if let Some(angle) = options.angular_bound() {
        let c = (0.5 * (lo.0 + hi.0), 0.5 * (lo.1 + hi.1));
        let normal = |w: (f64, f64)| {
            let p = warp.xy(w);
            charted.normal_at(p.0, p.1)
        };
        let n0 = normal(c);
        let spread = |n: usize| {
            let h = 0.5 * extent / n as f64;
            [(-h, -h), (h, -h), (-h, h), (h, h)]
                .iter()
                .map(|&(dx, dy)| {
                    normal((c.0 + dx, c.1 + dy))
                        .dot(&n0)
                        .clamp(-1.0, 1.0)
                        .pacos()
                })
                .fold(0.0_f64, f64::max)
        };
        while spread(n) > angle && n < max_n {
            n = grow(n);
        }
    }
    let n = n.min(max_n);
    let cell = extent / n as f64;
    let rim_cell = cell * warp.rim_scale();
    CapLattice {
        n,
        cell,
        lo,
        warp,
        rim_cell,
    }
}

/// The cell (arc length) a cap about `cut` will be meshed at, estimated
/// before its rings are known from its row's disc (assembly samples the
/// row to it in space, so the cap's boundary is as fine as its lattice).
pub(in super::super) fn estimated_cap_cell(
    descriptor: &ChartDescriptor,
    surface: &SurfaceGeometry,
    cut: &Cut,
    options: &TessellationOptions,
) -> Option<f64> {
    let chart = cap_chart(descriptor, surface, &cut.cap)?;
    let rim = popkern_core::geometry::chart::ops::evaluate(descriptor, surface, cut.uv(cut.start));
    let (x, y) = chart.local_at(cut.uv(cut.start), &rim);
    let r = x.phypot(y);
    if !(r.is_finite() && r > 0.0) {
        return None;
    }
    let charted = ChartedSurface::new(surface, chart, ((-r, r), (-r, r)));
    Some(cap_lattice(&charted, ((-r, -r), (r, r)), r, options).cell)
}

/// The band's lattice spacing along each cut's row, in space, estimated
/// before the band exists from the face's box (the face's rings `built`,
/// the along span at most a period, extended across to the singular rows
/// of the caps `reached`): rows are sampled to the finer of it and the
/// cap's cell, since a row segment has no edge to split and the band's
/// lattice encroaching a long one leaves it uncovered (a rod's wall
/// passing a sphere's poles: once the cap's lattice coarsened to its
/// own sagitta, its row segments outran the band's cells and the band
/// declined every round).
pub(in super::super) fn band_row_cells(
    (descriptor, surface): (&ChartDescriptor, &SurfaceGeometry),
    (cuts, reached): (&[Cut], &[bool]),
    built: &[super::super::BoundaryRing],
    options: &TessellationOptions,
) -> Vec<Option<f64>> {
    let (mut lo, mut hi) = (
        (f64::INFINITY, f64::INFINITY),
        (f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for &(u, v) in built.iter().flat_map(|r| r.uv.iter()) {
        lo = (lo.0.min(u), lo.1.min(v));
        hi = (hi.0.max(u), hi.1.max(v));
    }
    let mut extend = |axis: Axis, at: f64| match axis {
        Axis::V => (lo.1, hi.1) = (lo.1.min(at), hi.1.max(at)),
        Axis::U => (lo.0, hi.0) = (lo.0.min(at), hi.0.max(at)),
    };
    for (cut, _) in cuts.iter().zip(reached).filter(|(_, r)| **r) {
        extend(cut.cap.singular.axis, cut.cap.singular.at);
    }
    let Some(cut) = cuts.first() else {
        return Vec::new();
    };
    // The along span at most a period (consecutively unwrapped rings run
    // past it).
    let along_span = |a: f64, b: f64| match cut.period {
        Some(p) => (b - a).min(p),
        None => b - a,
    };
    let (box_u, box_v) = match cut.cap.singular.axis {
        Axis::V => ((lo.0, lo.0 + along_span(lo.0, hi.0)), (lo.1, hi.1)),
        Axis::U => ((lo.0, hi.0), (lo.1, lo.1 + along_span(lo.1, hi.1))),
    };
    if !(box_u.0 < box_u.1 && box_v.0 < box_v.1) {
        return vec![None; cuts.len()];
    }
    let (_, n_u, n_v) = super::super::super::boundary_first_lattice::face_stations(
        (surface, descriptor),
        (box_u, box_v),
        None,
        options,
    );
    cuts.iter()
        .map(|cut| {
            let (span, n) = match cut.cap.singular.axis {
                Axis::V => (box_u.1 - box_u.0, n_u),
                Axis::U => (box_v.1 - box_v.0, n_v),
            };
            let at = cut.uv(cut.start);
            let (pu, pv) = surface.first_partials(at.0, at.1);
            let speed = match cut.cap.singular.axis {
                Axis::V => pu.length(),
                Axis::U => pv.length(),
            };
            let cell = speed * span / n as f64;
            (cell.is_finite() && cell > 0.0).then_some(cell)
        })
        .collect()
}
