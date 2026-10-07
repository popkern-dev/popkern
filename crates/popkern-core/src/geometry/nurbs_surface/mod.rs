//! NURBS surface evaluation via tensor-product de Boor algorithm.

use crate::geometry::nurbs_ops::{Homogeneous, MAX_DE_BOOR_DEGREE, de_boor_1d_local, find_span};
use crate::geometry::surface::SurfaceEval;
use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

mod degree_and_split;
mod eval_derivs;
mod iso_curve;
mod knot_ops;
mod projection_seed_cache;
pub use projection_seed_cache::ProjectionSeedField;

/// A Non-Uniform Rational B-Spline surface.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NurbsSurface {
    /// Polynomial degree in the u direction.
    pub degree_u: usize,
    /// Polynomial degree in the v direction.
    pub degree_v: usize,
    /// Knot vector in the u direction.
    pub knots_u: Vec<f64>,
    /// Knot vector in the v direction.
    pub knots_v: Vec<f64>,
    /// Control point grid, stored u-major: the flat index of grid
    /// position `(i_u, i_v)` is `i_u * num_v + i_v` (v varies fastest).
    /// Length = num_u * num_v.
    pub control_points: Vec<Point3>,
    /// Number of control points in the u direction.
    pub num_u: usize,
    /// Number of control points in the v direction.
    pub num_v: usize,
    /// Weights, one per control point.
    pub weights: Vec<f64>,
}

impl NurbsSurface {
    /// Knot-span-midpoint + chart-border projection seeds for
    /// multi-sheet closest-point queries (PK-446): every span
    /// midpoint pair plus the chart borders — a border minimum (a rim
    /// point ON the v=1 line of a coiled wall) is farther from every
    /// span MIDPOINT than from the adjacent turn's surface, so
    /// midpoint-only ranking walked the wrong sheet. Returns `None`
    /// for few-span charts (≤5 stations per axis), where a generic
    /// multi-start already covers the basin structure.
    ///
    /// Subsampling is budgeted on the TOTAL pair count (48² = 2304),
    /// not per axis (PK-461): a 26-turn thread-tool wall has ~200
    /// v-spans against 1–3 u-spans, and the old 48-per-axis cap
    /// crushed the v-axis to half-turn seed gaps while using a
    /// fraction of the pair budget. With half-turn gaps the coarse
    /// distance ranking prefers the ADJACENT TURN's same-phase seed
    /// (one pitch away) over the right turn's phase-offset seed (a
    /// cross-coil chord), the Gauss–Newton refinement converges in
    /// the wrong turn's basin, and every consumer of the projection
    /// misreads an on-surface point as ~a pitch off the wall — the
    /// curve-surface validator declined a fully valid 25-turn thread
    /// serve on a phantom 0.747 deviation (the 24-vs-25-turn boundary
    /// was seed-phase stride resonance, not geometry).
    pub fn projection_seeds(&self) -> Option<Vec<(f64, f64)>> {
        let axis_mids = |knots: &[f64]| -> (Vec<f64>, Option<(f64, f64)>) {
            let mut uniq: Vec<f64> = Vec::new();
            for &k in knots {
                if uniq.last().is_none_or(|&last| k > last) {
                    uniq.push(k);
                }
            }
            let mids: Vec<f64> = uniq.windows(2).map(|w| 0.5 * (w[0] + w[1])).collect();
            let ends = match (uniq.first(), uniq.last()) {
                (Some(&lo), Some(&hi)) => Some((lo, hi)),
                _ => None,
            };
            (mids, ends)
        };
        let subsample = |mids: Vec<f64>, cap: usize| -> Vec<f64> {
            if mids.len() > cap {
                let step = mids.len() as f64 / cap as f64;
                (0..cap).map(|i| mids[(i as f64 * step) as usize]).collect()
            } else {
                mids
            }
        };
        const PAIR_BUDGET: usize = 48 * 48;
        let (mids_u, ends_u) = axis_mids(&self.knots_u);
        let (mids_v, ends_v) = axis_mids(&self.knots_v);
        // Give each axis at least the historical 48, then spend the
        // remaining pair budget on the denser axis.
        let cap_u = (PAIR_BUDGET / mids_v.len().clamp(1, 48)).max(48);
        let cap_v = (PAIR_BUDGET / mids_u.len().clamp(1, 48)).max(48);
        let mut us = subsample(mids_u, cap_u);
        let mut vs = subsample(mids_v, cap_v);
        if let Some((lo, hi)) = ends_u {
            us.insert(0, lo);
            us.push(hi);
        }
        if let Some((lo, hi)) = ends_v {
            vs.insert(0, lo);
            vs.push(hi);
        }
        if us.len() <= 5 && vs.len() <= 5 {
            return None;
        }
        let mut seeds = Vec::with_capacity(us.len() * vs.len());
        for &u in &us {
            for &v in &vs {
                seeds.push((u, v));
            }
        }
        Some(seeds)
    }

    /// Access control point at grid position (i_u, i_v) in the flat array.
    pub(super) fn cp(&self, i_u: usize, i_v: usize) -> &Point3 {
        &self.control_points[i_u * self.num_v + i_v]
    }

    /// Access weight at grid position (i_u, i_v) in the flat array.
    pub(super) fn wt(&self, i_u: usize, i_v: usize) -> f64 {
        self.weights[i_u * self.num_v + i_v]
    }

    /// The same surface with the parameter axes swapped: `S'(u, v) =
    /// S(v, u)`, exactly (degrees, knot vectors and the control grid
    /// transpose; no refit). Transposition reverses the chart
    /// orientation — `∂u' × ∂v' = −(∂u × ∂v)` — so a face that keeps
    /// its material side must flip its stored orientation flag with
    /// it. Introduced for PK-456: surfaces of revolution author their
    /// closed (angle) direction in `u` like every other closed chart
    /// in the kernel.
    pub fn transposed(&self) -> NurbsSurface {
        let mut control_points = Vec::with_capacity(self.control_points.len());
        let mut weights = Vec::with_capacity(self.weights.len());
        for j_v in 0..self.num_v {
            for i_u in 0..self.num_u {
                control_points.push(*self.cp(i_u, j_v));
                weights.push(self.wt(i_u, j_v));
            }
        }
        NurbsSurface {
            degree_u: self.degree_v,
            degree_v: self.degree_u,
            knots_u: self.knots_v.clone(),
            knots_v: self.knots_u.clone(),
            control_points,
            num_u: self.num_v,
            num_v: self.num_u,
            weights,
        }
    }

    /// Evaluate the surface at (u, v) using tensor-product de Boor in
    /// homogeneous coordinates.
    ///
    /// Span-local: a de Boor evaluation at (u, v) only ever reads the
    /// `(degree_u + 1) × (degree_v + 1)` control window around the two
    /// knot spans, so that is all this builds — one `find_span` per
    /// direction, `degree_u + 1` row recursions, one column recursion.
    /// Cost is O((p+1)·(q+1)) regardless of grid size, always
    /// allocation-free (the scratch is bounded by the degree cap, not
    /// by `num_u`/`num_v`), and bit-identical to evaluating the full
    /// grid because [`de_boor_1d_local`] runs the same lerps on the
    /// same values in the same order (PK-281; the earlier full-grid
    /// form evaluated every row and discarded all but `degree_u + 1`
    /// intermediates, which is what made `de_boor_1d` dominate the
    /// G2 cube profile at 74/100 samples).
    fn eval_homogeneous(&self, u: f64, v: f64) -> Homogeneous {
        let zero = Homogeneous {
            wx: 0.0,
            wy: 0.0,
            wz: 0.0,
            w: 0.0,
        };
        let span_u = find_span(self.degree_u, &self.knots_u, self.num_u, u);
        let span_v = find_span(self.degree_v, &self.knots_v, self.num_v, v);
        let mut intermediate = [zero; MAX_DE_BOOR_DEGREE];
        let mut row = [zero; MAX_DE_BOOR_DEGREE];
        for (j_u, intermediate_slot) in intermediate.iter_mut().enumerate().take(self.degree_u + 1)
        {
            let i_u = span_u - self.degree_u + j_u;
            for (j_v, row_slot) in row.iter_mut().enumerate().take(self.degree_v + 1) {
                let i_v = span_v - self.degree_v + j_v;
                *row_slot = Homogeneous::from_point_weight(self.cp(i_u, i_v), self.wt(i_u, i_v));
            }
            *intermediate_slot = de_boor_1d_local(
                self.degree_v,
                &self.knots_v,
                span_v,
                &row[..=self.degree_v],
                v,
            );
        }
        de_boor_1d_local(
            self.degree_u,
            &self.knots_u,
            span_u,
            &intermediate[..=self.degree_u],
            u,
        )
    }
}

impl SurfaceEval for NurbsSurface {
    fn point_at(&self, u: f64, v: f64) -> Point3 {
        self.eval_homogeneous(u, v).to_point3()
    }

    // First partials are exact rational derivatives via the fused
    // basis-contraction evaluation (PK-281) — previously 1e-7 central
    // differences, which cost two full evaluations each and carried
    // ~1e-9 truncation error into every Gauss–Newton projection.
    fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        self.point_and_first_partials(u, v).1
    }

    fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        self.point_and_first_partials(u, v).2
    }

    fn normal_at(&self, u: f64, v: f64) -> Vec3 {
        let (_, du, dv) = self.point_and_first_partials(u, v);
        du.cross(&dv)
    }

    fn parameter_range(&self) -> ((f64, f64), (f64, f64)) {
        let p = self.degree_u;
        let n_u = self.num_u - 1;
        let q = self.degree_v;
        let n_v = self.num_v - 1;
        (
            (self.knots_u[p], self.knots_u[n_u + 1]),
            (self.knots_v[q], self.knots_v[n_v + 1]),
        )
    }

    // Second partials are exact rational derivatives via the fused
    // two-derivative basis contraction (PK-285) — previously 1e-5
    // central differences of point_at (3–5 evaluations each), whose
    // truncation error fed every curvature and G2 metric on NURBS.
    fn second_partial_uu(&self, u: f64, v: f64) -> Vec3 {
        self.point_and_two_partials(u, v).3
    }

    fn second_partial_uv(&self, u: f64, v: f64) -> Vec3 {
        self.point_and_two_partials(u, v).4
    }

    fn second_partial_vv(&self, u: f64, v: f64) -> Vec3 {
        self.point_and_two_partials(u, v).5
    }

    fn closest_point(&self, point: &Point3) -> Option<(f64, f64)> {
        // KNOT-SPAN-AWARE seeding: the generic 5×5 start grid cannot
        // resolve a multi-sheet chart — on a 6.5-turn helical tube
        // wall every start landed in the wrong turn's basin and the
        // projection for a point ON the rim answered a border
        // parameter 6e-2 of surface away (PK-446, the phantom
        // watertightness gap). Seed every knot-span-midpoint pair
        // (subsampled to 48/axis), coarse-select by distance, and
        // Newton-refine only the best few.
        use crate::geometry::surface::{closest_point_from_seeds, closest_point_gauss_newton};
        match self.projection_seeds() {
            // Few spans — the generic multi-start already covers the
            // basin structure; keep the tuned path.
            None => closest_point_gauss_newton(self, point),
            Some(seeds) => closest_point_from_seeds(self, point, &seeds, (false, false), 6),
        }
    }

    fn bounding_box(
        &self,
        u_range: (f64, f64),
        v_range: (f64, f64),
    ) -> crate::geometry::intersection::surface_subdivision::Aabb {
        use crate::geometry::surface::bounding_box_by_sampling;
        bounding_box_by_sampling(self, u_range, v_range)
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
#[path = "projection_seeds_tests.rs"]
mod projection_seeds_tests;
#[cfg(test)]
#[path = "transpose_tests.rs"]
mod transpose_tests;

#[cfg(test)]
mod eval_span_local_tests;

#[cfg(test)]
mod iso_curve_tests;
