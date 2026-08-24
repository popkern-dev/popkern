//! Offset classification oracle (Rossignac–Requicha semantics made
//! executable, design doc §9 / PK-323).
//!
//! The mathematical offset of a solid `S` by distance `d` is
//! `{ p : dist(p, S) ≤ d }` (dilation), `{ p : dist(p, ∁S) ≥ d }`
//! written via the signed distance as `sdf(p) ≤ −d` (erosion), and a
//! thickened wall occupies the band `−t ≤ sdf(p) ≤ 0`. The oracle
//! samples a grid, evaluates the source body's mesh SDF at each point,
//! and compares the predicted membership against actual containment in
//! the served result mesh. It is kernel-independent — both sides are
//! judged purely from tessellations — so it catches
//! wrong-shape-served-clean defects that no structural validator can
//! see, and it validates trimming correctness.
//!
//! Every excluded or boundary-skipped sample is counted in the report:
//! silent truncation of the sample space would read as "covered
//! everything" when it didn't.

use popkern_core::geometry::Point3;

use crate::point_query::PointMeshIndex;

/// Which offset predicate the served result is expected to satisfy.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OffsetSemantics {
    /// Minkowski dilation by `distance` > 0: a point belongs to the
    /// result iff its signed distance to the source is ≤ `distance`.
    Grow {
        /// Offset distance (> 0).
        distance: f64,
    },
    /// Erosion by `distance` > 0: a point belongs to the result iff
    /// its signed distance to the source is ≤ `−distance`. Exact for
    /// convex bodies; for non-convex bodies erosion is still the
    /// correct predicate (sharp inward joins coincide with it).
    Shrink {
        /// Inset distance (> 0).
        distance: f64,
    },
    /// Wall material band `−inward ≤ sdf ≤ outward`: a shell wall is
    /// `{ inward: t, outward: 0 }` (away from its openings), a
    /// symmetric thicken of the same skin is
    /// `{ inward: t/2, outward: t/2 }`.
    Wall {
        /// Wall depth on the material side (≥ 0).
        inward: f64,
        /// Wall growth on the outside (≥ 0).
        outward: f64,
    },
}

/// Sampling parameters for the classification sweep.
#[derive(Debug, Clone)]
pub struct OffsetOracleOptions {
    /// Grid samples per axis over the inflated joint bounding box.
    pub resolution: usize,
    /// Half-width of the skip band around every membership boundary
    /// (both the predicted iso-surface and the served mesh itself):
    /// classification within tessellation chord error of a boundary is
    /// noise, not signal. Skipped samples are counted, never silent.
    pub boundary_band: f64,
    /// Cap on individually recorded disagreement points (the total
    /// count is always exact).
    pub max_recorded_disagreements: usize,
}

impl Default for OffsetOracleOptions {
    fn default() -> Self {
        Self {
            resolution: 40,
            boundary_band: 0.02,
            max_recorded_disagreements: 16,
        }
    }
}

/// One grid point whose predicted and actual membership disagree.
#[derive(Debug, Clone)]
pub struct OffsetDisagreement {
    /// Sample position.
    pub point: Point3,
    /// Signed distance from the sample to the source body.
    pub source_signed_distance: f64,
    /// Membership the offset predicate expected (`true` = inside the
    /// result); the served body reported the opposite.
    pub expected_inside: bool,
}

/// Outcome of a classification sweep.
#[derive(Debug, Clone)]
pub struct OffsetOracleReport {
    /// Grid samples visited.
    pub total_samples: usize,
    /// Samples where prediction and served membership agree.
    pub agreements: usize,
    /// Samples where they disagree (exact count).
    pub disagreement_count: usize,
    /// Samples skipped inside the boundary band.
    pub boundary_skipped: usize,
    /// Samples excluded by the caller's predicate (e.g. shell-opening
    /// neighbourhoods where the wall predicate does not apply).
    pub excluded: usize,
    /// Up to `max_recorded_disagreements` witness points.
    pub disagreements: Vec<OffsetDisagreement>,
}

impl OffsetOracleReport {
    /// Whether every classified sample agreed.
    pub fn is_clean(&self) -> bool {
        self.disagreement_count == 0
    }

    /// Fraction of visited samples that were actually classified
    /// (neither boundary-skipped nor excluded). A low value means the
    /// sweep barely tested anything — gate on this too.
    pub fn classified_fraction(&self) -> f64 {
        if self.total_samples == 0 {
            return 0.0;
        }
        (self.agreements + self.disagreement_count) as f64 / self.total_samples as f64
    }
}

/// Classify a served offset/shell result against the source body.
///
/// `exclude` removes samples where the declared semantics do not apply
/// (shell openings, sharp-join wedge rings); exclusions are counted in
/// the report so a test can bound them.
pub fn classify_offset_agreement(
    source: &PointMeshIndex,
    result: &PointMeshIndex,
    semantics: OffsetSemantics,
    options: &OffsetOracleOptions,
    exclude: impl Fn(Point3) -> bool,
) -> OffsetOracleReport {
    let mut report = OffsetOracleReport {
        total_samples: 0,
        agreements: 0,
        disagreement_count: 0,
        boundary_skipped: 0,
        excluded: 0,
        disagreements: Vec::new(),
    };
    let n = options.resolution.max(2);
    let band = options.boundary_band;
    let (lo, hi) = joint_bounds(source, result, band * 2.0 + growth_margin(semantics));
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let p = Point3::new(
                    lo.x + (hi.x - lo.x) * i as f64 / (n - 1) as f64,
                    lo.y + (hi.y - lo.y) * j as f64 / (n - 1) as f64,
                    lo.z + (hi.z - lo.z) * k as f64 / (n - 1) as f64,
                );
                report.total_samples += 1;
                if exclude(p) {
                    report.excluded += 1;
                    continue;
                }
                let sdf = source.signed_distance(p);
                let (expected, boundary_margin) = match semantics {
                    OffsetSemantics::Grow { distance } => (sdf <= distance, (sdf - distance).abs()),
                    OffsetSemantics::Shrink { distance } => {
                        (sdf <= -distance, (sdf + distance).abs())
                    }
                    OffsetSemantics::Wall { inward, outward } => (
                        (-inward..=outward).contains(&sdf),
                        (sdf - outward).abs().min((sdf + inward).abs()),
                    ),
                };
                if boundary_margin <= band || result.distance(p) <= band {
                    report.boundary_skipped += 1;
                    continue;
                }
                let actual = result.contains(p);
                if actual == expected {
                    report.agreements += 1;
                } else {
                    report.disagreement_count += 1;
                    if report.disagreements.len() < options.max_recorded_disagreements {
                        report.disagreements.push(OffsetDisagreement {
                            point: p,
                            source_signed_distance: sdf,
                            expected_inside: expected,
                        });
                    }
                }
            }
        }
    }
    report
}

/// Volume of the Minkowski dilation `{ p : sdf(p) ≤ distance }`
/// integrated on a grid with sub-cell fractional occupancy (the SDF
/// linearizes the boundary within a cell, giving O(h²) accuracy).
/// `distance` may be negative (erosion volume).
pub fn grid_offset_volume(source: &PointMeshIndex, distance: f64, resolution: usize) -> f64 {
    let n = resolution.max(2);
    let (blo, bhi) = source.bounds();
    let pad = distance.max(0.0) * 1.05 + 1e-9;
    let lo = Point3::new(blo.x - pad, blo.y - pad, blo.z - pad);
    let hi = Point3::new(bhi.x + pad, bhi.y + pad, bhi.z + pad);
    let hx = (hi.x - lo.x) / n as f64;
    let hy = (hi.y - lo.y) / n as f64;
    let hz = (hi.z - lo.z) / n as f64;
    let h = hx.max(hy).max(hz);
    let cell = hx * hy * hz;
    let mut volume = 0.0_f64;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let p = Point3::new(
                    lo.x + hx * (i as f64 + 0.5),
                    lo.y + hy * (j as f64 + 0.5),
                    lo.z + hz * (k as f64 + 0.5),
                );
                let sdf = source.signed_distance(p);
                let frac = (0.5 - (sdf - distance) / h).clamp(0.0, 1.0);
                volume += cell * frac;
            }
        }
    }
    volume
}

/// Steiner/tube formula for the volume of the Minkowski dilation:
/// `V + A·d + M·d² + (4π/3)·(χ/2)·d³`, where `M` is the total
/// mean-curvature integral `∫H dA + Σ (θₑ/2)·ℓₑ` over smooth faces and
/// exterior-angle edge terms, and `χ` is the boundary's Euler
/// characteristic (2 for sphere-like bodies, 0 for a torus — the
/// Gauss–Bonnet form of the cubic coefficient). Exact for convex
/// bodies at any `d ≥ 0`, and for non-convex bodies while `d` stays
/// below the reach (no self-intersection of the offset).
pub fn steiner_offset_volume(
    volume: f64,
    surface_area: f64,
    mean_curvature_integral: f64,
    euler_characteristic: f64,
    distance: f64,
) -> f64 {
    volume
        + surface_area * distance
        + mean_curvature_integral * distance * distance
        + 4.0 * std::f64::consts::PI / 3.0 * (euler_characteristic / 2.0) * distance.powi(3)
}

fn growth_margin(semantics: OffsetSemantics) -> f64 {
    match semantics {
        OffsetSemantics::Grow { distance } => distance * 0.25,
        OffsetSemantics::Shrink { .. } | OffsetSemantics::Wall { .. } => 0.0,
    }
}

fn joint_bounds(a: &PointMeshIndex, b: &PointMeshIndex, pad: f64) -> (Point3, Point3) {
    let (alo, ahi) = a.bounds();
    let (blo, bhi) = b.bounds();
    (
        Point3::new(
            alo.x.min(blo.x) - pad,
            alo.y.min(blo.y) - pad,
            alo.z.min(blo.z) - pad,
        ),
        Point3::new(
            ahi.x.max(bhi.x) + pad,
            ahi.y.max(bhi.y) + pad,
            ahi.z.max(bhi.z) + pad,
        ),
    )
}

#[cfg(test)]
#[path = "offset_oracle_tests.rs"]
mod tests;
