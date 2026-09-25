//! Ring-level pairwise crossing detection for K-P bi-nonic cap rings.
//!
//! The per-patch fold gate ([`GregoryPatch::self_intersects_with_tolerance`](
//! super::GregoryPatch::self_intersects_with_tolerance)) cannot see two
//! *different* ring patches crossing each other: an optimizer trajectory can
//! keep every patch individually fold-free while driving adjacent patches
//! through one another (the PK-283 class, first observed on
//! `pyramid4_uniform` under a deeper LM trajectory — every per-patch accept
//! gate passed, `validate_all` flagged the pair downstream). This module is
//! the ring-level, crossing-honest complement: it scans every unordered
//! patch pair for interior contact beyond the ring's legitimate shared
//! contact loci.
//!
//! # Ring convention
//!
//! The same convention as the `build_binonic_ring*` builders and the fillet
//! promotion ladder's spoke table: spoke `k` glues patch `k`'s `u = 1` side
//! to patch `(k + 1) % n`'s `v = 0` side (reversed), and every patch's
//! `(u, v) = (1, 0)` corner is the shared central apex.
//!
//! # Legitimate contact vs crossing
//!
//! Contact **along the shared curve is legitimate; interior overlap is not**
//! (the PK-252 extended-Boundary lesson). The gate's strictness follows the
//! pair's topological relation, because that is what decides which
//! downstream oracle exists:
//!
//! - **Apex-only pairs** (no shared spoke): an intersection sample is
//!   exempt only within the exemption band (`10 × tolerance` in 3D) of the
//!   apex point. Deliberately strict: pairs without a shared edge are
//!   exactly the class `validate_all`'s non-adjacent face check flags
//!   downstream (the PK-275 escape), so the gate mirrors that verdict at
//!   acceptance time.
//! - **Spoke-sharing pairs** (ring-consecutive; both directions are
//!   checked mod `n`, so `n = 2` and `n = 3` rings — where every pair is
//!   consecutive — get every shared spoke): a sample is exempt within the
//!   scale-relative macro band ([`MACRO_DEPTH_FACTOR`] × the pair's joint
//!   AABB diagonal) of a shared spoke or the apex. `validate_all` never
//!   examines adjacent faces, and measured serving rings legitimately
//!   carry spoke-hover branches well above solver noise
//!   (pyramid3_uniform's validate_all-clean, OCC-green ring: 2.38e-5 on
//!   one spoke and 1.886e-3 — reaching `v = 0.076` on the VZero side — on
//!   another), so only macro-scale penetration is defect evidence here.
//!
//! Any sample beyond these exemptions is an interior crossing and produces
//! a [`BinonicRingCrossing`] report. A `CoincidentRegion` result between
//! two distinct ring patches is always a crossing (patches may share a
//! boundary curve — measure zero — never an area).
//!
//! Known limitation (fixture-less as of 2026-08-16, see the decline-class
//! registry): a spoke-neighbour crossing below the macro band is exempt.
//! No fixture reaches a defective instance of that class; `validate_all`
//! cannot see it either (adjacent-face exemption), so a tighter gate
//! there would be unmeasurable containment — the measured serving
//! envelope (1.886e-3) already sits within ~1.5× of the macro cap.
//!
//! # Cost structure
//!
//! Full NURBS×NURBS surface intersection per pair is far too expensive for
//! a routine acceptance gate (measured: the A.6 uniform fixture set went
//! 2.47 s → 395 s when every pair ran the solver). Each pair therefore
//! first attempts a recursive **separation certificate**: exact de
//! Casteljau subdivision of both control grids with an interval bound —
//! a sub-patch box pair is *clean* when the boxes don't overlap (padded
//! by the strict band) or when every potential close approach in their
//! overlap provably lies within the pair's exemption depth of the shared
//! contact locus (`center distance + half diagonal` bound). Only pairs
//! the certificate cannot clear (real crossings, or pathological
//! geometry that exhausts the subdivision budget) fall through to the
//! full surface-surface solver. Clean rings certify in microseconds.
//!
//! # Failure containment
//!
//! The surface-surface solver failing (`Err`) or returning a non-finite
//! curve is treated as "no evidence of crossing" (fail-open), consistent
//! with the per-patch subdivision check: this is an *accept gate*, so a
//! spurious decline would silently degrade serving quality, while a missed
//! crossing is still caught by `validate_all` downstream.

use super::RectangularGregoryBiNonic;
use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::intersect_surfaces;
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;
use ring_pair_crossing_certificate::{distance_to_spoke, pair_certified_clean};

#[path = "ring_pair_crossing_bounded.rs"]
mod ring_pair_crossing_bounded;
#[path = "ring_pair_crossing_certificate.rs"]
mod ring_pair_crossing_certificate;
#[path = "ring_pair_crossing_setup.rs"]
mod ring_pair_crossing_setup;
#[path = "ring_pair_crossing_witness.rs"]
mod ring_pair_crossing_witness;
pub use ring_pair_crossing_bounded::{BoundedRingCrossing, binonic_ring_crossing_bounded};
use ring_pair_crossing_setup::{PairSetup, pair_setup, patch_aabbs, spoke_polylines};

/// How two patches of a bi-nonic cap ring relate topologically.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RingPairRelation {
    /// The pair shares at least one spoke curve (ring-consecutive).
    SpokeNeighbors,
    /// The pair shares only the central apex point.
    ApexOnly,
}

/// A detected interior crossing between two patches of a bi-nonic ring.
///
/// Produced by [`binonic_ring_crossing`]; carries the worst offending
/// sample of the first crossing pair in deterministic scan order, so
/// accept-gate declines are reproducible and localizable.
#[derive(Debug, Clone)]
pub struct BinonicRingCrossing {
    /// Ring index of the first patch of the crossing pair (scan order).
    pub patch_a: usize,
    /// Ring index of the second patch of the crossing pair.
    pub patch_b: usize,
    /// Topological relation of the pair within the ring.
    pub relation: RingPairRelation,
    /// 3D location of the worst offending intersection sample.
    pub point: Point3,
    /// UV of the sample on patch `patch_a`.
    pub uv_a: (f64, f64),
    /// UV of the sample on patch `patch_b`.
    pub uv_b: (f64, f64),
    /// Distance from the sample to the nearest legitimate contact locus
    /// (shared spoke curve or apex). `f64::INFINITY` for coincident
    /// regions, which have no legitimate extent.
    pub separation: f64,
    /// The exemption depth the sample exceeded: the `10 × tolerance`
    /// band for strict (apex-only) pairs, the scale-relative macro cap
    /// for spoke-neighbour pairs and the feasible-start screen.
    pub exemption_band: f64,
    /// Whether the solver reported the pair coincident over a region
    /// (always a crossing between distinct ring patches).
    pub coincident_region: bool,
}

impl BinonicRingCrossing {
    /// Stable diagnostic string for decline traces and structured facts.
    pub fn diagnostic_label(&self) -> String {
        format!(
            "ring_crossing pair=({},{}) relation={:?} coincident_region={} separation={:.3e} band={:.3e} point=({:.6},{:.6},{:.6}) uv_a=({:.4},{:.4}) uv_b=({:.4},{:.4})",
            self.patch_a,
            self.patch_b,
            self.relation,
            self.coincident_region,
            self.separation,
            self.exemption_band,
            self.point.x,
            self.point.y,
            self.point.z,
            self.uv_a.0,
            self.uv_a.1,
            self.uv_b.0,
            self.uv_b.1
        )
    }
}

/// Samples per spoke polyline used for distance-to-spoke estimates.
const SPOKE_POLYLINE_SAMPLES: usize = 129;
/// Samples taken along each reported intersection curve.
const CURVE_SAMPLES: usize = 17;
/// Exemption band multiplier over the caller's linear tolerance.
const EXEMPTION_BAND_FACTOR: f64 = 10.0;
/// Spoke-neighbour exemption depth as a fraction of the pair's joint
/// AABB diagonal. Measured anchors (2026-08-16): legitimate spoke-hover
/// on the serving, validate_all-clean pyramid3_uniform ring reaches
/// 1.886e-3 absolute on a 0.463 joint diagonal — 4.1e-3 relative,
/// 2.4× below this cap — while a defective neighbour crossing (a sheet
/// sweeping through the neighbouring sector, the DOF class the
/// spoke-row optimizer moves) measures ≥ 5e-2 relative on the
/// synthetic reproducer. Re-measure the serving envelope (run the
/// fillet corpus with `POPKERN_TRACE=kp-cap` and read the decline
/// labels) before tightening.
pub const MACRO_DEPTH_FACTOR: f64 = 1.0e-2;

/// Scan a bi-nonic cap ring for interior crossings between patch pairs.
///
/// `tolerance` is the linear tolerance of the acceptance gate (the same
/// value the caller uses for the per-patch fold gate); the legitimate
/// shared-contact exemption band is `10 × tolerance`. See the module docs
/// for the ring convention, exemption semantics, and failure containment.
///
/// Returns the first crossing pair in deterministic scan order
/// (lexicographic over `(patch_a, patch_b)`, `patch_a < patch_b`) with the
/// pair's worst sample, or `None` when the ring is crossing-free. Rings
/// with fewer than two patches are trivially crossing-free.
///
/// Cost: one exact NURBS conversion per patch plus one surface-surface
/// intersection per AABB-overlapping pair, run once per acceptance — not
/// per optimizer step.
pub fn binonic_ring_crossing(
    ring: &[RectangularGregoryBiNonic],
    tolerance: f64,
) -> Option<BinonicRingCrossing> {
    ring_crossing_scan(ring, tolerance, false)
}

/// Coarse feasible-start variant of [`binonic_ring_crossing`] for gating
/// *raw* (pre-optimization) rings: only macro-scale crossings flag — any
/// sample within the scale-relative macro band
/// ([`MACRO_DEPTH_FACTOR`] × the pair's joint AABB diagonal) is exempt
/// regardless of UV position. Measured basis: legitimate pre-LM rings
/// carry spoke-hover contact an order deeper and wider than converged
/// rings (pyramid3_uniform's initial guide ring: 1.886e-3 at
/// `v_b = 0.076`, later polished to 2.38e-5 by the optimizer and served
/// validate_all-clean), so the full accept-gate semantics would decline
/// starts the optimizer routinely repairs.
pub fn binonic_ring_feasible_start_crossing(
    ring: &[RectangularGregoryBiNonic],
    tolerance: f64,
) -> Option<BinonicRingCrossing> {
    ring_crossing_scan(ring, tolerance, true)
}

/// A cheap, sound-positive crossing witness (PK-547): a crossing of the
/// ring proven at a refined intersection point beyond the pair's
/// exemption depth, found from a coarse tessellation in milliseconds.
/// A crossing the full gate needs the surface-surface solver for (27 s
/// on the pyramid3_uniform welded-retry ring) is found here in
/// milliseconds. `Some` is a crossing the full gate also declines.
/// `None` proves nothing, and [`binonic_ring_crossing`] stays the
/// acceptance authority. Optimizers use it as a per-step feasibility
/// screen.
pub fn binonic_ring_crossing_witness(
    ring: &[RectangularGregoryBiNonic],
    tolerance: f64,
) -> Option<BinonicRingCrossing> {
    ring_pair_crossing_witness::ring_crossing_witness(ring, tolerance)
}

fn ring_crossing_scan(
    ring: &[RectangularGregoryBiNonic],
    tolerance: f64,
    macro_only: bool,
) -> Option<BinonicRingCrossing> {
    let n = ring.len();
    if n < 2 {
        return None;
    }
    let band = EXEMPTION_BAND_FACTOR * tolerance;

    let aabbs = patch_aabbs(ring);

    let nurbs: Vec<SurfaceGeometry> = ring
        .iter()
        .map(|patch| SurfaceGeometry::Nurbs(patch.to_nurbs_exact()))
        .collect();

    let spokes = spoke_polylines(ring);

    let ssi_tol = Tolerance::default().with_linear(tolerance.max(1e-6));

    for a in 0..n {
        for b in (a + 1)..n {
            if let (Some(bb_a), Some(bb_b)) = (&aabbs[a], &aabbs[b])
                && !bb_a.overlaps(bb_b, band)
            {
                continue;
            }
            let PairSetup {
                relation,
                shared_spokes,
                apex,
                macro_depth,
                exempt_depth,
            } = pair_setup(ring, &aabbs, &spokes, (a, b), band, macro_only);
            // Fast path: prove the pair clean by subdivision before
            // paying for the surface-surface solver (see module docs).
            if pair_certified_clean(
                &ring[a].grid,
                &ring[b].grid,
                &shared_spokes,
                apex,
                band,
                exempt_depth,
            ) {
                continue;
            }

            if let Some(crossing) = pair_crossing(
                &nurbs[a],
                &nurbs[b],
                (a, b),
                relation,
                &shared_spokes,
                apex,
                band,
                macro_depth,
                macro_only,
                &ssi_tol,
            ) {
                return Some(crossing);
            }
        }
    }
    None
}

/// One shared spoke of a ring-consecutive pair.
pub(super) struct SharedSpoke<'r> {
    /// Precomputed samples of the spoke owner's `u = 1` boundary iso.
    polyline: &'r [Point3],
    /// The patch whose `u = 1` side is this spoke (for refinement).
    patch: &'r RectangularGregoryBiNonic,
}

/// Test one patch pair for interior crossing. Returns the pair's worst
/// non-exempt sample, or `None` when every intersection sample lies
/// within an exemption arm (see the module docs).
#[allow(clippy::too_many_arguments)]
fn pair_crossing(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    pair: (usize, usize),
    relation: RingPairRelation,
    shared_spokes: &[SharedSpoke<'_>],
    apex: Point3,
    band: f64,
    macro_depth: f64,
    macro_only: bool,
    ssi_tol: &Tolerance,
) -> Option<BinonicRingCrossing> {
    let (SurfaceGeometry::Nurbs(na), SurfaceGeometry::Nurbs(nb)) = (surf_a, surf_b) else {
        return None;
    };
    let (u_a, v_a) = nurbs_uv_range(na);
    let (u_b, v_b) = nurbs_uv_range(nb);
    let Ok(result) = intersect_surfaces(surf_a, u_a, v_a, surf_b, u_b, v_b, ssi_tol) else {
        // Solver breakdown: fail-open (see module docs).
        return None;
    };

    let mut worst: Option<BinonicRingCrossing> = None;
    for hit in &result.intersections {
        match hit {
            crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection::CoincidentRegion {
                u1_range,
                v1_range,
                u2_range,
                v2_range,
            } => {
                // Distinct ring patches may share a curve, never an area.
                let uv_a = (
                    0.5 * (u1_range.0 + u1_range.1),
                    0.5 * (v1_range.0 + v1_range.1),
                );
                let uv_b = (
                    0.5 * (u2_range.0 + u2_range.1),
                    0.5 * (v2_range.0 + v2_range.1),
                );
                use crate::geometry::surface::SurfaceEval;
                return Some(BinonicRingCrossing {
                    patch_a: pair.0,
                    patch_b: pair.1,
                    relation,
                    point: surf_a.point_at(uv_a.0, uv_a.1),
                    uv_a,
                    uv_b,
                    separation: f64::INFINITY,
                    exemption_band: band,
                    coincident_region: true,
                });
            }
            crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection::Curve(
                cos,
            ) => {
                let (t0, t1) = cos.curve_3d.parameter_range();
                let (a0, a1) = cos.curve_uv1.parameter_range();
                let (b0, b1) = cos.curve_uv2.parameter_range();
                if ![t0, t1, a0, a1, b0, b1].iter().all(|v| v.is_finite()) {
                    // Non-finite curve: fail-open (see module docs).
                    continue;
                }
                for s in 0..CURVE_SAMPLES {
                    let alpha = s as f64 / (CURVE_SAMPLES - 1) as f64;
                    let point = cos.curve_3d.point_at(t0 + (t1 - t0) * alpha);
                    if !point.x.is_finite() || !point.y.is_finite() || !point.z.is_finite() {
                        continue;
                    }
                    let mut separation = point.distance_to(&apex);
                    for spoke in shared_spokes {
                        separation =
                            separation.min(distance_to_spoke(point, spoke.polyline, spoke.patch));
                    }
                    // Strictness by relation (see module docs): apex-only
                    // pairs get the validator-consistent band; spoke
                    // neighbours (and every pair on the feasible-start
                    // screen) get the measured macro cap.
                    let exempt_depth =
                        if macro_only || relation == RingPairRelation::SpokeNeighbors {
                            macro_depth.max(band)
                        } else {
                            band
                        };
                    if separation <= exempt_depth {
                        continue;
                    }
                    let uv_a_pt = cos.curve_uv1.point_at(a0 + (a1 - a0) * alpha);
                    let uv_b_pt = cos.curve_uv2.point_at(b0 + (b1 - b0) * alpha);
                    if worst.as_ref().is_none_or(|w| separation > w.separation) {
                        worst = Some(BinonicRingCrossing {
                            patch_a: pair.0,
                            patch_b: pair.1,
                            relation,
                            point,
                            uv_a: (uv_a_pt.x, uv_a_pt.y),
                            uv_b: (uv_b_pt.x, uv_b_pt.y),
                            separation,
                            exemption_band: exempt_depth,
                            coincident_region: false,
                        });
                    }
                }
            }
        }
    }
    worst
}

/// Full parameter range of a NURBS surface from its knot vectors.
fn nurbs_uv_range(
    nurbs: &crate::geometry::nurbs_surface::NurbsSurface,
) -> ((f64, f64), (f64, f64)) {
    let knots_u = &nurbs.knots_u;
    let knots_v = &nurbs.knots_v;
    let du = nurbs.degree_u;
    let dv = nurbs.degree_v;
    (
        (knots_u[du], knots_u[knots_u.len() - du - 1]),
        (knots_v[dv], knots_v[knots_v.len() - dv - 1]),
    )
}
