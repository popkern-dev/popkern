//! Make adjacent exact-boundary cap sectors G¹ across their shared
//! spoke.
//!
//! # The defect
//!
//! `apply_derivative_rows` constrains each wedge against its *outer*
//! boundaries — the fillet arc and the face bridge — so every sector
//! meets its supports well. Nothing ever constrains the two *inner*
//! spokes, which are the boundaries a sector shares with its
//! neighbours. The sectors therefore meet each other at G⁰ only, and
//! measured across the direct-apex matrix all twelve committed caps
//! exceed the crate's 5e-2 mid-seam G¹ bar, the best at 0.80 (~53°).
//!
//! # The control net
//!
//! `coons_controls` lays a sector out as `controls[i * num_v + j]`:
//!
//! ```text
//!   j = 0            bottom  — fillet arc      (OUTER RIM)
//!   i = 0            left    — face bridge     (OUTER RIM)
//!   i = num_u - 1    right   — fillet spoke    (SHARED)
//!   j = num_v - 1    top     — face spoke      (SHARED)
//!   (num_u-1, num_v-1)       — the pole
//! ```
//!
//! So one sector's *right column* and its neighbour's *top row* are
//! the same physical spoke, each running rim → pole.
//!
//! # The condition
//!
//! For a shared tangent plane across the spoke, each spoke control
//! point must be collinear with — and midway between — the interior
//! control points on either side of it. The spoke itself is a trim
//! curve and must not move, so the correction is a common offset
//! applied to the two interior points:
//!
//! ```text
//!   offset = S(k) - (A_int(k) + B_int(k)) / 2
//!   A_int(k) += offset;   B_int(k) += offset
//! ```
//!
//! # What must not move
//!
//! Both ends of the spoke are excluded. `k = 0` sits on the outer rim,
//! and moving its interior neighbour breaks the exact-boundary
//! contract that the cap *contains* the rim curve. The far end is the
//! pole, whose interior neighbour lies on the sector's *other* shared
//! spoke — moving it would break that spoke's G⁰. Only the strictly
//! interior span is corrected.
//!
//! # Status
//!
//! Default **on**, after measurement rather than on principle:
//!
//! ```text
//!   cap        off        on
//!   1     8.475e-1   1.490e-8
//!   2     8.490e-1   1.426e-5
//!   3     9.082e-1   6.925e-1
//! ```
//!
//! Two caps go from ~53° of mid-seam misalignment to machine zero.
//! The third stayed creased until 2026-07-26: its sectors solve to
//! non-square nets that index-wise pairing cannot see, which
//! `harmonize_mismatched_spoke_nets` now repairs by exact knot
//! insertion before pairing (n=3 mixed apex: cos_min 0.00214 ->
//! 0.999994, pinned by `cross_spoke_correction_fixes_both_low_valences`).
//!
//! Enabled against the suites that could falsify it: direct-apex
//! matrix 20/20, polar commercial gates 11/11, fillet corpus 18/22
//! with 18/18 OCC-clean — all unchanged. `POPKERN_CAP_SPOKE_G1=0`
//! restores the previous behaviour for A/B.

use popkern_core::geometry::Point3;

use super::solve::ExactBoundarySolvedSurfaces;

/// Whether the cross-spoke correction runs.
pub(crate) fn cross_spoke_g1_enabled() -> bool {
    std::env::var("POPKERN_CAP_SPOKE_G1")
        .map(|v| v != "0")
        .unwrap_or(true)
}

/// A spoke shared by two solved sectors.
struct SharedSpoke {
    /// Index of the sector contributing its right column.
    column_sector: usize,
    /// Index of the sector contributing its top row.
    row_sector: usize,
    /// Whether the row runs opposite to the column.
    reversed: bool,
    /// Number of control points along the spoke.
    len: usize,
}

/// What the correction managed on one corner.
pub(crate) struct CrossSpokeOutcome {
    /// Control-point pairs corrected.
    pub corrected: usize,
    /// Shared spokes matched and corrected.
    pub spokes: usize,
    /// Neighbours that ARE geometrically adjacent but whose control
    /// nets have incompatible sizes, so no index-wise correspondence
    /// exists to correct through.
    ///
    /// Non-square nets — the n = 3 class, measured `(49, 41)` and
    /// `(56, 48)` — are first harmonized by exact knot insertion
    /// (`harmonize_mismatched_spoke_nets`); anything still counted
    /// here is genuinely incompatible (degree mismatch, or refined
    /// nets that do not coincide along the spoke). Counted rather
    /// than skipped, because a silent skip is what let the original
    /// crease survive unmeasured.
    pub dimension_mismatched: usize,
}

/// Enforce cross-spoke G¹ across every shared spoke.
pub(crate) fn enforce_cross_spoke_g1(
    solved: &mut ExactBoundarySolvedSurfaces,
    tolerance: f64,
) -> CrossSpokeOutcome {
    harmonize_mismatched_spoke_nets(solved, tolerance);
    let (spokes, dimension_mismatched) = find_shared_spokes(solved, tolerance);
    let spoke_count = spokes.len();
    let mut corrected = 0;
    for spoke in spokes {
        corrected += correct_spoke(solved, &spoke);
    }
    CrossSpokeOutcome {
        corrected,
        spokes: spoke_count,
        dimension_mismatched,
    }
}

/// Bring dimension-mismatched neighbouring nets to a common spoke
/// knot vector so the index-wise correction can serve them.
///
/// This is the n = 3 residual's mechanism: its sectors solve to
/// non-square nets — measured `(49, 41)` and `(56, 48)` — so a
/// 41-point column cannot be index-matched against a 49-point row
/// and the crease survived the correction. Knot insertion is exact
/// (the refined surface is the same surface), so harmonizing the two
/// spoke-direction knot vectors to their union creates the missing
/// index correspondence without moving any geometry. Pairs that
/// still cannot harmonize — degree mismatch, incompatible parameter
/// ranges, or refined nets that still do not coincide along the
/// spoke — are left for `find_shared_spokes` to count in
/// `dimension_mismatched`, never silently skipped.
fn harmonize_mismatched_spoke_nets(solved: &mut ExactBoundarySolvedSurfaces, tolerance: f64) {
    let n = solved.surfaces.len();
    for a in 0..n {
        for b in 0..n {
            if a == b {
                continue;
            }
            let sa = &solved.surfaces[a].surface;
            let sb = &solved.surfaces[b].surface;
            if sa.num_v == sb.num_u || sa.num_v < 3 || sb.num_u < 3 {
                continue;
            }
            if sa.degree_v != sb.degree_u {
                continue;
            }
            if !spoke_endpoints_touch(&solved.surfaces[a], &solved.surfaces[b], tolerance) {
                continue;
            }
            // Column side runs along v; row side along u. Establish
            // orientation from the endpoints before mapping knots.
            let (Some(a0), Some(a1)) = (
                column_cp(sa, 0, 0),
                column_cp(sa, 0, sa.num_v.saturating_sub(1)),
            ) else {
                continue;
            };
            let (Some(b0), _) = (row_cp(sb, 0, 0), ()) else {
                continue;
            };
            let reversed = !(a0.distance_to(&b0) <= tolerance);
            let _ = a1;
            let (va0, va1) = knot_range(&sa.knots_v);
            let (ub0, ub1) = knot_range(&sb.knots_u);
            if (va1 - va0).abs() < 1e-12 || (ub1 - ub0).abs() < 1e-12 {
                continue;
            }
            // Interior knot multiset of each side, mapped into the
            // other side's parameter range (mirrored when reversed).
            let map_b_to_a = |t: f64| {
                let f = (t - ub0) / (ub1 - ub0);
                let f = if reversed { 1.0 - f } else { f };
                va0 + f * (va1 - va0)
            };
            let map_a_to_b = |t: f64| {
                let f = (t - va0) / (va1 - va0);
                let f = if reversed { 1.0 - f } else { f };
                ub0 + f * (ub1 - ub0)
            };
            let missing_in_a =
                missing_interior_knots(&sa.knots_v, sa.degree_v, &sb.knots_u, sb.degree_u, |t| {
                    map_b_to_a(t)
                });
            let missing_in_b =
                missing_interior_knots(&sb.knots_u, sb.degree_u, &sa.knots_v, sa.degree_v, |t| {
                    map_a_to_b(t)
                });
            let refined_a = if missing_in_a.is_empty() {
                sa.clone()
            } else {
                match sa.insert_knots_v(&missing_in_a) {
                    Ok(s) => s,
                    Err(_) => continue,
                }
            };
            let refined_b = if missing_in_b.is_empty() {
                sb.clone()
            } else {
                match sb.insert_knots_u(&missing_in_b) {
                    Ok(s) => s,
                    Err(_) => continue,
                }
            };
            if refined_a.num_v != refined_b.num_u {
                continue;
            }
            // Only commit when the refined nets actually coincide
            // along the spoke — otherwise the pair is genuinely
            // incompatible and must keep counting as mismatched.
            let column: Vec<Point3> = (0..refined_a.num_v)
                .filter_map(|j| column_cp(&refined_a, 0, j))
                .collect();
            let row: Vec<Point3> = (0..refined_b.num_u)
                .filter_map(|i| row_cp(&refined_b, 0, i))
                .collect();
            let forward = column
                .iter()
                .zip(row.iter())
                .all(|(p, q)| p.distance_to(q) <= tolerance);
            let backward = column
                .iter()
                .zip(row.iter().rev())
                .all(|(p, q)| p.distance_to(q) <= tolerance);
            if forward || backward {
                solved.surfaces[a].surface = refined_a;
                solved.surfaces[b].surface = refined_b;
            }
        }
    }
}

/// First/last knot value (the parameter range of a clamped vector).
fn knot_range(knots: &[f64]) -> (f64, f64) {
    (
        knots.first().copied().unwrap_or(0.0),
        knots.last().copied().unwrap_or(1.0),
    )
}

/// Interior knots of `other` (mapped into `own`'s range) that `own`
/// is missing, respecting multiplicity up to the other's count.
fn missing_interior_knots(
    own: &[f64],
    own_degree: usize,
    other: &[f64],
    other_degree: usize,
    map: impl Fn(f64) -> f64,
) -> Vec<f64> {
    const KNOT_EPS: f64 = 1e-9;
    let interior = |knots: &[f64], degree: usize| -> Vec<f64> {
        let end = knots.len().saturating_sub(degree + 1);
        knots
            .get(degree + 1..end.max(degree + 1))
            .map(|s| s.to_vec())
            .unwrap_or_default()
    };
    let own_interior = interior(own, own_degree);
    let mut missing = Vec::new();
    for &t in &interior(other, other_degree) {
        let mapped = map(t);
        let mult_other = 1 + missing
            .iter()
            .filter(|&&m: &&f64| (m - mapped).abs() <= KNOT_EPS)
            .count();
        let mult_own = own_interior
            .iter()
            .filter(|&&k| (k - mapped).abs() <= KNOT_EPS)
            .count();
        if mult_own < mult_other {
            missing.push(mapped);
        }
    }
    missing.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    missing
}

fn column_cp(
    surface: &popkern_core::geometry::NurbsSurface,
    offset_from_edge: usize,
    j: usize,
) -> Option<Point3> {
    let i = surface.num_u.checked_sub(1 + offset_from_edge)?;
    surface.control_points.get(i * surface.num_v + j).copied()
}

fn row_cp(
    surface: &popkern_core::geometry::NurbsSurface,
    offset_from_edge: usize,
    i: usize,
) -> Option<Point3> {
    let j = surface.num_v.checked_sub(1 + offset_from_edge)?;
    surface.control_points.get(i * surface.num_v + j).copied()
}

/// Pair sectors whose right column and top row trace the same spoke.
fn find_shared_spokes(
    solved: &ExactBoundarySolvedSurfaces,
    tolerance: f64,
) -> (Vec<SharedSpoke>, usize) {
    let mut spokes = Vec::new();
    let mut dimension_mismatched = 0;
    for (a, sa) in solved.surfaces.iter().enumerate() {
        for (b, sb) in solved.surfaces.iter().enumerate() {
            if a == b {
                continue;
            }
            let (na, nb) = (sa.surface.num_v, sb.surface.num_u);
            if na < 3 || nb < 3 {
                continue;
            }
            if na != nb {
                // Geometrically adjacent sectors can still have
                // incompatible net sizes; count it so the gap is
                // visible rather than silently unserved.
                if spoke_endpoints_touch(sa, sb, tolerance) {
                    dimension_mismatched += 1;
                }
                continue;
            }
            let column: Vec<Point3> = (0..na)
                .filter_map(|j| column_cp(&sa.surface, 0, j))
                .collect();
            let row: Vec<Point3> = (0..nb).filter_map(|i| row_cp(&sb.surface, 0, i)).collect();
            if column.len() != na || row.len() != nb {
                continue;
            }
            let forward = column
                .iter()
                .zip(row.iter())
                .all(|(p, q)| p.distance_to(q) <= tolerance);
            let backward = column
                .iter()
                .zip(row.iter().rev())
                .all(|(p, q)| p.distance_to(q) <= tolerance);
            if forward || backward {
                spokes.push(SharedSpoke {
                    column_sector: a,
                    row_sector: b,
                    reversed: !forward,
                    len: na,
                });
            }
        }
    }
    (spokes, dimension_mismatched)
}

/// Whether two sectors' spoke ends coincide, i.e. they really are
/// neighbours even when their nets cannot be index-matched.
fn spoke_endpoints_touch(
    a: &super::solve::ExactBoundarySectorSurface,
    b: &super::solve::ExactBoundarySectorSurface,
    tolerance: f64,
) -> bool {
    let (Some(a0), Some(a1)) = (
        column_cp(&a.surface, 0, 0),
        column_cp(&a.surface, 0, a.surface.num_v.saturating_sub(1)),
    ) else {
        return false;
    };
    let (Some(b0), Some(b1)) = (
        row_cp(&b.surface, 0, 0),
        row_cp(&b.surface, 0, b.surface.num_u.saturating_sub(1)),
    ) else {
        return false;
    };
    (a0.distance_to(&b0) <= tolerance && a1.distance_to(&b1) <= tolerance)
        || (a0.distance_to(&b1) <= tolerance && a1.distance_to(&b0) <= tolerance)
}

/// Apply the midpoint correction along one shared spoke.
fn correct_spoke(solved: &mut ExactBoundarySolvedSurfaces, spoke: &SharedSpoke) -> usize {
    let mut corrected = 0;
    // Strictly interior only: k = 0 is on the outer rim and the far
    // end is the pole, whose neighbour lies on the other shared spoke.
    for k in 1..spoke.len.saturating_sub(1) {
        let k_row = if spoke.reversed { spoke.len - 1 - k } else { k };
        let (Some(shared), Some(a_int), Some(b_int)) = (
            column_cp(&solved.surfaces[spoke.column_sector].surface, 0, k),
            column_cp(&solved.surfaces[spoke.column_sector].surface, 1, k),
            row_cp(&solved.surfaces[spoke.row_sector].surface, 1, k_row),
        ) else {
            continue;
        };
        let midpoint = Point3::new(
            0.5 * (a_int.x + b_int.x),
            0.5 * (a_int.y + b_int.y),
            0.5 * (a_int.z + b_int.z),
        );
        let offset = shared - midpoint;
        set_column_cp(
            &mut solved.surfaces[spoke.column_sector].surface,
            1,
            k,
            a_int + offset,
        );
        set_row_cp(
            &mut solved.surfaces[spoke.row_sector].surface,
            1,
            k_row,
            b_int + offset,
        );
        corrected += 1;
    }
    corrected
}

fn set_column_cp(
    surface: &mut popkern_core::geometry::NurbsSurface,
    offset_from_edge: usize,
    j: usize,
    value: Point3,
) {
    if let Some(i) = surface.num_u.checked_sub(1 + offset_from_edge) {
        let idx = i * surface.num_v + j;
        if let Some(slot) = surface.control_points.get_mut(idx) {
            *slot = value;
        }
    }
}

fn set_row_cp(
    surface: &mut popkern_core::geometry::NurbsSurface,
    offset_from_edge: usize,
    i: usize,
    value: Point3,
) {
    if let Some(j) = surface.num_v.checked_sub(1 + offset_from_edge) {
        let idx = i * surface.num_v + j;
        if let Some(slot) = surface.control_points.get_mut(idx) {
            *slot = value;
        }
    }
}
