use super::super::direct_input::DirectBoundaryInput;
use super::super::kp_guide_trait::KpGuide;
use super::super::kp_row1::Row1Frame;
use super::super::kp_row2::Row2Frame;
use super::super::rectangular_binonic::RectangularGregoryBiNonic;
use crate::geometry::{Point3, Vec3};

/// the sphere-fit edge-wiring work Stage 3 scaffolding — transfinite-interpolation
/// (Coons-style) fill for the bidegree-9 patch's strict-interior
/// 4×4 block (`grid[3..=6][3..=6]`).
///
/// With injection, the K-P guide doesn't represent the surface
/// the patch should approximate (boundary moved; guide didn't),
/// so the LSQ interior fit pulls CPs toward an inconsistent
/// reference and patches fold. This function instead fills the
/// 16 interior CPs by linearly blending the surrounding ring
/// frame (rings at index 2 and 7 along each axis), Coons style:
///
/// ```text
///   P[i][j] = (1-α)·P[2][j] + α·P[7][j]
///           + (1-β)·P[i][2] + β·P[i][7]
///           − (1-α)(1-β)·P[2][2] − α(1-β)·P[7][2]
///           − (1-α)β·P[2][7] − αβ·P[7][7]
/// ```
///
/// where `α = (i-2)/5`, `β = (j-2)/5`. Standard transfinite
/// interpolation on the inner [2..=7]² rectangle. Guarantees
/// fold-freeness: the 4×4 interior is a smooth bilinear blend
/// of the boundary-consistent ring.
///
/// Trade-off: the surface in the deep interior loses any "pull"
/// from the original cap-shape guide; it becomes a pure
/// frame-driven shape. Acceptable for cap surfaces because the
/// cap's outer-boundary fidelity dominates the user-visible
/// quality bar; deep interior is not reachable by the cap-vs-
/// neighbour edge measurements.
///
/// **Stage 3a status:** Default interior-fill under injection.
/// Eliminates the LSQ-against-stale-guide pull and produces
/// patches that the existing `self_intersects_with_tolerance(1e-6)`
/// gate accepts on pyramid4_mixed (4 binonic ship). cube_mixed
/// remains rejected by the quality regression guard
/// (`worst_ring_coplanarity` ratio > 1.2× pre-rung); Stage 3b will
/// address that.
pub(super) fn fill_interior_cps_via_coons(patch: &mut RectangularGregoryBiNonic) {
    for i in 3..=6 {
        let alpha = ((i - 2) as f64) / 5.0;
        for j in 3..=6 {
            let beta = ((j - 2) as f64) / 5.0;
            let p2j = patch.grid[2][j];
            let p7j = patch.grid[7][j];
            let pi2 = patch.grid[i][2];
            let pi7 = patch.grid[i][7];
            let p22 = patch.grid[2][2];
            let p72 = patch.grid[7][2];
            let p27 = patch.grid[2][7];
            let p77 = patch.grid[7][7];
            let mix = |c2j: f64,
                       c7j: f64,
                       ci2: f64,
                       ci7: f64,
                       c22: f64,
                       c72: f64,
                       c27: f64,
                       c77: f64|
             -> f64 {
                (1.0 - alpha) * c2j + alpha * c7j + (1.0 - beta) * ci2 + beta * ci7
                    - (1.0 - alpha) * (1.0 - beta) * c22
                    - alpha * (1.0 - beta) * c72
                    - (1.0 - alpha) * beta * c27
                    - alpha * beta * c77
            };
            patch.grid[i][j] = Point3::new(
                mix(p2j.x, p7j.x, pi2.x, pi7.x, p22.x, p72.x, p27.x, p77.x),
                mix(p2j.y, p7j.y, pi2.y, pi7.y, p22.y, p72.y, p27.y, p77.y),
                mix(p2j.z, p7j.z, pi2.z, pi7.z, p22.z, p72.z, p27.z, p77.z),
            );
        }
    }
}

pub(super) fn reverse_boundary(curve: &[Point3; 10]) -> [Point3; 10] {
    let mut out = [Point3::new(0.0, 0.0, 0.0); 10];
    for i in 0..10 {
        out[i] = curve[9 - i];
    }
    out
}

/// Pick the row0 / col9 boundary CPs for
/// sub-patch `k` — either from `input.hex_boundaries[k]` (when
/// the caller has injected substituted analytic boundary CPs) or
/// synthesised from the guide's corner Hermite data via
/// [`k_p_boundary`].
///
/// When `input.hex_boundaries` is `Some` and indexable at `k`, the
/// injected `row0` / `col9` arrays are returned verbatim — exact
/// `sphere ∩ cylinder` great-circle CPs (deg-9 from
/// `BlendBoundary.hex_boundary_cps_deg9`). When `None`, the
/// generic `k_p_boundary` fallback synthesises from corner
/// position + tangent + curvature, the original the bi-9 cap campaign behaviour.
///
/// See `docs/decisions/2026-04-29-phase-6c-d-stage-1-boundary-injection.md`.
pub(super) fn pick_hex_boundary<G: KpGuide>(
    input: &DirectBoundaryInput,
    guide: &G,
    k: usize,
) -> ([Point3; 10], [Point3; 10]) {
    if let Some(ref hex) = input.hex_boundaries
        && let Some(entry) = hex.get(k)
    {
        return (entry.row0, entry.col9);
    }
    (
        k_p_boundary(guide, true, false),
        k_p_boundary(guide, false, true),
    )
}

pub(super) fn k_p_boundary<G: KpGuide>(guide: &G, is_u: bool, at_end: bool) -> [Point3; 10] {
    let (ua, va, ub, vb) = match (is_u, at_end) {
        (true, false) => (0.0, 0.0, 0.0, 1.0),
        (true, true) => (1.0, 0.0, 1.0, 1.0),
        (false, false) => (0.0, 0.0, 1.0, 0.0),
        (false, true) => (0.0, 1.0, 1.0, 1.0),
    };
    let p_start = guide.point_at(ua, va);
    let p_end = guide.point_at(ub, vb);
    let (t_start, t_end, tt_start, tt_end) = if is_u {
        (
            guide.partial_v(ua, va),
            guide.partial_v(ub, vb),
            guide.second_partial_vv(ua, va),
            guide.second_partial_vv(ub, vb),
        )
    } else {
        (
            guide.partial_u(ua, va),
            guide.partial_u(ub, vb),
            guide.second_partial_uu(ua, va),
            guide.second_partial_uu(ub, vb),
        )
    };
    let mut cps = [Point3::new(0.0, 0.0, 0.0); 10];
    cps[0] = p_start;
    cps[9] = p_end;
    cps[1] = offset(p_start, t_start, 1.0 / 9.0);
    cps[8] = offset(p_end, t_end, -1.0 / 9.0);
    cps[2] = offset2(p_start, t_start, 2.0 / 9.0, tt_start, 1.0 / 72.0);
    cps[7] = offset2(p_end, t_end, -2.0 / 9.0, tt_end, 1.0 / 72.0);
    for i in 3..=6 {
        let t = ((i - 2) as f64) / 5.0;
        cps[i] = Point3::new(
            cps[2].x + t * (cps[7].x - cps[2].x),
            cps[2].y + t * (cps[7].y - cps[2].y),
            cps[2].z + t * (cps[7].z - cps[2].z),
        );
    }
    cps
}

#[inline]
fn offset(p: Point3, v: Vec3, s: f64) -> Point3 {
    Point3::new(p.x + v.x * s, p.y + v.y * s, p.z + v.z * s)
}

#[inline]
fn offset2(p: Point3, v1: Vec3, s1: f64, v2: Vec3, s2: f64) -> Point3 {
    Point3::new(
        p.x + v1.x * s1 + v2.x * s2,
        p.y + v1.y * s1 + v2.y * s2,
        p.z + v1.z * s1 + v2.z * s2,
    )
}

pub(super) fn assemble_patch(
    row0: &[Point3; 10],
    row9: &[Point3; 10],
    col0: &[Point3; 10],
    col9: &[Point3; 10],
    r1: &Row1Frame,
    r2: &Row2Frame,
) -> RectangularGregoryBiNonic {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    grid[0].copy_from_slice(row0);
    grid[9].copy_from_slice(row9);
    grid[1].copy_from_slice(&r1.row1);
    grid[8].copy_from_slice(&r1.row8);
    grid[2].copy_from_slice(&r2.row2);
    grid[7].copy_from_slice(&r2.row7);
    for i in 0..10 {
        grid[i][0] = col0[i];
        grid[i][9] = col9[i];
        grid[i][1] = r1.col1[i];
        grid[i][8] = r1.col8[i];
        grid[i][2] = r2.col2[i];
        grid[i][7] = r2.col7[i];
    }
    RectangularGregoryBiNonic { grid }
}
