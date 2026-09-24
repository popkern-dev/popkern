//! Dual cross-jet corner substrate for the setback cap fan (PK-519).
//!
//! The KP20 fan's corner piece carries TWO rim sides meeting at the
//! setback corner `c`, on two different fillet strips. Near `c`, one
//! side's depth-2 (curvature) row crosses the other side's depth-0/1
//! rows, and the shared cells cannot serve both strips' transverse
//! cross-curvature fields (measured 0.96 against the 5e-2 gate on the
//! pentagon flagship; the full falsification inventory is in
//! `docs/decisions/2026-09-23-setback-corner-dual-jet-substrate.md`).
//!
//! This module owns the substrate plumbing:
//! - **Promotion**: after the rim jet reconcile, every fan piece with
//!   a registered row rim side AND a registered column rim side is
//!   promoted to [`RectangularGregoryBiNonicDualJet`] with the shared
//!   corner activated (claims initialize to the grid — bit-identical
//!   until the curvature pass writes a claim).
//! - **Claim-aware row reading**: a side's effective boundary rows
//!   substitute its own claims into activated cells — exactly what
//!   that boundary's jets evaluate to under the squared-weight blend.
//! - **Claim writes**: the curvature pass's per-side end corrections
//!   route into the side's own claims instead of a shared constant.

use popkern_core::geometry::gregory::{
    DualJetCell, GregoryPatch, RectangularGregoryBiNonicDualJet, dual_jet_cell_of,
};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{Point3, Vec3};

use super::bezier9_eval;
use super::setback_rim_reconcile::distance_to_polylines;

/// Rim-band proximity shared with the jet and curvature reconciles.
pub(super) const RIM_NEAR: f64 = 2.0e-3;

/// Grid coordinate of side `side`'s depth-`d` row at position `k`
/// (side 0 = row `i = 0`, 1 = row `i = 9`, 2 = column `j = 0`,
/// 3 = column `j = 9`; `k` runs along the side).
pub(super) fn coord_for(side: usize, d: usize, k: usize) -> (usize, usize) {
    match side {
        0 => (d, k),
        1 => (9 - d, k),
        2 => (k, d),
        _ => (k, 9 - d),
    }
}

/// Whether a side is a u-constant boundary (grid row). Row sides
/// realize the `P_u` claim slot of an activated cell; column sides
/// realize `P_v`.
pub(super) fn side_is_row(side: usize) -> bool {
    side < 2
}

/// The patch corner index (dual-jet CCW convention) at the given end
/// of a side's rows: `at_k9 = false` is the `k = 0` end.
pub(super) fn corner_for_side_end(side: usize, at_k9: bool) -> usize {
    match (side, at_k9) {
        (0, false) => 0, // (0,0)
        (0, true) => 1,  // (0,9)
        (1, false) => 3, // (9,0)
        (1, true) => 2,  // (9,9)
        (2, false) => 0,
        (2, true) => 3,
        (3, false) => 1,
        _ => 2,
    }
}

/// The claim cell classes a side's depth-2 row addresses at a dual
/// end, ordered by distance from the corner: 1 step (the other
/// side's tangency band) then 2 steps (the shared curvature cell).
pub(super) fn dual_end_cells(side: usize) -> [DualJetCell; 2] {
    if side_is_row(side) {
        [DualJetCell::Cell21, DualJetCell::Cell22]
    } else {
        [DualJetCell::Cell12, DualJetCell::Cell22]
    }
}

/// The activated corner at the given end of a side, when the side's
/// patch is a dual-jet piece with that corner active.
pub(super) fn active_corner_at_end(
    patches: &[GregoryPatch],
    pi: usize,
    side: usize,
    at_k9: bool,
) -> Option<usize> {
    let corner = corner_for_side_end(side, at_k9);
    match &patches[pi] {
        GregoryPatch::RectangularBiNonicDualJet(r) if r.corners[corner].is_some() => Some(corner),
        _ => None,
    }
}

/// A side's effective depth-`d` row: grid cells, with activated cells
/// replaced by THIS side's claims (what the side's boundary jets
/// evaluate to under the squared-weight blend).
pub(super) fn row_for_side(
    patches: &[GregoryPatch],
    pi: usize,
    side: usize,
    d: usize,
) -> Option<[Point3; 10]> {
    let grid = patches[pi].binonic_grid()?;
    let dual: Option<&RectangularGregoryBiNonicDualJet> = match &patches[pi] {
        GregoryPatch::RectangularBiNonicDualJet(r) => Some(r),
        _ => None,
    };
    Some(core::array::from_fn(|k| {
        let (i, j) = coord_for(side, d, k);
        if let Some(r) = dual
            && let Some((corner, cell)) = dual_jet_cell_of(i, j)
            && let Some(block) = r.corners[corner].as_ref()
        {
            let pair = match cell {
                DualJetCell::Cell11 => block.cell_11,
                DualJetCell::Cell12 => block.cell_12,
                DualJetCell::Cell21 => block.cell_21,
                DualJetCell::Cell22 => block.cell_22,
            };
            return if side_is_row(side) { pair.0 } else { pair.1 };
        }
        grid[i][j]
    }))
}

/// Add `delta` to a side's own claim of `cell` at `corner` on patch
/// `pi`. A non-dual patch or inactive corner is a no-op (the caller
/// routes those writes through the grid path).
pub(super) fn add_claim_delta(
    patches: &mut [GregoryPatch],
    pi: usize,
    side: usize,
    corner: usize,
    cell: DualJetCell,
    delta: Vec3,
) {
    let GregoryPatch::RectangularBiNonicDualJet(r) = &mut patches[pi] else {
        return;
    };
    let Some(block) = r.corners[corner].as_mut() else {
        return;
    };
    let pair = match cell {
        DualJetCell::Cell11 => &mut block.cell_11,
        DualJetCell::Cell12 => &mut block.cell_12,
        DualJetCell::Cell21 => &mut block.cell_21,
        DualJetCell::Cell22 => &mut block.cell_22,
    };
    let claim = if side_is_row(side) {
        &mut pair.0
    } else {
        &mut pair.1
    };
    *claim = Point3::new(claim.x + delta.x, claim.y + delta.y, claim.z + delta.z);
}

/// Register rim sides: near a support along the whole row AND on the
/// setback border loop (the same classifier as the jet reconcile,
/// extended to read dual-jet pieces).
pub(super) fn rim_sides(
    patches: &[GregoryPatch],
    candidates: &[SurfaceGeometry],
    border_loop: &[Vec<Point3>],
    rim_off_loop_allowance: f64,
) -> Vec<(usize, usize)> {
    let mut sides: Vec<(usize, usize)> = Vec::new();
    for (pi, patch) in patches.iter().enumerate() {
        let Some(grid) = patch.binonic_grid() else {
            continue;
        };
        for side in 0..4 {
            let row: [Point3; 10] = core::array::from_fn(|k| {
                let (i, j) = coord_for(side, 0, k);
                grid[i][j]
            });
            let on_support = (0..=16).all(|i| {
                let pt = bezier9_eval(&row, i as f64 / 16.0);
                candidates.iter().any(|sfc| {
                    project_point_to_uv(sfc, &pt)
                        .map(|(u, v)| sfc.point_at(u, v).distance_to(&pt) < RIM_NEAR)
                        .unwrap_or(false)
                })
            });
            if !on_support {
                continue;
            }
            let worst_off_loop = (0..=16)
                .map(|i| distance_to_polylines(bezier9_eval(&row, i as f64 / 16.0), border_loop))
                .fold(0.0_f64, f64::max);
            if worst_off_loop > rim_off_loop_allowance.max(RIM_NEAR) {
                continue;
            }
            sides.push((pi, side));
        }
    }
    sides
}

/// Promote every fan piece carrying a row rim side AND a column rim
/// side to a dual-jet patch with the shared corner activated.
/// Returns the number of promoted pieces. Bit-identical evaluation:
/// claims initialize from the grid; only the curvature pass writes
/// them apart.
pub(super) fn promote_dual_jet_corners(
    patches: &mut [GregoryPatch],
    candidates: &[SurfaceGeometry],
    border_loop: &[Vec<Point3>],
    rim_off_loop_allowance: f64,
) -> usize {
    let sides = rim_sides(patches, candidates, border_loop, rim_off_loop_allowance);
    let mut per_patch: std::collections::BTreeMap<usize, Vec<usize>> =
        std::collections::BTreeMap::new();
    for (pi, side) in sides {
        per_patch.entry(pi).or_default().push(side);
    }
    let mut promoted = 0usize;
    for (pi, patch_sides) in per_patch {
        let rows: Vec<usize> = patch_sides.iter().copied().filter(|&s| s < 2).collect();
        let cols: Vec<usize> = patch_sides.iter().copied().filter(|&s| s >= 2).collect();
        // The fan's corner piece carries exactly one rim row side and
        // one rim column side. Any other census is not the corner
        // class this substrate addresses — leave those pieces plain
        // (the shared-constant fallback semantics still apply there).
        if rows.len() != 1 || cols.len() != 1 {
            continue;
        }
        let corner = shared_corner(rows[0], cols[0]);
        let plain = match &patches[pi] {
            GregoryPatch::RectangularBiNonic(bp) => (**bp).clone(),
            _ => continue,
        };
        let mut dual = RectangularGregoryBiNonicDualJet::from_plain(&plain);
        dual.activate_corner(corner);
        patches[pi] = GregoryPatch::RectangularBiNonicDualJet(Box::new(dual));
        promoted += 1;
        if popkern_core::trace::trace_enabled("setback") {
            eprintln!(
                "  [dual-jet] promoted patch {pi} corner {corner} (row side {} / col side {})",
                rows[0], cols[0]
            );
        }
    }
    promoted
}

/// The patch corner where a row side and a column side meet.
fn shared_corner(row_side: usize, col_side: usize) -> usize {
    match (row_side, col_side) {
        (0, 2) => 0,
        (0, 3) => 1,
        (1, 3) => 2,
        _ => 3,
    }
}

/// Grid position of a side's depth-2 row cell at distance `step`
/// (1 or 2) from the corner end — used by the write path to assert
/// the claim cells line up with [`coord_for`].
#[cfg(test)]
pub(super) fn dual_end_cell_position(side: usize, at_k9: bool, step: usize) -> (usize, usize) {
    let corner = corner_for_side_end(side, at_k9);
    let cells = dual_end_cells(side);
    popkern_core::geometry::gregory::dual_jet_cell_position(corner, cells[step - 1])
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The claim-cell positions must be exactly the side's own d2-row
    /// cells nearest each dual end, for every side/end combination.
    #[test]
    fn dual_end_cells_line_up_with_row_coordinates() {
        for side in 0..4usize {
            for at_k9 in [false, true] {
                for step in [1usize, 2] {
                    let k = if at_k9 { 9 - step } else { step };
                    let expected = coord_for(side, 2, k);
                    let got = dual_end_cell_position(side, at_k9, step);
                    assert_eq!(
                        got, expected,
                        "side {side} at_k9 {at_k9} step {step}: cell {got:?} vs row {expected:?}"
                    );
                }
            }
        }
    }
}
