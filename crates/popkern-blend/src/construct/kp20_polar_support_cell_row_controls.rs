//! Row-control helpers for constructive support-cell re-authoring.

use std::collections::{HashMap, HashSet};

use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::FaceHandle;

use super::super::BoundarySide;
use super::super::chain_fit::ChainRowKey;
use super::super::chain_ribbon_field::row_derivative_coefficients;
use super::CellReauthorOutcome;

/// Rotate a row's two end legs into the local derivation span
/// `span(C'(end), R(end))`, preserving magnitude and in-plane
/// direction, and write the repaired transverse controls back.
pub(super) fn repair_row_end_legs(
    patches: &mut HashMap<FaceHandle, RectangularGregoryBiNonic>,
    key: ChainRowKey,
    row0: &[Point3; 10],
    field_segment: &[Vec3; 9],
    chained_rows: &HashSet<ChainRowKey>,
    outcome: &mut CellReauthorOutcome,
) {
    let derivative = row_derivative_coefficients(row0);
    for (corner_k, tangent, transverse) in [
        (0usize, derivative[0], field_segment[0]),
        (9, derivative[8], field_segment[8]),
    ] {
        // A chained transverse row's control is chain-authored
        // (exact G0); it must never move. Its leg is made
        // span-compatible on the field side instead (the ribbon
        // fit is constrained to contain it).
        if chained_rows.contains(&ChainRowKey {
            face: key.face,
            side: transverse_row_side(key.side, corner_k),
        }) {
            continue;
        }
        let Some(plane_normal) = tangent.cross(&transverse).normalized() else {
            continue;
        };
        let Some(patch) = patches.get(&key.face) else {
            continue;
        };
        let row1 = second_row_controls(patch, key.side);
        let leg = row1[corner_k] - row0[corner_k];
        let magnitude = leg.dot(&leg).sqrt();
        if magnitude <= 1.0e-300 {
            continue;
        }
        let out_of_plane = leg.dot(&plane_normal);
        let ratio = out_of_plane.abs() / magnitude;
        if ratio <= 1.0e-9 {
            continue;
        }
        outcome.max_leg_out_of_plane_ratio = outcome.max_leg_out_of_plane_ratio.max(ratio);
        let Some(direction) = (leg - plane_normal * out_of_plane).normalized() else {
            continue;
        };
        let (i, j) = second_row_grid_index(key.side, corner_k);
        if let Some(patch) = patches.get_mut(&key.face) {
            patch.grid[i][j] = row0[corner_k] + direction * magnitude;
            outcome.repaired_leg_count += 1;
        }
    }
}

/// Discrete Coons fill of the interior controls (indices 2..=7)
/// from the four derived band rows/columns at indices 1 and 8.
pub(super) fn refill_interior_coons(patch: &mut RectangularGregoryBiNonic) {
    let grid = patch.grid;
    let weighted = |terms: &[(f64, Point3)]| -> Point3 {
        let mut out = Point3::new(0.0, 0.0, 0.0);
        for &(weight, point) in terms {
            out = Point3::new(
                out.x + weight * point.x,
                out.y + weight * point.y,
                out.z + weight * point.z,
            );
        }
        out
    };
    for i in 2..8 {
        let a = (i as f64 - 1.0) / 7.0;
        for j in 2..8 {
            let b = (j as f64 - 1.0) / 7.0;
            patch.grid[i][j] = weighted(&[
                (1.0 - a, grid[1][j]),
                (a, grid[8][j]),
                (1.0 - b, grid[i][1]),
                (b, grid[i][8]),
                (-(1.0 - a) * (1.0 - b), grid[1][1]),
                (-(1.0 - a) * b, grid[1][8]),
                (-a * (1.0 - b), grid[8][1]),
                (-a * b, grid[8][8]),
            ]);
        }
    }
}

/// The boundary side owning the transverse control at a row
/// corner: the perpendicular boundary through that corner.
pub(crate) fn transverse_row_side(side: BoundarySide, corner_k: usize) -> BoundarySide {
    match (side, corner_k) {
        (BoundarySide::UMin | BoundarySide::UMax, 0) => BoundarySide::VMin,
        (BoundarySide::UMin | BoundarySide::UMax, _) => BoundarySide::VMax,
        (BoundarySide::VMin | BoundarySide::VMax, 0) => BoundarySide::UMin,
        (BoundarySide::VMin | BoundarySide::VMax, _) => BoundarySide::UMax,
    }
}

pub(crate) fn second_row_controls(
    patch: &RectangularGregoryBiNonic,
    side: BoundarySide,
) -> [Point3; 10] {
    let mut out = [Point3::new(0.0, 0.0, 0.0); 10];
    for (index, point) in out.iter_mut().enumerate() {
        *point = match side {
            BoundarySide::UMin => patch.grid[1][index],
            BoundarySide::UMax => patch.grid[8][index],
            BoundarySide::VMin => patch.grid[index][1],
            BoundarySide::VMax => patch.grid[index][8],
        };
    }
    out
}

pub(crate) fn second_row_grid_index(side: BoundarySide, k: usize) -> (usize, usize) {
    match side {
        BoundarySide::UMin => (1, k),
        BoundarySide::UMax => (8, k),
        BoundarySide::VMin => (k, 1),
        BoundarySide::VMax => (k, 8),
    }
}
