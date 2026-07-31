//! Welded polar-wedge spoke pairing (PK-219).
//!
//! A welded sector's face side collapses to the corner point, so
//! BOTH of its spoke lines are columns: A's u=0 column coincides
//! with its neighbour B's right column. A collapsed line (the pole
//! row) would match everything, so candidates must have real extent.

use popkern_core::geometry::Point3;

use super::{SharedSpoke, SpokePairing, column_cp};
use crate::apex_construction::exact_boundary_cap::solve::ExactBoundarySolvedSurfaces;

pub(super) fn find_welded_shared_spokes(
    solved: &ExactBoundarySolvedSurfaces,
    tolerance: f64,
    spokes: &mut Vec<SharedSpoke>,
) {
    for (a, sa) in solved.surfaces.iter().enumerate() {
        for (b, sb) in solved.surfaces.iter().enumerate() {
            if a == b {
                continue;
            }
            let (na, nb) = (sa.surface.num_v, sb.surface.num_v);
            if na < 3 || na != nb {
                continue;
            }
            let left: Vec<Point3> = (0..na)
                .filter_map(|j| sa.surface.control_points.get(j).copied())
                .collect();
            let right: Vec<Point3> = (0..nb)
                .filter_map(|j| column_cp(&sb.surface, 0, j))
                .collect();
            if left.len() != na || right.len() != nb {
                continue;
            }
            if line_is_collapsed(&left, tolerance) || line_is_collapsed(&right, tolerance) {
                continue;
            }
            let forward = left
                .iter()
                .zip(right.iter())
                .all(|(p, q)| p.distance_to(q) <= tolerance);
            let backward = left
                .iter()
                .zip(right.iter().rev())
                .all(|(p, q)| p.distance_to(q) <= tolerance);
            if forward || backward {
                spokes.push(SharedSpoke {
                    column_sector: b,
                    row_sector: a,
                    reversed: !forward,
                    len: na,
                    pairing: SpokePairing::LeftColumnToColumn,
                });
            }
        }
    }
}

/// True when the line has no real extent (every point within
/// tolerance of the first — e.g. the collapsed pole row).
fn line_is_collapsed(line: &[Point3], tolerance: f64) -> bool {
    let Some(first) = line.first() else {
        return true;
    };
    line.iter().all(|p| p.distance_to(first) <= tolerance)
}
