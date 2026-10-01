//! The fit's negative-control corruption (row 878): split from
//! `setback_fit_fan` for the file budget.

use popkern_core::geometry::NurbsSurface;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// The negative-control corruption behind
/// [`crate::dispatch_overrides::with_fit_cap_rim_tilted`]: lift the
/// control row next to the `u = 1` rim along the unit surface normal
/// by a fifth of the row spacing. Only interior columns move — columns
/// 0 and 1 carry the spoke `v = 0` and its cross tangent, the last two
/// the `v = 1` rim and its — so exactly one thing breaks: the `u = 1`
/// rim's cross tangent (position intact).
pub(super) fn tilt_rim_row(sector: &mut NurbsSurface) {
    let (m_u, m_v) = (sector.num_u, sector.num_v);
    if m_u < 2 || m_v < 5 {
        return;
    }
    for j in 2..m_v - 2 {
        let (rim, inner) = ((m_u - 1) * m_v + j, (m_u - 2) * m_v + j);
        let spacing = sector.control_points[rim].distance_to(&sector.control_points[inner]);
        let v = j as f64 / (m_v - 1) as f64;
        let Some(normal) = sector.normal_at(1.0, v).normalized() else {
            continue;
        };
        sector.control_points[inner] = sector.control_points[inner] + normal * (0.2 * spacing);
    }
}

/// The served face after [`tilt_rim_row`] moved its plain net: a NURBS
/// face takes the net; a dual-jet face its grid and the u-claims of the
/// dual cells (the u = 1 rim reads them), the v-claims kept.
pub(super) fn tilted_face(face: &mut SurfaceGeometry, sector: &NurbsSurface) {
    use popkern_core::geometry::gregory::GregoryPatch;
    match face {
        SurfaceGeometry::Nurbs(n) => *n = sector.clone(),
        SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonicDualJet(patch)) => {
            let m = sector.num_v;
            for (i, row) in patch.grid.iter_mut().enumerate() {
                for (j, cell) in row.iter_mut().enumerate() {
                    *cell = sector.control_points[i * m + j];
                }
            }
            for block in patch.corners.iter_mut().flatten() {
                block.cell_11.0 = patch.grid[8][8];
                block.cell_12.0 = patch.grid[8][7];
                block.cell_21.0 = patch.grid[7][8];
                block.cell_22.0 = patch.grid[7][7];
            }
        }
        _ => {}
    }
}
