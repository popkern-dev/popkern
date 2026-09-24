//! Gregory-patch control-net transform helper.
//!
//! Extracted from `brep_transforms.rs` to stay under the 500-line
//! structural cap. Every Gregory family transforms by mapping each
//! control point; parameterization is knot/UV-space and invariant.

use crate::geometry::Transform3;
use crate::geometry::gregory::{
    GregoryPatch, RectangularGregory, RectangularGregoryBiNonic,
    RectangularGregoryBiNonicDualTwist, RectangularGregoryBiquinticFull,
    RectangularGregoryBiseptic, RectangularGregoryDeg5, RectangularGregoryTwoStripe,
    TriangularGregory,
};

/// Apply a rigid transform to a Gregory patch by transforming every
/// control point (boundary + inner candidates).
pub(super) fn transform_gregory(g: &GregoryPatch, t: &Transform3) -> GregoryPatch {
    match g {
        GregoryPatch::Rectangular(r) => {
            let mut boundary = r.boundary_cps;
            for row in boundary.iter_mut() {
                for cell in row.iter_mut() {
                    *cell = t.apply_point(cell);
                }
            }
            let inner_cps = r
                .inner_cps
                .map(|(p_u, p_v)| (t.apply_point(&p_u), t.apply_point(&p_v)));
            GregoryPatch::Rectangular(Box::new(RectangularGregory {
                boundary_cps: boundary,
                inner_cps,
            }))
        }
        GregoryPatch::RectangularDeg5(r) => {
            let mut boundary = r.boundary_cps;
            for row in boundary.iter_mut() {
                for cell in row.iter_mut() {
                    *cell = t.apply_point(cell);
                }
            }
            let inner_cps = r
                .inner_cps
                .map(|(p_u, p_v)| (t.apply_point(&p_u), t.apply_point(&p_v)));
            GregoryPatch::RectangularDeg5(Box::new(RectangularGregoryDeg5 {
                boundary_cps: boundary,
                inner_cps,
            }))
        }
        GregoryPatch::RectangularBiquinticFull(r) => {
            let mut grid = r.grid;
            for row in grid.iter_mut() {
                for cell in row.iter_mut() {
                    *cell = t.apply_point(cell);
                }
            }
            let corner_inner_cps = r
                .corner_inner_cps
                .map(|(p_u, p_v)| (t.apply_point(&p_u), t.apply_point(&p_v)));
            GregoryPatch::RectangularBiquinticFull(Box::new(RectangularGregoryBiquinticFull {
                grid,
                corner_inner_cps,
            }))
        }
        GregoryPatch::RectangularBiseptic(r) => {
            let mut grid = r.grid;
            for row in grid.iter_mut() {
                for cell in row.iter_mut() {
                    *cell = t.apply_point(cell);
                }
            }
            let corner_inner_cps = r
                .corner_inner_cps
                .map(|(p_u, p_v)| (t.apply_point(&p_u), t.apply_point(&p_v)));
            GregoryPatch::RectangularBiseptic(Box::new(RectangularGregoryBiseptic {
                grid,
                corner_inner_cps,
            }))
        }
        GregoryPatch::RectangularTwoStripe(r) => {
            let mut grid = r.grid;
            for row in grid.iter_mut() {
                for cell in row.iter_mut() {
                    *cell = t.apply_point(cell);
                }
            }
            GregoryPatch::RectangularTwoStripe(Box::new(RectangularGregoryTwoStripe { grid }))
        }
        GregoryPatch::RectangularBiNonic(r) => {
            let mut grid = r.grid;
            for row in grid.iter_mut() {
                for cell in row.iter_mut() {
                    *cell = t.apply_point(cell);
                }
            }
            GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid }))
        }
        GregoryPatch::RectangularBiNonicDualTwist(r) => {
            let mut grid = r.grid;
            for row in grid.iter_mut() {
                for cell in row.iter_mut() {
                    *cell = t.apply_point(cell);
                }
            }
            let inner_cps = r
                .inner_cps
                .map(|(p_u, p_v)| (t.apply_point(&p_u), t.apply_point(&p_v)));
            GregoryPatch::RectangularBiNonicDualTwist(Box::new(
                RectangularGregoryBiNonicDualTwist { grid, inner_cps },
            ))
        }
        GregoryPatch::RectangularBiNonicDualJet(r) => {
            let mut grid = r.grid;
            for row in grid.iter_mut() {
                for cell in row.iter_mut() {
                    *cell = t.apply_point(cell);
                }
            }
            let corners = r.corners.map(|corner| {
                corner.map(|block| crate::geometry::gregory::DualJetCorner {
                    cell_12: (
                        t.apply_point(&block.cell_12.0),
                        t.apply_point(&block.cell_12.1),
                    ),
                    cell_21: (
                        t.apply_point(&block.cell_21.0),
                        t.apply_point(&block.cell_21.1),
                    ),
                    cell_22: (
                        t.apply_point(&block.cell_22.0),
                        t.apply_point(&block.cell_22.1),
                    ),
                })
            });
            GregoryPatch::RectangularBiNonicDualJet(Box::new(
                crate::geometry::gregory::RectangularGregoryBiNonicDualJet { grid, corners },
            ))
        }
        GregoryPatch::Triangular(tr) => {
            let boundary_cps = tr.boundary_cps.map(|p| t.apply_point(&p));
            let inner_cps = tr.inner_cps.map(|p| t.apply_point(&p));
            GregoryPatch::Triangular(Box::new(TriangularGregory {
                boundary_cps,
                inner_cps,
            }))
        }
    }
}
