//! The boundary surrogate carries the dual-jet patch's full 2-jet on
//! its boundary, with all four corners activated and distinct claims.

use super::*;
use crate::geometry::Point3;
use crate::geometry::Vec3;

fn all_corners_dual_jet() -> RectangularGregoryBiNonicDualJet {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, p) in row.iter_mut().enumerate() {
            let (x, y) = (i as f64 / 9.0, j as f64 / 9.0);
            *p = Point3::new(
                x + 0.04 * y,
                y - 0.03 * x,
                0.25 * (1.9 * x + 0.2).sin() * (2.3 * y - 0.4).cos(),
            );
        }
    }
    let mut patch =
        RectangularGregoryBiNonicDualJet::from_plain(&RectangularGregoryBiNonic { grid });
    let bump = |p: Point3, d: f64| Point3::new(p.x + d, p.y - 0.6 * d, p.z + 1.7 * d);
    for corner in 0..4 {
        patch.activate_corner(corner);
        let k = corner as f64;
        if let Some(block) = patch.corners[corner].as_mut() {
            block.cell_11 = (
                bump(block.cell_11.0, 0.008 + 0.001 * k),
                bump(block.cell_11.1, -0.010),
            );
            block.cell_12 = (
                bump(block.cell_12.0, 0.011),
                bump(block.cell_12.1, -0.007 * k),
            );
            block.cell_21 = (bump(block.cell_21.0, -0.009), bump(block.cell_21.1, 0.013));
            block.cell_22 = (
                bump(block.cell_22.0, 0.005 * k),
                bump(block.cell_22.1, -0.012),
            );
        }
    }
    patch
}

fn jet(p: Point3, parts: [Vec3; 5]) -> [f64; 18] {
    let mut out = [
        p.x, p.y, p.z, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ];
    for (k, d) in parts.iter().enumerate() {
        out[3 + 3 * k] = d.x;
        out[4 + 3 * k] = d.y;
        out[5 + 3 * k] = d.z;
    }
    out
}

fn dual_jet_2jet(patch: &RectangularGregoryBiNonicDualJet, u: f64, v: f64) -> [f64; 18] {
    jet(
        patch.point_at(u, v),
        [
            patch.partial_u(u, v),
            patch.partial_v(u, v),
            patch.second_partial_uu(u, v),
            patch.second_partial_uv(u, v),
            patch.second_partial_vv(u, v),
        ],
    )
}

fn plain_2jet(patch: &RectangularGregoryBiNonic, u: f64, v: f64) -> [f64; 18] {
    jet(
        patch.point_at(u, v),
        [
            patch.partial_u(u, v),
            patch.partial_v(u, v),
            patch.second_partial_uu(u, v),
            patch.second_partial_uv(u, v),
            patch.second_partial_vv(u, v),
        ],
    )
}

fn worst_boundary_mismatch(
    patch: &RectangularGregoryBiNonicDualJet,
    surrogate: &RectangularGregoryBiNonic,
    boundary: PatchBoundary,
) -> f64 {
    let mut worst = 0.0_f64;
    // The open boundary plus its corner limits: at the corner point
    // itself the twist is direction-dependent (see
    // `boundary_surrogate`), so the evaluator's corner convention is
    // not the along-boundary limit.
    let near = [1e-7, 1.0 - 1e-7];
    let open = (1..64).map(|k| k as f64 / 64.0);
    for t in open.chain(near) {
        let (u, v) = match boundary {
            PatchBoundary::U0 => (0.0, t),
            PatchBoundary::U1 => (1.0, t),
            PatchBoundary::V0 => (t, 0.0),
            PatchBoundary::V1 => (t, 1.0),
        };
        let (a, b) = (dual_jet_2jet(patch, u, v), plain_2jet(surrogate, u, v));
        for (x, y) in a.iter().zip(b.iter()) {
            worst = worst.max((x - y).abs());
        }
    }
    worst
}

#[test]
fn the_surrogate_carries_the_patch_2jet_on_every_boundary() {
    let patch = all_corners_dual_jet();
    for boundary in [
        PatchBoundary::U0,
        PatchBoundary::U1,
        PatchBoundary::V0,
        PatchBoundary::V1,
    ] {
        let surrogate = patch.boundary_surrogate(boundary);
        let worst = worst_boundary_mismatch(&patch, &surrogate, boundary);
        assert!(worst < 1e-9, "{boundary:?}: 2-jet mismatch {worst:e}");
    }
}

/// Negative control: the other boundary's surrogate reads the wrong
/// claims at the shared corners, and the comparison sees it.
#[test]
fn the_crossing_boundary_surrogate_is_rejected() {
    let patch = all_corners_dual_jet();
    for (boundary, wrong) in [
        (PatchBoundary::U0, PatchBoundary::V0),
        (PatchBoundary::V1, PatchBoundary::U1),
    ] {
        let surrogate = patch.boundary_surrogate(wrong);
        let worst = worst_boundary_mismatch(&patch, &surrogate, boundary);
        assert!(worst > 1e-3, "{boundary:?} vs {wrong:?}: {worst:e}");
    }
}

#[test]
fn chart_lines_name_only_exact_boundaries() {
    assert_eq!(
        PatchBoundary::of_chart_line(true, 0.0),
        Some(PatchBoundary::U0)
    );
    assert_eq!(
        PatchBoundary::of_chart_line(true, 1.0),
        Some(PatchBoundary::U1)
    );
    assert_eq!(
        PatchBoundary::of_chart_line(false, 0.0),
        Some(PatchBoundary::V0)
    );
    assert_eq!(
        PatchBoundary::of_chart_line(false, 1.0),
        Some(PatchBoundary::V1)
    );
    assert_eq!(PatchBoundary::of_chart_line(true, 1e-15), None);
    assert_eq!(PatchBoundary::of_chart_line(false, 0.5), None);
}
