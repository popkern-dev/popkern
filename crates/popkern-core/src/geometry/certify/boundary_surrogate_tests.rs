//! Exact patch-boundary detection for the dual-jet surrogate: only a
//! pcurve lying on a boundary in real arithmetic qualifies.

use super::*;
use crate::geometry::Point3;
use crate::geometry::gregory::RectangularGregoryBiNonic;
use crate::geometry::nurbs_curve::NurbsCurve;

fn nurbs_pcurve(points: &[(f64, f64)], weights: Vec<f64>) -> Pcurve {
    let n = points.len();
    let mut knots = vec![0.0; 2];
    knots.extend((1..n - 1).map(|k| k as f64 / (n - 1) as f64));
    knots.extend([1.0, 1.0]);
    Pcurve::from_nurbs(NurbsCurve {
        degree: 1,
        knots,
        control_points: points
            .iter()
            .map(|&(u, v)| Point3::new(u, v, 0.0))
            .collect(),
        weights,
    })
}

#[test]
fn only_exact_boundary_pcurves_name_a_boundary() {
    let line = |a, b| Pcurve::linear(a, b, 0.0, 1.0);
    assert_eq!(
        pcurve_patch_boundary(&line((0.0, 0.1), (0.0, 0.9))),
        Some(PatchBoundary::U0)
    );
    assert_eq!(
        pcurve_patch_boundary(&line((0.2, 1.0), (0.7, 1.0))),
        Some(PatchBoundary::V1)
    );
    // Off the boundary by one ulp-scale step, or crossing the patch.
    assert_eq!(
        pcurve_patch_boundary(&line((1e-16, 0.1), (1e-16, 0.9))),
        None
    );
    assert_eq!(pcurve_patch_boundary(&line((0.0, 0.1), (1e-9, 0.9))), None);
    assert_eq!(pcurve_patch_boundary(&line((0.5, 0.1), (0.5, 0.9))), None);

    let on_u1 = [(1.0, 0.0), (1.0, 0.4), (1.0, 1.0)];
    assert_eq!(
        pcurve_patch_boundary(&nurbs_pcurve(&on_u1, vec![1.0, 2.0, 1.0])),
        Some(PatchBoundary::U1)
    );
    // A non-positive weight breaks the convex-combination argument.
    assert_eq!(
        pcurve_patch_boundary(&nurbs_pcurve(&on_u1, vec![1.0, -0.5, 1.0])),
        None
    );
    let bent = [(1.0, 0.0), (1.0 - 1e-12, 0.4), (1.0, 1.0)];
    assert_eq!(
        pcurve_patch_boundary(&nurbs_pcurve(&bent, vec![1.0; 3])),
        None
    );
}

#[test]
fn the_surrogate_applies_only_to_activated_dual_jet_patches() {
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    for (i, row) in grid.iter_mut().enumerate() {
        for (j, p) in row.iter_mut().enumerate() {
            *p = Point3::new(i as f64 / 9.0, j as f64 / 9.0, 0.1 * (i * j) as f64 / 81.0);
        }
    }
    let plain = RectangularGregoryBiNonic { grid };
    let boundary = Pcurve::linear((0.0, 0.0), (0.0, 1.0), 0.0, 1.0);
    let inactive = RectangularGregoryBiNonicDualJet::from_plain(&plain);
    let surface = |p: &RectangularGregoryBiNonicDualJet| {
        SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonicDualJet(Box::new(p.clone())))
    };
    assert!(boundary_surrogate_surface(&surface(&inactive), &boundary).is_none());
    let mut active = inactive.clone();
    active.activate_corner(0);
    assert!(boundary_surrogate_surface(&surface(&active), &boundary).is_some());
    let interior = Pcurve::linear((0.3, 0.0), (0.3, 1.0), 0.0, 1.0);
    assert!(boundary_surrogate_surface(&surface(&active), &interior).is_none());
    let plain_surface = SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonic(Box::new(plain)));
    assert!(boundary_surrogate_surface(&plain_surface, &boundary).is_none());
}
