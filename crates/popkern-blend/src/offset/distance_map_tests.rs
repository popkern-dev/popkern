//! Distance-map self-fold detection gates (PK-327 tier 3): a
//! C-section sheet with converging tips folds GLOBALLY under an
//! into-the-C offset (every local curvature radius stays above the
//! distance, so tier 1 passes) and must report fold regions at the
//! tips; the away-from-the-C offset of the same sheet must scan
//! clean.

use super::distance_map::scan_offset_self_fold;
use super::offset_surface_measured;
use popkern_core::geometry::NurbsSurface;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, surface::SurfaceGeometry};
use popkern_core::tolerance::Tolerance;

/// Extruded C-section: bottom leg near z = 0 from x = 1 to the bend at
/// x ≈ −0.55, top leg returning to a tip at (1, 0.2) — tip gap 0.2,
/// bend radius ≈ 0.3, every local curvature radius above the offset.
fn c_sheet() -> SurfaceGeometry {
    let profile = [
        (1.0, 0.0),
        (0.2, 0.0),
        (-0.35, 0.05),
        (-0.55, 0.3),
        (-0.35, 0.55),
        (0.2, 0.6),
        (1.0, 0.2),
    ];
    // Storage is u-major (index = i_u * num_v + i_v).
    let mut control_points = Vec::with_capacity(14);
    for (x, z) in profile {
        for y in [0.0, 1.0] {
            control_points.push(Point3::new(x, y, z));
        }
    }
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 2,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points,
        num_u: 7,
        num_v: 2,
        weights: vec![1.0; 14],
    })
}

/// Offset the sheet toward / away from the C interior, decided by
/// measurement (the bottom-leg offset moving +z is into the C).
fn fitted_offset(source: &SurfaceGeometry, d: f64, into_c: bool) -> SurfaceGeometry {
    for inward in [false, true] {
        let offset = offset_surface_measured(
            source,
            d,
            inward,
            (0.0, 1.0),
            (0.0, 1.0),
            &Tolerance::default(),
        )
        .expect("offset fit")
        .surface;
        let s = source.point_at(0.05, 0.5);
        let o = offset.point_at(0.05, 0.5);
        if (o.z > s.z) == into_c {
            return offset;
        }
    }
    unreachable!("one of the two directions moves the bottom leg up");
}

#[test]
fn converging_tips_fold_globally_and_are_detected() {
    let source = c_sheet();
    let d = 0.12;
    let offset = fitted_offset(&source, d, true);
    let report = scan_offset_self_fold(&source, &offset, d, (0.0, 1.0), (0.0, 1.0));
    assert!(
        !report.regions.is_empty(),
        "tip fold must be detected: {report:?}"
    );
    // The fold localizes at the tips (u near 0 or 1, where the legs
    // end at x = 1).
    let near_tip = report.regions.iter().any(|r| {
        let mid = (r.r_range.0 + r.r_range.1) * 0.5;
        !(0.25..=0.75).contains(&mid)
    });
    assert!(
        near_tip,
        "fold regions must localize at the tips: {report:?}"
    );
    assert!(!report.truncated, "budget must suffice: {report:?}");
}

#[test]
fn away_from_the_c_offset_scans_clean() {
    let source = c_sheet();
    let d = 0.12;
    let offset = fitted_offset(&source, d, false);
    let report = scan_offset_self_fold(&source, &offset, d, (0.0, 1.0), (0.0, 1.0));
    assert!(
        report.regions.is_empty(),
        "outward offset must scan clean: {report:?}"
    );
}
