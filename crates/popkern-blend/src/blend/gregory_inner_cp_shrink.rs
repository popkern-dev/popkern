//! Adaptive shrinkage of Chiyokura-Kimura inner-CP pairs for
//! narrow pie-slice Gregory sub-patches.
//!
//! At high valence (N ≥ 5 vertex blends), the classical raw
//! parallelogram construction `P21 = b_prev[m-1] + b_curr[1] −
//! corner` places inner CPs that overshoot the sub-patch's
//! angular sector, producing self-intersecting surfaces. This
//! module linearly shrinks each inner CP toward its corner
//! until the resulting patch is fold-free.
//!
//! **Tradeoff:** shrinking breaks the exact zero-twist
//! constraint at the corner, slightly degrading G1 at the
//! OUTER corners of the sub-patch. The spoke-G1 solver at
//! the central vertex plus downstream solver rungs partially
//! recover quality. Outer-corner G1 enforcement is later solver work.
//!
//! **Limitation:** shrinkage CANNOT fix folds that originate
//! from the boundary curves themselves. At genuinely degenerate
//! pie slices (very high N and/or very tight radii), even
//! fully-collapsed inner CPs (`factor = 0`) produce folded
//! boundary-driven surfaces. Those cases remain `no-Gregory` in
//! the pipeline's polar-NURBS fallback.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{GregoryPatch, RectangularGregory, RectangularGregoryDeg5};

pub(super) fn shrink_inner_cps_if_folded_cubic(
    boundaries: &[[Point3; 4]; 4],
    corners: &[Point3; 4],
    inner_cps: [(Point3, Point3); 4],
) -> [(Point3, Point3); 4] {
    let mut current = inner_cps;
    let factors = [1.0, 0.7, 0.5, 0.35, 0.2];
    for &f in &factors {
        let attempt = scale_inner_cps(corners, &inner_cps, f);
        if let Ok(patch) = RectangularGregory::from_boundary_beziers(*boundaries, attempt, 1e-3) {
            let gp = GregoryPatch::Rectangular(Box::new(patch));
            if !gp.self_intersects_with_tolerance(1e-6) {
                return attempt;
            }
        }
        current = attempt;
    }
    current
}

pub(super) fn shrink_inner_cps_if_folded_deg5(
    boundaries: &[[Point3; 6]; 4],
    corners: &[Point3; 4],
    inner_cps: [(Point3, Point3); 4],
) -> [(Point3, Point3); 4] {
    let mut current = inner_cps;
    let factors = [1.0, 0.7, 0.5, 0.35, 0.2, 0.1, 0.05];
    for &f in &factors {
        let attempt = scale_inner_cps(corners, &inner_cps, f);
        if let Ok(patch) = RectangularGregoryDeg5::from_boundary_beziers(*boundaries, attempt, 1e-3)
        {
            let gp = GregoryPatch::RectangularDeg5(Box::new(patch));
            if !gp.self_intersects_with_tolerance(1e-6) {
                return attempt;
            }
        }
        current = attempt;
    }
    current
}

/// Move each inner CP (both `P_u` and `P_v` of the pair) toward
/// its corresponding corner by `(1 − factor)`: `factor = 1` = no
/// shrink, `factor = 0` = collapse to corner.
fn scale_inner_cps(
    corners: &[Point3; 4],
    inner_cps: &[(Point3, Point3); 4],
    factor: f64,
) -> [(Point3, Point3); 4] {
    let mut out = *inner_cps;
    for (i, slot) in out.iter_mut().enumerate() {
        let c = corners[i];
        let (pu, pv) = inner_cps[i];
        slot.0 = Point3::new(
            c.x + factor * (pu.x - c.x),
            c.y + factor * (pu.y - c.y),
            c.z + factor * (pu.z - c.z),
        );
        slot.1 = Point3::new(
            c.x + factor * (pv.x - c.x),
            c.y + factor * (pv.y - c.y),
            c.z + factor * (pv.z - c.z),
        );
    }
    out
}
