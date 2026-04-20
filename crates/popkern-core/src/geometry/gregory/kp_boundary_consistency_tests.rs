//! Tests for [`super::check_boundary_consistency`].

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::gregory::kp_guide_build::{GuideCornerData, build_guide_from_hermite};

fn simple_guide() -> KpGuidePatch {
    let corners = [
        GuideCornerData {
            position: Point3::new(0.0, 0.0, 0.0),
            tangent_u: Vec3::new(1.0, 0.1, 0.2),
            tangent_v: Vec3::new(0.05, 1.1, -0.15),
            twist: Vec3::new(0.2, -0.1, 0.3),
        },
        GuideCornerData {
            position: Point3::new(1.0, 0.1, 0.3),
            tangent_u: Vec3::new(1.1, -0.2, 0.25),
            tangent_v: Vec3::new(0.1, 1.0, 0.2),
            twist: Vec3::new(-0.15, 0.2, -0.1),
        },
        GuideCornerData {
            position: Point3::new(1.15, 1.2, 0.1),
            tangent_u: Vec3::new(0.9, 0.15, -0.1),
            tangent_v: Vec3::new(-0.1, 0.95, 0.25),
            twist: Vec3::new(0.1, 0.15, -0.2),
        },
        GuideCornerData {
            position: Point3::new(0.05, 1.1, -0.05),
            tangent_u: Vec3::new(1.05, -0.1, 0.15),
            tangent_v: Vec3::new(0.0, 1.05, -0.1),
            twist: Vec3::new(-0.05, -0.2, 0.15),
        },
    ];
    build_guide_from_hermite(corners)
}

/// Build a K-P G2-consistent boundary (same shape as the tests
/// in kp_row2_corner_tests).
fn consistent_boundary(guide: &KpGuidePatch, is_u_boundary: bool, at_end: bool) -> [Point3; 10] {
    let (ua, va, ub, vb) = match (is_u_boundary, at_end) {
        (true, false) => (0.0, 0.0, 0.0, 1.0),
        (true, true) => (1.0, 0.0, 1.0, 1.0),
        (false, false) => (0.0, 0.0, 1.0, 0.0),
        (false, true) => (0.0, 1.0, 1.0, 1.0),
    };
    let p_start = guide.point_at(ua, va);
    let p_end = guide.point_at(ub, vb);
    let (t_start, t_end, tt_start, tt_end) = if is_u_boundary {
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
    cps[1] = off(p_start, t_start, 1.0 / 9.0);
    cps[8] = off(p_end, t_end, -1.0 / 9.0);
    cps[2] = off2(p_start, t_start, 2.0 / 9.0, tt_start, 1.0 / 72.0);
    cps[7] = off2(p_end, t_end, -2.0 / 9.0, tt_end, 1.0 / 72.0);
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

fn off(p: Point3, v: Vec3, s: f64) -> Point3 {
    Point3::new(p.x + v.x * s, p.y + v.y * s, p.z + v.z * s)
}
fn off2(p: Point3, v1: Vec3, s1: f64, v2: Vec3, s2: f64) -> Point3 {
    Point3::new(
        p.x + v1.x * s1 + v2.x * s2,
        p.y + v1.y * s1 + v2.y * s2,
        p.z + v1.z * s1 + v2.z * s2,
    )
}

#[test]
fn k_p_consistent_boundaries_pass_at_machine_precision() {
    let guide = simple_guide();
    let row0 = consistent_boundary(&guide, true, false);
    let row9 = consistent_boundary(&guide, true, true);
    let col0 = consistent_boundary(&guide, false, false);
    let col9 = consistent_boundary(&guide, false, true);
    let report = check_boundary_consistency(&guide, &row0, &row9, &col0, &col9);
    assert!(
        report.is_consistent(1e-9),
        "consistent boundaries should have errors < 1e-9: {report:?}"
    );
    assert!(report.max_position_error < 1e-12);
    assert!(report.max_tangent_error < 1e-10);
    assert!(report.max_second_derivative_error < 1e-8);
}

#[test]
fn position_inconsistency_flagged() {
    let guide = simple_guide();
    let mut row0 = consistent_boundary(&guide, true, false);
    let row9 = consistent_boundary(&guide, true, true);
    let col0 = consistent_boundary(&guide, false, false);
    let col9 = consistent_boundary(&guide, false, true);
    // Inject a position error at the (0, 0) corner.
    row0[0].z += 0.5;
    let report = check_boundary_consistency(&guide, &row0, &row9, &col0, &col9);
    assert!(report.max_position_error > 0.1);
    assert!(!report.is_consistent(1e-3));
}

#[test]
fn tangent_inconsistency_flagged() {
    let guide = simple_guide();
    let mut row0 = consistent_boundary(&guide, true, false);
    let row9 = consistent_boundary(&guide, true, true);
    let col0 = consistent_boundary(&guide, false, false);
    let col9 = consistent_boundary(&guide, false, true);
    // Nudge P[1] so the boundary's tangent at P[0] differs from
    // guide's G_v(0, 0).
    row0[1].x += 0.2;
    let report = check_boundary_consistency(&guide, &row0, &row9, &col0, &col9);
    assert!(report.max_tangent_error > 1.0); // 9 × 0.2 = 1.8 magnitude
    assert!(!report.is_consistent(1e-3));
}

#[test]
fn curvature_inconsistency_flagged() {
    let guide = simple_guide();
    let row0 = consistent_boundary(&guide, true, false);
    let row9 = consistent_boundary(&guide, true, true);
    let mut col0 = consistent_boundary(&guide, false, false);
    let col9 = consistent_boundary(&guide, false, true);
    // Nudge P[2] of col0 → second-derivative mismatch at (0, 0).
    col0[2].y += 0.1;
    let report = check_boundary_consistency(&guide, &row0, &row9, &col0, &col9);
    assert!(report.max_second_derivative_error > 1.0); // 72 × 0.1
    assert!(!report.is_consistent(1e-3));
}

#[test]
fn is_consistent_matches_manual_threshold() {
    let guide = simple_guide();
    let row0 = consistent_boundary(&guide, true, false);
    let row9 = consistent_boundary(&guide, true, true);
    let col0 = consistent_boundary(&guide, false, false);
    let col9 = consistent_boundary(&guide, false, true);
    let report = check_boundary_consistency(&guide, &row0, &row9, &col0, &col9);
    assert!(report.is_consistent(1e-6));
    // Very tight tolerance — curvature errors are in 1e-10 range,
    // so still consistent.
    assert!(report.is_consistent(1e-7));
}

#[test]
fn per_corner_breakdown_localises_inconsistency() {
    // If only ONE corner's boundary is inconsistent, only that
    // corner's per-corner errors should be elevated.
    let guide = simple_guide();
    let row0 = consistent_boundary(&guide, true, false);
    let row9 = consistent_boundary(&guide, true, true);
    let mut col0 = consistent_boundary(&guide, false, false);
    let col9 = consistent_boundary(&guide, false, true);
    col0[1].z += 0.3; // Breaks tangent at corner 0 only.
    let report = check_boundary_consistency(&guide, &row0, &row9, &col0, &col9);
    assert!(report.per_corner[0].tangent_u_error > 1.0);
    // Other corners should stay clean.
    assert!(report.per_corner[2].tangent_u_error < 1e-6);
    assert!(report.per_corner[2].tangent_v_error < 1e-6);
}
