//! Tests for `chamfer_corner_2d_by_angle` — the §5.7
//! distance+angle chamfer mode. Covers the classic line-line
//! path (regression-guarding law-of-sines equivalence) and the
//! curved-support extensions (line-arc, arc-arc).

use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

// ---------------------------------------------------------------------------
// Line-Line: regression guard — the new ray-intersection path
// must still agree with the historical law-of-sines shortcut.
// ---------------------------------------------------------------------------

/// 90° corner with `angle_rad = 3π/4` gives a symmetric 0.1 ×
/// 0.1 chamfer on a `sheet_quad`. Valid topology + counts.
#[test]
fn chamfer_by_angle_symmetric_on_quad_90deg() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0)
        .map(|(h, _)| h)
        .unwrap();

    let result = chamfer_corner_2d_by_angle(&sheet, corner, 0.1, 3.0 * std::f64::consts::FRAC_PI_4)
        .expect("line-line 90° symmetric chamfer should succeed");
    assert_eq!(result.vertices.len(), 5);
    assert_eq!(result.edges.len(), 5);
    result.validate().unwrap();

    // Symmetric: the chamfer line goes from (0.1, 0, 0) to
    // (0, 0.1, 0). Verify both tangent-point vertices exist.
    let has_t1 = result.vertices.values().any(|v| {
        (v.position.x - 0.1).abs() < 1e-9 && v.position.y.abs() < 1e-9 && v.position.z.abs() < 1e-9
    });
    let has_t2 = result.vertices.values().any(|v| {
        v.position.x.abs() < 1e-9 && (v.position.y - 0.1).abs() < 1e-9 && v.position.z.abs() < 1e-9
    });
    assert!(has_t1, "missing tangent point at (0.1, 0, 0)");
    assert!(has_t2, "missing tangent point at (0, 0.1, 0)");
}

/// `angle_rad` at or below the corner angle is rejected as
/// `InvalidTopology` (the chamfer ray wouldn't reach `edge_b`
/// forward from `t1`).
#[test]
fn chamfer_by_angle_rejects_angle_below_corner_angle() {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, tol()).unwrap();
    let corner = sheet.vertices.keys().next().unwrap();
    // Corner angle on a quad is π/2. Choose angle_rad = π/4
    // which is below the corner angle — should reject.
    let err =
        chamfer_corner_2d_by_angle(&sheet, corner, 0.1, std::f64::consts::FRAC_PI_4).unwrap_err();
    assert!(matches!(err, BlendError::PartialFailure { .. }));
}

// ---------------------------------------------------------------------------
// Line-Arc (v0 corner of sheet_arched_triangle)
// ---------------------------------------------------------------------------

/// Distance+angle chamfer at the v0 line-arc corner of a
/// `sheet_arched_triangle`. `angle_rad = 3π/4` on a 90° corner
/// is the symmetric-looking chamfer; the curved-support path
/// computes `t2` via ray-arc intersection rather than law of
/// sines, and the result must still validate.
#[test]
fn chamfer_by_angle_on_line_arc_corner() {
    let r = 1.0_f64;
    let sheet = BRep::sheet_arched_triangle(r, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| (v.position.x - (-r)).abs() < 1e-9 && v.position.y.abs() < 1e-9)
        .map(|(h, _)| h)
        .unwrap();

    let result = chamfer_corner_2d_by_angle(&sheet, corner, 0.1, 3.0 * std::f64::consts::FRAC_PI_4)
        .expect("line-arc by-angle chamfer should succeed");
    assert_eq!(result.vertices.len(), 4);
    assert_eq!(result.edges.len(), 4);
    assert_eq!(result.faces.len(), 2);
    let euler =
        result.vertices.len() as i64 - result.edges.len() as i64 + result.faces.len() as i64;
    assert_eq!(euler, 2);
    result.validate().unwrap();
}

/// The ray-arc intersection must leave the arc edge still on
/// the original circle: after clipping, every arc edge's
/// midpoint must sit on `|p| = r`.
#[test]
fn chamfer_by_angle_clipped_arc_stays_on_circle() {
    use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
    let r = 1.0_f64;
    let sheet = BRep::sheet_arched_triangle(r, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| (v.position.x - (-r)).abs() < 1e-9 && v.position.y.abs() < 1e-9)
        .map(|(h, _)| h)
        .unwrap();
    let result =
        chamfer_corner_2d_by_angle(&sheet, corner, 0.1, 3.0 * std::f64::consts::FRAC_PI_4).unwrap();

    for edge in result.edges.values() {
        let curve = edge.curve.as_ref().unwrap();
        if matches!(&curve.geometry, CurveGeometry::Arc { .. }) {
            let t_mid = 0.5 * (curve.t_start + curve.t_end);
            let p = curve.geometry.point_at(t_mid);
            let d = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
            assert!(
                (d - r).abs() < 1e-9,
                "clipped arc midpoint off unit circle: |{p:?}| = {d}"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Arc-Arc (apex of sheet_biased_arched_triangle)
// ---------------------------------------------------------------------------

/// Distance+angle chamfer at the arc-arc apex of a
/// `sheet_biased_arched_triangle`. With radii `(1.5, 1.0)` the
/// apex corner angle measured between the two arc tangents is
/// very obtuse (~168.6°) — the arcs bulge outward hard and leave
/// a nearly-straight-through corner. A valid chamfer therefore
/// needs `angle_rad` in `(168.6°, 180°)`, with a small `dist_a`
/// that still fits the narrow wedge. This exercises the
/// ray-arc intersection path on *both* supports and regression-
/// guards the arc-arc by-angle case.
#[test]
fn chamfer_by_angle_on_arc_arc_apex() {
    let sheet = BRep::sheet_biased_arched_triangle(1.5, 1.0, 0.0, tol()).unwrap();
    let apex = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x.abs() < 1e-9 && (v.position.y - 1.5).abs() < 1e-9)
        .map(|(h, _)| h)
        .unwrap();

    let angle = 178.0_f64.to_radians();
    let result = chamfer_corner_2d_by_angle(&sheet, apex, 0.01, angle)
        .expect("arc-arc by-angle chamfer should succeed on the biased apex");
    assert_eq!(result.vertices.len(), 4);
    assert_eq!(result.edges.len(), 4);
    result.validate().unwrap();
}

/// On the same arc-arc apex, an `angle_rad` well below the
/// corner angle (the chamfer ray would go away from `edge_b`)
/// must be rejected.
#[test]
fn chamfer_by_angle_rejects_below_corner_angle_arc_arc() {
    let sheet = BRep::sheet_biased_arched_triangle(1.5, 1.0, 0.0, tol()).unwrap();
    let apex = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x.abs() < 1e-9 && (v.position.y - 1.5).abs() < 1e-9)
        .map(|(h, _)| h)
        .unwrap();
    let err = chamfer_corner_2d_by_angle(&sheet, apex, 0.01, 30.0_f64.to_radians()).unwrap_err();
    assert!(matches!(err, BlendError::PartialFailure { .. }));
}
