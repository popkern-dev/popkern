//! `chamfer_2d` curved-boundary tests split out of
//! `chamfer_2d_tests.rs` so each test file stays under the 500-
//! line structural cap.

use super::*;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Chamfer the `v0` corner of a `sheet_arched_triangle` — one
/// adjacent edge is the straight diameter, the other is a
/// quarter-circle arc. Expect:
/// - The corner vertex is deleted.
/// - Two new tangent-point vertices are inserted.
/// - The diameter edge is shortened.
/// - The arc edge is clipped (its `start_angle` / `end_angle`
///   updates to match the new tangent point on the same circle).
/// - A new straight-line chamfer edge spans the two tangent points.
/// - Full BRep validation passes.
#[test]
fn chamfer_corner_2d_on_line_arc_corner_of_arched_triangle() {
    let r = 1.0_f64;
    let sheet = BRep::sheet_arched_triangle(r, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| (v.position.x - (-r)).abs() < 1e-9 && v.position.y.abs() < 1e-9)
        .map(|(h, _)| h)
        .unwrap();

    let result =
        chamfer_corner_2d(&sheet, corner, 0.1, 0.1).expect("line-arc chamfer should succeed");

    assert_eq!(result.vertices.len(), 4);
    assert_eq!(result.edges.len(), 4);
    assert_eq!(result.faces.len(), 2);
    let euler =
        result.vertices.len() as i64 - result.edges.len() as i64 + result.faces.len() as i64;
    assert_eq!(euler, 2);
    result.validate().unwrap();
}

/// The chamfer's new edge is a straight `Line`; the diameter
/// edge is still a `Line` (clipped); the arc edge is still an
/// `Arc` (clipped) — total counts: 2 Arcs + 2 Lines.
#[test]
fn chamfer_corner_2d_on_arched_triangle_preserves_arc_curve_type() {
    let r = 1.0_f64;
    let sheet = BRep::sheet_arched_triangle(r, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| (v.position.x - (-r)).abs() < 1e-9 && v.position.y.abs() < 1e-9)
        .map(|(h, _)| h)
        .unwrap();
    let result = chamfer_corner_2d(&sheet, corner, 0.15, 0.15).unwrap();

    let arc_count = result
        .edges
        .values()
        .filter(|e| {
            matches!(
                e.curve.as_ref().map(|c| &c.geometry),
                Some(CurveGeometry::Arc { .. })
            )
        })
        .count();
    let line_count = result
        .edges
        .values()
        .filter(|e| {
            matches!(
                e.curve.as_ref().map(|c| &c.geometry),
                Some(CurveGeometry::Line { .. })
            )
        })
        .count();
    assert_eq!(arc_count, 2, "expected 2 arc edges, got {arc_count}");
    assert_eq!(line_count, 2, "expected 2 line edges, got {line_count}");
}

/// The clipped arc edge's curve still evaluates to points on the
/// original circle.
#[test]
fn chamfer_corner_2d_clipped_arc_stays_on_circle() {
    use popkern_core::geometry::curve::CurveEval;
    let r = 1.0_f64;
    let sheet = BRep::sheet_arched_triangle(r, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| (v.position.x - (-r)).abs() < 1e-9 && v.position.y.abs() < 1e-9)
        .map(|(h, _)| h)
        .unwrap();
    let result = chamfer_corner_2d(&sheet, corner, 0.1, 0.1).unwrap();

    for edge in result.edges.values() {
        let curve = edge.curve.as_ref().unwrap();
        if matches!(&curve.geometry, CurveGeometry::Arc { .. }) {
            let t_mid = 0.5 * (curve.t_start + curve.t_end);
            let p = curve.geometry.point_at(t_mid);
            let d_from_origin = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
            assert!(
                (d_from_origin - r).abs() < 1e-9,
                "clipped arc midpoint not on circle r={r}: |{p:?}| = {d_from_origin}"
            );
        }
    }
}

/// Oversized chamfer distance on an arc-adjacent corner returns
/// a proper "too large" error. Arc length of the v0-v1 arc is
/// `r · π/2 ≈ 1.571` for `r = 1.0`.
#[test]
fn chamfer_corner_2d_rejects_oversize_arc_distance() {
    let sheet = BRep::sheet_arched_triangle(1.0, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| (v.position.x - (-1.0)).abs() < 1e-9 && v.position.y.abs() < 1e-9)
        .map(|(h, _)| h)
        .unwrap();
    let err = chamfer_corner_2d(&sheet, corner, 2.0, 0.1).unwrap_err();
    assert!(matches!(err, BlendError::PartialFailure { .. }));
}

// ---------------------------------------------------------------------------
// Arc-Arc corner (apex of sheet_biased_arched_triangle)
// ---------------------------------------------------------------------------

/// Chamfer the arc-arc apex of a `sheet_biased_arched_triangle`
/// with distinct radii (so the apex is a real corner). The
/// chamfer path is already curve-type-agnostic via
/// `tangent_at_corner_outgoing` / `point_at_arc_length_from_corner`,
/// so this primarily regression-guards that the Arc-Arc combo
/// flows through without Line-only assertions or panics.
#[test]
fn chamfer_corner_2d_on_arc_arc_apex_of_biased() {
    let sheet = BRep::sheet_biased_arched_triangle(1.5, 1.0, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x.abs() < 1e-9 && (v.position.y - 1.5).abs() < 1e-9)
        .map(|(h, _)| h)
        .unwrap();

    let result = chamfer_corner_2d(&sheet, corner, 0.1, 0.1)
        .expect("arc-arc chamfer should succeed on biased apex");

    assert_eq!(result.vertices.len(), 4);
    assert_eq!(result.edges.len(), 4);
    assert_eq!(result.faces.len(), 2);
    let euler =
        result.vertices.len() as i64 - result.edges.len() as i64 + result.faces.len() as i64;
    assert_eq!(euler, 2);
    result.validate().unwrap();
}

/// After arc-arc chamfer, the new edge is a straight `Line`,
/// the two support arcs are still `Arc`s (with clipped angles),
/// and the diameter stays a `Line` — total 2 `Arc`s + 2 `Line`s.
#[test]
fn chamfer_corner_2d_arc_arc_curve_census() {
    let sheet = BRep::sheet_biased_arched_triangle(1.5, 1.0, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x.abs() < 1e-9 && (v.position.y - 1.5).abs() < 1e-9)
        .map(|(h, _)| h)
        .unwrap();
    let result = chamfer_corner_2d(&sheet, corner, 0.15, 0.15).unwrap();

    let arc_count = result
        .edges
        .values()
        .filter(|e| {
            matches!(
                e.curve.as_ref().map(|c| &c.geometry),
                Some(CurveGeometry::Arc { .. })
            )
        })
        .count();
    let line_count = result
        .edges
        .values()
        .filter(|e| {
            matches!(
                e.curve.as_ref().map(|c| &c.geometry),
                Some(CurveGeometry::Line { .. })
            )
        })
        .count();
    assert_eq!(arc_count, 2, "expected 2 arc edges, got {arc_count}");
    assert_eq!(line_count, 2, "expected 2 line edges, got {line_count}");
}
