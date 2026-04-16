//! `fillet_2d` curved-boundary tests split out of
//! `fillet_2d_tests.rs` so each test file stays under the 500-
//! line structural cap and the Line-Line and Line-Arc code
//! paths exercise independent test vehicles.

use super::*;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Fillet the `v0` corner of a `sheet_arched_triangle` — one
/// adjacent edge is the straight diameter, the other is a
/// quarter-circle arc. Expect:
/// - The corner vertex is deleted.
/// - Two tangent-point vertices are inserted.
/// - The diameter edge is shortened (still a `Line`).
/// - The arc edge is clipped (still an `Arc`, but with new
///   `start_angle`/`end_angle` on the same circle).
/// - A fresh fillet `Arc` edge spans the two tangent points,
///   tangent to both supports.
/// - Full BRep validation passes, `V − E + F = 2`.
#[test]
fn fillet_corner_2d_on_line_arc_corner_of_arched_triangle() {
    let r = 1.0_f64;
    let sheet = BRep::sheet_arched_triangle(r, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| (v.position.x - (-r)).abs() < 1e-9 && v.position.y.abs() < 1e-9)
        .map(|(h, _)| h)
        .unwrap();

    let result = fillet_corner_2d(&sheet, corner, 0.1).expect("line-arc fillet should succeed");

    assert_eq!(result.vertices.len(), 4);
    assert_eq!(result.edges.len(), 4);
    assert_eq!(result.faces.len(), 2);
    let euler =
        result.vertices.len() as i64 - result.edges.len() as i64 + result.faces.len() as i64;
    assert_eq!(euler, 2);
    result.validate().unwrap();
}

/// After filleting a Line-Arc corner the curve-kind census is:
/// - 2 `Arc` edges (the clipped original + the new fillet arc)
/// - 2 `Line` edges (the shortened diameter + the other arc
///   edge… wait no, the other arc on the sheet_arched_triangle is
///   *also* an `Arc`). Actual census: 3 Arcs + 1 Line.
#[test]
fn fillet_corner_2d_on_arched_triangle_preserves_arc_curve_types() {
    let r = 1.0_f64;
    let sheet = BRep::sheet_arched_triangle(r, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| (v.position.x - (-r)).abs() < 1e-9 && v.position.y.abs() < 1e-9)
        .map(|(h, _)| h)
        .unwrap();
    let result = fillet_corner_2d(&sheet, corner, 0.15).unwrap();

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
    assert_eq!(arc_count, 3, "expected 3 arc edges, got {arc_count}");
    assert_eq!(line_count, 1, "expected 1 line edge, got {line_count}");
}

/// The clipped half of the original sheet arc still has every
/// evaluated point on the unit circle. This regression-guards
/// `rebuild_edge_curve_clipped` against frame / wrap drift.
#[test]
fn fillet_corner_2d_clipped_arc_stays_on_circle() {
    let r = 1.0_f64;
    let sheet = BRep::sheet_arched_triangle(r, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| (v.position.x - (-r)).abs() < 1e-9 && v.position.y.abs() < 1e-9)
        .map(|(h, _)| h)
        .unwrap();
    let result = fillet_corner_2d(&sheet, corner, 0.1).unwrap();

    // Sample every arc edge at its midpoint; the new fillet arc
    // has radius 0.1 so it's *not* on the unit circle, but the
    // two remaining halves of the original arcs should be.
    for edge in result.edges.values() {
        let curve = edge.curve.as_ref().unwrap();
        if let CurveGeometry::Arc { radius, .. } = &curve.geometry {
            if (*radius - r).abs() > 1e-9 {
                // This is the fillet arc, not the original circle.
                continue;
            }
            let t_mid = 0.5 * (curve.t_start + curve.t_end);
            let p = curve.geometry.point_at(t_mid);
            let d_from_origin = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
            assert!(
                (d_from_origin - r).abs() < 1e-9,
                "clipped original arc midpoint not on unit circle: |{p:?}| = {d_from_origin}"
            );
        }
    }
}

/// The new fillet arc's radius and center should match the
/// theoretical values for a line-arc fillet. On the
/// `sheet_arched_triangle` v0 corner with fillet radius 0.1:
///
/// - Interior half-plane is `y ≥ 0`; line tangent at corner is
///   `(1, 0, 0)`, arc tangent at corner is `(0, 1, 0)`.
/// - Interior normal to the line = `(0, 1, 0)`.
/// - Offset line at `y = 0.1`, offset circle (internal tangency)
///   of radius `R − r = 0.9` around origin.
/// - Fillet center `F = (−√(1 − 2r), r, 0)` ≈ `(−0.894427…, 0.1, 0)`.
#[test]
fn fillet_corner_2d_on_arched_triangle_fillet_center_matches_theory() {
    let r_arc = 1.0_f64;
    let r_fillet = 0.1_f64;
    let sheet = BRep::sheet_arched_triangle(r_arc, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| (v.position.x - (-r_arc)).abs() < 1e-9 && v.position.y.abs() < 1e-9)
        .map(|(h, _)| h)
        .unwrap();
    let result = fillet_corner_2d(&sheet, corner, r_fillet).unwrap();

    // Find the fillet arc edge (the one with radius 0.1).
    let fillet_arc = result
        .edges
        .values()
        .find_map(|e| match e.curve.as_ref().map(|c| &c.geometry) {
            Some(CurveGeometry::Arc { center, radius, .. })
                if (*radius - r_fillet).abs() < 1e-9 =>
            {
                Some(*center)
            }
            _ => None,
        })
        .expect("missing fillet arc edge with radius 0.1");

    let expected_x = -(1.0 - 2.0 * r_fillet).sqrt();
    let expected_y = r_fillet;
    assert!(
        (fillet_arc.x - expected_x).abs() < 1e-9,
        "fillet center x: got {}, expected {}",
        fillet_arc.x,
        expected_x
    );
    assert!(
        (fillet_arc.y - expected_y).abs() < 1e-9,
        "fillet center y: got {}, expected {}",
        fillet_arc.y,
        expected_y
    );
    assert!(
        fillet_arc.z.abs() < 1e-9,
        "fillet center z: got {}, expected 0",
        fillet_arc.z
    );
}

/// Oversized fillet radius on an arc-adjacent corner returns a
/// proper "too large" error wrapped in `PartialFailure` ⇒
/// `InvalidTopology`.
///
/// A fillet radius of 1.0 on the unit-circle arched triangle is
/// geometrically impossible: the required internal-tangency
/// offset radius would be `R − r = 0`, and the solver rejects
/// it as "too large for internal tangency".
#[test]
fn fillet_corner_2d_rejects_oversize_on_arched_triangle() {
    let sheet = BRep::sheet_arched_triangle(1.0, 0.0, tol()).unwrap();
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| (v.position.x - (-1.0)).abs() < 1e-9 && v.position.y.abs() < 1e-9)
        .map(|(h, _)| h)
        .unwrap();
    let err = fillet_corner_2d(&sheet, corner, 1.0).unwrap_err();
    match err {
        BlendError::PartialFailure { failures, .. } => {
            let (_, f) = &failures[0];
            assert!(
                matches!(
                    f,
                    crate::BlendFailure::InvalidTopology { reason, .. }
                        if reason.contains("too large")
                ),
                "expected 'too large' InvalidTopology, got {f:?}"
            );
        }
        other => panic!("expected PartialFailure, got {other:?}"),
    }
}

// ---------------------------------------------------------------------------
// Arc-Arc corner (v1 apex of sheet_biased_arched_triangle)
// ---------------------------------------------------------------------------

/// Find the apex vertex `(0, 1.5, 0)` of a
/// `sheet_biased_arched_triangle`, where the two different-radius
/// arcs meet.
fn apex_of_biased(sheet: &BRep) -> popkern_core::topology::VertexHandle {
    sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x.abs() < 1e-9 && (v.position.y - 1.5).abs() < 1e-9)
        .map(|(h, _)| h)
        .unwrap()
}

/// Fillet the arc-arc apex of a `sheet_biased_arched_triangle`.
/// Expect 4 vertices, 4 edges, 2 faces (one corner replaced by
/// tangent-points + arc edge, Euler identity preserved).
#[test]
fn fillet_corner_2d_on_arc_arc_apex_of_biased() {
    let sheet = BRep::sheet_biased_arched_triangle(1.5, 1.0, 0.0, tol()).unwrap();
    let corner = apex_of_biased(&sheet);

    let result = fillet_corner_2d(&sheet, corner, 0.1).expect("arc-arc fillet should succeed");

    assert_eq!(result.vertices.len(), 4);
    assert_eq!(result.edges.len(), 4);
    assert_eq!(result.faces.len(), 2);
    let euler =
        result.vertices.len() as i64 - result.edges.len() as i64 + result.faces.len() as i64;
    assert_eq!(euler, 2);
    result.validate().unwrap();
}

/// Curve-kind census after arc-arc fillet: 3 arcs (2 clipped
/// originals + 1 fillet arc) and 1 line (the diameter).
#[test]
fn fillet_corner_2d_arc_arc_curve_census() {
    let sheet = BRep::sheet_biased_arched_triangle(1.5, 1.0, 0.0, tol()).unwrap();
    let corner = apex_of_biased(&sheet);
    let result = fillet_corner_2d(&sheet, corner, 0.1).unwrap();

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
    assert_eq!(arc_count, 3, "expected 3 arc edges (2 clipped + 1 fillet)");
    assert_eq!(line_count, 1, "expected 1 line edge (diameter)");
}

/// The fillet arc's centre must be at distance `R_a − r_f` from
/// `arc_a_center` and `R_b − r_f` from `arc_b_center` (internal
/// tangency on both arcs, which is correct for the convex apex
/// of the biased arched triangle).
#[test]
fn fillet_corner_2d_arc_arc_center_on_offset_circles() {
    let r_a = 1.5_f64;
    let r_b = 1.0_f64;
    let r_fillet = 0.1_f64;
    let sheet = BRep::sheet_biased_arched_triangle(r_a, r_b, 0.0, tol()).unwrap();
    let corner = apex_of_biased(&sheet);
    let result = fillet_corner_2d(&sheet, corner, r_fillet).unwrap();

    let (fillet_center, fillet_radius) = result
        .edges
        .values()
        .find_map(|e| match e.curve.as_ref().map(|c| &c.geometry) {
            Some(CurveGeometry::Arc { center, radius, .. })
                if (*radius - r_fillet).abs() < 1e-9 =>
            {
                Some((*center, *radius))
            }
            _ => None,
        })
        .expect("missing fillet arc");

    // Centres of the two original arcs (must match the constructor).
    let (center_a_x, center_a_y, center_b_x, center_b_y) = {
        let chord_len = (1.0_f64 + 1.5 * 1.5).sqrt();
        let offset_a = (r_a * r_a - 0.25 * (1.0 + 1.5 * 1.5)).sqrt();
        let offset_b = (r_b * r_b - 0.25 * (1.0 + 1.5 * 1.5)).sqrt();
        (
            -0.5 + offset_a * 1.5 / chord_len,
            1.5 * 0.5 + offset_a * (-1.0) / chord_len,
            0.5 + offset_b * (-1.5) / chord_len,
            1.5 * 0.5 + offset_b * (-1.0) / chord_len,
        )
    };

    let d_to_a =
        ((fillet_center.x - center_a_x).powi(2) + (fillet_center.y - center_a_y).powi(2)).sqrt();
    let d_to_b =
        ((fillet_center.x - center_b_x).powi(2) + (fillet_center.y - center_b_y).powi(2)).sqrt();

    assert!(
        (d_to_a - (r_a - fillet_radius)).abs() < 1e-9,
        "|F - C_a| = {d_to_a}, expected {}",
        r_a - fillet_radius
    );
    assert!(
        (d_to_b - (r_b - fillet_radius)).abs() < 1e-9,
        "|F - C_b| = {d_to_b}, expected {}",
        r_b - fillet_radius
    );
}

/// Fillet radius larger than either support arc's radius makes
/// internal tangency impossible — the offset circle would have
/// radius `R − r ≤ 0`. Expected: rejected with a "too large"
/// `InvalidTopology`.
#[test]
fn fillet_corner_2d_rejects_oversize_on_biased_apex() {
    let sheet = BRep::sheet_biased_arched_triangle(1.5, 1.0, 0.0, tol()).unwrap();
    let corner = apex_of_biased(&sheet);
    // r_b = 1.0, so r_fillet = 1.1 gives internal offset −0.1 on
    // arc b — infeasible.
    let err = fillet_corner_2d(&sheet, corner, 1.1).unwrap_err();
    match err {
        BlendError::PartialFailure { failures, .. } => {
            let (_, f) = &failures[0];
            assert!(
                matches!(
                    f,
                    crate::BlendFailure::InvalidTopology { reason, .. }
                        if reason.contains("too large")
                ),
                "expected 'too large' InvalidTopology, got {f:?}"
            );
        }
        other => panic!("expected PartialFailure, got {other:?}"),
    }
}
