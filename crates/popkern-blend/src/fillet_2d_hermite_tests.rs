//! Hermite corner blends (PK-508): a sheet square's corner replaced by
//! a G1 / G2 / G3 polynomial blend at given setbacks and speeds.

use super::*;
use crate::curve_join::{JoinTolerances, analyze_join};
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::tolerance::Tolerance;

fn square_corner() -> (BRep, VertexHandle) {
    let sheet = BRep::sheet_quad(0.0, 1.0, 0.0, 1.0, 0.0, Tolerance::default()).expect("sheet");
    let corner = sheet
        .vertices
        .iter()
        .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0)
        .map(|(h, _)| h)
        .expect("corner");
    (sheet, corner)
}

fn blend(order: u8) -> HermiteCorner2d {
    HermiteCorner2d {
        order,
        setback_a: 0.2,
        setback_b: 0.3,
        speed_a: 0.3,
        speed_b: 0.4,
    }
}

const TOL: JoinTolerances = JoinTolerances {
    gap: 1e-9,
    angle: 1e-9,
    curvature: 1e-9,
    curvature_derivative: 1e-7,
};

/// The blend edge of `out` (its one NURBS edge) and the continuity it
/// reaches with each straight support at its ends.
fn continuity_at_both_ends(out: &BRep) -> (usize, i32, i32) {
    let curve = out
        .edges
        .values()
        .find_map(|e| match e.curve.as_ref().map(|c| &c.geometry) {
            Some(CurveGeometry::Nurbs(n)) => Some(n.clone()),
            _ => None,
        })
        .expect("the blend edge");
    let geometry = CurveGeometry::Nurbs(curve.clone());
    let (start, end) = (JoinJet::of(&geometry, 0.0), JoinJet::of(&geometry, 1.0));
    // A straight support arriving at a contact along the blend's own
    // tangent (supports on the axes): its jet is the unit tangent and
    // zero higher derivatives.
    let line = |at: Point3, d: Vec3| JoinJet {
        point: at,
        d: [
            d.normalized().expect("tangent"),
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 0.0),
        ],
    };
    let into = analyze_join(&line(start.point, start.d[0]), &start, &TOL).continuity;
    let out_of = analyze_join(&end, &line(end.point, end.d[0]), &TOL).continuity;
    (curve.degree, into, out_of)
}

#[test]
fn each_order_matches_its_supports_through_that_order() {
    let (sheet, corner) = square_corner();
    for order in 1u8..=3 {
        let out = hermite_corner_2d(&sheet, corner, &blend(order)).expect("blend");
        assert!(out.validate().is_ok());
        assert_eq!(
            out.vertices.len() as i64 - out.edges.len() as i64 + out.faces.len() as i64,
            2
        );
        let (degree, a, b) = continuity_at_both_ends(&out);
        assert_eq!(degree, 2 * usize::from(order) + 1);
        assert!(
            a >= i32::from(order) && b >= i32::from(order),
            "G{order}: {a}, {b}"
        );
    }
}

#[test]
fn the_contacts_sit_at_their_setbacks_and_leave_at_their_speeds() {
    let (sheet, corner) = square_corner();
    let out = hermite_corner_2d(&sheet, corner, &blend(2)).expect("blend");
    let curve = out
        .edges
        .values()
        .find_map(|e| match e.curve.as_ref().map(|c| &c.geometry) {
            Some(CurveGeometry::Nurbs(n)) => Some(n.clone()),
            _ => None,
        })
        .expect("blend");
    let d0 = curve.derivatives(0.0, 1);
    let d1 = curve.derivatives(1.0, 1);
    let (p0, p1) = (d0[0], d1[0]);
    // One contact 0.2 along one axis, the other 0.3 along the other.
    let mut setbacks = [p0.x.max(p0.y), p1.x.max(p1.y)];
    setbacks.sort_by(f64::total_cmp);
    assert!(
        (setbacks[0] - 0.2).abs() < 1e-9 && (setbacks[1] - 0.3).abs() < 1e-9,
        "{setbacks:?}"
    );
    let mut speeds = [d0[1].length(), d1[1].length()];
    speeds.sort_by(f64::total_cmp);
    assert!(
        (speeds[0] - 0.3).abs() < 1e-9 && (speeds[1] - 0.4).abs() < 1e-9,
        "{speeds:?}"
    );
}

#[test]
fn an_overdriven_blend_declines_before_any_surgery() {
    let (sheet, corner) = square_corner();
    // Speeds far beyond the setbacks loop the blend back over itself.
    let wild = HermiteCorner2d {
        speed_a: 40.0,
        speed_b: 40.0,
        ..blend(3)
    };
    let err = hermite_corner_2d(&sheet, corner, &wild).expect_err("declines");
    let text = format!("{err:?}");
    assert!(
        text.contains("regular") || text.contains("crosses"),
        "{text}"
    );
}

#[test]
fn bad_parameters_decline_by_name() {
    let (sheet, corner) = square_corner();
    for (b, name) in [
        (
            HermiteCorner2d {
                order: 4,
                ..blend(1)
            },
            "order",
        ),
        (
            HermiteCorner2d {
                speed_a: 0.0,
                ..blend(1)
            },
            "speed_a",
        ),
    ] {
        let err = hermite_corner_2d(&sheet, corner, &b).expect_err("declines");
        assert!(format!("{err:?}").contains(name), "{err:?}");
    }
    let long = HermiteCorner2d {
        setback_b: 2.0,
        ..blend(1)
    };
    assert!(
        format!(
            "{:?}",
            hermite_corner_2d(&sheet, corner, &long).expect_err("too long")
        )
        .contains("reaches past")
    );
}
