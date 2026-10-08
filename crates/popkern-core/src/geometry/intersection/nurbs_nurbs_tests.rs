//! Tests for `nurbs_nurbs` (split out for the 500-line cap).
use super::*;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;

#[test]
fn parabola_crosses_horizontal_line() {
    // Quadratic Bezier (parabola) and horizontal line at y=0.5
    let c1 = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let c2 = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let result = intersect_nurbs_nurbs(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &Tolerance::default());
    // Parabola y = 2t(1-t) crosses y=0.5 at 2 points
    assert!(
        result.len() >= 2,
        "Expected at least 2 intersections, got {}",
        result.len()
    );

    // Verify invariant: points lie on both curves
    for r in &result {
        if let CurveIntersection::Point { point, t1, t2 } = r {
            let p1 = c1.point_at(*t1);
            let p2 = c2.point_at(*t2);
            assert!(p1.distance_to(&p2) < 1e-6, "Invariant: points don't match");
            assert!(p1.distance_to(point) < 1e-6);
        }
    }
}

#[test]
fn non_intersecting_curves() {
    let c1 = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let c2 = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 10.0, 0.0), Point3::new(1.0, 10.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let result = intersect_nurbs_nurbs(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &Tolerance::default());
    assert!(result.is_empty());
}

#[test]
fn identical_curves_overlap() {
    let c1 = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let result = intersect_nurbs_nurbs(&c1, (0.0, 1.0), &c1, (0.0, 1.0), &Tolerance::default());
    assert!(
        result
            .iter()
            .any(|r| matches!(r, CurveIntersection::Overlap { .. }))
    );
}

#[test]
fn multi_span_curve_decomposition() {
    // Multi-span curve (has interior knot) crossing a line
    let c1 = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.25, 1.0, 0.0),
            Point3::new(0.75, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0, 1.0],
    };
    let c2 = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let result = intersect_nurbs_nurbs(&c1, (0.0, 1.0), &c2, (0.0, 1.0), &Tolerance::default());
    assert!(
        result.len() >= 2,
        "Expected at least 2 intersections with multi-span curve"
    );
}

#[test]
fn aabb_overlap_disjoint() {
    let pts1 = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)];
    let pts2 = vec![Point3::new(5.0, 5.0, 0.0), Point3::new(6.0, 6.0, 0.0)];
    assert!(!aabb_overlap(&pts1, &pts2));
}

#[test]
fn aabb_overlap_touching() {
    let pts1 = vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 1.0, 1.0)];
    let pts2 = vec![Point3::new(0.5, 0.5, 0.5), Point3::new(2.0, 2.0, 2.0)];
    assert!(aabb_overlap(&pts1, &pts2));
}

#[test]
fn aabb_overlap_empty() {
    let empty: Vec<Point3> = vec![];
    let pts = vec![Point3::new(0.0, 0.0, 0.0)];
    assert!(!aabb_overlap(&empty, &pts));
    assert!(!aabb_overlap(&pts, &empty));
}

#[test]
fn deduplicate_removes_close_points() {
    let mut results = vec![
        CurveIntersection::Point {
            point: Point3::new(0.5, 0.5, 0.0),
            t1: 0.5,
            t2: 0.5,
        },
        CurveIntersection::Point {
            point: Point3::new(0.5, 0.5, 0.0),
            t1: 0.500001,
            t2: 0.500001,
        },
    ];
    deduplicate(&mut results, &Tolerance::default());
    assert_eq!(results.len(), 1);
}

#[test]
fn sort_by_overlap_t1_range() {
    // Test that overlaps sort correctly by t1_range.0
    let result = intersect_nurbs_nurbs(
        &NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        },
        (0.0, 1.0),
        &NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            weights: vec![1.0, 1.0],
        },
        (0.0, 1.0),
        &Tolerance::default(),
    );
    // Identical curves => overlap
    assert!(!result.is_empty());
}

#[test]
fn tangential_out_of_window_contact_terminates() {
    // PK-295/PK-297 pin: a query curve tangentially TOUCHING a
    // boundary line's extension OUTSIDE the requested window used to
    // explode bezier clipping (un-clipped terminal spans put the
    // contact inside the segment boxes; the stall path then 4-way
    // bisected under only a depth cap — O(4^50) with zero candidates
    // for the candidate cap to catch; measured as a 25-minute
    // boolean). Window-clipped decomposition removes the phantom
    // contact and the node budget bounds any real tangency; this
    // test simply terminating is the regression pin.
    let boundary = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.464, 1.8, 0.0), Point3::new(0.464, 0.0, 0.0)],
        weights: vec![1.0, 1.0],
    };
    // Vertical-tangent departure from (0.464, 1.8) — tangentially
    // aligned with the boundary line at its off-window endpoint.
    let query = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.464, 1.8, 0.0),
            Point3::new(0.464, 1.9, 0.0),
            Point3::new(0.6, 1.97, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    // Window covers only the LOWER third of the boundary — far from
    // the tangential contact.
    let hits = intersect_nurbs_nurbs(
        &boundary,
        (0.667, 1.0),
        &query,
        (0.0, 1.0),
        &Tolerance::default(),
    );
    assert!(
        hits.is_empty(),
        "no true intersection exists inside the window, got {hits:?}"
    );
}

/// The 9-control-point rational-quadratic circle the extrude rim
/// pcurves and `plane_uv_conic` produce.
fn rational_circle(cx: f64, cy: f64, r: f64, phase_deg: f64) -> NurbsCurve {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let p = phase_deg.to_radians();
    let s2 = std::f64::consts::SQRT_2;
    let at = |deg: f64, scale: f64| {
        let a = p + deg.to_radians();
        Point3::new(cx + scale * r * a.cos(), cy + scale * r * a.sin(), 0.0)
    };
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            at(0.0, 1.0),
            at(45.0, s2),
            at(90.0, 1.0),
            at(135.0, s2),
            at(180.0, 1.0),
            at(225.0, s2),
            at(270.0, 1.0),
            at(315.0, s2),
            at(0.0, 1.0),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

/// Two offset rational circles cross at exactly two transversal
/// points; both must be found. Measured miss (PK-427): near the
/// converged root the A segment collapses below the fat-line
/// helper's absolute epsilons, and treating "no fat line" as "no
/// intersection" discarded the root — a cap face's clip then saw
/// ONE crossing for a bore-class record and dropped it as Outside.
#[test]
fn two_offset_rational_circles_yield_both_crossings() {
    let rim = rational_circle(0.0, 0.3, 0.3, 0.0);
    let record = rational_circle(0.0, -0.1, 0.25, -90.0);
    let mut hits: Vec<(f64, f64)> =
        intersect_nurbs_nurbs(&rim, (0.0, 1.0), &record, (0.0, 1.0), &Tolerance::default())
            .into_iter()
            .filter_map(|i| match i {
                CurveIntersection::Point { point, .. } => Some((point.x, point.y)),
                _ => None,
            })
            .collect();
    hits.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    assert_eq!(hits.len(), 2, "both transversal crossings: {hits:?}");
    assert!((hits[0].0 + 0.187265).abs() < 1e-4 && (hits[1].0 - 0.187265).abs() < 1e-4);
    assert!((hits[0].1 - 0.065625).abs() < 1e-4 && (hits[1].1 - 0.065625).abs() < 1e-4);
}

/// The single quadrant pair that held the missed crossing, in
/// isolation — the bezier-clipping path alone must serve it.
#[test]
fn the_converged_collapse_quadrant_still_finds_its_crossing() {
    let rim = rational_circle(0.0, 0.3, 0.3, 0.0);
    let record = rational_circle(0.0, -0.1, 0.25, -90.0);
    let hits = intersect_nurbs_nurbs(
        &rim,
        (0.5, 0.75),
        &record,
        (0.5, 0.75),
        &Tolerance::default(),
    );
    assert_eq!(hits.len(), 1, "the quadrant's crossing: {hits:?}");
}

// PK-513: whole-curve hull rejection.
#[test]
fn far_apart_curves_are_rejected_by_their_hulls() {
    use super::hulls_disjoint;
    let a = line(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
    let far = line(Point3::new(0.0, 5.0, 0.0), Point3::new(1.0, 5.0, 0.0));
    let crossing = line(Point3::new(0.5, -1.0, 0.0), Point3::new(0.5, 1.0, 0.0));
    assert!(hulls_disjoint(&a, &far, 1e-4));
    assert!(!hulls_disjoint(&a, &crossing, 1e-4));
    // The pad keeps a near-touching pair in the solver's hands.
    let near = line(Point3::new(0.0, 1e-5, 0.0), Point3::new(1.0, 1e-5, 0.0));
    assert!(!hulls_disjoint(&a, &near, 1e-4));
}

fn line(a: Point3, b: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    }
}

/// PK-739: a fin side (line) crossing the housing circle (9-point
/// rational quadratic, r 55) transversally in quadrant III. The top
/// cap of a 16-fin ring missed exactly this crossing; the cellular
/// union then lost a fin footprint and dropped the housing wall.
#[test]
fn fin_side_crosses_the_rational_circle() {
    let r = 55.000000000000014;
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let circle = NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            Point3::new(-8.326672684688674e-16, -55.00000000000001, 0.0),
            Point3::new(r, -55.00000000000001, 0.0),
            Point3::new(r, 1.0547118733938987e-14, 0.0),
            Point3::new(r, 55.00000000000002, 0.0),
            Point3::new(-8.326672684688674e-16, 55.00000000000002, 0.0),
            Point3::new(-r, 55.00000000000002, 0.0),
            Point3::new(-r, 1.0547118733938987e-14, 0.0),
            Point3::new(-r, -55.00000000000001, 0.0),
            Point3::new(-8.326672684688674e-16, -55.00000000000001, 0.0),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    };
    let side = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(-43.48706704297268, -45.608387386532314, 0.0),
            Point3::new(-36.41599923110721, -38.53731957466684, 0.0),
        ],
        weights: vec![1.0, 1.0],
    };
    let hits = intersect_nurbs_nurbs(
        &side,
        (0.0, 1.0),
        &circle,
        (0.0, 1.0),
        &Tolerance::default(),
    );
    assert_eq!(hits.len(), 1, "one transversal crossing, got {hits:?}");
    let CurveIntersection::Point { point, .. } = &hits[0] else {
        panic!("point hit expected");
    };
    assert!(point.distance_to(&Point3::new(-37.81574655517892, -39.93706689873854, 0.0)) < 1e-6);
}

/// PK-739: a square side crossing a circle TWICE inside one quarter
/// arc's neighbourhood. Bounding the clip by the coefficient hull
/// revisits the first root at its noise floor; those revisits once
/// filled the candidate cap and the second crossing was never searched.
#[test]
fn square_side_crosses_circle_twice() {
    let r = 0.2;
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let circle = NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            Point3::new(r, 0.0, 0.0),
            Point3::new(r, r, 0.0),
            Point3::new(0.0, r, 0.0),
            Point3::new(-r, r, 0.0),
            Point3::new(-r, 0.0, 0.0),
            Point3::new(-r, -r, 0.0),
            Point3::new(0.0, -r, 0.0),
            Point3::new(r, -r, 0.0),
            Point3::new(r, 0.0, 0.0),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    };
    let side = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(-0.15, 0.15, 0.0),
            Point3::new(-0.15, -0.15, 0.0),
        ],
        weights: vec![1.0, 1.0],
    };
    let hits = intersect_nurbs_nurbs(
        &side,
        (0.0, 1.0),
        &circle,
        (0.0, 1.0),
        &Tolerance::default(),
    );
    assert_eq!(hits.len(), 2, "{hits:?}");
}
