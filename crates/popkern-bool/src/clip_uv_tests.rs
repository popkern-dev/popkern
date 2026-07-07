// Structural lint marker: this file is test code, included via #[cfg(test)] in clip.rs.
use popkern_core::{
    geometry::{Point3, curve_on_surface::CurveOnSurface, nurbs_curve::NurbsCurve},
    tolerance::Tolerance,
    topology::BRep,
};
use std::f64::consts::{FRAC_1_SQRT_2, PI, TAU};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Build a degree-2 rational NURBS circle in the XY plane.
fn nurbs_circle(center: Point3, radius: f64) -> NurbsCurve {
    let w = FRAC_1_SQRT_2;
    let r = radius;
    let cx = center.x;
    let cy = center.y;
    let cz = center.z;
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            Point3::new(cx + r, cy, cz),
            Point3::new(cx + r, cy + r, cz),
            Point3::new(cx, cy + r, cz),
            Point3::new(cx - r, cy + r, cz),
            Point3::new(cx - r, cy, cz),
            Point3::new(cx - r, cy - r, cz),
            Point3::new(cx, cy - r, cz),
            Point3::new(cx + r, cy - r, cz),
            Point3::new(cx + r, cy, cz),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

fn nurbs_line(p0: Point3, p1: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![p0, p1],
        weights: vec![1.0, 1.0],
    }
}

/// Build a CurveOnSurface with a degree-2 circle for 3D and UV curves.
fn cos_circle(center_3d: Point3, center_uv: Point3, radius: f64) -> CurveOnSurface {
    CurveOnSurface {
        curve_3d: nurbs_circle(center_3d, radius),
        curve_uv1: nurbs_circle(center_uv, radius),
        curve_uv2: nurbs_circle(center_uv, radius),
        is_closed: true,
    }
}

// ── clip_curve_to_faces: UV-space clipping for non-linear curves ─────────

#[test]
fn clip_uv_space_circle_inside_both_faces_passes_through() {
    // Both faces are unit cubes with pcurves — circle at (0.5, 0.5) with
    // r=0.1 is entirely inside the [0,1] UV domain of both faces.
    let cube_a = BRep::unit_cube(tol()).expect("cube_a");
    let face_a = cube_a.faces.keys().next().expect("face_a");
    let cube_b = BRep::unit_cube(tol()).expect("cube_b");
    let face_b = cube_b.faces.keys().next().expect("face_b");

    let curve = cos_circle(Point3::new(0.5, 0.5, 0.0), Point3::new(0.5, 0.5, 0.0), 0.1);

    let result = super::clip_curve_to_faces(&cube_a, face_a, &cube_b, face_b, &curve);
    assert!(
        result.is_some(),
        "closed circle inside both faces should pass through"
    );
    let r = result.expect("verified");
    assert!(r.is_closed, "should remain closed");
    assert_eq!(r.curve_3d.degree, 2, "should remain degree 2");
}

#[test]
fn clip_uv_space_circle_outside_face_returns_none() {
    // Circle at (5.0, 5.0) with r=0.1 — far outside the [0,1] UV domain.
    let cube_a = BRep::unit_cube(tol()).expect("cube_a");
    let face_a = cube_a.faces.keys().next().expect("face_a");
    let cube_b = BRep::unit_cube(tol()).expect("cube_b");
    let face_b = cube_b.faces.keys().next().expect("face_b");

    let curve = CurveOnSurface {
        curve_3d: nurbs_circle(Point3::new(5.0, 5.0, 0.0), 0.1),
        curve_uv1: nurbs_circle(Point3::new(5.0, 5.0, 0.0), 0.1),
        curve_uv2: nurbs_circle(Point3::new(5.0, 5.0, 0.0), 0.1),
        is_closed: false, // not closed, so 0 crossings => None
    };

    let result = super::clip_curve_to_faces(&cube_a, face_a, &cube_b, face_b, &curve);
    assert!(
        result.is_none(),
        "non-closed curve outside both faces should return None"
    );
}

#[test]
fn clip_uv_space_non_closed_outside_returns_none() {
    // Non-closed circle far outside the face UV domain — 0 crossings => None.
    let cube_a = BRep::unit_cube(tol()).expect("cube_a");
    let face_a = cube_a.faces.keys().next().expect("face_a");
    let cube_b = BRep::unit_cube(tol()).expect("cube_b");
    let face_b = cube_b.faces.keys().next().expect("face_b");

    let curve = CurveOnSurface {
        curve_3d: nurbs_circle(Point3::new(10.0, 10.0, 0.0), 0.1),
        curve_uv1: nurbs_circle(Point3::new(10.0, 10.0, 0.0), 0.1),
        curve_uv2: nurbs_circle(Point3::new(10.0, 10.0, 0.0), 0.1),
        is_closed: false,
    };

    let result = super::clip_curve_to_faces(&cube_a, face_a, &cube_b, face_b, &curve);
    assert!(
        result.is_none(),
        "non-closed curve outside both faces should return None"
    );
}

#[test]
fn clip_uv_space_no_pcurves_falls_back_to_unclipped() {
    // Build faces without pcurves (manual Euler ops).
    let mut brep_a = BRep::new(tol());
    let mvfs_a = brep_a.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs_a");
    let mev1 = brep_a
        .mev_at_vertex(mvfs_a.vertex, mvfs_a.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev1");
    let mev2 = brep_a
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .expect("mev2");
    let mev3 = brep_a
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .expect("mev3");
    let mef_a = brep_a
        .mef(mev3.half_edge_twin, mev1.half_edge)
        .expect("mef_a");

    let mut brep_b = BRep::new(tol());
    let mvfs_b = brep_b.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs_b");
    let mev4 = brep_b
        .mev_at_vertex(mvfs_b.vertex, mvfs_b.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev4");
    let mev5 = brep_b
        .mev(mev4.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .expect("mev5");
    let mev6 = brep_b
        .mev(mev5.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .expect("mev6");
    let mef_b = brep_b
        .mef(mev6.half_edge_twin, mev4.half_edge)
        .expect("mef_b");

    let curve = cos_circle(Point3::new(0.5, 0.5, 0.0), Point3::new(0.5, 0.5, 0.0), 0.3);

    let result = super::clip_curve_to_faces(&brep_a, mef_a.face, &brep_b, mef_b.face, &curve);
    assert!(
        result.is_some(),
        "non-linear curve with no pcurves should fall back to unclipped"
    );
    let r = result.expect("verified");
    assert!(r.is_closed, "should keep is_closed flag");
    assert_eq!(r.curve_3d.degree, 2, "should keep degree");
}

// ── uv_range_from_crossings unit tests ───────────────────────────────────

#[test]
fn uv_range_from_crossings_single_crossing_returns_none() {
    let dummy_he = popkern_core::topology::HalfEdgeHandle::default();
    let crossings = vec![super::UvCrossing {
        t_curve: 0.5,
        boundary_he: dummy_he,
        t_boundary: 0.0,
    }];
    let curve = nurbs_line(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
    let result = super::uv_range_from_crossings(&crossings, &curve, false);
    assert!(result.is_none(), "single crossing should return None");
}

#[test]
fn uv_range_from_crossings_closed_no_crossings_returns_full() {
    use popkern_core::geometry::curve::CurveEval;
    let crossings = vec![];
    let curve = nurbs_circle(Point3::new(0.5, 0.5, 0.0), 0.1);
    let result = super::uv_range_from_crossings(&crossings, &curve, true);
    assert!(
        result.is_some(),
        "closed curve with 0 crossings => full range"
    );
    let (t_enter, t_exit) = result.expect("verified");
    let (t_min, t_max) = curve.parameter_range();
    assert!((t_enter - t_min).abs() < 1e-12);
    assert!((t_exit - t_max).abs() < 1e-12);
}

#[test]
fn uv_range_from_crossings_two_crossings_returns_range() {
    let dummy_he = popkern_core::topology::HalfEdgeHandle::default();
    let crossings = vec![
        super::UvCrossing {
            t_curve: 0.2,
            boundary_he: dummy_he,
            t_boundary: 0.0,
        },
        super::UvCrossing {
            t_curve: 0.7,
            boundary_he: dummy_he,
            t_boundary: 0.0,
        },
    ];
    let curve = nurbs_line(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
    let result = super::uv_range_from_crossings(&crossings, &curve, false);
    assert!(result.is_some());
    let (t_enter, t_exit) = result.expect("verified");
    assert!((t_enter - 0.2).abs() < 1e-12);
    assert!((t_exit - 0.7).abs() < 1e-12);
}

#[test]
fn uv_range_from_crossings_open_no_crossings_returns_none() {
    let crossings: Vec<super::UvCrossing> = vec![];
    let curve = nurbs_line(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
    let result = super::uv_range_from_crossings(&crossings, &curve, false);
    assert!(
        result.is_none(),
        "open curve with 0 crossings should return None"
    );
}

// ── is_near_seam tests ──────────────────────────────────────────────────

#[test]
fn is_near_seam_returns_true_for_small_u() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.1, 0.5, 0.0), Point3::new(0.5, 0.5, 0.0)],
        weights: vec![1.0, 1.0],
    };
    assert!(super::is_near_seam(&curve), "u=0.1 is below pi/4");
}

#[test]
fn is_near_seam_returns_true_for_large_u() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(6.0, 0.5, 0.0), Point3::new(PI, 0.5, 0.0)],
        weights: vec![1.0, 1.0],
    };
    assert!(super::is_near_seam(&curve), "u=6.0 is above 7*pi/4");
}

#[test]
fn is_near_seam_returns_false_for_mid_range() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(PI, 0.5, 0.0), Point3::new(PI + 0.5, 0.5, 0.0)],
        weights: vec![1.0, 1.0],
    };
    assert!(!super::is_near_seam(&curve), "u around pi is far from seam");
}

// ── shift_curve_u tests ─────────────────────────────────────────────────

#[test]
fn shift_curve_u_wraps_to_period() {
    let curve = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.1, 0.5, 0.0), Point3::new(6.0, 0.5, 0.0)],
        weights: vec![1.0, 1.0],
    };
    let shifted = super::shift_curve_u(&curve, PI);
    // 0.1 + pi should stay as-is (< 2pi).
    assert!(
        (shifted.control_points[0].x - (0.1 + PI)).abs() < 1e-12,
        "0.1 + pi should not wrap"
    );
    // 6.0 + pi wraps: (6.0 + pi) mod 2pi
    let expected = (6.0 + PI).rem_euclid(TAU);
    assert!(
        (shifted.control_points[1].x - expected).abs() < 1e-12,
        "6.0 + pi should wrap"
    );
    // Y coordinates unchanged.
    assert!((shifted.control_points[0].y - 0.5).abs() < 1e-12);
}

#[test]
fn shift_curve_u_preserves_degree_and_knots() {
    let curve = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
            Point3::new(3.0, 0.0, 0.0),
        ],
        weights: vec![1.0, 1.0, 1.0],
    };
    let shifted = super::shift_curve_u(&curve, 0.5);
    assert_eq!(shifted.degree, 2);
    assert_eq!(shifted.knots, curve.knots);
    assert_eq!(shifted.weights, curve.weights);
}

// ── find_uv_crossings with seam-crossing boundary ───────────────────────

/// Helper: build a quad face and assign pcurves to its half-edges.
fn build_face_with_pcurves(
    uv_corners: [(Point3, Point3); 4],
) -> (BRep, popkern_core::topology::FaceHandle) {
    let tol = Tolerance::default();
    let mut brep = BRep::new(tol);
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).expect("mvfs");
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .expect("mev1");
    let mev2 = brep
        .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
        .expect("mev2");
    let mev3 = brep
        .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
        .expect("mev3");
    let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).expect("mef");

    let face_data = &brep.faces[mef.face];
    let loop_data = brep.loops.get(face_data.outer_loop).expect("loop");
    let start = loop_data.half_edge;
    let mut cursor = start;
    for (i, (p0, p1)) in uv_corners.iter().enumerate() {
        let pcurve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![*p0, *p1],
            weights: vec![1.0, 1.0],
        };
        brep.half_edges[cursor].pcurve = Some(pcurve.into());
        cursor = brep.half_edges[cursor].next;
        if cursor == start && i < 3 {
            break;
        }
    }
    (brep, mef.face)
}

#[test]
fn find_uv_crossings_near_seam_finds_crossings() {
    // Boundary spans the seam: u=[5.5, 0.8]. The seam-shift dance
    // is gated on the face carrying a PERIODIC surface (the raw
    // u-value thresholds fired spuriously on planar faces whose UV
    // coordinates are plain lengths), so this fixture must attach
    // the cylinder its pcurves imply.
    let (mut brep, face) = build_face_with_pcurves([
        (Point3::new(5.5, 0.0, 0.0), Point3::new(0.8, 0.0, 0.0)),
        (Point3::new(0.8, 0.0, 0.0), Point3::new(0.8, 1.0, 0.0)),
        (Point3::new(0.8, 1.0, 0.0), Point3::new(5.5, 1.0, 0.0)),
        (Point3::new(5.5, 1.0, 0.0), Point3::new(5.5, 0.0, 0.0)),
    ]);
    brep.faces[face].surface = Some(popkern_core::geometry::SurfaceOnFace {
        geometry: popkern_core::geometry::surface::SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        },
        orientation: true,
    });

    // Curve near u=0.1 crossing the face in v direction.
    let curve_uv = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.1, -0.5, 0.0), Point3::new(0.1, 1.5, 0.0)],
        weights: vec![1.0, 1.0],
    };

    let crossings = super::find_uv_crossings(&brep, face, &curve_uv, &tol());
    assert!(crossings.is_some(), "should produce crossings");
    let crossings = crossings.expect("verified");
    assert!(
        crossings.len() >= 2,
        "curve crossing seam-spanning face should have >= 2 crossings, got {}",
        crossings.len()
    );
}

#[test]
fn find_uv_crossings_far_from_seam_unaffected() {
    // Boundary in [1,2]x[0,1] — no seam.
    let (brep, face) = build_face_with_pcurves([
        (Point3::new(1.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)),
        (Point3::new(2.0, 0.0, 0.0), Point3::new(2.0, 1.0, 0.0)),
        (Point3::new(2.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
        (Point3::new(1.0, 1.0, 0.0), Point3::new(1.0, 0.0, 0.0)),
    ]);

    let curve_uv = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(1.5, -0.5, 0.0), Point3::new(1.5, 1.5, 0.0)],
        weights: vec![1.0, 1.0],
    };

    let crossings = super::find_uv_crossings(&brep, face, &curve_uv, &tol());
    assert!(crossings.is_some());
    let crossings = crossings.expect("verified");
    assert!(
        crossings.len() >= 2,
        "curve through mid-range face should have >= 2 crossings, got {}",
        crossings.len()
    );
}
