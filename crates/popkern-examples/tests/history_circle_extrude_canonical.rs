//! PK-387 — every extrude of a full-circle profile serves ONE periodic
//! cylinder wall with single-edge caps, the same body `popkern-part`'s
//! wire sweep serves.
//!
//! The sketch builder splits a full circle into two arcs (its edge
//! graph avoids self-loop edges), so `BRep::extrude` of the compiled
//! face made two co-cylindrical laterals and two-edge caps. The exact
//! periodic wire sweep is now the route for full-circle profiles from
//! BOTH sketch dialects and for every end condition and second
//! direction; draft is the one class with no sweep arm, and that
//! boundary is recorded in the replay report rather than silent.

mod circle_extrude_support;

use circle_extrude_support::*;
use popkern_core::naming::RoleCategory;
use popkern_core::tolerance::Tolerance;
use popkern_history::*;
use popkern_part::prelude::*;
use popkern_sketch::Sketch;

#[test]
fn constraint_circle_blind_extrude_is_one_periodic_cylinder() {
    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: SketchPayload::Constraint(circle_sketch_json(0.3, -0.2, 0.35)),
        },
        "circle",
    );
    let dist = add_param(&mut tree, "depth", 0.8);
    tree.add_feature(
        Feature::Extrude {
            face: name(1, RoleCategory::Generated, 0),
            direction: None,
            end_condition: ExtrudeEnd::Blind { distance: dist },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "boss",
    );
    let body = served(&mut tree);
    assert_eq!(surface_census(&body), (2, 1, 0), "one cylinder wall");
    assert_eq!(
        one_edge_planar_loops(&body),
        2,
        "both caps are single closed edges"
    );
    let part = Profile::on(Plane::XY)
        .circle(0.35)
        .at(0.3, -0.2)
        .extrude(0.8)
        .expect("part");
    assert_same_body("constraint circle blind", &body, part.brep());
}

#[test]
fn analytic_circle_second_direction_is_one_periodic_cylinder() {
    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: analytic_circle(0.0, 0.0, 0.25),
        },
        "circle",
    );
    let up = add_param(&mut tree, "up", 0.3);
    let down = add_param(&mut tree, "down", 0.2);
    tree.add_feature(
        Feature::Extrude {
            face: name(1, RoleCategory::Generated, 0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: up },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: Some(ExtrudeSecondDir {
                direction: Vec3Desc::Literal([0.0, 0.0, -1.0]),
                end_condition: ExtrudeEnd::Blind { distance: down },
                draft_angle: None,
            }),
        },
        "both_ways",
    );
    let body = served(&mut tree);
    assert_eq!(
        surface_census(&body),
        (2, 1, 0),
        "one wall spanning both directions"
    );
    let part = Profile::on(Plane::xy_at(-0.2))
        .circle(0.25)
        .extrude(0.5)
        .expect("part");
    assert_same_body("second direction", &body, part.brep());
}

#[test]
fn analytic_circle_mid_plane_spans_both_halves() {
    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: analytic_circle(0.0, 0.0, 0.25),
        },
        "circle",
    );
    let half = add_param(&mut tree, "half", 0.3);
    tree.add_feature(
        Feature::Extrude {
            face: name(1, RoleCategory::Generated, 0),
            direction: None,
            end_condition: ExtrudeEnd::MidPlane {
                half_distance: half,
            },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "mid",
    );
    let body = served(&mut tree);
    assert_eq!(surface_census(&body), (2, 1, 0));
    let part = Profile::on(Plane::xy_at(-0.3))
        .circle(0.25)
        .extrude(0.6)
        .expect("part");
    assert_same_body("mid-plane", &body, part.brep());
}

#[test]
fn constraint_rect_mid_plane_grows_both_halves_on_the_generic_path() {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(-0.5, -0.4);
    let p1 = sk.add_point(0.5, -0.4);
    let p2 = sk.add_point(0.5, 0.4);
    let p3 = sk.add_point(-0.5, 0.4);
    sk.add_line(p0, p1);
    sk.add_line(p1, p2);
    sk.add_line(p2, p3);
    sk.add_line(p3, p0);
    sk.solve().expect("rect solves");
    let json = serde_json::to_string(&sk).expect("serialize");
    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: SketchPayload::Constraint(json),
        },
        "rect",
    );
    let half = add_param(&mut tree, "half", 0.3);
    tree.add_feature(
        Feature::Extrude {
            face: name(1, RoleCategory::Generated, 0),
            direction: None,
            end_condition: ExtrudeEnd::MidPlane {
                half_distance: half,
            },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "mid",
    );
    let body = served(&mut tree);
    let zs: Vec<f64> = body.vertices.values().map(|v| v.position.z).collect();
    let (lo, hi) = zs
        .iter()
        .fold((f64::MAX, f64::MIN), |(lo, hi), z| (lo.min(*z), hi.max(*z)));
    assert!(
        (lo + 0.3).abs() < 1e-9 && (hi - 0.3).abs() < 1e-9,
        "z span [{lo}, {hi}]"
    );
    let v = volume(&body);
    assert!((v - 1.0 * 0.8 * 0.6).abs() < 1e-9, "volume {v}");
}

#[test]
fn draft_on_a_circle_serves_one_cone() {
    // PK-395: the exact route drafts its single periodic lateral with
    // the same operator the generic path uses — one exact Cone, where
    // the generic face extrude served two co-conical halves.
    use popkern_core::geometry::surface::SurfaceGeometry;
    for (label, payload) in [
        ("analytic", analytic_circle(0.0, 0.0, 0.3)),
        (
            "constraint",
            SketchPayload::Constraint(circle_sketch_json(0.0, 0.0, 0.3)),
        ),
    ] {
        let mut tree = FeatureTree::new();
        tree.add_feature(
            Feature::Sketch {
                plane: SketchPlaneDesc::XY,
                sketch_data: payload,
            },
            "circle",
        );
        let dist = add_param(&mut tree, "depth", 0.5);
        let draft = add_param(&mut tree, "draft", 0.1);
        let extrude = tree.add_feature(
            Feature::Extrude {
                face: name(1, RoleCategory::Generated, 0),
                direction: None,
                end_condition: ExtrudeEnd::Blind { distance: dist },
                mode: OperationMode::NewBody,
                draft_angle: Some(draft),
                second_direction: None,
            },
            "drafted",
        );
        let result = tree.replay(Tolerance::default()).expect("replay");
        assert!(
            result.report.failures().next().is_none(),
            "{label}: {:?}",
            result.report
        );
        let report = result.report.feature(extrude).expect("report");
        assert!(
            report.fallbacks.is_empty(),
            "{label}: the exact drafted arm serves; no fallback: {:?}",
            report.fallbacks
        );
        let body = tree.current_brep().expect("body");
        assert!(
            body.validate_all().is_empty(),
            "{label}: {:?}",
            body.validate_all()
        );
        assert_eq!(
            surface_census(&body),
            (2, 0, 1),
            "{label}: one Cone lateral"
        );
        assert_eq!(
            one_edge_planar_loops(&body),
            2,
            "{label}: single closed-edge caps"
        );
        // The lateral is an exact cone whose half-angle is the draft
        // angle, and the cap radii differ by h·tan(angle).
        let cone = body
            .faces
            .values()
            .find_map(|f| match f.surface.as_ref().map(|s| &s.geometry) {
                Some(SurfaceGeometry::Cone { half_angle, .. }) => Some(*half_angle),
                _ => None,
            })
            .expect("cone lateral");
        assert!(
            (cone.abs() - 0.1).abs() < 1e-9,
            "{label}: cone half-angle {cone} vs draft 0.1"
        );
        let radius_at = |z: f64| -> f64 {
            body.vertices
                .values()
                .filter(|v| (v.position.z - z).abs() < 1e-9)
                .map(|v| (v.position.x * v.position.x + v.position.y * v.position.y).sqrt())
                .fold(0.0, f64::max)
        };
        let (r0, r1) = (radius_at(0.0), radius_at(0.5));
        assert!((r0 - 0.3).abs() < 1e-9, "{label}: base radius {r0}");
        assert!(
            ((r0 - r1).abs() - 0.5 * 0.1f64.tan()).abs() < 1e-9,
            "{label}: |Δr| {} vs h·tan(angle) {}",
            (r0 - r1).abs(),
            0.5 * 0.1f64.tan()
        );
        let v = volume(&body);
        let frustum = std::f64::consts::PI * 0.5 * (r0 * r0 + r0 * r1 + r1 * r1) / 3.0;
        // The default tessellation inscribes the cone (~0.8 % volume
        // deficit); a cylinder would read 18 % off, so 1 % is decisive.
        assert!(
            (v - frustum).abs() < 1e-2 * frustum,
            "{label}: volume {v} vs frustum {frustum}"
        );
    }
}
