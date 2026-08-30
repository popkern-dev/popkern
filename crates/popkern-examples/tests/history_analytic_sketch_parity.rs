//! PK-364 acceptance: a `Feature::Sketch` in the ANALYTIC dialect
//! (shape + placement + parameter bindings, compiled exactly on
//! replay — no solver) extruded through the feature tree replays to
//! the same body as the equivalent `popkern-part` `Profile` call,
//! survives `.pkn` save/load, and keeps the exactness contract: the
//! circle's walls are true `Cylinder`s at the exact centre and radius
//! (design invariant 6), not NURBS.
//!
//! Parity is asserted through numeric gates — face/edge/vertex
//! census, sorted vertex positions, body volume, surface-kind census
//! — never visual inspection. Where the two pipelines legitimately
//! differ in seam census (the sweep serves one seamed cylinder wall,
//! the extrude operator two half walls from the two half-arc edges),
//! the gate compares geometry, not handle counts.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_history::*;
use popkern_part::prelude::*;
use popkern_sketch::{AnalyticProfile, PathSegment, ProfileShape};

fn name(op: u64, category: RoleCategory, index: u32) -> EntityName {
    EntityName {
        operation: OperationId(op),
        role: EntityRole { category, index },
        ancestors: Vec::new(),
    }
}

fn add_param(tree: &mut FeatureTree, label: &str, value: f64, unit: Unit) -> ParameterId {
    tree.parameters_mut()
        .add_global(label, value, unit)
        .expect("add param")
}

fn sorted_vertices(brep: &BRep) -> Vec<[f64; 3]> {
    let mut pts: Vec<[f64; 3]> = brep
        .vertices
        .values()
        .map(|v| [v.position.x, v.position.y, v.position.z])
        .collect();
    pts.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    pts
}

fn volume(brep: &BRep) -> f64 {
    popkern_measure::body_volume(brep, &Tolerance::default()).expect("volume")
}

/// Full census + vertex + volume parity.
fn assert_same_body(label: &str, a: &BRep, b: &BRep, tol: f64) {
    assert_eq!(a.faces.len(), b.faces.len(), "{label}: face count");
    assert_eq!(a.edges.len(), b.edges.len(), "{label}: edge count");
    assert_eq!(a.vertices.len(), b.vertices.len(), "{label}: vertex count");
    for (p, q) in sorted_vertices(a).iter().zip(sorted_vertices(b).iter()) {
        for k in 0..3 {
            assert!((p[k] - q[k]).abs() < tol, "{label}: vertex {p:?} vs {q:?}");
        }
    }
    // Volumes of curved bodies come through a fine tessellation on
    // each side; the two pipelines chart the same cylinder differently,
    // so the budget is the mesh's, not the vertex gate's.
    let (va, vb) = (volume(a), volume(b));
    assert!(
        (va - vb).abs() < 1e-6 * va.abs().max(1.0),
        "{label}: volume {va} vs {vb}"
    );
}

fn surface_census(brep: &BRep) -> (usize, usize, usize) {
    let (mut planes, mut cylinders, mut other) = (0, 0, 0);
    for f in brep.faces.values() {
        match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Plane { .. }) => planes += 1,
            Some(SurfaceGeometry::Cylinder { .. }) => cylinders += 1,
            _ => other += 1,
        }
    }
    (planes, cylinders, other)
}

/// Analytic sketch on XY + blind extrude along +Z, saved, loaded, and
/// replayed under the full validation gate. Returns the loaded tree
/// and the served body.
fn replay_extruded(
    profile: AnalyticProfile<ScalarBinding>,
    depth: f64,
    with_tree: impl FnOnce(&mut FeatureTree),
) -> (FeatureTree, ReplayResult, Option<BRep>) {
    let mut tree = FeatureTree::new();
    with_tree(&mut tree);
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: SketchPayload::Analytic(profile),
        },
        "AnalyticSketch",
    );
    let dist = add_param(&mut tree, "_depth", depth, Unit::Meters);
    tree.add_feature(
        Feature::Extrude {
            face: name(1, RoleCategory::Generated, 0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude",
    );
    let bytes = tree.save().expect("save");
    let mut loaded = FeatureTree::load(&bytes).expect("load");
    let result = loaded.replay(Tolerance::default()).expect("replay runs");
    let body = loaded.current_brep().cloned();
    (loaded, result, body)
}

fn assert_replay_clean(result: &ReplayResult) {
    let failures: Vec<_> = result
        .feature_results
        .iter()
        .filter(|(_, r)| !matches!(r, FeatureResult::Success { .. }))
        .collect();
    assert!(failures.is_empty(), "replay failures: {failures:?}");
}

fn lit(v: f64) -> ScalarBinding {
    ScalarBinding::Literal(v)
}

#[test]
fn rect_and_polygon_profiles_replay_to_part_parity() {
    let rect = AnalyticProfile {
        shape: ProfileShape::Rect {
            width: lit(2.0),
            height: lit(1.0),
        },
        center: [lit(0.5), lit(-0.25)],
    };
    let (_, result, body) = replay_extruded(rect, 0.5, |_| {});
    assert_replay_clean(&result);
    let body = body.expect("body");
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
    let part = Profile::on(Plane::XY)
        .rect(2.0, 1.0)
        .at(0.5, -0.25)
        .extrude(0.5)
        .expect("part rect");
    assert_same_body("rect", &body, part.brep(), 1e-9);
    assert_eq!(surface_census(&body), (6, 0, 0));

    let tri = AnalyticProfile {
        shape: ProfileShape::Polygon {
            points: vec![
                [lit(0.0), lit(0.0)],
                [lit(4.0), lit(0.0)],
                [lit(0.0), lit(3.0)],
            ],
        },
        center: [lit(0.0), lit(0.0)],
    };
    let (_, result, body) = replay_extruded(tri, 2.0, |_| {});
    assert_replay_clean(&result);
    let body = body.expect("body");
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
    let part = Profile::on(Plane::XY)
        .polygon(&[(0.0, 0.0), (4.0, 0.0), (0.0, 3.0)])
        .extrude(2.0)
        .expect("part polygon");
    assert_same_body("polygon", &body, part.brep(), 1e-9);
}

/// The exactness contract end to end: a circle bound to a PARAMETER
/// replays to true `Cylinder` walls at the exact centre and radius;
/// the volume matches the sweep-built part; a parameter edit replays
/// to the new radius.
#[test]
fn circle_profile_replays_to_exact_cylinder_and_tracks_its_parameter() {
    let circle = |r: ScalarBinding| AnalyticProfile {
        shape: ProfileShape::Circle { radius: r },
        center: [lit(1.0), lit(0.0)],
    };
    let mut tree = FeatureTree::new();
    let r = add_param(&mut tree, "radius", 0.25, Unit::Meters);
    let (mut loaded, result, body) = replay_extruded(circle(r.into()), 0.5, |t| {
        *t = tree.clone();
    });
    assert_replay_clean(&result);
    let body = body.expect("body");
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
    // One seamed cylinder wall (the exact wire sweep `popkern-part`
    // uses — PK-382 unified the analytic + blind path; the compiled
    // face's two half arcs no longer leak into the solid), one axis,
    // exact radius.
    assert_eq!(surface_census(&body), (2, 1, 0));
    for f in body.faces.values() {
        if let Some(SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        }) = f.surface.as_ref().map(|s| &s.geometry)
        {
            assert!((origin.x - 1.0).abs() < 1e-12 && origin.y.abs() < 1e-12);
            assert!((axis.z.abs() - 1.0).abs() < 1e-12);
            assert!((radius - 0.25).abs() < 1e-12);
        }
    }
    let part = Profile::on(Plane::XY)
        .circle(0.25)
        .at(1.0, 0.0)
        .extrude(0.5)
        .expect("part circle");
    assert_eq!(
        surface_census(part.brep()),
        (2, 1, 0),
        "sweep serves one seamed wall"
    );
    assert_eq!(
        (body.faces.len(), body.edges.len(), body.vertices.len()),
        (
            part.brep().faces.len(),
            part.brep().edges.len(),
            part.brep().vertices.len()
        ),
        "tree and direct build share one topology"
    );
    let (vh, vp) = (volume(&body), volume(part.brep()));
    assert!((vh - vp).abs() < 1e-6 * vp, "volume {vh} vs {vp}");
    let exact = std::f64::consts::PI * 0.25 * 0.25 * 0.5;
    assert!(
        (vh - exact).abs() < 5e-3 * exact,
        "history volume {vh} vs π r² h {exact}"
    );

    // Parameter edit: the loaded document's `radius` (same id — the
    // table's keys are document state) drives the profile.
    loaded
        .set_parameter_value(r, 0.4, Unit::Meters)
        .expect("set radius");
    let result = loaded.replay(Tolerance::default()).expect("replay");
    assert_replay_clean(&result);
    let body = loaded.current_brep().expect("body");
    assert!(body.validate_all().is_empty());
    let exact = std::f64::consts::PI * 0.4 * 0.4 * 0.5;
    let v = volume(body);
    assert!(
        (v - exact).abs() < 5e-3 * exact,
        "edited volume {v} vs {exact}"
    );
    let radii: Vec<f64> = body
        .faces
        .values()
        .filter_map(|f| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Cylinder { radius, .. }) => Some(*radius),
            _ => None,
        })
        .collect();
    assert_eq!(radii.len(), 1, "one seamed wall after the edit");
    assert!(radii.iter().all(|r| (r - 0.4).abs() < 1e-12));
}

/// A mixed line/arc outline (the part API's `path` class): the history
/// body and the part body agree in full census — one exact arc edge
/// per arc on both sides — and the arc wall is a true `Cylinder`.
#[test]
fn mixed_path_profile_replays_to_part_parity_with_a_true_cylinder_wall() {
    let path = AnalyticProfile {
        shape: ProfileShape::Path {
            start: [lit(-30.0), lit(0.0)],
            segments: vec![
                PathSegment::Line {
                    to: [lit(30.0), lit(0.0)],
                },
                PathSegment::Arc {
                    via: [lit(0.0), lit(18.0)],
                    to: [lit(-30.0), lit(0.0)],
                },
            ],
        },
        center: [lit(0.0), lit(0.0)],
    };
    let (_, result, body) = replay_extruded(path, 20.0, |_| {});
    assert_replay_clean(&result);
    let body = body.expect("body");
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
    let part = Profile::on(Plane::XY)
        .path(-30.0, 0.0)
        .line_to(30.0, 0.0)
        .arc_to((0.0, 18.0), (-30.0, 0.0))
        .close()
        .extrude(20.0)
        .expect("part path");
    assert_same_body("mixed path", &body, part.brep(), 1e-9);
    assert_eq!(surface_census(&body), (3, 1, 0));
    assert_eq!(surface_census(part.brep()), (3, 1, 0));
}

/// An uncompilable description declines transactionally with the
/// defect named, and the tree serves nothing rather than garbage.
#[test]
fn degenerate_profile_declines_transactionally() {
    let bad = AnalyticProfile {
        shape: ProfileShape::Circle { radius: lit(0.0) },
        center: [lit(0.0), lit(0.0)],
    };
    let (_, result, body) = replay_extruded(bad, 1.0, |_| {});
    let failed = result.feature_results.iter().any(|(_, r)| {
        matches!(r, FeatureResult::Failed { error, .. } if error.contains("circle radius must be positive"))
    });
    assert!(
        failed,
        "expected a named decline, got {:?}",
        result.feature_results
    );
    assert!(body.is_none() || body.as_ref().map(|b| b.faces.is_empty()).unwrap_or(true));
}

/// E3 / PK-371: the CONSTRAINT dialect reaches the same exactness —
/// a solved constraint sketch with an arc extrudes to a body whose
/// arc wall is a true `Cylinder`, not NURBS, and the full circle (the
/// case that replayed as two identical half-arcs) serves watertight.
#[test]
fn constraint_sketch_arc_extrudes_to_a_true_cylinder() {
    use popkern_sketch::{Constraint, Sketch};
    // "D" outline: line (-30,0)→(30,0) + the arc through (0,18) back.
    let mut sk = Sketch::new();
    let p0 = sk.add_point(-30.0, 0.0);
    let p1 = sk.add_point(30.0, 0.0);
    sk.add_line(p0, p1);
    let cy = (18.0_f64.powi(2) - 30.0_f64.powi(2)) / 36.0;
    let r = (30.0_f64.powi(2) + cy * cy).sqrt();
    let a0 = (0.0 - cy).atan2(30.0);
    sk.add_arc(0.0, cy, r, a0, std::f64::consts::PI - a0);
    sk.add_constraint(Constraint::Fixed(p0));
    sk.add_constraint(Constraint::Fixed(p1));
    sk.solve().expect("D outline solves");
    let json = serde_json::to_string(&sk).expect("serialize");

    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: SketchPayload::Constraint(json),
        },
        "ConstraintSketch",
    );
    let dist = add_param(&mut tree, "_depth", 20.0, Unit::Meters);
    tree.add_feature(
        Feature::Extrude {
            face: name(1, RoleCategory::Generated, 0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude",
    );
    let result = tree.replay(Tolerance::default()).expect("replay");
    assert_replay_clean(&result);
    let body = tree.current_brep().expect("body");
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
    assert_eq!(
        surface_census(body),
        (3, 1, 0),
        "arc wall must be a Cylinder"
    );
    let part = Profile::on(Plane::XY)
        .path(-30.0, 0.0)
        .line_to(30.0, 0.0)
        .arc_to((0.0, 18.0), (-30.0, 0.0))
        .close()
        .extrude(20.0)
        .expect("part");
    assert_same_body("constraint D", body, part.brep(), 1e-6);

    // Full circle through the constraint dialect.
    let mut sk = Sketch::new();
    sk.add_arc(1.0, 0.0, 0.25, 0.0, std::f64::consts::TAU);
    sk.solve().expect("circle solves");
    let json = serde_json::to_string(&sk).expect("serialize");
    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: SketchPayload::Constraint(json),
        },
        "Circle",
    );
    let dist = add_param(&mut tree, "_depth", 0.5, Unit::Meters);
    tree.add_feature(
        Feature::Extrude {
            face: name(1, RoleCategory::Generated, 0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude",
    );
    let result = tree.replay(Tolerance::default()).expect("replay");
    assert_replay_clean(&result);
    let body = tree.current_brep().expect("body");
    assert!(body.validate_all().is_empty(), "{:?}", body.validate_all());
    // One periodic cylinder wall (PK-387): the constraint dialect's full
    // circle takes the exact wire sweep like the analytic one.
    assert_eq!(surface_census(body), (2, 1, 0));
    let part = Profile::on(Plane::XY)
        .circle(0.25)
        .at(1.0, 0.0)
        .extrude(0.5)
        .expect("part circle");
    assert_same_body("constraint circle", body, part.brep(), 1e-9);
    let exact = std::f64::consts::PI * 0.25 * 0.25 * 0.5;
    let v = volume(body);
    assert!(
        (v - exact).abs() < 5e-3 * exact,
        "constraint circle volume {v} vs {exact}"
    );
}
