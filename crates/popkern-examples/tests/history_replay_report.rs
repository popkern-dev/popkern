//! PK-365 acceptance: the replay report channel. A replayed feature
//! tree returns per-feature evidence beside the BRep — a replayed
//! fillet carries its measured continuity record (equal to a direct
//! measurement of the served body), chamfers carry a crease record,
//! suppressed and failed features say why, and a partial replay keeps
//! the untouched records on the tree.
//!
//! Every claim is a numeric or structural gate on the report payload,
//! never a visual check.

use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_history::*;
use popkern_sketch::{Constraint, Sketch};

fn rectangle_sketch_json() -> String {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0);
    let p2 = sk.add_point(1.0, 1.0);
    let p3 = sk.add_point(0.0, 1.0);
    sk.add_line(p0, p1);
    sk.add_line(p1, p2);
    sk.add_line(p2, p3);
    sk.add_line(p3, p0);
    for p in [p0, p1, p2, p3] {
        sk.add_constraint(Constraint::Fixed(p));
    }
    sk.solve().expect("rectangle should solve");
    serde_json::to_string(&sk).expect("rectangle should serialize")
}

fn sketch_face_name(op_id: u64) -> EntityName {
    EntityName {
        operation: OperationId(op_id),
        role: EntityRole {
            category: RoleCategory::Generated,
            index: 0,
        },
        ancestors: Vec::new(),
    }
}

fn add_param(tree: &mut FeatureTree, name: &str, value: f64) -> ParameterId {
    tree.parameters_mut()
        .add_global(name, value, Unit::Meters)
        .expect("add param")
}

/// Unit cube via sketch + extrude, replayed. Returns the tree, the
/// cube body, and the extrude feature id.
fn cube_tree() -> (FeatureTree, BRep, FeatureId) {
    let mut tree = FeatureTree::new();
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: SketchPayload::Constraint(rectangle_sketch_json()),
        },
        "Sketch1",
    );
    let dist = add_param(&mut tree, "_extrude_distance", 1.0);
    let extrude = tree.add_feature(
        Feature::Extrude {
            face: sketch_face_name(1),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: dist },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Extrude1",
    );
    tree.replay(Tolerance::default()).expect("cube replays");
    let brep = tree.current_brep().expect("cube body").clone();
    assert_eq!(brep.faces.len(), 6);
    (tree, brep, extrude)
}

fn vertical_edge_names(brep: &BRep) -> Vec<EntityName> {
    let mut names = Vec::new();
    for (handle, edge) in &brep.edges {
        let tangent = popkern_measure::edge_tangent(brep, handle).expect("edge tangent");
        if tangent.z.abs() > 0.99 {
            let nid = edge.name.expect("cube edges are named");
            names.push(brep.name_pool.get(nid).expect("name pool").clone());
        }
    }
    assert_eq!(names.len(), 4, "a cube has four vertical edges");
    names
}

fn new_faces(pre: &BRep, post: &BRep) -> Vec<FaceHandle> {
    let before: std::collections::BTreeSet<_> = pre.faces.keys().collect();
    post.faces.keys().filter(|f| !before.contains(f)).collect()
}

fn direct_evidence(
    pre: &BRep,
    post: &BRep,
    requested: popkern_blend::ContinuityLevel,
) -> popkern_blend::BlendContinuityEvidence {
    popkern_blend::measure_blend_continuity_evidence(
        post,
        &new_faces(pre, post),
        requested,
        &popkern_blend::ContinuityGates::commercial(),
        popkern_history::replay_report::EVIDENCE_SAMPLES_PER_SEAM,
    )
}

fn add_fillet(
    tree: &mut FeatureTree,
    cube: &BRep,
    recipe: FilletRecipeDesc,
    radius: f64,
) -> FeatureId {
    let edges = vertical_edge_names(cube);
    let r = add_param(tree, "_r", radius);
    tree.add_feature(
        Feature::Fillet {
            edges,
            radius: r,
            recipe,
        },
        "Fillet1",
    )
}

fn blend_of(report: &ReplayReport, id: FeatureId) -> BlendEvidenceRecord {
    let feature = report.feature(id).expect("feature reported");
    assert_eq!(feature.outcome, ReplayOutcome::Served, "{feature:?}");
    feature
        .blend
        .clone()
        .expect("blend feature carries evidence")
}

/// A replayed fillet's report carries the continuity measured over
/// the faces it created — identical to measuring the served body
/// directly — and the report is returned beside the BRep AND kept on
/// the tree.
#[test]
fn replayed_fillet_carries_measured_continuity_evidence() {
    let (mut tree, cube, _) = cube_tree();
    let fillet = add_fillet(&mut tree, &cube, FilletRecipeDesc::default(), 0.1);
    let result = tree.replay(Tolerance::default()).expect("replay");
    let body = tree.current_brep().expect("body").clone();

    let labels: Vec<(&str, &str)> = result
        .report
        .features
        .iter()
        .map(|r| (r.label.as_str(), r.operation.as_str()))
        .collect();
    assert_eq!(
        labels,
        [
            ("Sketch1", "sketch"),
            ("Extrude1", "extrude"),
            ("Fillet1", "fillet")
        ],
        "timeline order, one record per feature"
    );
    assert!(
        result.report.features.iter().all(|r| r.elapsed.is_some()),
        "executed features carry wall time"
    );
    assert!(result.report.is_clean(), "{:?}", result.report);
    assert_eq!(tree.last_report(), &result.report);

    let record = blend_of(&result.report, fillet);
    assert_eq!(record.edges_requested, 4);
    assert_eq!(
        record.requested_continuity,
        popkern_blend::ContinuityLevel::G1
    );
    assert!(!record.downgraded());
    let reported = record.continuity.clone().expect("fillets always measure");
    let direct = direct_evidence(&cube, &body, popkern_blend::ContinuityLevel::G1);
    assert_eq!(reported.achieved, popkern_blend::ContinuityLevel::G1);
    assert_eq!(reported, direct, "the report IS the measurement");
    assert!(reported.measured() && !reported.downgraded());
}

/// The G2 recipe reports measured G2; a G1 body asked to be G2 would
/// read as a downgrade — the census a GUI shows next to the fillet.
#[test]
fn g2_recipe_reports_measured_g2_and_no_downgrade() {
    let (mut tree, cube, _) = cube_tree();
    let recipe = FilletRecipeDesc {
        continuity: ContinuityKind::G2,
        ..Default::default()
    };
    let fillet = add_fillet(&mut tree, &cube, recipe, 0.1);
    let result = tree.replay(Tolerance::default()).expect("replay");
    let record = blend_of(&result.report, fillet);
    assert_eq!(
        record.requested_continuity,
        popkern_blend::ContinuityLevel::G2
    );
    let evidence = record.continuity.as_ref().expect("measured");
    assert_eq!(
        evidence.achieved,
        popkern_blend::ContinuityLevel::G2,
        "{evidence:?}"
    );
    assert_eq!(result.report.downgraded_blends().count(), 0);
}

/// Chamfers record the blend (edge count, corner patches) but no
/// continuity — their boundary is a designed crease.
#[test]
fn chamfer_reports_a_crease_record_without_continuity() {
    let (mut tree, cube, _) = cube_tree();
    let edges = vertical_edge_names(&cube);
    let d = add_param(&mut tree, "_d", 0.1);
    let chamfer = tree.add_feature(
        Feature::Chamfer {
            edges,
            distance: d,
            kind: ChamferKindDesc::Equal,
        },
        "Chamfer1",
    );
    let result = tree.replay(Tolerance::default()).expect("replay");
    let record = blend_of(&result.report, chamfer);
    assert_eq!(record.edges_requested, 4);
    assert_eq!(
        record.requested_continuity,
        popkern_blend::ContinuityLevel::G0
    );
    assert!(record.continuity.is_none());
    assert!(!record.downgraded());
    assert!(result.report.is_clean());
}

/// A declined feature reports `Failed` with the structured error's
/// text and the body stays at its pre-feature state; a suppressed
/// feature reports why it was skipped.
#[test]
fn failed_and_suppressed_features_report_their_outcome() {
    let (mut tree, cube, _) = cube_tree();
    // A radius larger than the cube cannot serve.
    let fillet = add_fillet(&mut tree, &cube, FilletRecipeDesc::default(), 10.0);
    let result = tree.replay(Tolerance::default()).expect("replay");
    let report = result.report.feature(fillet).expect("reported");
    match &report.outcome {
        ReplayOutcome::Failed { error } => assert!(!error.is_empty()),
        other => panic!("expected Failed, got {other:?}"),
    }
    assert!(report.blend.is_none(), "a declined blend has no evidence");
    assert_eq!(result.report.failures().count(), 1);
    assert!(!result.report.is_clean());
    assert_eq!(tree.current_brep().expect("body").faces.len(), 6);

    tree.suppress(fillet);
    let result = tree.replay(Tolerance::default()).expect("replay");
    let report = result.report.feature(fillet).expect("reported");
    assert_eq!(
        report.outcome,
        ReplayOutcome::Suppressed {
            reason: SuppressionReason::Explicit
        }
    );
    assert!(report.elapsed.is_none());
    assert!(result.report.is_clean(), "suppression is not a finding");
}

/// A partial replay reports only what it executed; the tree's report
/// keeps the earlier features' records in timeline order.
#[test]
fn partial_replay_keeps_untouched_records_on_the_tree() {
    let (mut tree, cube, _) = cube_tree();
    let fillet = add_fillet(&mut tree, &cube, FilletRecipeDesc::default(), 0.1);
    tree.replay(Tolerance::default()).expect("replay");
    assert_eq!(tree.last_report().features.len(), 3);

    let partial = tree
        .replay_from(fillet, Tolerance::default())
        .expect("replay_from");
    let executed: Vec<&str> = partial
        .report
        .features
        .iter()
        .map(|r| r.label.as_str())
        .collect();
    assert_eq!(
        executed,
        ["Fillet1"],
        "only the resumed range is re-reported"
    );

    let kept: Vec<&str> = tree
        .last_report()
        .features
        .iter()
        .map(|r| r.label.as_str())
        .collect();
    assert_eq!(kept, ["Sketch1", "Extrude1", "Fillet1"]);
    assert_eq!(
        tree.last_report().feature(fillet),
        partial.report.feature(fillet)
    );
}

/// The report is the daemon/MCP payload: it round-trips through JSON
/// with the measured evidence intact.
#[test]
fn replay_report_round_trips_through_json() {
    let (mut tree, cube, _) = cube_tree();
    add_fillet(&mut tree, &cube, FilletRecipeDesc::default(), 0.1);
    tree.replay(Tolerance::default()).expect("replay");
    let json = serde_json::to_string(tree.last_report()).expect("serialize");
    let back: ReplayReport = serde_json::from_str(&json).expect("deserialize");
    assert_eq!(&back, tree.last_report());
    assert!(json.contains("\"achieved\":\"G1\""), "{json}");
}

/// The `.pkn` document is JSON; without correctly-rounded float
/// parsing a parameter value shifted by one ULP on load (found by the
/// report round trip above — seam residuals came back one ULP off).
/// Pin bit-exact document round trips on values with long shortest
/// representations.
#[test]
fn pkn_document_round_trips_parameter_values_bit_exactly() {
    let mut tree = FeatureTree::new();
    let awkward = [
        3.302913498259841e-13,
        6.489005608163026e-14,
        0.1 + 0.2,
        1.0 / 3.0,
        123456.78901234568,
    ];
    let ids: Vec<ParameterId> = awkward
        .iter()
        .enumerate()
        .map(|(i, v)| add_param(&mut tree, &format!("_p{i}"), *v))
        .collect();
    let bytes = tree.save().expect("save");
    let loaded = FeatureTree::load(&bytes).expect("load");
    for (id, expected) in ids.iter().zip(awkward) {
        let got = loaded.parameters().get(*id).expect("param survives").value;
        assert_eq!(
            got.to_bits(),
            expected.to_bits(),
            "parameter {expected:e} must load bit-exactly, got {got:e}"
        );
    }
}
