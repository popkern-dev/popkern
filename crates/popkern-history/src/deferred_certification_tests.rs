//! Deferred seam certification: a deferred replay serves on sampled
//! evidence, and certifying every deferred seam reproduces the inline
//! record exactly.

use popkern_blend::{BlendContinuityEvidence, ContinuityCertification};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sketch::{Constraint, Sketch};

use super::DeferredCertificationError;
use crate::FeatureTree;
use crate::feature::{ExtrudeEnd, Feature, OperationMode, SketchPlaneDesc, Vec3Desc};
use crate::feature_blend::{ContinuityKind, FilletRecipeDesc};
use crate::node::FeatureId;
use crate::unit::Unit;

fn entity_name(op_id: u64, category: RoleCategory, index: u32) -> EntityName {
    EntityName {
        operation: OperationId(op_id),
        role: EntityRole { category, index },
        ancestors: Vec::new(),
    }
}

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
    sk.add_constraint(Constraint::Fixed(p0));
    sk.add_constraint(Constraint::Fixed(p1));
    sk.add_constraint(Constraint::Fixed(p2));
    sk.add_constraint(Constraint::Fixed(p3));
    serde_json::to_string(&sk).expect("sketch json")
}

/// A unit box with one vertical edge rounded at G2.
fn g2_box(policy: ContinuityCertification) -> (FeatureTree, FeatureId) {
    let mut tree = FeatureTree::new();
    tree.set_continuity_certification(policy);
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: crate::SketchPayload::Constraint(rectangle_sketch_json()),
        },
        "sketch",
    );
    let height = tree
        .parameters_mut()
        .add_global("height", 1.0, Unit::Meters)
        .expect("height param");
    let _extrude = tree.add_feature(
        Feature::Extrude {
            face: entity_name(1, RoleCategory::Generated, 0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: height },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "extrude",
    );
    tree.replay(Tolerance::default()).expect("base replay");
    let edge = sharp_edge_name(
        &tree
            .brep_at(2, Tolerance::default())
            .expect("extrude state"),
    );
    let r = tree
        .parameters_mut()
        .add_global("r", 0.2, Unit::Meters)
        .expect("radius param");
    let fillet = tree.add_feature(
        Feature::Fillet {
            edges: vec![edge],
            radius: r,
            recipe: FilletRecipeDesc {
                continuity: ContinuityKind::G2,
                ..Default::default()
            },
        },
        "round",
    );
    tree.replay(Tolerance::default()).expect("replay");
    (tree, fillet)
}

/// A named edge between two non-parallel planes (a box edge).
fn sharp_edge_name(state: &BRep) -> EntityName {
    state
        .name_pool
        .names()
        .map(|(_, n)| n)
        .find(|n| {
            crate::resolve::resolve_edge(n, state).is_ok_and(|eh| {
                let ed = &state.edges[eh];
                let hed = &state.half_edges[ed.half_edge];
                let twin = &state.half_edges[hed.twin];
                let normal = |lh| {
                    let f = state.loops.get(lh).map(|l| l.face)?;
                    match &state.faces.get(f)?.surface.as_ref()?.geometry {
                        SurfaceGeometry::Plane { normal, .. } => Some(*normal),
                        _ => None,
                    }
                };
                ed.curve.is_some()
                    && matches!(
                        (normal(hed.loop_), normal(twin.loop_)),
                        (Some(a), Some(b)) if a.dot(&b).abs() < 0.9
                    )
            })
        })
        .cloned()
        .expect("a named box edge")
}

fn continuity(tree: &FeatureTree, feature: FeatureId) -> BlendContinuityEvidence {
    tree.last_report()
        .feature(feature)
        .and_then(|row| row.blend.as_ref())
        .and_then(|blend| blend.continuity.clone())
        .expect("fillet continuity evidence")
}

#[test]
fn inline_is_the_default_and_defers_nothing() {
    let (tree, fillet) = g2_box(ContinuityCertification::Inline);
    assert_eq!(
        FeatureTree::new().continuity_certification(),
        ContinuityCertification::Inline
    );
    let evidence = continuity(&tree, fillet);
    assert_eq!(tree.deferred_seam_count(), 0);
    assert_eq!(evidence.deferred_seams, 0);
    assert!(evidence.certified_seams > 0, "{evidence:?}");
    assert!(tree.deferred.bodies.is_empty());
}

#[test]
fn deferred_replay_serves_on_sampled_evidence_without_a_certified_claim() {
    let (tree, fillet) = g2_box(ContinuityCertification::Deferred);
    let evidence = continuity(&tree, fillet);
    assert!(evidence.deferred_seams > 0, "{evidence:?}");
    assert_eq!(tree.deferred_seam_count(), evidence.deferred_seams);
    assert_eq!(evidence.certified_seams, 0);
    assert_eq!(
        evidence.certified_achieved, None,
        "no certified claim while deferred"
    );
    assert!(
        evidence
            .seams
            .iter()
            .all(|s| s.certification_deferred == (s.gated && s.certificate.is_none())),
        "deferred exactly on the gated seams, none certified"
    );
    assert!(tree.deferred.bodies.contains_key(&fillet));
}

#[test]
fn certifying_every_deferred_seam_reproduces_the_inline_record() {
    let (inline, fillet) = g2_box(ContinuityCertification::Inline);
    let (mut deferred, deferred_fillet) = g2_box(ContinuityCertification::Deferred);
    let mut remaining = deferred.deferred_seam_count();
    let mut steps = 0;
    while let Some(progress) = deferred.certify_next_deferred_seam().expect("certify") {
        assert_eq!(progress.feature, deferred_fillet);
        assert_eq!(progress.remaining, remaining - 1, "one seam per call");
        remaining = progress.remaining;
        steps += 1;
    }
    assert!(steps > 0);
    assert_eq!(deferred.deferred_seam_count(), 0);
    assert_eq!(
        continuity(&deferred, deferred_fillet),
        continuity(&inline, fillet),
        "deferred-then-certified must equal inline, certificate for certificate"
    );
    assert!(
        deferred.deferred.bodies.is_empty(),
        "a fully certified feature releases its BRep"
    );
}

#[test]
fn a_replay_replaces_kept_bodies_with_its_own_evidence() {
    let (mut tree, fillet) = g2_box(ContinuityCertification::Deferred);
    tree.set_continuity_certification(ContinuityCertification::Inline);
    tree.replay(Tolerance::default()).expect("replay");
    assert_eq!(tree.deferred_seam_count(), 0);
    assert!(
        !tree.deferred.bodies.contains_key(&fillet),
        "an inline replay leaves no deferred BRep behind"
    );
    assert_eq!(tree.certify_next_deferred_seam(), Ok(None));
}

#[test]
fn a_missing_measured_body_is_reported_not_certified_elsewhere() {
    let (mut tree, fillet) = g2_box(ContinuityCertification::Deferred);
    tree.deferred.bodies.clear();
    assert_eq!(
        tree.certify_next_deferred_seam(),
        Err(DeferredCertificationError::BodyUnavailable {
            feature: fillet,
            label: "round".into(),
        })
    );
}
