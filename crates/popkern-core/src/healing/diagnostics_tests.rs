//! Additional tests for healing diagnostics types.
//!
//! Covers: DefectRecord without measured_value, HealingReport::into_defects,
//! HealingReport::defects, RemediationAction variants, EntityRef variants,
//! DefectLocation variants, Remediation struct, for_vertex/for_edge/for_face
//! entity filtering.

use super::*;
use crate::geometry::Point3;

fn sample_defect(severity: Severity, category: DefectCategory) -> DefectRecord {
    DefectRecord {
        category,
        severity,
        description: format!("Test {category} defect"),
        entities: Vec::new(),
        location: None,
        measured_value: None,
        threshold: None,
        remediation: None,
    }
}

#[test]
fn defect_record_without_measured_value_display() {
    let d = DefectRecord {
        category: DefectCategory::Topological,
        severity: Severity::Error,
        description: "Broken twin".into(),
        entities: Vec::new(),
        location: None,
        measured_value: None,
        threshold: None,
        remediation: None,
    };
    let s = d.to_string();
    assert!(s.contains("Broken twin"));
    // Without measured/threshold, the parenthetical suffix should be absent.
    assert!(!s.contains("measured"));
    assert!(!s.contains("threshold"));
}

#[test]
fn defect_record_with_only_measured_no_threshold() {
    let d = DefectRecord {
        category: DefectCategory::Geometric,
        severity: Severity::Warning,
        description: "Partial data".into(),
        entities: Vec::new(),
        location: None,
        measured_value: Some(0.1),
        threshold: None,
        remediation: None,
    };
    let s = d.to_string();
    assert!(s.contains("Partial data"));
    // Both must be Some for the parenthetical to appear.
    assert!(!s.contains("measured"));
}

#[test]
fn healing_report_into_defects() {
    let mut r = HealingReport::new();
    r.push(sample_defect(Severity::Warning, DefectCategory::Geometric));
    let defects = r.into_defects();
    assert_eq!(defects.len(), 1);
    assert_eq!(defects[0].severity, Severity::Warning);
}

#[test]
fn healing_report_defects_slice() {
    let mut r = HealingReport::new();
    r.push(sample_defect(Severity::Error, DefectCategory::Topological));
    assert_eq!(r.defects().len(), 1);
    assert_eq!(r.defects()[0].category, DefectCategory::Topological);
}

#[test]
fn healing_report_defects_slice_empty() {
    let r = HealingReport::new();
    assert!(r.defects().is_empty());
}

#[test]
fn remediation_action_variants_match() {
    let variants = [
        RemediationAction::AbsorbTolerance,
        RemediationAction::AdjustGeometry,
        RemediationAction::RemoveEntity,
        RemediationAction::SplitEntity,
        RemediationAction::SimplifySurface,
        RemediationAction::EnforceToleranceHierarchy,
        RemediationAction::RestorePeriodicTopology,
        RemediationAction::RestoreSingularityTopology,
        RemediationAction::RecomputeEdgeCurve,
        RemediationAction::ReparameterizeSurface,
        RemediationAction::RecomputePCurve,
    ];
    // All variants should be distinguishable.
    assert!(matches!(variants[0], RemediationAction::AbsorbTolerance));
    assert!(matches!(variants[1], RemediationAction::AdjustGeometry));
    assert!(matches!(variants[10], RemediationAction::RecomputePCurve));
}

#[test]
fn remediation_struct_construction() {
    let r = Remediation {
        action: RemediationAction::AdjustGeometry,
        parameter: Some("max_displacement".into()),
        value: Some(0.01),
        rationale: "Close a gap".into(),
    };
    assert_eq!(r.action, RemediationAction::AdjustGeometry);
    assert_eq!(r.parameter.as_deref(), Some("max_displacement"));
    assert!((r.value.unwrap_or(0.0) - 0.01).abs() < 1e-15);
}

#[test]
fn entity_ref_vertex_variant() {
    use slotmap::SlotMap;
    let mut sm: SlotMap<VertexHandle, ()> = SlotMap::with_key();
    let h = sm.insert(());
    let e = EntityRef::Vertex(h);
    assert!(matches!(e, EntityRef::Vertex(_)));
}

#[test]
fn entity_ref_edge_variant() {
    use slotmap::SlotMap;
    let mut sm: SlotMap<EdgeHandle, ()> = SlotMap::with_key();
    let h = sm.insert(());
    let e = EntityRef::Edge(h);
    assert!(matches!(e, EntityRef::Edge(_)));
}

#[test]
fn entity_ref_halfedge_variant() {
    use slotmap::SlotMap;
    let mut sm: SlotMap<HalfEdgeHandle, ()> = SlotMap::with_key();
    let h = sm.insert(());
    let e = EntityRef::HalfEdge(h);
    assert!(matches!(e, EntityRef::HalfEdge(_)));
}

#[test]
fn entity_ref_face_variant() {
    use slotmap::SlotMap;
    let mut sm: SlotMap<FaceHandle, ()> = SlotMap::with_key();
    let h = sm.insert(());
    let e = EntityRef::Face(h);
    assert!(matches!(e, EntityRef::Face(_)));
}

#[test]
fn entity_ref_loop_variant() {
    use slotmap::SlotMap;
    let mut sm: SlotMap<LoopHandle, ()> = SlotMap::with_key();
    let h = sm.insert(());
    let e = EntityRef::Loop(h);
    assert!(matches!(e, EntityRef::Loop(_)));
}

#[test]
fn entity_ref_shell_variant() {
    use slotmap::SlotMap;
    let mut sm: SlotMap<ShellHandle, ()> = SlotMap::with_key();
    let h = sm.insert(());
    let e = EntityRef::Shell(h);
    assert!(matches!(e, EntityRef::Shell(_)));
}

#[test]
fn entity_ref_edge_pair_variant() {
    use slotmap::SlotMap;
    let mut sm: SlotMap<EdgeHandle, ()> = SlotMap::with_key();
    let h1 = sm.insert(());
    let h2 = sm.insert(());
    let e = EntityRef::EdgePair(h1, h2);
    assert!(matches!(e, EntityRef::EdgePair(_, _)));
}

#[test]
fn entity_ref_face_pair_variant() {
    use slotmap::SlotMap;
    let mut sm: SlotMap<FaceHandle, ()> = SlotMap::with_key();
    let h1 = sm.insert(());
    let h2 = sm.insert(());
    let e = EntityRef::FacePair(h1, h2);
    assert!(matches!(e, EntityRef::FacePair(_, _)));
}

#[test]
fn defect_location_point() {
    let loc = DefectLocation::Point(Point3::new(1.0, 2.0, 3.0));
    assert!(matches!(loc, DefectLocation::Point(_)));
}

#[test]
fn defect_location_parameter_point() {
    let loc = DefectLocation::ParameterPoint { u: 0.5, v: 0.5 };
    assert!(matches!(loc, DefectLocation::ParameterPoint { .. }));
}

#[test]
fn defect_location_parameter_region() {
    let loc = DefectLocation::ParameterRegion {
        u_range: (0.0, 1.0),
        v_range: (0.0, 1.0),
    };
    assert!(matches!(loc, DefectLocation::ParameterRegion { .. }));
}

#[test]
fn defect_location_bounding_box() {
    use crate::geometry::intersection::surface_subdivision::Aabb;
    let aabb = Aabb {
        min: Point3::new(0.0, 0.0, 0.0),
        max: Point3::new(1.0, 1.0, 1.0),
    };
    let loc = DefectLocation::BoundingBox(aabb);
    assert!(matches!(loc, DefectLocation::BoundingBox(_)));
}

#[test]
fn defect_category_display_all_variants() {
    assert_eq!(DefectCategory::Topological.to_string(), "topological");
    assert_eq!(DefectCategory::Geometric.to_string(), "geometric");
    assert_eq!(DefectCategory::Degenerate.to_string(), "degenerate");
    assert_eq!(DefectCategory::Tolerance.to_string(), "tolerance");
    assert_eq!(DefectCategory::Import.to_string(), "import");
}

#[test]
fn severity_display_all_variants() {
    assert_eq!(Severity::Info.to_string(), "info");
    assert_eq!(Severity::Warning.to_string(), "warning");
    assert_eq!(Severity::Error.to_string(), "error");
    assert_eq!(Severity::Critical.to_string(), "critical");
}

#[test]
fn worst_severity_returns_critical() {
    let mut report = HealingReport::new();
    report.push(sample_defect(Severity::Info, DefectCategory::Geometric));
    report.push(sample_defect(Severity::Critical, DefectCategory::Tolerance));
    report.push(sample_defect(Severity::Warning, DefectCategory::Import));
    assert_eq!(report.worst_severity(), Some(Severity::Critical));
}

#[test]
fn for_vertex_finds_matching_defects() {
    use slotmap::SlotMap;
    let mut sm: SlotMap<VertexHandle, ()> = SlotMap::with_key();
    let v1 = sm.insert(());
    let v2 = sm.insert(());

    let mut report = HealingReport::new();
    report.push(DefectRecord {
        category: DefectCategory::Topological,
        severity: Severity::Error,
        description: "Bad vertex".into(),
        entities: vec![EntityRef::Vertex(v1)],
        location: None,
        measured_value: None,
        threshold: None,
        remediation: None,
    });
    report.push(sample_defect(Severity::Warning, DefectCategory::Geometric));

    assert_eq!(report.for_vertex(v1).len(), 1);
    assert_eq!(report.for_vertex(v2).len(), 0);
}

#[test]
fn for_edge_finds_matching_defects() {
    use slotmap::SlotMap;
    let mut sm: SlotMap<EdgeHandle, ()> = SlotMap::with_key();
    let e1 = sm.insert(());
    let e2 = sm.insert(());

    let mut report = HealingReport::new();
    report.push(DefectRecord {
        category: DefectCategory::Geometric,
        severity: Severity::Warning,
        description: "Edge gap".into(),
        entities: vec![EntityRef::Edge(e1)],
        location: None,
        measured_value: None,
        threshold: None,
        remediation: None,
    });

    assert_eq!(report.for_edge(e1).len(), 1);
    assert_eq!(report.for_edge(e2).len(), 0);
}

#[test]
fn for_face_finds_matching_defects() {
    use slotmap::SlotMap;
    let mut sm: SlotMap<FaceHandle, ()> = SlotMap::with_key();
    let f1 = sm.insert(());
    let f2 = sm.insert(());

    let mut report = HealingReport::new();
    report.push(DefectRecord {
        category: DefectCategory::Degenerate,
        severity: Severity::Info,
        description: "Small face".into(),
        entities: vec![EntityRef::Face(f1)],
        location: None,
        measured_value: None,
        threshold: None,
        remediation: None,
    });

    assert_eq!(report.for_face(f1).len(), 1);
    assert_eq!(report.for_face(f2).len(), 0);
}

#[test]
fn healing_report_iter() {
    let mut report = HealingReport::new();
    report.push(sample_defect(Severity::Info, DefectCategory::Geometric));
    report.push(sample_defect(
        Severity::Warning,
        DefectCategory::Topological,
    ));
    let collected: Vec<_> = report.iter().collect();
    assert_eq!(collected.len(), 2);
}

#[test]
fn healing_summary_default_is_all_zeros() {
    let s = HealingSummary::default();
    assert_eq!(s.total, 0);
    assert_eq!(s.critical, 0);
    assert_eq!(s.error, 0);
    assert_eq!(s.warning, 0);
    assert_eq!(s.info, 0);
}

#[test]
fn defect_record_with_remediation_display() {
    let d = DefectRecord {
        category: DefectCategory::Geometric,
        severity: Severity::Warning,
        description: "Gap found".into(),
        entities: Vec::new(),
        location: None,
        measured_value: Some(0.01),
        threshold: Some(0.001),
        remediation: Some(Remediation {
            action: RemediationAction::AdjustGeometry,
            parameter: None,
            value: None,
            rationale: "Close the gap".into(),
        }),
    };
    let s = d.to_string();
    assert!(s.contains("Gap found"));
    assert!(s.contains("0.010000"));
    assert!(s.contains("0.001000"));
}
