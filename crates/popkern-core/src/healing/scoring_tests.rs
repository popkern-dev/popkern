use super::*;
use crate::healing::diagnostics::{
    DefectCategory, DefectRecord, EntityRef, HealingReport, Severity,
};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

fn defect_with_entity(severity: Severity, entity: EntityRef) -> DefectRecord {
    DefectRecord {
        category: DefectCategory::Geometric,
        severity,
        description: format!("test defect ({severity})"),
        entities: vec![entity],
        location: None,
        measured_value: None,
        threshold: None,
        remediation: None,
    }
}

fn defect_no_entity(severity: Severity) -> DefectRecord {
    DefectRecord {
        category: DefectCategory::Geometric,
        severity,
        description: format!("test defect ({severity})"),
        entities: Vec::new(),
        location: None,
        measured_value: None,
        threshold: None,
        remediation: None,
    }
}

#[test]
fn valid_cube_gets_grade_a() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("unit cube creation failed");
    let report = crate::healing::analyze(&brep, &Tolerance::default());
    let profile = ScoringProfile::manufacturing();
    let score = score_model(&report, &profile);

    assert_eq!(
        score.grade,
        QualityGrade::A,
        "valid cube should get grade A, got {} (composite={:.3})",
        score.grade,
        score.composite,
    );
    assert!(score.composite >= 0.95);
}

#[test]
fn report_with_critical_defects_gets_grade_f() {
    let mut report = HealingReport::new();
    let entity = EntityRef::Face(Default::default());

    // Push enough critical defects to drive the score well below 0.40.
    for _ in 0..3 {
        report.push(defect_with_entity(Severity::Critical, entity.clone()));
    }

    let profile = ScoringProfile::manufacturing();
    let score = score_model(&report, &profile);

    assert_eq!(
        score.grade,
        QualityGrade::F,
        "3 critical defects should get grade F, got {} (composite={:.3})",
        score.grade,
        score.composite,
    );
    assert!(score.composite < 0.40);
}

#[test]
fn empty_report_is_perfect() {
    let report = HealingReport::new();
    let profile = ScoringProfile::manufacturing();
    let score = score_model(&report, &profile);

    assert_eq!(score.grade, QualityGrade::A);
    assert!((score.composite - 1.0).abs() < f64::EPSILON);
    assert!(score.entity_scores.is_empty());
    assert!(score.worst_entity.is_none());
}

#[test]
fn score_delta_positive_means_improvement() {
    let mut before_report = HealingReport::new();
    before_report.push(defect_with_entity(
        Severity::Error,
        EntityRef::Edge(Default::default()),
    ));

    let after_report = HealingReport::new();

    let profile = ScoringProfile::manufacturing();
    let before = score_model(&before_report, &profile);
    let after = score_model(&after_report, &profile);

    let delta = score_delta(&before, &after);
    assert!(
        delta > 0.0,
        "fixing defects should give positive delta, got {delta}"
    );
}

#[test]
fn score_delta_negative_means_regression() {
    let before_report = HealingReport::new();
    let mut after_report = HealingReport::new();
    after_report.push(defect_with_entity(
        Severity::Critical,
        EntityRef::Edge(Default::default()),
    ));

    let profile = ScoringProfile::manufacturing();
    let before = score_model(&before_report, &profile);
    let after = score_model(&after_report, &profile);

    let delta = score_delta(&before, &after);
    assert!(
        delta < 0.0,
        "adding defects should give negative delta, got {delta}"
    );
}

#[test]
fn grade_boundaries() {
    // Test each grade boundary with precisely crafted reports.
    let profile = ScoringProfile::manufacturing();

    // A single info defect on one entity: penalty = 0.05 * 0.05 = 0.0025.
    // Score = 0.9975 -> grade A.
    let mut info_report = HealingReport::new();
    info_report.push(defect_with_entity(
        Severity::Info,
        EntityRef::Vertex(Default::default()),
    ));
    let score = score_model(&info_report, &profile);
    assert_eq!(score.grade, QualityGrade::A);

    // A single warning: penalty = 0.2 * 0.2 = 0.04. Score = 0.96 -> A.
    let mut warn_report = HealingReport::new();
    warn_report.push(defect_with_entity(
        Severity::Warning,
        EntityRef::Vertex(Default::default()),
    ));
    let score = score_model(&warn_report, &profile);
    assert_eq!(score.grade, QualityGrade::A);

    // A single error: penalty = 0.5 * 0.5 = 0.25. Score = 0.75 -> C.
    let mut err_report = HealingReport::new();
    err_report.push(defect_with_entity(
        Severity::Error,
        EntityRef::Vertex(Default::default()),
    ));
    let score = score_model(&err_report, &profile);
    assert_eq!(score.grade, QualityGrade::C);

    // A single critical: penalty = 1.0 * 1.0 = 1.0. Score = 0.0 -> F.
    let mut crit_report = HealingReport::new();
    crit_report.push(defect_with_entity(
        Severity::Critical,
        EntityRef::Vertex(Default::default()),
    ));
    let score = score_model(&crit_report, &profile);
    assert_eq!(score.grade, QualityGrade::F);
}

#[test]
fn entity_scores_sorted_ascending() {
    let mut report = HealingReport::new();
    // Entity A: 1 critical -> score 0.0
    report.push(defect_with_entity(
        Severity::Critical,
        EntityRef::Edge(Default::default()),
    ));
    // Entity B: 1 info -> high score
    report.push(defect_with_entity(
        Severity::Info,
        EntityRef::Vertex(Default::default()),
    ));

    let profile = ScoringProfile::manufacturing();
    let score = score_model(&report, &profile);

    assert!(score.entity_scores.len() >= 2);
    for window in score.entity_scores.windows(2) {
        assert!(
            window[0].score <= window[1].score,
            "entity scores should be sorted ascending"
        );
    }
}

#[test]
fn worst_entity_has_lowest_score() {
    let mut report = HealingReport::new();
    report.push(defect_with_entity(
        Severity::Critical,
        EntityRef::Face(Default::default()),
    ));
    report.push(defect_with_entity(
        Severity::Info,
        EntityRef::Vertex(Default::default()),
    ));

    let profile = ScoringProfile::manufacturing();
    let score = score_model(&report, &profile);

    let worst = score
        .worst_entity
        .as_ref()
        .expect("should have worst entity");
    assert!(
        worst.score <= score.mean_score,
        "worst entity score should be <= mean"
    );
}

#[test]
fn profiles_have_distinct_names() {
    let sim = ScoringProfile::simulation();
    let mfg = ScoringProfile::manufacturing();
    let viz = ScoringProfile::visualization();

    assert_ne!(sim.name, mfg.name);
    assert_ne!(mfg.name, viz.name);
    assert_ne!(sim.name, viz.name);
}

#[test]
fn defects_without_entity_refs_still_penalize() {
    let mut report = HealingReport::new();
    report.push(defect_no_entity(Severity::Critical));

    let profile = ScoringProfile::manufacturing();
    let score = score_model(&report, &profile);

    assert!(
        score.composite < 1.0,
        "defects without entities should still penalize composite"
    );
}

#[test]
fn score_clamps_to_zero_not_negative() {
    let mut report = HealingReport::new();
    let entity = EntityRef::Edge(Default::default());
    // 10 critical defects on the same entity.
    for _ in 0..10 {
        report.push(defect_with_entity(Severity::Critical, entity.clone()));
    }

    let profile = ScoringProfile::manufacturing();
    let score = score_model(&report, &profile);

    assert!(
        score.composite >= 0.0,
        "score should never go below 0.0, got {}",
        score.composite,
    );
    for es in &score.entity_scores {
        assert!(es.score >= 0.0);
    }
}

#[test]
fn quality_grade_ordering() {
    assert!(QualityGrade::F < QualityGrade::D);
    assert!(QualityGrade::D < QualityGrade::C);
    assert!(QualityGrade::C < QualityGrade::B);
    assert!(QualityGrade::B < QualityGrade::A);
}

#[test]
fn quality_score_display() {
    let score = QualityScore {
        composite: 0.85,
        grade: QualityGrade::B,
        worst_entity: None,
        mean_score: 0.85,
        p5_score: 0.70,
        entity_scores: Vec::new(),
    };
    let s = score.to_string();
    assert!(s.contains("Grade B"));
    assert!(s.contains("0.850"));
}

#[test]
fn quality_grade_display() {
    assert_eq!(QualityGrade::A.to_string(), "A");
    assert_eq!(QualityGrade::F.to_string(), "F");
}
