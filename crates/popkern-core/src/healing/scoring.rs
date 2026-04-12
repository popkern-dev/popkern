//! Per-entity quality scoring and model-wide grading.
//!
//! Maps defect severity to normalized \[0,1\] scores, aggregates into
//! composite scores, and classifies models into letter grades.

use crate::healing::diagnostics::{EntityRef, HealingReport, Severity};
use std::collections::HashMap;
use std::fmt;

/// Quality grade classification.
///
/// Grades range from A (best) to F (worst) and are derived from the
/// composite quality score.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum QualityGrade {
    /// Failing quality (composite < 0.40).
    F,
    /// Poor quality (composite 0.40..0.60).
    D,
    /// Acceptable quality (composite 0.60..0.80).
    C,
    /// Good quality (composite 0.80..0.95).
    B,
    /// Excellent quality (composite >= 0.95).
    A,
}

impl fmt::Display for QualityGrade {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            QualityGrade::A => write!(f, "A"),
            QualityGrade::B => write!(f, "B"),
            QualityGrade::C => write!(f, "C"),
            QualityGrade::D => write!(f, "D"),
            QualityGrade::F => write!(f, "F"),
        }
    }
}

/// Per-entity quality score.
///
/// Tracks the normalized score and defect count for a single
/// topological entity referenced by defect records.
#[derive(Debug, Clone)]
pub struct EntityScore {
    /// The entity being scored.
    pub entity: EntityRef,
    /// Normalized quality score: 0.0 (worst) to 1.0 (perfect).
    pub score: f64,
    /// Number of defects affecting this entity.
    pub defect_count: usize,
}

/// Model-wide quality assessment.
///
/// Aggregates per-entity scores into composite metrics and a
/// letter grade classification.
#[derive(Debug, Clone)]
pub struct QualityScore {
    /// Composite quality score: 0.0 to 1.0.
    pub composite: f64,
    /// Letter grade derived from the composite score.
    pub grade: QualityGrade,
    /// The entity with the lowest score, if any entities were scored.
    pub worst_entity: Option<EntityScore>,
    /// Mean of all entity scores.
    pub mean_score: f64,
    /// 5th percentile score (worst 5% of entities).
    pub p5_score: f64,
    /// All per-entity scores, sorted ascending by score.
    pub entity_scores: Vec<EntityScore>,
}

impl fmt::Display for QualityScore {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Grade {}: composite={:.3}, mean={:.3}, p5={:.3}, entities={}",
            self.grade,
            self.composite,
            self.mean_score,
            self.p5_score,
            self.entity_scores.len(),
        )
    }
}

/// Scoring weights for different use cases.
///
/// Controls how much each defect severity level penalizes an entity's
/// quality score. Different profiles emphasize different aspects of
/// model quality.
#[derive(Debug, Clone)]
pub struct ScoringProfile {
    /// Human-readable profile name.
    pub name: String,
    /// Weight applied to critical-severity defects.
    pub critical_weight: f64,
    /// Weight applied to error-severity defects.
    pub error_weight: f64,
    /// Weight applied to warning-severity defects.
    pub warning_weight: f64,
    /// Weight applied to info-severity defects.
    pub info_weight: f64,
}

impl ScoringProfile {
    /// Profile tuned for simulation — heavy weight on geometric accuracy.
    ///
    /// Simulation meshes require tight geometry; even warnings have
    /// noticeable impact on solver convergence.
    pub fn simulation() -> Self {
        Self {
            name: "simulation".into(),
            critical_weight: 1.0,
            error_weight: 0.7,
            warning_weight: 0.4,
            info_weight: 0.1,
        }
    }

    /// Profile tuned for manufacturing — heavy weight on tolerance/degenerate issues.
    ///
    /// Manufacturing processes need valid topology and clean tolerances;
    /// critical defects are absolute show-stoppers.
    pub fn manufacturing() -> Self {
        Self {
            name: "manufacturing".into(),
            critical_weight: 1.0,
            error_weight: 0.5,
            warning_weight: 0.2,
            info_weight: 0.05,
        }
    }

    /// Profile tuned for visualization — heavy weight on continuity/appearance.
    ///
    /// Visual quality degrades gracefully; only critical issues truly
    /// break rendering, and warnings are mostly cosmetic.
    pub fn visualization() -> Self {
        Self {
            name: "visualization".into(),
            critical_weight: 1.0,
            error_weight: 0.3,
            warning_weight: 0.1,
            info_weight: 0.02,
        }
    }
}

/// Penalty value for a given severity level.
///
/// These are the base penalty values before profile weighting is applied.
fn severity_penalty(severity: Severity) -> f64 {
    match severity {
        Severity::Critical => 1.0,
        Severity::Error => 0.5,
        Severity::Warning => 0.2,
        Severity::Info => 0.05,
    }
}

/// Weight multiplier from the profile for a given severity.
fn profile_weight(profile: &ScoringProfile, severity: Severity) -> f64 {
    match severity {
        Severity::Critical => profile.critical_weight,
        Severity::Error => profile.error_weight,
        Severity::Warning => profile.warning_weight,
        Severity::Info => profile.info_weight,
    }
}

/// Derive a letter grade from a composite score.
fn grade_from_score(score: f64) -> QualityGrade {
    if score >= 0.95 {
        QualityGrade::A
    } else if score >= 0.80 {
        QualityGrade::B
    } else if score >= 0.60 {
        QualityGrade::C
    } else if score >= 0.40 {
        QualityGrade::D
    } else {
        QualityGrade::F
    }
}

/// Compute a stable string key for an entity reference.
///
/// Used to group defects by entity in a HashMap. The key format is
/// deterministic but not intended for display.
fn entity_key(entity: &EntityRef) -> String {
    match entity {
        EntityRef::Vertex(h) => format!("v:{h:?}"),
        EntityRef::Edge(h) => format!("e:{h:?}"),
        EntityRef::HalfEdge(h) => format!("he:{h:?}"),
        EntityRef::Face(h) => format!("f:{h:?}"),
        EntityRef::Loop(h) => format!("l:{h:?}"),
        EntityRef::Shell(h) => format!("s:{h:?}"),
        EntityRef::EdgePair(a, b) => format!("ep:{a:?},{b:?}"),
        EntityRef::FacePair(a, b) => format!("fp:{a:?},{b:?}"),
    }
}

/// Compute quality scores from a healing report.
///
/// Builds a map from entities to defects, computes per-entity scores
/// using the profile weights, then aggregates into composite metrics.
///
/// # Algorithm
///
/// 1. Build entity -> defect severity list map.
/// 2. Per entity: score = 1.0 - sum(severity_penalty * profile_weight). Clamped to \[0, 1\].
/// 3. Composite = mean of entity scores (or 1.0 if no entities).
/// 4. Grade from composite: A >= 0.95, B >= 0.80, C >= 0.60, D >= 0.40, F < 0.40.
pub fn score_model(report: &HealingReport, profile: &ScoringProfile) -> QualityScore {
    // Step 1: group defect severities by entity.
    let mut entity_defects: HashMap<String, (EntityRef, Vec<Severity>)> = HashMap::new();

    for defect in report.iter() {
        for entity in &defect.entities {
            let key = entity_key(entity);
            entity_defects
                .entry(key)
                .or_insert_with(|| (entity.clone(), Vec::new()))
                .1
                .push(defect.severity);
        }
    }

    // Step 2: compute per-entity scores.
    let mut entity_scores: Vec<EntityScore> = entity_defects
        .into_values()
        .map(|(entity, severities)| {
            let total_penalty: f64 = severities
                .iter()
                .map(|s| severity_penalty(*s) * profile_weight(profile, *s))
                .sum();
            let score = (1.0 - total_penalty).clamp(0.0, 1.0);
            EntityScore {
                entity,
                score,
                defect_count: severities.len(),
            }
        })
        .collect();

    // Sort ascending by score so worst entities come first.
    entity_scores.sort_by(|a, b| {
        a.score
            .partial_cmp(&b.score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Step 3: aggregate.
    if entity_scores.is_empty() {
        // No entities referenced by defects. If the report has defects
        // with no entity refs, we still penalize based on overall severity.
        let overall_penalty: f64 = report
            .iter()
            .map(|d| severity_penalty(d.severity) * profile_weight(profile, d.severity))
            .sum();
        let composite = (1.0 - overall_penalty).clamp(0.0, 1.0);
        return QualityScore {
            composite,
            grade: grade_from_score(composite),
            worst_entity: None,
            mean_score: composite,
            p5_score: composite,
            entity_scores: Vec::new(),
        };
    }

    let sum: f64 = entity_scores.iter().map(|e| e.score).sum();
    let mean_score = sum / entity_scores.len() as f64;

    // 5th percentile: index at 5% of the sorted array (worst end).
    let p5_idx = ((entity_scores.len() as f64 * 0.05).ceil() as usize).max(1);
    let p5_score = entity_scores[..p5_idx].iter().map(|e| e.score).sum::<f64>() / p5_idx as f64;

    let worst_entity = entity_scores.first().cloned();
    let composite = mean_score;
    let grade = grade_from_score(composite);

    QualityScore {
        composite,
        grade,
        worst_entity,
        mean_score,
        p5_score,
        entity_scores,
    }
}

/// Compare scores before and after repair.
///
/// Returns the difference `after.composite - before.composite`.
/// Positive values indicate improvement, negative values indicate regression.
pub fn score_delta(before: &QualityScore, after: &QualityScore) -> f64 {
    after.composite - before.composite
}

#[cfg(test)]
#[path = "scoring_tests.rs"]
mod tests;
