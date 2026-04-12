//! Synthetic test models with planted defects for regression testing.
//!
//! Each test model is paired with a ground-truth defect manifest
//! describing the expected defects by type, severity, and location.
//! The oracle runs the healing analyzer against each model and compares
//! the results to the manifest.

use crate::healing::diagnostics::{DefectCategory, Severity};
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// A planted defect in a test model.
///
/// Describes an intentionally introduced defect so the oracle can
/// verify that the analyzer detects it.
#[derive(Debug, Clone)]
pub struct PlantedDefect {
    /// Category of the planted defect.
    pub category: DefectCategory,
    /// Expected severity.
    pub severity: Severity,
    /// Human-readable description of what was planted.
    pub description: String,
    /// Description of the entity that was corrupted.
    pub entity_description: String,
}

/// A test model with known defects.
///
/// Bundles a BRep model together with its ground-truth defect manifest.
/// The oracle compares analyzer output against the manifest.
#[derive(Debug, Clone)]
pub struct TestModel {
    /// Model name for reporting.
    pub name: String,
    /// The BRep model (possibly corrupted).
    pub brep: BRep,
    /// Expected defects planted in this model.
    pub planted_defects: Vec<PlantedDefect>,
}

/// Result of running the oracle on a single model.
///
/// Indicates whether the analyzer found all planted defects and
/// whether any false positives were reported.
#[derive(Debug, Clone)]
pub struct OracleResult {
    /// Name of the test model.
    pub model_name: String,
    /// Number of planted defects.
    pub planted: usize,
    /// Number of planted defects found by the analyzer.
    pub found: usize,
    /// Planted defects that the analyzer missed.
    pub missed: Vec<PlantedDefect>,
    /// Number of defects reported that don't match any planted defect.
    pub false_positives: usize,
    /// Whether this model passed the oracle check.
    pub passed: bool,
}

// ---------------------------------------------------------------------------
// Test model factories
// ---------------------------------------------------------------------------

/// Create a clean unit cube (no defects expected).
///
/// If cube construction fails, returns an empty BRep as a fallback
/// (which itself may produce topological defects, but the test model
/// is marked as expecting zero defects from a valid cube).
pub fn clean_cube() -> TestModel {
    let tol = Tolerance::default();
    let brep = BRep::unit_cube(tol).unwrap_or_else(|_| BRep::new(tol));
    TestModel {
        name: "clean_cube".into(),
        brep,
        planted_defects: Vec::new(),
    }
}

/// Create a cube with a vertex displaced to introduce a gap.
///
/// Moves the first vertex by a large offset so that the loop closure
/// check and/or edge-surface consistency check detects a geometric
/// deviation.
pub fn cube_with_gap() -> TestModel {
    let tol = Tolerance::default();
    let mut brep = BRep::unit_cube(tol).unwrap_or_else(|_| BRep::new(tol));

    // Displace the first vertex to create a gap.
    if let Some((_, vertex)) = brep.vertices.iter_mut().next() {
        vertex.position.x += 0.1;
    }

    TestModel {
        name: "cube_with_gap".into(),
        brep,
        planted_defects: vec![PlantedDefect {
            category: DefectCategory::Topological,
            severity: Severity::Error,
            description: "Vertex displaced to break loop closure".into(),
            entity_description: "first vertex offset by 0.1 in X".into(),
        }],
    }
}

/// Create a cube with a tolerance hierarchy violation.
///
/// Sets a vertex tolerance larger than its incident edge tolerance,
/// violating the vertex <= edge invariant.
pub fn cube_with_tolerance_violation() -> TestModel {
    let tol = Tolerance::default();
    let mut brep = BRep::unit_cube(tol).unwrap_or_else(|_| BRep::new(tol));

    // Set a vertex tolerance larger than the edge tolerance.
    if let Some((_, vertex)) = brep.vertices.iter_mut().next() {
        vertex.tolerance = Some(1.0);
    }
    // Set a small edge tolerance on the first edge.
    if let Some((_, edge)) = brep.edges.iter_mut().next() {
        edge.tolerance = Some(0.001);
    }

    TestModel {
        name: "cube_with_tolerance_violation".into(),
        brep,
        planted_defects: vec![PlantedDefect {
            category: DefectCategory::Tolerance,
            severity: Severity::Warning,
            description: "Vertex tolerance exceeds edge tolerance".into(),
            entity_description: "vertex tol=1.0 > edge tol=0.001".into(),
        }],
    }
}

/// Create a cube with missing edge geometry.
///
/// Removes the curve geometry from the first edge that has one.
pub fn cube_with_missing_geometry() -> TestModel {
    let tol = Tolerance::default();
    let mut brep = BRep::unit_cube(tol).unwrap_or_else(|_| BRep::new(tol));

    // Remove curve from the first edge.
    if let Some((_, edge)) = brep.edges.iter_mut().next() {
        edge.curve = None;
    }

    TestModel {
        name: "cube_with_missing_geometry".into(),
        brep,
        planted_defects: vec![PlantedDefect {
            category: DefectCategory::Degenerate,
            severity: Severity::Info,
            description: "Edge curve geometry removed".into(),
            entity_description: "first edge has curve set to None".into(),
        }],
    }
}

/// Create a cube with a small/degenerate edge.
///
/// Moves two adjacent vertices very close together so that the
/// edge connecting them is shorter than the linear tolerance.
pub fn cube_with_small_edge() -> TestModel {
    let tol = Tolerance::default();
    let mut brep = BRep::unit_cube(tol).unwrap_or_else(|_| BRep::new(tol));

    // Find the first edge and collapse its endpoints.
    let first_edge_he = brep.edges.iter().next().map(|(_, e)| e.half_edge);
    if let Some(he_handle) = first_edge_he
        && let Some(he) = brep.half_edges.get(he_handle)
    {
        let twin = he.twin;
        let origin = he.origin;
        if let Some(twin_he) = brep.half_edges.get(twin) {
            let twin_origin = twin_he.origin;
            // Get origin position and move twin_origin very close.
            if let Some(v0) = brep.vertices.get(origin) {
                let target = v0.position;
                if let Some(v1) = brep.vertices.get_mut(twin_origin) {
                    // Move to within 1e-8 of the other vertex.
                    v1.position = target;
                    v1.position.x += 1e-8;
                }
            }
        }
    }

    TestModel {
        name: "cube_with_small_edge".into(),
        brep,
        planted_defects: vec![PlantedDefect {
            category: DefectCategory::Degenerate,
            severity: Severity::Warning,
            description: "Edge endpoints collapsed to near-coincident positions".into(),
            entity_description: "first edge length ≈ 1e-8".into(),
        }],
    }
}

/// Get all test models.
///
/// Returns one clean model and one model for each defect type.
pub fn all_test_models() -> Vec<TestModel> {
    vec![
        clean_cube(),
        cube_with_gap(),
        cube_with_tolerance_violation(),
        cube_with_missing_geometry(),
        cube_with_small_edge(),
    ]
}

// ---------------------------------------------------------------------------
// Oracle runner
// ---------------------------------------------------------------------------

/// Run the oracle: analyze each model and compare against planted defects.
///
/// For each test model:
/// 1. Runs `healing::analyze(&model.brep, &tolerance)`.
/// 2. For each planted defect, checks if a matching `DefectRecord` exists
///    (matched by category).
/// 3. Counts found, missed, and false positives.
/// 4. A model passes if all planted defects are found. For the clean model,
///    it additionally must have zero defects.
pub fn run_oracle(tolerance: &Tolerance) -> Vec<OracleResult> {
    let models = all_test_models();
    let mut results = Vec::with_capacity(models.len());

    for model in &models {
        let report = crate::healing::analyze(&model.brep, tolerance);
        let detected_categories: Vec<DefectCategory> = report.iter().map(|d| d.category).collect();

        let mut found = 0usize;
        let mut missed = Vec::new();

        // Track which detected defects are accounted for.
        let mut accounted = vec![false; detected_categories.len()];

        for planted in &model.planted_defects {
            // Find the first unaccounted detected defect matching category.
            let matched = accounted
                .iter()
                .enumerate()
                .position(|(i, &used)| !used && detected_categories[i] == planted.category);

            if let Some(idx) = matched {
                accounted[idx] = true;
                found += 1;
            } else {
                missed.push(planted.clone());
            }
        }

        let false_positives = accounted.iter().filter(|&&used| !used).count();

        let passed = if model.planted_defects.is_empty() {
            // Clean model: must have zero defects.
            report.is_empty()
        } else {
            // Corrupted model: all planted defects must be found.
            missed.is_empty()
        };

        results.push(OracleResult {
            model_name: model.name.clone(),
            planted: model.planted_defects.len(),
            found,
            missed,
            false_positives,
            passed,
        });
    }

    results
}

#[cfg(test)]
#[path = "test_oracle_tests.rs"]
mod tests;
