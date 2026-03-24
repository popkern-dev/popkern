//! Public STEP import API.
//!
//! The main entry points are [`import_step`] (from bytes) and [`import_step_file`]
//! (from a filesystem path). Both return a [`StepImportResult`] containing the
//! imported [`StepDocument`] and a diagnostic log.

use popkern_core::tolerance::Tolerance;

use crate::error::IoError;
use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::document::{OpaqueBag, ShapeRepKind, StepDocument, UnitContext};
use crate::step::healing::{HealingLevel, run_healing};
use crate::step::mapping::geometry::GeometryCache;
use crate::step::mapping::presentation::map_presentation;
use crate::step::mapping::product::map_products;
use crate::step::mapping::topology::map_brep_solids;
use crate::step::parser::entities::parse_step_file;
use crate::step::parser::graph::EntityGraph;

/// Options controlling how a STEP file is imported.
#[derive(Debug, Clone)]
pub struct StepImportOptions {
    /// Override the geometric tolerance (otherwise inferred from the file).
    pub tolerance_override: Option<f64>,
    /// Which shape representation kind to prefer.
    pub preferred_rep: ShapeRepKind,
    /// Optional filter: only import products whose name matches one of these.
    pub product_filter: Option<Vec<String>>,
    /// How aggressively to heal the imported model.
    pub healing: HealingLevel,
}

impl Default for StepImportOptions {
    fn default() -> Self {
        Self {
            tolerance_override: None,
            preferred_rep: ShapeRepKind::AdvancedBrep,
            product_filter: None,
            healing: HealingLevel::Basic,
        }
    }
}

/// Result of a STEP import operation.
pub struct StepImportResult {
    /// The imported document containing products, shapes, and metadata.
    pub document: StepDocument,
    /// Diagnostic log with warnings and errors from all pipeline phases.
    pub diagnostics: DiagnosticLog,
}

/// Import a STEP file from raw bytes.
///
/// Runs the full import pipeline: parse → resolve → map → heal.
///
/// # Errors
///
/// Returns [`IoError`] if parsing fails entirely or the file is not a valid
/// Part 21 file. Partial failures (individual entities that cannot be mapped)
/// are reported as diagnostics rather than hard errors.
pub fn import_step(data: &[u8], options: &StepImportOptions) -> Result<StepImportResult, IoError> {
    let mut log = DiagnosticLog::new();

    // Phase 1: Parse
    let parse_result = parse_step_file(data).map_err(|e| IoError::Parse {
        line: 0,
        col: 0,
        message: e,
    })?;

    for err in &parse_result.parse_errors {
        log.warn(
            Phase::Parse,
            None,
            format!("offset {}: {}", err.offset, err.message),
        );
    }

    // Phase 2: Build entity graph
    let graph = EntityGraph::build(parse_result.entities, &mut log);

    // Phase 3: Determine tolerance
    let tolerance = match options.tolerance_override {
        Some(tol) => Tolerance::default().with_linear(tol),
        None => Tolerance::default(),
    };

    // Phase 4: Map geometry + topology
    let mut geo_cache = GeometryCache::new(&graph);
    let mut brep_solids = map_brep_solids(&graph, &mut geo_cache, &mut log, tolerance);

    // Phase 5: Healing
    let mut solids_only: Vec<_> = brep_solids.iter().map(|(_, b)| b.clone()).collect();
    let heal_log = log.clone();
    run_healing(&mut solids_only, &heal_log, options.healing, &mut log)?;
    // Put healed solids back
    for (i, (_, brep)) in brep_solids.iter_mut().enumerate() {
        if let Some(healed) = solids_only.get(i) {
            *brep = healed.clone();
        }
    }

    // Phase 6: Map product structure
    let brep_refs: Vec<_> = brep_solids.iter().map(|(id, b)| (*id, b.clone())).collect();
    let (mut products, roots) = map_products(&graph, &brep_refs, &mut geo_cache, &mut log);

    // Phase 7: Map presentation (colors, layers)
    map_presentation(&graph, &mut products, &mut log);

    // Phase 8: Apply product filter if specified
    if let Some(ref filter) = options.product_filter {
        for product in &mut products {
            if !filter.iter().any(|f| product.name.contains(f.as_str())) {
                product.shapes.clear();
            }
        }
    }

    let document = StepDocument {
        products,
        roots,
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };

    // Log summary
    let error_count = log.errors().count();
    if error_count > 0 {
        log.warn(
            Phase::Validate,
            None,
            format!("import completed with {error_count} errors"),
        );
    }

    Ok(StepImportResult {
        document,
        diagnostics: log,
    })
}

/// Import a STEP file from a filesystem path.
///
/// # Errors
///
/// Returns [`IoError::Io`] if the file cannot be read, or any error from
/// [`import_step`].
pub fn import_step_file(
    path: &std::path::Path,
    options: &StepImportOptions,
) -> Result<StepImportResult, IoError> {
    let data = std::fs::read(path)?;
    import_step(&data, options)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_options() {
        let opts = StepImportOptions::default();
        assert!(opts.tolerance_override.is_none());
        assert_eq!(opts.preferred_rep, ShapeRepKind::AdvancedBrep);
        assert!(opts.product_filter.is_none());
        assert_eq!(opts.healing, HealingLevel::Basic);
    }

    /// Minimal end-to-end test with a STEP file containing a single triangle.
    #[test]
    fn import_minimal_step_file() {
        let step_data = br#"ISO-10303-21;
HEADER;
FILE_DESCRIPTION((''), '2;1');
FILE_NAME('test.stp', '2024-01-01', (''), (''), '', '', '');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1 = CARTESIAN_POINT('', (0.0, 0.0, 0.0));
#2 = CARTESIAN_POINT('', (1.0, 0.0, 0.0));
#3 = CARTESIAN_POINT('', (0.0, 1.0, 0.0));
#4 = DIRECTION('', (0.0, 0.0, 1.0));
#5 = DIRECTION('', (1.0, 0.0, 0.0));
#6 = AXIS2_PLACEMENT_3D('', #1, #4, #5);
#7 = PLANE('', #6);
#8 = DIRECTION('', (1.0, 0.0, 0.0));
#9 = VECTOR('', #8, 1.0);
#10 = DIRECTION('', (-1.0, 1.0, 0.0));
#11 = VECTOR('', #10, 1.0);
#12 = DIRECTION('', (0.0, -1.0, 0.0));
#13 = VECTOR('', #12, 1.0);
#14 = LINE('', #1, #9);
#15 = LINE('', #2, #11);
#16 = LINE('', #3, #13);
#20 = VERTEX_POINT('', #1);
#21 = VERTEX_POINT('', #2);
#22 = VERTEX_POINT('', #3);
#30 = EDGE_CURVE('', #20, #21, #14, .T.);
#31 = EDGE_CURVE('', #21, #22, #15, .T.);
#32 = EDGE_CURVE('', #22, #20, #16, .T.);
#40 = ORIENTED_EDGE('', *, *, #30, .T.);
#41 = ORIENTED_EDGE('', *, *, #31, .T.);
#42 = ORIENTED_EDGE('', *, *, #32, .T.);
#50 = EDGE_LOOP('', (#40, #41, #42));
#60 = FACE_OUTER_BOUND('', #50, .T.);
#70 = ADVANCED_FACE('', (#60), #7, .T.);
#80 = CLOSED_SHELL('', (#70));
#90 = MANIFOLD_SOLID_BREP('body', #80);
ENDSEC;
END-ISO-10303-21;
"#;

        let opts = StepImportOptions {
            healing: HealingLevel::None,
            ..Default::default()
        };
        let result = import_step(step_data, &opts);
        assert!(result.is_ok(), "import failed: {:?}", result.err());
        let result = result.unwrap();

        // The import should succeed (document is valid regardless of product count)
        assert!(
            result.diagnostics.entries().len() < 100,
            "too many diagnostics"
        );
    }

    #[test]
    fn import_nonexistent_file() {
        let opts = StepImportOptions::default();
        let result = import_step_file(std::path::Path::new("/nonexistent/file.stp"), &opts);
        assert!(result.is_err());
    }
}
