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

    // Parse
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

    // Build entity graph
    let graph = EntityGraph::build(parse_result.entities, &mut log);

    // The file's units: geometry is mapped in them, angles converted at
    // the source and lengths converted once the document is mapped.
    let units = crate::step::mapping::units::resolve_file_units(&graph, &mut log);

    // Determine tolerance from override (millimetres), STEP uncertainty
    // (file units), or default
    let default_linear = Tolerance::default().linear() / units.mm_per_length;
    // Geometry maps at the strict default (file pcurves and edges must
    // meet it to be used as written; otherwise they are rebuilt); the
    // file's uncertainty, when looser, is the precision its writer
    // guarantees and is honoured only where it RESCUES data: a vertex
    // within it of its curve end is that end, and an edge end gap
    // within it is recorded as that edge's tolerance (PK-639). Used as
    // the mapping tolerance it accepted OCC's 2e-3 pcurves over exact
    // rebuilt ones (offset_bracket: seven seams G0 -> open).
    let (tolerance, allowance) = match options.tolerance_override {
        Some(tol) => {
            let t = tol / units.mm_per_length;
            (Tolerance::default().with_linear(t), t)
        }
        None => {
            let file = crate::step::import_uncertainty::extract_step_uncertainty(
                &graph,
                units.mm_per_length,
            );
            let allowance = file.map_or(default_linear, |u| u.max(default_linear));
            if allowance > default_linear {
                log.warn(
                    Phase::MapTopology,
                    None,
                    format!(
                        "trim allowance widened to the file's uncertainty {:.2e} mm",
                        allowance * units.mm_per_length
                    ),
                );
            }
            (Tolerance::default().with_linear(default_linear), allowance)
        }
    };

    // Map geometry + topology
    let mut geo_cache = GeometryCache::new(&graph);
    geo_cache.set_radians_per_angle(units.radians_per_angle);
    geo_cache.set_trim_allowance(allowance);
    let mut brep_solids = map_brep_solids(&graph, &mut geo_cache, &mut log, tolerance);

    // Healing
    let mut solids_only: Vec<_> = brep_solids.iter().map(|(_, b)| b.clone()).collect();
    let heal_log = log.clone();
    run_healing(&mut solids_only, &heal_log, options.healing, &mut log)?;
    // Put healed solids back and measure edge deviations for precision tracking
    for (i, (id, brep)) in brep_solids.iter_mut().enumerate() {
        if let Some(healed) = solids_only.get(i) {
            *brep = healed.clone();
        }
        // Faces whose half-edges carry no pcurve (writers that omit
        // them, or whose pcurves the mapper could not use) tessellate
        // their boundaries face by face, and neighbouring faces' samples
        // never meet: measured, the plain imported unit cube carried 24
        // open mesh edges. Synthesize the missing pcurves by projecting
        // each edge onto its face, so every face samples a shared edge
        // the same way (PK-636); recorded, never silent.
        let missing = brep
            .half_edges
            .values()
            .filter(|h| h.pcurve.is_none())
            .count();
        if missing > 0 {
            let synthesized = brep.recompute_pcurves_where_missing();
            let left = brep
                .half_edges
                .values()
                .filter(|h| h.pcurve.is_none())
                .count();
            log.warn(
                Phase::MapTopology,
                Some(*id),
                format!(
                    "synthesized {} missing pcurve(s) by projection ({left} left{})",
                    missing - left,
                    synthesized
                        .err()
                        .map(|e| format!("; {e}"))
                        .unwrap_or_default()
                ),
            );
        }
        // Deviations are measured once every half-edge has its pcurve:
        // an edge without one measures zero, so measuring first left
        // the synthesized edges' file imprecision untracked and part
        // validation refused them (measured: hose-fitting, vaccase,
        // weldment bodies off by 6e-5 to 4e-4).
        crate::step::import_precision::measure_import_deviations(
            brep,
            default_linear,
            allowance,
            &mut log,
        );
        // Import is a serving boundary for UNTRUSTED external data:
        // the stored orientation flag is derived state (loop winding
        // × surface), so a file whose `same_sense` contradicts its
        // own loop winding (dishonest exporter or mapping defect)
        // must not ship that lie into downstream operators. Derive
        // the honest value and record the repair (PK-245).
        let flipped = brep.reconcile_orientation_flags();
        if flipped > 0 {
            log.warn(
                Phase::MapTopology,
                Some(*id),
                format!(
                    "reconciled {flipped} orientation flag(s) contradicting loop winding \
                     (file same_sense was dishonest)"
                ),
            );
        }
    }

    // Map product structure
    let brep_refs: Vec<_> = brep_solids.iter().map(|(id, b)| (*id, b.clone())).collect();
    let (mut products, mut roots) = map_products(&graph, &brep_refs, &mut geo_cache, &mut log);
    // A file with solids but no product structure (some writers emit
    // bare geometry): one product holding every solid, so the
    // geometry is never silently dropped (PK-635).
    if products
        .iter()
        .all(|p| p.shapes.iter().all(|s| s.solids.is_empty()))
        && !brep_refs.is_empty()
    {
        log.warn(
            Phase::MapProduct,
            None,
            "no product carries the file's solids: importing them as one product",
        );
        products.push(crate::step::document::StepProductNode {
            name: "geometry".into(),
            id: "geometry".into(),
            children: Vec::new(),
            shapes: vec![std::sync::Arc::new(crate::step::document::ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: brep_refs.iter().map(|(_, b)| b.clone()).collect(),
                context: crate::step::document::RepresentationContext::default(),
                non_brep_items: Vec::new(),
            })],
            metadata: Default::default(),
        });
        roots = vec![products.len() - 1];
    }

    crate::step::mapping::units::convert_lengths_to_mm(
        &mut products,
        units.mm_per_length,
        &mut log,
    );

    // Map presentation (colors, layers)
    map_presentation(&graph, &mut products, &mut log);

    // Apply product filter if specified
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
#[path = "import_malformed_tests.rs"]
mod malformed_tests;

#[cfg(test)]
mod import_precision_tests {
    use super::*;
    use popkern_core::geometry::pcurve::Pcurve;
    use popkern_core::precision::EntityHandle;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    #[test]
    fn import_deviation_tracking_expands_precision_envelope() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let edge = brep.edges.keys().next().unwrap();
        let half_edge = brep.edges[edge].half_edge;
        brep.half_edges[half_edge].pcurve =
            Some(Pcurve::linear((0.25, 0.25), (1.25, 0.25), 0.0, 1.0));
        let dev = brep.measure_edge_deviation(edge);
        assert!(dev > brep.max_imprecision);
        let mut log = DiagnosticLog::new();

        let record_above = brep.tolerance.linear();
        crate::step::import_precision::measure_import_deviations(
            &mut brep,
            record_above,
            record_above,
            &mut log,
        );

        let tracked = brep.edge_tolerance(edge);
        assert!(tracked >= dev);
        assert!(brep.max_imprecision > tracked);
        assert!(brep.validate_precision_fast().is_ok());
        assert!(
            log.warnings()
                .any(|entry| { entry.message.contains("import precision envelope") })
        );
        assert!(matches!(
            brep.precision_summary().worst_entity,
            Some((EntityHandle::Edge(e), _)) if e == edge
        ));
    }
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
