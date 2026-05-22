use std::sync::Arc;

use popkern_core::topology::{
    BRep, BrepEdgeSurfaceCompatibilityOptions, BrepExportQualityOptions, BrepExportReadiness,
    BrepTopologyOutputPolicy, BrepTrimRepairOptions, EdgeHandle, FaceHandle, HalfEdgeHandle,
    analyze_brep_edge_surface_compatibility, analyze_brep_export_quality,
};
use popkern_io::{
    HealingLevel, OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind,
    StepDocument, StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step,
    import_step,
};

#[derive(Debug, Clone)]
pub struct StepRoundTripValidationOptions {
    pub import: StepImportOptions,
    pub export: StepExportOptions,
    pub quality: BrepExportQualityOptions,
}

impl Default for StepRoundTripValidationOptions {
    fn default() -> Self {
        Self {
            import: StepImportOptions {
                healing: HealingLevel::None,
                ..Default::default()
            },
            export: StepExportOptions::default(),
            quality: BrepExportQualityOptions {
                topology_output: BrepTopologyOutputPolicy::OpenSheet,
                trim: BrepTrimRepairOptions {
                    endpoint_tolerance: 1e-6,
                    sampled_tolerance: 1e-6,
                    samples: 5,
                    validate_output: true,
                },
                require_edge_curves: true,
                allow_collapsed_trim_edge_curves: true,
                max_pcurve_deviation: Some(1e-6),
                max_edge_surface_deviation: Some(1e-6),
                validate_brep: true,
            },
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct StepRoundTripReport {
    pub label: String,
    pub bytes: usize,
    pub expected_bodies: usize,
    pub imported_bodies: usize,
    pub body_reports: Vec<StepRoundTripBodyReport>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepStrictAcceptanceStatus {
    Accepted,
    Declined,
}

#[derive(Debug, Clone, PartialEq)]
pub struct StepStrictAcceptanceReport {
    pub label: String,
    pub status: StepStrictAcceptanceStatus,
    pub failure_reason: Option<String>,
    pub pre_export_body_reports: Vec<StepRoundTripBodyReport>,
    pub round_trip_report: Option<StepRoundTripReport>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct StepRoundTripBodyReport {
    pub index: usize,
    pub vertices: usize,
    pub edges: usize,
    pub faces: usize,
    pub boundary_edges: usize,
    pub readiness: BrepExportReadiness,
    pub topology_blockers: usize,
    pub trim_blockers: usize,
    pub edge_curve_blockers: usize,
    pub pcurve_blockers: usize,
    pub max_pcurve_deviation: f64,
    pub pcurve_evidence: Vec<StepRoundTripEdgePcurveEvidence>,
    pub edge_surface_blockers: usize,
    pub max_edge_surface_deviation: f64,
    pub edge_surface_evidence: Vec<StepRoundTripEdgeSurfaceEvidence>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct StepRoundTripEdgePcurveEvidence {
    pub edge: EdgeHandle,
    pub max_deviation: f64,
    pub half_edges: Vec<StepRoundTripHalfEdgePcurveEvidence>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct StepRoundTripHalfEdgePcurveEvidence {
    pub half_edge: HalfEdgeHandle,
    pub face: FaceHandle,
    pub max_deviation: f64,
    pub same_fraction_deviation: f64,
    pub trace_deviation: f64,
    pub worst_fraction: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct StepRoundTripEdgeSurfaceEvidence {
    pub edge: EdgeHandle,
    pub half_edge: HalfEdgeHandle,
    pub face: FaceHandle,
    pub max_deviation: f64,
    pub max_curve_parameter: f64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StepRoundTripFailure {
    pub label: String,
    pub reason: String,
}

impl std::fmt::Display for StepRoundTripFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.label, self.reason)
    }
}

impl std::error::Error for StepRoundTripFailure {}

pub fn export_and_validate_breps(
    label: &str,
    breps: &[BRep],
    options: &StepRoundTripValidationOptions,
) -> Result<(Vec<u8>, StepRoundTripReport), StepRoundTripFailure> {
    let bytes = export_step(&step_document(label, breps.to_vec()), &options.export)
        .map_err(|err| failure(label, format!("STEP export failed: {err}")))?;
    let report = validate_step_bytes(label, &bytes, breps.len(), options)?;
    Ok((bytes, report))
}

pub fn export_validate_and_import_breps(
    label: &str,
    breps: &[BRep],
    options: &StepRoundTripValidationOptions,
) -> Result<(Vec<u8>, StepRoundTripReport, Vec<BRep>), StepRoundTripFailure> {
    let bytes = export_step(&step_document(label, breps.to_vec()), &options.export)
        .map_err(|err| failure(label, format!("STEP export failed: {err}")))?;
    let imported = import_step(&bytes, &options.import)
        .map_err(|err| failure(label, format!("STEP import failed: {err}")))?;
    let report = validate_imported_document(label, bytes.len(), breps.len(), &imported, options)?;
    let bodies = imported_bodies(&imported.document)
        .into_iter()
        .cloned()
        .collect();
    Ok((bytes, report, bodies))
}

pub fn validate_strict_step_acceptance(
    label: &str,
    breps: &[BRep],
    options: &StepRoundTripValidationOptions,
) -> StepStrictAcceptanceReport {
    let pre_export_body_reports = breps
        .iter()
        .enumerate()
        .map(|(index, body)| step_body_report(index, body, options))
        .collect();
    match export_validate_and_import_breps(label, breps, options) {
        Ok((_bytes, report, _bodies)) => StepStrictAcceptanceReport {
            label: label.to_string(),
            status: StepStrictAcceptanceStatus::Accepted,
            failure_reason: None,
            pre_export_body_reports,
            round_trip_report: Some(report),
        },
        Err(err) => StepStrictAcceptanceReport {
            label: label.to_string(),
            status: StepStrictAcceptanceStatus::Declined,
            failure_reason: Some(err.reason),
            pre_export_body_reports,
            round_trip_report: None,
        },
    }
}

pub fn inspect_brep_for_step(
    index: usize,
    body: &BRep,
    options: &StepRoundTripValidationOptions,
) -> StepRoundTripBodyReport {
    step_body_report(index, body, options)
}

pub fn validate_step_bytes(
    label: &str,
    bytes: &[u8],
    expected_bodies: usize,
    options: &StepRoundTripValidationOptions,
) -> Result<StepRoundTripReport, StepRoundTripFailure> {
    let imported = import_step(bytes, &options.import)
        .map_err(|err| failure(label, format!("STEP import failed: {err}")))?;
    validate_imported_document(label, bytes.len(), expected_bodies, &imported, options)
}

fn validate_imported_document(
    label: &str,
    byte_count: usize,
    expected_bodies: usize,
    imported: &popkern_io::StepImportResult,
    options: &StepRoundTripValidationOptions,
) -> Result<StepRoundTripReport, StepRoundTripFailure> {
    let errors: Vec<_> = imported
        .diagnostics
        .errors()
        .map(|entry| format!("{:?}: {}", entry.phase, entry.message))
        .collect();
    if !errors.is_empty() {
        return Err(failure(
            label,
            format!(
                "STEP import diagnostics contain errors: {}",
                errors.join("; ")
            ),
        ));
    }

    let bodies = imported_bodies(&imported.document);
    if bodies.len() != expected_bodies {
        return Err(failure(
            label,
            format!(
                "imported body count {} does not match expected {expected_bodies}",
                bodies.len()
            ),
        ));
    }

    let mut body_reports = Vec::with_capacity(bodies.len());
    for (index, body) in bodies.into_iter().enumerate() {
        body.validate().map_err(|err| {
            failure(
                label,
                format!("STEP import body {index} failed validation: {err}"),
            )
        })?;
        let quality = analyze_brep_export_quality(body, options.quality);
        if !quality.is_ready() {
            return Err(failure(
                label,
                format!("body {index} export-quality blocked: {quality:?}"),
            ));
        }
        body_reports.push(step_body_report(index, body, options));
    }

    Ok(StepRoundTripReport {
        label: label.into(),
        bytes: byte_count,
        expected_bodies,
        imported_bodies: body_reports.len(),
        body_reports,
    })
}

fn step_body_report(
    index: usize,
    body: &BRep,
    options: &StepRoundTripValidationOptions,
) -> StepRoundTripBodyReport {
    let quality = analyze_brep_export_quality(body, options.quality);
    StepRoundTripBodyReport {
        index,
        vertices: body.vertices.len(),
        edges: body.edges.len(),
        faces: body.faces.len(),
        boundary_edges: body.boundary_edges.len(),
        readiness: quality.readiness,
        topology_blockers: quality.topology_blockers,
        trim_blockers: quality.trim_blockers,
        edge_curve_blockers: quality.edge_curve_blockers,
        pcurve_blockers: quality.pcurve_blockers,
        max_pcurve_deviation: quality.max_pcurve_deviation,
        pcurve_evidence: pcurve_evidence(body, options.quality.max_pcurve_deviation),
        edge_surface_blockers: quality.edge_surface_blockers,
        max_edge_surface_deviation: quality.max_edge_surface_deviation,
        edge_surface_evidence: edge_surface_evidence(body, options.quality),
    }
}

fn pcurve_evidence(body: &BRep, tolerance: Option<f64>) -> Vec<StepRoundTripEdgePcurveEvidence> {
    let mut evidence: Vec<_> = body
        .edges
        .keys()
        .filter_map(|edge| body.measure_edge_pcurve_deviation(edge))
        .filter(|report| {
            tolerance
                .map(|tolerance| report.max_deviation > tolerance)
                .unwrap_or(report.max_deviation > 0.0)
        })
        .map(|report| StepRoundTripEdgePcurveEvidence {
            edge: report.edge,
            max_deviation: report.max_deviation,
            half_edges: report
                .half_edges
                .into_iter()
                .map(|half_edge| StepRoundTripHalfEdgePcurveEvidence {
                    half_edge: half_edge.half_edge,
                    face: half_edge.face,
                    max_deviation: half_edge.max_deviation,
                    same_fraction_deviation: half_edge.max_same_fraction_deviation,
                    trace_deviation: half_edge.max_trace_deviation,
                    worst_fraction: half_edge.worst_fraction,
                })
                .collect(),
        })
        .collect();
    evidence.sort_by(|a, b| b.max_deviation.total_cmp(&a.max_deviation));
    evidence
}

fn edge_surface_evidence(
    body: &BRep,
    quality: BrepExportQualityOptions,
) -> Vec<StepRoundTripEdgeSurfaceEvidence> {
    let report = analyze_brep_edge_surface_compatibility(
        body,
        BrepEdgeSurfaceCompatibilityOptions {
            tolerance: quality.max_edge_surface_deviation.unwrap_or(f64::INFINITY),
            max_anchor_motion: 0.0,
            samples: quality.trim.samples,
            require_all_compatible: false,
            validate_output: false,
        },
    );
    let mut evidence: Vec<_> = report
        .findings
        .into_iter()
        .map(|finding| StepRoundTripEdgeSurfaceEvidence {
            edge: finding.edge,
            half_edge: finding.half_edge,
            face: finding.face,
            max_deviation: finding.max_deviation,
            max_curve_parameter: finding.max_curve_parameter,
        })
        .collect();
    evidence.sort_by(|a, b| b.max_deviation.total_cmp(&a.max_deviation));
    evidence
}

pub fn imported_bodies(document: &StepDocument) -> Vec<&BRep> {
    document
        .products
        .iter()
        .flat_map(|product| &product.shapes)
        .flat_map(|shape| &shape.solids)
        .collect()
}

fn step_document(label: &str, breps: Vec<BRep>) -> StepDocument {
    StepDocument {
        products: vec![StepProductNode {
            name: label.into(),
            id: label.into(),
            children: Vec::new(),
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::Other("SHAPE_REPRESENTATION".into()),
                solids: breps,
                context: RepresentationContext::default(),
                non_brep_items: Vec::new(),
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    }
}

fn failure(label: &str, reason: String) -> StepRoundTripFailure {
    StepRoundTripFailure {
        label: label.into(),
        reason,
    }
}
