use popkern_core::error::ValidationError;
use popkern_core::geometry::CurveEval;
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle};

use crate::step_validation::StepRoundTripBodyReport;

#[derive(Debug, Clone, PartialEq)]
pub struct StepRoundTripFailure {
    pub label: String,
    pub reason: String,
    pub imported_validation: Option<Box<StepImportedValidationFailure>>,
    pub export_quality: Option<Box<StepRoundTripBodyReport>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct StepImportedValidationFailure {
    pub body_index: usize,
    pub kind: StepImportedValidationFailureKind,
}

#[derive(Debug, Clone, PartialEq)]
pub enum StepImportedValidationFailureKind {
    ShellGap(StepImportedShellGapEvidence),
    Other { message: String },
}

#[derive(Debug, Clone, PartialEq)]
pub struct StepImportedShellGapEvidence {
    pub edge: EdgeHandle,
    pub gap: f64,
    pub half_edge: HalfEdgeHandle,
    pub twin: HalfEdgeHandle,
    pub half_edge_origin: String,
    pub twin_origin: String,
    pub curve_kind: String,
    pub curve_t_start: f64,
    pub curve_t_end: f64,
    pub curve_start: String,
    pub curve_end: String,
    pub natural_start_gap: f64,
    pub natural_end_gap: f64,
    pub natural_gap: f64,
    pub swapped_start_gap: f64,
    pub swapped_end_gap: f64,
    pub swapped_gap: f64,
}

impl std::fmt::Display for StepRoundTripFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.label, self.reason)
    }
}

impl std::error::Error for StepRoundTripFailure {}

pub(super) fn imported_validation_failure(
    body_index: usize,
    body: &BRep,
    err: &ValidationError,
) -> StepImportedValidationFailure {
    StepImportedValidationFailure {
        body_index,
        kind: match err {
            ValidationError::ShellGap { edge, gap } => StepImportedValidationFailureKind::ShellGap(
                imported_shell_gap_evidence(body, *edge, *gap),
            ),
            _ => StepImportedValidationFailureKind::Other {
                message: err.to_string(),
            },
        },
    }
}

fn imported_shell_gap_evidence(
    body: &BRep,
    edge: EdgeHandle,
    gap: f64,
) -> StepImportedShellGapEvidence {
    let edge_data = body.edges.get(edge);
    let half_edge = edge_data
        .map(|edge_data| edge_data.half_edge)
        .unwrap_or_default();
    let twin = body
        .half_edges
        .get(half_edge)
        .map(|half_edge| half_edge.twin)
        .unwrap_or_default();
    let half_edge_origin = body
        .half_edges
        .get(half_edge)
        .and_then(|half_edge| body.vertices.get(half_edge.origin))
        .map(|vertex| vertex.position);
    let twin_origin = body
        .half_edges
        .get(twin)
        .and_then(|half_edge| body.vertices.get(half_edge.origin))
        .map(|vertex| vertex.position);
    let curve = edge_data.and_then(|edge_data| edge_data.curve.as_ref());
    let curve_start = curve.map(|curve| curve.geometry.point_at(curve.t_start));
    let curve_end = curve.map(|curve| curve.geometry.point_at(curve.t_end));

    let natural_start_gap = distance_or_infinity(half_edge_origin, curve_start);
    let natural_end_gap = distance_or_infinity(twin_origin, curve_end);
    let swapped_start_gap = distance_or_infinity(half_edge_origin, curve_end);
    let swapped_end_gap = distance_or_infinity(twin_origin, curve_start);

    StepImportedShellGapEvidence {
        edge,
        gap,
        half_edge,
        twin,
        half_edge_origin: point_label(half_edge_origin),
        twin_origin: point_label(twin_origin),
        curve_kind: curve
            .map(|curve| curve_kind(&curve.geometry))
            .unwrap_or("missing_curve")
            .to_string(),
        curve_t_start: curve.map(|curve| curve.t_start).unwrap_or(f64::NAN),
        curve_t_end: curve.map(|curve| curve.t_end).unwrap_or(f64::NAN),
        curve_start: point_label(curve_start),
        curve_end: point_label(curve_end),
        natural_start_gap,
        natural_end_gap,
        natural_gap: natural_start_gap.max(natural_end_gap),
        swapped_start_gap,
        swapped_end_gap,
        swapped_gap: swapped_start_gap.max(swapped_end_gap),
    }
}

fn distance_or_infinity(
    a: Option<popkern_core::geometry::Point3>,
    b: Option<popkern_core::geometry::Point3>,
) -> f64 {
    match (a, b) {
        (Some(a), Some(b)) => a.distance_to(&b),
        _ => f64::INFINITY,
    }
}

fn point_label(point: Option<popkern_core::geometry::Point3>) -> String {
    point
        .map(|point| format!("{:.9},{:.9},{:.9}", point.x, point.y, point.z))
        .unwrap_or_else(|| "missing".to_string())
}

fn curve_kind(curve: &popkern_core::geometry::CurveGeometry) -> &'static str {
    match curve {
        popkern_core::geometry::CurveGeometry::Line { .. } => "line",
        popkern_core::geometry::CurveGeometry::Circle { .. } => "circle",
        popkern_core::geometry::CurveGeometry::Arc { .. } => "arc",
        popkern_core::geometry::CurveGeometry::Ellipse { .. } => "ellipse",
        popkern_core::geometry::CurveGeometry::Parabola { .. } => "parabola",
        popkern_core::geometry::CurveGeometry::Hyperbola { .. } => "hyperbola",
        popkern_core::geometry::CurveGeometry::Nurbs(_) => "nurbs",
    }
}

pub(super) fn failure(label: &str, reason: String) -> StepRoundTripFailure {
    StepRoundTripFailure {
        label: label.into(),
        reason,
        imported_validation: None,
        export_quality: None,
    }
}

pub(super) fn failure_with_imported_validation(
    label: &str,
    reason: String,
    imported_validation: StepImportedValidationFailure,
) -> StepRoundTripFailure {
    StepRoundTripFailure {
        label: label.into(),
        reason,
        imported_validation: Some(Box::new(imported_validation)),
        export_quality: None,
    }
}

pub(super) fn failure_with_export_quality(
    label: &str,
    reason: String,
    export_quality: StepRoundTripBodyReport,
) -> StepRoundTripFailure {
    StepRoundTripFailure {
        label: label.into(),
        reason,
        imported_validation: None,
        export_quality: Some(Box::new(export_quality)),
    }
}
