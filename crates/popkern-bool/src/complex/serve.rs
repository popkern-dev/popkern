//! Serving a boolean's fragment selection from the cell complex
//! (PK-976).

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;

use super::model::CellComplex;
use super::pipeline::{ComplexInput, classify_with_complex};
use super::split_plan::{SplitPlan, split_plan};
use crate::error::{
    BoolDiagnostics, BoolError, BoolOperation, ContactKind, NonManifoldContact, Stage,
};
use crate::merge::fragment_classify::FragmentSelection;

/// Which classifier decides a cellular boolean's fragments.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SelectionSource {
    /// Per-fragment point classification (the cellular classifier).
    Cellular,
    /// Region membership in the merged cell complex.
    Complex,
    /// The complex, falling back to the cellular classifier when the
    /// complex declines (a non-manifold verdict is final).
    ComplexFirst,
    /// The complex, a non-manifold boundary split into manifold
    /// pieces by its split plan (`NonManifoldPolicy::SplitToManifold`).
    ComplexSplit,
}

/// A point on an edge of the complex: its curve's midpoint, else the
/// chord's.
fn edge_point(c: &CellComplex, edge: usize) -> Point3 {
    let e = &c.edges[edge];
    match &e.curve {
        Some(cv) => cv.geometry.point_at(0.5 * (cv.t_start + cv.t_end)),
        None => {
            let (a, b) = (
                c.vertices[e.start.index()].position,
                c.vertices[e.end.index()].position,
            );
            a + (b - a) * 0.5
        }
    }
}

fn failure(error: BoolError, message: String) -> BoolDiagnostics {
    let mut diag = BoolDiagnostics {
        error: Some(error),
        failed_stage: Some(Stage::Merge),
        ..Default::default()
    };
    diag.add(crate::error::DiagnosticEntry {
        severity: crate::error::Severity::Error,
        phase: Stage::Merge,
        message,
        entities: Vec::new(),
        measured_value: None,
        location: None,
        remediation: None,
    });
    diag
}

/// Whether the complex may serve where the cellular classifier could
/// not (`POPKERN_COMPLEX=0` is the diagnostic kill switch).
pub(crate) fn complex_rescue_enabled() -> bool {
    std::env::var_os("POPKERN_COMPLEX").is_none_or(|v| v != "0")
}

/// Which classifier decides the primary cellular attempt.
/// `POPKERN_COMPLEX_PRIMARY=1` makes it the complex (plan Stage 7c:
/// complex primary, cellular as the fallback); the default is the
/// cellular classifier with the complex as the rescue.
pub(crate) fn primary_source() -> SelectionSource {
    if std::env::var_os("POPKERN_COMPLEX_PRIMARY").is_some_and(|v| v == "1") {
        SelectionSource::ComplexFirst
    } else {
        SelectionSource::Cellular
    }
}

/// Serve an unscoped cellular boolean with complex selections, under
/// the cellular acceptance rules (valid body, complete pcurves).
/// `Ok(None)`: the complex declined or its body failed acceptance —
/// the caller falls back. `Err`: the result is non-manifold.
pub(crate) fn rescue(
    a: &popkern_core::topology::BRep,
    b: &popkern_core::topology::BRep,
    intersection: &crate::intersect::IntersectionResult,
    op: BoolOperation,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<Option<(popkern_core::topology::BRep, crate::serve_region::Carried)>, BoolError> {
    rescue_from(a, b, intersection, op, tolerance, SelectionSource::Complex)
}

/// `rescue` with the non-manifold boundary split into manifold pieces
/// by the complex's split plan (`NonManifoldPolicy::SplitToManifold`).
/// The pieces still meet at their contacts: the caller records them.
pub(crate) fn rescue_split(
    a: &popkern_core::topology::BRep,
    b: &popkern_core::topology::BRep,
    intersection: &crate::intersect::IntersectionResult,
    op: BoolOperation,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<Option<(popkern_core::topology::BRep, crate::serve_region::Carried)>, BoolError> {
    rescue_from(
        a,
        b,
        intersection,
        op,
        tolerance,
        SelectionSource::ComplexSplit,
    )
}

fn rescue_from(
    a: &popkern_core::topology::BRep,
    b: &popkern_core::topology::BRep,
    intersection: &crate::intersect::IntersectionResult,
    op: BoolOperation,
    tolerance: &popkern_core::tolerance::Tolerance,
    source: SelectionSource,
) -> Result<Option<(popkern_core::topology::BRep, crate::serve_region::Carried)>, BoolError> {
    use crate::complex_shadow::{ComplexRescueOutcome, record_complex_rescue};
    let merged = crate::cellular::boolean::cellular_graph_merge(
        a,
        b,
        intersection,
        op,
        tolerance,
        None,
        source,
    );
    match merged {
        Ok((mut r, carried)) => {
            if let Err(e) = r.validate() {
                record_complex_rescue(op, ComplexRescueOutcome::Invalid(e.to_string()));
                return Ok(None);
            }
            let _ = r.recompute_pcurves_where_missing();
            if crate::operation::first_missing_pcurve_face(&r).is_some() {
                record_complex_rescue(
                    op,
                    ComplexRescueOutcome::Invalid("pcurve-incomplete".into()),
                );
                return Ok(None);
            }
            record_complex_rescue(op, ComplexRescueOutcome::Served);
            Ok(Some((r, carried)))
        }
        Err(diag) => match diag.into_error() {
            nm @ BoolError::NonManifoldResult { .. } => {
                record_complex_rescue(op, ComplexRescueOutcome::NonManifold);
                Err(nm)
            }
            other => {
                record_complex_rescue(op, ComplexRescueOutcome::Declined(other.to_string()));
                Ok(None)
            }
        },
    }
}

/// The complex's fragment selections for `op`, or why it cannot serve:
/// `NonManifoldResult` naming the contacts when the selected boundary
/// is non-manifold (unless `split` asks for the split plan, which then
/// comes with the selections), a merge failure carrying the complex
/// artifact when the complex declined.
pub(crate) fn complex_selections(
    input: &ComplexInput<'_>,
    op: BoolOperation,
    split: bool,
) -> Result<(Vec<FragmentSelection>, Option<SplitPlan>), BoolDiagnostics> {
    match classify_with_complex(input, op) {
        Ok(verdict) => {
            let s = &verdict.selection;
            if s.nonmanifold_edges.is_empty() && s.nonmanifold_vertices.is_empty() {
                return Ok((verdict.selection.selections, None));
            }
            let c = &verdict.lift.complex;
            let plan =
                split.then(|| split_plan(c, &verdict.sides, &s.kept_face, &s.selected_shell));
            let declined = match plan {
                Some(Ok(plan)) => return Ok((verdict.selection.selections, Some(plan))),
                Some(Err(decline)) => format!("; split declined: {decline:?}"),
                None => String::new(),
            };
            let mut contacts: Vec<NonManifoldContact> = s
                .nonmanifold_edges
                .iter()
                .map(|e| NonManifoldContact {
                    kind: ContactKind::Edge,
                    location: edge_point(c, e.index()),
                })
                .collect();
            contacts.extend(s.nonmanifold_vertices.iter().map(|v| NonManifoldContact {
                kind: ContactKind::Vertex,
                location: c.vertices[v.index()].position,
            }));
            let message = format!(
                "complex: selected boundary is non-manifold ({} edge(s), {} vertex(es)){declined}",
                s.nonmanifold_edges.len(),
                s.nonmanifold_vertices.len()
            );
            Err(failure(BoolError::NonManifoldResult { contacts }, message))
        }
        Err(decline) => {
            let message = format!("complex declined: {decline:?}");
            if popkern_core::trace::trace_enabled("complex") {
                eprintln!("[complex] {message}");
            }
            Err(failure(
                BoolError::MergeFailed {
                    stage: "Complex".to_string(),
                    reason: message.clone(),
                },
                message,
            ))
        }
    }
}
