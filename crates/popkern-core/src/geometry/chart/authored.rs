//! Collapsed rows a construction builds, stated through the service
//! (chart-service.md §6c, W2e).
//!
//! The service owns the question "where is this surface singular, and of
//! which kind". A construction that builds a collapsed row (a polar wedge
//! meeting its fan's hub, a corner cap's apex patch, a ruled fan) answers
//! it at build time by stating the row here; the statement is checked
//! against the descriptor's measured classification, so a construction
//! whose row is not the one it meant (not collapsed, or an apex where it
//! built a smooth hub) is caught where it is built, with the measurement,
//! rather than read wrong by every consumer after it. Consumers read the
//! row from the descriptor like any other.

use super::descriptor::{Axis, ChartDescriptor, SingularKind, SingularRow};
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

/// How a stated collapsed row disagrees with the measured chart.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum AuthoredRowMismatch {
    /// No singular row lies on the stated line.
    #[error("no collapsed row at {axis:?} = {at}: the descriptor measures rows {measured:?}")]
    Missing {
        /// The stated axis.
        axis: Axis,
        /// The stated parameter.
        at: f64,
        /// The rows measured, as (axis, parameter, kind).
        measured: Vec<(Axis, f64, SingularKind)>,
    },
    /// The row is there but of the other kind.
    #[error("the row at {axis:?} = {at} was built {stated:?} and measures {measured:?}")]
    KindDiffers {
        /// The stated axis.
        axis: Axis,
        /// The stated parameter.
        at: f64,
        /// The kind the construction meant.
        stated: SingularKind,
        /// The kind the descriptor measures.
        measured: SingularKind,
    },
}

/// State a collapsed row a construction built on `surface` (the line
/// `axis = at`, of kind `kind`, or of whichever kind its geometry gives
/// when `kind` is `None`: a ruled fan's hub is smooth on a flat fan and an
/// apex otherwise): the descriptor's measured row there, or how the
/// measurement disagrees. A row matches within a billionth of the
/// surface's parameter span along that axis (rows lie on domain
/// boundaries, so the comparison is of two domain ends).
pub fn authored_row(
    surface: &SurfaceGeometry,
    tolerance: &Tolerance,
    (axis, at): (Axis, f64),
    kind: Option<SingularKind>,
) -> Result<SingularRow, AuthoredRowMismatch> {
    let descriptor = ChartDescriptor::of(surface, tolerance);
    let ((u0, u1), (v0, v1)) = descriptor.domain;
    let span = match axis {
        Axis::U => u1 - u0,
        Axis::V => v1 - v0,
    };
    let eps = if span.is_finite() {
        1e-9 * span.abs().max(1.0)
    } else {
        1e-9
    };
    let found = descriptor
        .singular_rows
        .iter()
        .find(|r| r.axis == axis && (r.at - at).abs() <= eps);
    match found {
        None => Err(AuthoredRowMismatch::Missing {
            axis,
            at,
            measured: descriptor
                .singular_rows
                .iter()
                .map(|r| (r.axis, r.at, r.kind))
                .collect(),
        }),
        Some(row) if kind.is_some_and(|k| k != row.kind) => Err(AuthoredRowMismatch::KindDiffers {
            axis,
            at,
            stated: kind.unwrap_or(row.kind),
            measured: row.kind,
        }),
        Some(row) => Ok(*row),
    }
}

#[cfg(test)]
#[path = "authored_tests.rs"]
mod tests;
