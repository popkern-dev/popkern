//! Hard boundary-row constraints for topology-owned cap/support samples.

use super::{BoundaryRowSamples, BoundaryRowTargetSource, CapSupportBoundarySolveReport, row_fit};

pub(super) fn hard_constraints_for_row(
    samples: &BoundaryRowSamples,
    report: &mut CapSupportBoundarySolveReport,
) -> Vec<row_fit::BoundaryRowFitConstraint> {
    let mut constraints: Vec<row_fit::BoundaryRowFitConstraint> = Vec::new();
    for source_sample in &samples.source_samples {
        if source_sample.target_source == BoundaryRowTargetSource::SupportPcurve {
            continue;
        }
        if let Some(existing) = constraints
            .iter()
            .find(|existing| (existing.parameter - source_sample.parameter).abs() <= 1.0e-9)
        {
            if existing.target.distance_to(&source_sample.target) > 1.0e-8 {
                report.hard_constraint_conflicts += 1;
            }
            continue;
        }
        constraints.push(row_fit::BoundaryRowFitConstraint {
            parameter: source_sample.parameter,
            target: source_sample.target,
        });
    }
    if constraints.len() > 10 {
        report.hard_constraint_overfull_rows += 1;
        return Vec::new();
    }
    constraints
}
