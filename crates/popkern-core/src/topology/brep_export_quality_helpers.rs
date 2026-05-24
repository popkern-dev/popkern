use crate::topology::{
    BRep, BrepEdgeSurfaceCompatibilityOptions, BrepTrimRepairOptions,
    analyze_brep_edge_surface_compatibility,
};

pub(super) fn edge_surface_fidelity_blockers(
    brep: &BRep,
    tolerance: Option<f64>,
    trim: BrepTrimRepairOptions,
) -> (usize, f64) {
    let report = analyze_brep_edge_surface_compatibility(
        brep,
        BrepEdgeSurfaceCompatibilityOptions {
            tolerance: tolerance.unwrap_or(f64::INFINITY),
            max_anchor_motion: 0.0,
            samples: trim.samples,
            require_all_compatible: false,
            validate_output: false,
        },
    );
    let blockers = tolerance
        .map(|_| report.incompatible_half_edges)
        .unwrap_or_default();
    (blockers, report.max_deviation)
}

pub(super) fn pcurve_fidelity_blockers(brep: &BRep, tolerance: Option<f64>) -> (usize, f64) {
    let mut blockers = 0usize;
    let mut max_deviation = 0.0_f64;
    for edge in brep.edges.keys() {
        let Some(report) = brep.measure_edge_pcurve_deviation(edge) else {
            continue;
        };
        max_deviation = max_deviation.max(report.max_deviation);
        if let Some(tolerance) = tolerance
            && report.max_deviation > tolerance
        {
            blockers += 1;
        }
    }
    (blockers, max_deviation)
}
