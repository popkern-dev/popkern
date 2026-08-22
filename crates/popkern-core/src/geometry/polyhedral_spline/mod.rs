//! Polyhedral spline construction for ACM TOMS Algorithm 1032 style control nets.
//!
//! This module owns popkern's implementation model: callers provide an already
//! ordered local control stencil for a supported polyhedral layout, and the
//! constructor returns Bernstein-Bezier patches with explicit family metadata.
//! Mesh/topology adapters sit above this layer.
#![allow(missing_docs)]

mod adapter;
mod atlas;
mod brep_export;
mod catmull_clark;
mod complex;
mod complex_weld;
mod complex_weld_report;
mod complex_weld_scope;
mod construct;
mod coverage;
mod diagnostics;
mod mask;
mod net;
mod patch_subdivision;
mod readiness;
mod tables;
mod types;
mod validate;

pub use adapter::{
    ControlNetSpline, ControlNetSplineOrigin, construct_from_control_net,
    construct_from_control_net_with_origins,
};
pub use construct::{
    construct_extraordinary_star, construct_ngon, construct_polar, construct_regular_quad,
    construct_t0, construct_t1, construct_t2, construct_two_triangles_two_quads,
};
pub use net::{ControlFace, ControlNet, ControlVertex};

/// True when the ACM ToMS 1032 mask tables were embedded at build
/// time (the optional, deliberately untracked package was present —
/// see `PolyhedralSplineError::ToMsTablesUnavailable`). Table-backed
/// construction fails with that structured error when this is false;
/// tests exercising PnS at runtime should skip WITH A NOTICE instead
/// of failing on checkouts without the package.
pub const fn toms_tables_available() -> bool {
    cfg!(toms_package)
}
pub use readiness::{
    PatchDomain, PatchEdgeRole, PatchEdgeUse, PolyhedralSplineBrepReadinessReport,
    validate_brep_readiness,
};
pub use types::{
    PatchEdge, PolyhedralPatch, PolyhedralSpline, PolyhedralSplineError, PolyhedralSplineFamily,
};
pub use validate::{
    ContinuityReport, PatchRegularityReport, PolyhedralSplineValidationReport,
    validate_c0_edge_match, validate_finite_patch, validate_finite_spline,
    validate_patch_interior_regularity, validate_polyhedral_spline,
};

// The PnS test modules embed ACM ToMS fixtures and mask tables; they
// gate on the optional package the same way `tables.rs` does. The
// notice module keeps the skip VISIBLE in test output when the
// package is absent (fresh clones, CI) — the OCC-oracle
// skip-with-notice precedent.
#[cfg(all(test, not(toms_package)))]
mod toms_package_notice {
    #[test]
    fn toms_package_absent_notice() {
        eprintln!(
            "[toms] ACM ToMS 1032 package absent — PnS table tests skipped              (tables.rs returns ToMsTablesUnavailable; see PK-314/PK-315)"
        );
    }
}
#[cfg(all(test, toms_package))]
#[path = "brep_export_tests.rs"]
mod brep_export_tests;
#[cfg(test)]
#[path = "catmull_clark_tests.rs"]
mod catmull_clark_tests;
#[cfg(all(test, toms_package))]
#[path = "complex_tests.rs"]
mod complex_tests;
#[cfg(all(test, toms_package))]
#[path = "coverage_tests.rs"]
mod coverage_tests;
#[cfg(all(test, toms_package))]
#[path = "diagnostics_tests.rs"]
mod diagnostics_tests;
#[cfg(all(test, toms_package))]
#[path = "readiness_tests.rs"]
mod readiness_tests;
#[cfg(all(test, toms_package))]
#[path = "tests.rs"]
mod tests;
pub use atlas::{
    PolyhedralPatchAtlas, PolyhedralSeam, SeamContinuityReport, SeamOrientation, build_patch_atlas,
    validate_spline_atlas,
};
pub use brep_export::{
    PolyhedralBrepExport, PolyhedralBrepExportOptions, PolyhedralBrepFaceProvenance,
    build_polyhedral_spline_brep,
};
pub use complex::{
    PolyhedralComplexSewingReport, PolyhedralSplineRejection, PolyhedralSplineSurfaceComplex,
    PolyhedralSplineSurfaceComponent, analyze_surface_complex_sewing,
    build_control_net_surface_complex,
};
pub use complex_weld::weld_clean_surface_complex_boundaries;
pub use complex_weld_report::{
    PolyhedralBoundaryWeldReport, PolyhedralBoundaryWeldTiming, PolyhedralWeldedSurfaceComplex,
};
pub use complex_weld_scope::PolyhedralBoundaryWeldScope;
pub use coverage::{
    PolyhedralCandidateCoverage, PolyhedralControlNetCoverageReport, PolyhedralCoverageOptions,
    PolyhedralCoverageStatus, analyze_control_net_coverage,
};
pub use diagnostics::{
    PolyhedralSplineDebugArtifact, PolyhedralSplineWorkEstimate,
    build_polyhedral_spline_debug_artifact, estimate_polyhedral_spline_work,
};
