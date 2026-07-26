#![warn(missing_docs)]
//! popkern-fillet: fillets and chamfers. Entry points:
//! [`fillet_edges`], [`chamfer_edges`], [`vertex_blend_3d`].
pub mod analysis;
pub mod apex_construction;
pub mod blend;
pub mod blend_diagnostics;
pub(crate) mod blend_warnings;
pub mod cancellation;
pub mod cap;
pub mod cap_seam_continuity;
pub mod center_curve;
pub mod center_curve_helpers;
mod center_curve_station;
mod center_curve_sweep_channel;
pub mod center_curve_walker;
pub mod chamfer_2d;
pub mod chamfer_2d_by_angle;
pub(crate) mod chamfer_2d_ray;
pub mod chamfer_curve;
pub(crate) mod collapsed_skip;
pub mod construct;
/// Test-only hooks for gated work-in-progress construction arms.
pub mod construct_test_hooks;
pub mod contact;
pub(crate) mod continuity;
pub(crate) mod corner_2d_helpers;
pub mod curvature_comb;
pub(crate) mod cusp_subdivide;
pub(crate) mod degeneracy;
pub(crate) mod degeneracy_surface;
pub(crate) mod direct_apex_trace;
pub mod edges_2d;
pub mod error;
mod face_blend;
mod face_blend_attach;
mod face_blend_attach_ends;
mod face_blend_attach_geometry;
mod face_blend_attach_surgery;
mod face_blend_attach_through;
mod face_blend_parallel;
mod face_blend_sheet;
pub mod fillet_2d;
pub(crate) mod fillet_2d_plan;
pub(crate) mod fillet_2d_solvers;
pub mod fillet_surface;
pub mod fillet_surface_analytic;
#[cfg(test)]
#[path = "g2_acceptance_tests.rs"]
mod g2_acceptance_tests;
pub mod g2_cross_section;
pub(crate) mod geometry_compute;
pub(crate) mod geometry_compute_caps;
pub(crate) mod geometry_compute_cone;
pub(crate) mod geometry_compute_feasibility;
pub mod highlight_line;
pub mod history;
pub mod identify;
pub(crate) mod identify_helpers;
mod interference_seam;
pub(crate) mod loft_chart_extension;
pub mod modify;
pub mod modify_inplace;
pub(crate) mod modify_inplace_closed_sweep;
pub(crate) mod modify_inplace_geometry;
pub(crate) mod modify_inplace_sweep;
pub(crate) mod naming;
pub mod offset;
pub mod overlap;
mod overlap_strip_interference;
pub mod pcurve_fidelity;
pub mod precheck;
pub mod preview;
pub mod progress;
pub mod quality;
pub mod reflection_line;
pub mod remediation;
pub mod residual_scaffolding;
pub mod sensitivity;
pub mod stale_cap;
pub(crate) mod station_geometry;
mod surface_boundary_evidence;
pub mod surface_validate;
pub mod trim;
pub mod types;
#[cfg(test)]
#[path = "types_tests.rs"]
mod types_tests;
pub mod unfillet;
pub mod unfillet_closed;
pub(crate) mod unfillet_closed_helpers;
pub(crate) mod unfillet_corner_group;
pub(crate) mod unfillet_corner_group_geometry;
pub(crate) mod unfillet_helpers;
pub(crate) mod unfillet_surgery;
pub mod validate;
pub mod vertex_blend;
pub(crate) mod vertex_blend_topology;
pub mod vertex_classifier;
pub mod walker;
pub(crate) mod walker_adaptive;
pub mod zebra;
pub use cancellation::CancellationToken;
pub use chamfer_2d::{chamfer_corner_2d, chamfer_corner_2d_by_angle, chamfer_corners_2d};
pub use edges_2d::*;
pub use error::{
    BlendError, BlendFailure, BlendWarning, Remediation, RemediationAction, RemediationContext,
};
pub use fillet_2d::{fillet_corner_2d, fillet_corners_2d};
pub use history::{
    ChamferOperationRecord, FilletOperationRecord, chamfer_recorded, fillet_recorded,
};
pub use identify::{
    BlendInfo, BlendKind, IdentificationSource, blend_neighbours, blend_radius,
    blend_support_faces, identify_all_blends, identify_blend,
};
pub use modify::{modify_fillet_radius, modify_fillet_radius_via_record};
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};
pub use precheck::{
    BlendQuality, ComputeQualityInput, FilletFeasibility, FilletIssue, FilletRadiusBounds,
    feasible_radius_bounds, fillet_quality, precheck_fillet,
};
pub use preview::{
    CornerPatchPreview, PreviewEdgeGeometry, PreviewGeometry, preview_chamfer,
    preview_corner_patches, preview_fillet, preview_fillet_with_continuity,
};
pub use progress::*;
pub use stale_cap::{StaleCap, detect_stale_caps};
pub use types::*;
pub use unfillet::{UnfilletPlan, compute_unfillet_plan, unfillet_face, unfillet_via_record};
pub use vertex_blend::{vertex_blend_3d, vertex_blend_3d_variable};
/// Fillet the specified edges of a BRep solid. See [`BlendError`] for failure modes.
pub fn fillet_edges(
    brep: &BRep,
    specs: &[FilletSpec],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BlendResult, BlendError> {
    fillet_edges_internal(brep, specs, options, tolerance, tracking, None, None)
}
mod fillet_pipeline;
pub(crate) use fillet_pipeline::fillet_edges_internal;
#[path = "lib_chamfer.rs"]
mod lib_chamfer;
pub use lib_chamfer::chamfer_edges;
/// Fillet the specified edges and bundle the result with post-op
/// [`BlendQuality`] metrics. Same parameter validation and error
/// taxonomy as [`fillet_edges`]; quality comes from the new-face
/// delta via [`precheck::fillet_quality`].
#[path = "lib_convenience.rs"]
mod lib_convenience;
pub use face_blend::{FaceBlendAttach, FaceBlendResult, FaceBlendSpec, face_blend};
use geometry_compute::compute_chamfer_geometry;
pub use lib_convenience::{
    chamfer_edge, fillet_edge, fillet_edges_cancellable, fillet_edges_with_progress,
    fillet_edges_with_quality,
};
use naming::assign_fillet_names;
#[cfg(test)]
#[path = "corner_warning_tests.rs"]
mod corner_warning_tests;
#[cfg(test)]
#[path = "lib_tests.rs"]
mod tests;
