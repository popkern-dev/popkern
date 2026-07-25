//! Smooth solids from control meshes (polyhedral-net splines).
//!
//! A [`ControlNet`] — an arbitrary quad/n-gon cage of the kind a
//! subdivision modeller produces — is turned into a watertight BRep of
//! smooth spline patches, without the caller touching kernel
//! internals.
//!
//! ```text
//!   ControlNet  →  per-vertex / per-face polyhedral splines
//!               →  surface complex (one component per spline)
//!               →  boundary weld
//!               →  validated Part
//! ```
//!
//! The construction is the Karčiauskas–Peters polyhedral-net spline
//! family: `RegularQuad`, `ExtraordinaryStar`, `Ngon`, `Polar`, plus
//! the `T0`/`T1`/`T2` and `TwoTrianglesTwoQuads` transition
//! configurations. Every patch is a Bézier surface, so the result is
//! ordinary BRep geometry that booleans, fillets, tessellation, and
//! STEP export all consume.
//!
//! # Why this is a Part-level entry point
//!
//! The underlying machinery has existed in popkern-core for a long
//! time with no production caller — reachable only by assembling four
//! core calls in the right order and knowing which options each
//! needs. That is not an API a GUI can drive. This module is the
//! single call, and it returns a [`Part`], so a smooth body is
//! immediately filletable, boolean-able, and exportable like any
//! other.
//!
//! # Declining
//!
//! A control net can be unsuitable — non-manifold, a valence the
//! spline family does not cover, or a cage whose components cannot be
//! welded into a closed shell. Those are reported as
//! [`PartError::SmoothMesh`] carrying the per-spline rejections and
//! the weld report, not as a silently-degraded body.

use popkern_core::geometry::polyhedral_spline::{
    ControlNet, PolyhedralCoverageOptions, build_control_net_surface_complex,
    weld_clean_surface_complex_boundaries,
};
use popkern_core::tolerance::Tolerance;

use crate::error::PartError;
use crate::part::Part;

/// How a control net should be turned into a smooth body.
#[derive(Debug, Clone, Copy)]
pub struct SmoothMeshOptions {
    /// Samples per patch used when checking patch readiness and
    /// C⁰ agreement between neighbours.
    pub samples: usize,
    /// Positional agreement required between two patch boundaries
    /// before they are welded into one edge. This is a geometric
    /// tolerance in model units and is passed explicitly, never read
    /// from global state.
    pub weld_tolerance: f64,
    /// Whether an unwelded boundary edge is an error.
    ///
    /// `true` (the default) demands a closed body: any boundary edge
    /// left after welding fails the operation. Set `false` to accept
    /// an open sheet body, e.g. for a cage that is deliberately not
    /// closed.
    pub require_closed: bool,
}

impl Default for SmoothMeshOptions {
    fn default() -> Self {
        Self {
            samples: 8,
            weld_tolerance: 1e-6,
            require_closed: true,
        }
    }
}

/// What the construction produced, beyond the body itself.
#[derive(Debug, Clone)]
pub struct SmoothMeshReport {
    /// Control-net vertices consumed.
    pub source_vertices: usize,
    /// Control-net faces consumed.
    pub source_faces: usize,
    /// Splines the family produced from the net.
    pub candidate_splines: usize,
    /// Splines that became surface components.
    pub accepted_components: usize,
    /// Splines the family declined, with the reason for each.
    pub rejections: Vec<String>,
    /// Boundary pairs welded into shared edges.
    pub welded_pairs: usize,
    /// Boundary pairs the weld declined.
    pub declined_pairs: usize,
    /// Boundary edges still unshared after welding. Zero for a
    /// closed body.
    pub remaining_boundary_edges: usize,
    /// Largest positional gap closed by the weld.
    pub max_weld_gap: f64,
}

impl SmoothMeshReport {
    /// Whether every spline the net implied was realised.
    pub fn complete(&self) -> bool {
        self.rejections.is_empty() && self.declined_pairs == 0
    }
}

/// Build a smooth solid from a control net.
///
/// # Errors
///
/// [`PartError::SmoothMesh`] when the net cannot be turned into a
/// spline complex, when welding leaves the body open while
/// [`SmoothMeshOptions::require_closed`] is set, or when the welded
/// body fails kernel validation. The error carries the report so the
/// caller can see which splines or boundaries were responsible.
pub fn smooth_solid_from_control_net(
    net: &ControlNet,
    tolerance: Tolerance,
    options: SmoothMeshOptions,
) -> Result<(Part, SmoothMeshReport), PartError> {
    let coverage = PolyhedralCoverageOptions {
        tolerance,
        samples: options.samples,
        c0_tolerance: options.weld_tolerance,
        ..PolyhedralCoverageOptions::default()
    };

    let complex = build_control_net_surface_complex(net, coverage).map_err(|source| {
        PartError::SmoothMesh {
            stage: "spline construction",
            detail: source.to_string(),
        }
    })?;

    let rejections: Vec<String> = complex
        .rejections
        .iter()
        .map(|r| format!("{:?} ({:?}): {}", r.origin, r.family, r.reason))
        .collect();

    let welded = weld_clean_surface_complex_boundaries(&complex, coverage).map_err(|detail| {
        PartError::SmoothMesh {
            stage: "boundary weld",
            detail,
        }
    })?;

    let report = SmoothMeshReport {
        source_vertices: complex.source_vertices,
        source_faces: complex.source_faces,
        candidate_splines: complex.candidate_splines,
        accepted_components: complex.components.len(),
        rejections,
        welded_pairs: welded.report.welded_pairs,
        declined_pairs: welded.report.declined_pairs,
        remaining_boundary_edges: welded.report.remaining_boundary_edges,
        max_weld_gap: welded.report.max_weld_gap,
    };

    if options.require_closed && report.remaining_boundary_edges > 0 {
        return Err(PartError::SmoothMesh {
            stage: "closure",
            detail: format!(
                "{} boundary edge(s) remain unwelded after {} weld(s) (max gap {:.3e}); \
                 the cage did not close. Set require_closed=false to accept a sheet body.",
                report.remaining_boundary_edges, report.welded_pairs, report.max_weld_gap
            ),
        });
    }

    let part = Part::from_brep(welded.brep, tolerance, "smooth mesh")?;
    Ok((part, report))
}

/// Build a smooth solid from a Wavefront OBJ control cage.
///
/// Convenience over [`smooth_solid_from_control_net`] for the common
/// case where the cage came out of a mesh modeller. Only `v` and `f`
/// records are consulted; the cage is a *control* mesh, so its facets
/// are the spline's control polygons rather than the output surface.
///
/// # Errors
///
/// As [`smooth_solid_from_control_net`], plus
/// [`PartError::SmoothMesh`] with stage `"obj parse"` when the text is
/// not a control net this family can read.
pub fn smooth_solid_from_obj(
    obj_text: &str,
    tolerance: Tolerance,
    options: SmoothMeshOptions,
) -> Result<(Part, SmoothMeshReport), PartError> {
    let net = ControlNet::from_obj_str(obj_text).map_err(|source| PartError::SmoothMesh {
        stage: "obj parse",
        detail: source.to_string(),
    })?;
    smooth_solid_from_control_net(&net, tolerance, options)
}
