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
//!   (and, when the family does not cover the cage as given, one
//!    Catmull–Clark refinement step first — see below)
//! ```
//!
//! The construction is the Karčiauskas–Peters polyhedral-net spline
//! family: `RegularQuad`, `ExtraordinaryStar`, `Ngon`, `Polar`, plus
//! the `T0`/`T1`/`T2` and `TwoTrianglesTwoQuads` transition
//! configurations. Every patch is a Bézier surface, so the result is
//! ordinary BRep geometry that booleans, fillets, tessellation, and
//! STEP export all consume.
//!
//! # Cage normalization
//!
//! The spline family covers specific local configurations; a general
//! cage (an n-gon face whose corners are not valence 4, a vertex
//! whose neighborhood mixes quads and n-gons) may not decompose into
//! them at all. When the direct attempt fails to produce a closed
//! body, the cage is normalized by one Catmull–Clark refinement step
//! — after which every face is a quad and every interior vertex has
//! an all-quad neighborhood, the exact serving condition of the
//! vertex-origin constructions — and the pipeline retries. The
//! refined net represents the same Catmull–Clark limit surface the
//! cage implies, so this changes the realization, not the design
//! intent. Cages the family covers natively serve directly with
//! [`SmoothMeshReport::refinement_steps`] `== 0`. See
//! `docs/architecture/smooth-mesh-cage-normalization.md`.
//!
//! # Declining
//!
//! A control net can be unsuitable — non-manifold, a valence the
//! spline family does not cover even after refinement, or a cage
//! whose components cannot be welded into a closed shell. Those are
//! reported as [`PartError::SmoothMesh`] naming the pipeline stage,
//! the per-attempt failure trail, and the unserved cage elements —
//! not as a silently-degraded body.

use popkern_core::geometry::polyhedral_spline::{
    ControlNet, ControlNetSplineOrigin, PolyhedralCoverageOptions, PolyhedralSplineSurfaceComplex,
    build_control_net_surface_complex, weld_clean_surface_complex_boundaries,
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
    /// Catmull–Clark refinement steps allowed when the direct attempt
    /// does not produce a closed body (default 1 — one step already
    /// makes every interior vertex servable, and refinement preserves
    /// valences above the family's window, so more steps only help
    /// exotic mixed declines).
    pub max_refinement_steps: usize,
    /// Screen the served body for self-intersections (`validate_all`)
    /// and decline instead of serving a self-crossing body. On by
    /// default; cage-scale bodies keep the O(n²) check cheap.
    pub check_self_intersections: bool,
}

impl Default for SmoothMeshOptions {
    fn default() -> Self {
        Self {
            samples: 8,
            weld_tolerance: 1e-6,
            require_closed: true,
            max_refinement_steps: 1,
            check_self_intersections: true,
        }
    }
}

/// What the construction produced, beyond the body itself.
#[derive(Debug, Clone)]
pub struct SmoothMeshReport {
    /// Control-net vertices consumed (of the cage actually served —
    /// the refined cage when `refinement_steps > 0`).
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
    /// Catmull–Clark steps applied before the serving attempt
    /// succeeded. Zero when the cage served as given.
    pub refinement_steps: usize,
    /// Failure summary of each attempt that preceded the serving one
    /// (empty when the first attempt served).
    pub attempts: Vec<String>,
}

impl SmoothMeshReport {
    /// Whether every spline the net implied was realised.
    pub fn complete(&self) -> bool {
        self.rejections.is_empty() && self.declined_pairs == 0
    }
}

struct AttemptFailure {
    stage: &'static str,
    detail: String,
}

/// Build a smooth solid from a control net.
///
/// # Errors
///
/// [`PartError::SmoothMesh`] when the net cannot be turned into a
/// spline complex, when welding leaves the body open while
/// [`SmoothMeshOptions::require_closed`] is set, or when the welded
/// body fails kernel validation — after exhausting the allowed
/// Catmull–Clark normalization retries. The error names the failing
/// stage, the per-attempt trail, and the unserved cage elements.
pub fn smooth_solid_from_control_net(
    net: &ControlNet,
    tolerance: Tolerance,
    options: SmoothMeshOptions,
) -> Result<(Part, SmoothMeshReport), PartError> {
    let mut cage = net.clone();
    let mut attempts: Vec<String> = Vec::new();
    for step in 0..=options.max_refinement_steps {
        match attempt(&cage, tolerance, options) {
            Ok((part, mut report)) => {
                report.refinement_steps = step;
                report.attempts = attempts;
                return Ok((part, report));
            }
            Err(failure) => {
                attempts.push(format!(
                    "step {step}: {}: {}",
                    failure.stage, failure.detail
                ));
                if step == options.max_refinement_steps {
                    return Err(PartError::SmoothMesh {
                        stage: failure.stage,
                        detail: attempt_trail(&failure.detail, &attempts),
                    });
                }
                cage = cage
                    .catmull_clark_refined()
                    .map_err(|source| PartError::SmoothMesh {
                        stage: "cage refinement",
                        detail: source.to_string(),
                    })?;
            }
        }
    }
    unreachable!("the attempt loop always returns on its final step")
}

fn attempt_trail(final_detail: &str, attempts: &[String]) -> String {
    if attempts.len() <= 1 {
        return final_detail.to_string();
    }
    format!("{final_detail}; attempt trail: [{}]", attempts.join("; "))
}

fn attempt(
    cage: &ControlNet,
    tolerance: Tolerance,
    options: SmoothMeshOptions,
) -> Result<(Part, SmoothMeshReport), AttemptFailure> {
    let coverage = PolyhedralCoverageOptions {
        tolerance,
        samples: options.samples,
        c0_tolerance: options.weld_tolerance,
        ..PolyhedralCoverageOptions::default()
    };

    let complex =
        build_control_net_surface_complex(cage, coverage).map_err(|source| AttemptFailure {
            stage: "spline construction",
            detail: source.to_string(),
        })?;

    if complex.components.is_empty() {
        return Err(AttemptFailure {
            stage: "spline construction",
            detail: format!(
                "the spline family covers none of the cage; {}",
                unserved_census(cage, &complex)
            ),
        });
    }

    let rejections: Vec<String> = complex
        .rejections
        .iter()
        .map(|r| format!("{:?} ({:?}): {}", r.origin, r.family, r.reason))
        .collect();

    let welded = weld_clean_surface_complex_boundaries(&complex, coverage).map_err(|detail| {
        AttemptFailure {
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
        refinement_steps: 0,
        attempts: Vec::new(),
    };

    if options.require_closed && report.remaining_boundary_edges > 0 {
        return Err(AttemptFailure {
            stage: "closure",
            detail: format!(
                "{} boundary edge(s) remain unwelded after {} weld(s) (max gap {:.3e}); \
                 the cage did not close ({}). Set require_closed=false to accept a sheet body.",
                report.remaining_boundary_edges,
                report.welded_pairs,
                report.max_weld_gap,
                unserved_census(cage, &complex),
            ),
        });
    }

    let part = Part::from_brep(welded.brep, tolerance, "smooth mesh").map_err(|source| {
        AttemptFailure {
            stage: "validation",
            detail: source.to_string(),
        }
    })?;

    if options.check_self_intersections {
        let findings = part.brep().validate_all();
        if !findings.is_empty() {
            let summary: Vec<String> = findings.iter().take(4).map(|f| format!("{f}")).collect();
            return Err(AttemptFailure {
                stage: "self-intersection",
                detail: format!(
                    "the smooth body fails the full validation screen: {}",
                    summary.join("; ")
                ),
            });
        }
    }

    Ok((part, report))
}

/// Name the cage elements no spline serves: vertices that anchor no
/// vertex-origin spline and belong to no face-origin spline's face,
/// with their valence and incident face arities. This is what a
/// decline needs to be actionable — it points at the exact vertex
/// whose configuration the family does not cover.
fn unserved_census(cage: &ControlNet, complex: &PolyhedralSplineSurfaceComplex) -> String {
    let mut served_vertices = vec![false; cage.vertices.len()];
    let mut served_faces = vec![false; cage.faces.len()];
    for component in &complex.components {
        match component.origin {
            ControlNetSplineOrigin::Vertex { vertex } => {
                if let Some(slot) = served_vertices.get_mut(vertex) {
                    *slot = true;
                }
            }
            ControlNetSplineOrigin::Face { face } => {
                if let Some(slot) = served_faces.get_mut(face) {
                    *slot = true;
                }
            }
        }
    }
    for (face, served) in served_faces.iter().enumerate() {
        if *served {
            for &vertex in &cage.faces[face].vertices {
                served_vertices[vertex] = true;
            }
        }
    }

    let mut incident: Vec<Vec<usize>> = vec![Vec::new(); cage.vertices.len()];
    for (face, data) in cage.faces.iter().enumerate() {
        for &vertex in &data.vertices {
            incident[vertex].push(face);
        }
    }

    let unserved: Vec<usize> = (0..cage.vertices.len())
        .filter(|&v| !served_vertices[v])
        .collect();
    if unserved.is_empty() {
        return "every cage vertex anchors a spline".to_string();
    }
    let described: Vec<String> = unserved
        .iter()
        .take(8)
        .map(|&v| {
            let arities: Vec<usize> = incident[v]
                .iter()
                .map(|&f| cage.faces[f].vertices.len())
                .collect();
            format!(
                "vertex {v} (valence {}, faces {arities:?})",
                incident[v].len()
            )
        })
        .collect();
    format!(
        "{} of {} cage vertices unserved by the spline family: {}{}",
        unserved.len(),
        cage.vertices.len(),
        described.join(", "),
        if unserved.len() > described.len() {
            ", …"
        } else {
            ""
        }
    )
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
