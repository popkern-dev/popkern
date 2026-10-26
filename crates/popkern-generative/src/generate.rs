//! The generative feature: a design space BRep and a spec in, the
//! faceted design body and its report out.
//!
//! The design space is meshed as the analysis adapter meshes any served
//! body (popkern live's display tessellation), its faces addressed by
//! id, so the setup is built here rather than parsed: supports, loads
//! and kept faces are BRep face handles the caller resolved from stable
//! names. Lengths cross into the solver in millimetres (`mm_per_unit`)
//! and the body comes back in the design space's units.

use std::ops::ControlFlow;

use popkern_analysis::model::{FaceId, Provenance};
use popkern_analysis::report::Warning;
use popkern_analysis::setup::StructuralSetup;
use popkern_analysis::structural::from_setup::solve_generative_setup_observed;
use popkern_analysis::structural::generative::body::design_model;
use popkern_analysis::structural::generative::{Design, Progress};
use popkern_analysis::{AnalysisError, Degeneracy};
use popkern_analysis_brep::{
    FeatureTolerances, NamedModel, Resolution, model_from_breps_named, resolve_faces,
};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_mesh::TessellationOptions;

use crate::facet::{FacetError, FacetOptions, FacetStats, facet};
pub use crate::request::{GenerativeRequest, LoadGroup, SupportGroup};

/// The grid's shift against the design space, in cells, for the
/// fallback solve: far beyond the ulps PK-963 needs (a face an ulp off a
/// plane), small enough that a face on a cell face stays within a
/// thousandth of a cell (2 µm at 1.25 mm cells; a shift of a third of a
/// cell left the bracket's mounting face up to half a cell off).
const GRID_SHIFT: [f64; 3] = [1.37e-3, 1.13e-3, 1.71e-3];

/// Why a generative design declined.
#[derive(Debug, thiserror::Error)]
pub enum GenerativeError {
    /// A request value is out of range.
    #[error("generative design: {0}")]
    Request(String),
    /// The design space did not become an analysis model.
    #[error("the design space does not mesh for analysis: {0}")]
    Model(String),
    /// The solver declined (its typed error: setup, budget, cancel, ...).
    #[error(transparent)]
    Solver(#[from] AnalysisError),
    /// The design surface did not make a body.
    #[error(transparent)]
    Facet(#[from] FacetError),
    /// The design does not join its loads to its supports (too few
    /// updates, or a volume fraction too small for the span).
    #[error("the design is disconnected: {0}")]
    Disconnected(#[from] crate::components::Disconnected),
    /// The design is not solved and the host deferred solving
    /// ([`crate::cache::defer_with`]): the host asks for a solve.
    #[error("this design is not solved yet: solving is deferred on this host (re-solve to run it)")]
    NeedsSolve,
}

impl GenerativeError {
    /// A short machine-readable kind (the solver's own for its declines).
    pub fn kind(&self) -> &'static str {
        match self {
            GenerativeError::Request(_) => "InvalidRequest",
            GenerativeError::Model(_) => "Model",
            GenerativeError::Solver(e) => e.kind(),
            GenerativeError::Facet(_) => "DesignBody",
            GenerativeError::NeedsSolve => "NeedsSolve",
            GenerativeError::Disconnected(_) => "DesignDisconnected",
        }
    }
}

/// What the design came to.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct GenerativeReport {
    /// Compliance of the final design, N·mm.
    pub compliance_nmm: f64,
    /// Material volume over the design space's (the density field's).
    pub volume_fraction: f64,
    /// The body's volume over the design space's.
    pub body_fraction: f64,
    /// Updates run.
    pub iterations: usize,
    /// Whether the design settled before the cap.
    pub converged: bool,
    /// Non-discreteness, 0 for solid-and-void.
    pub grey_fraction: f64,
    /// Thinnest members, mm (10% quantile).
    pub member_thickness_mm: f64,
    /// Unsupported share when printed along the build direction.
    pub unsupported_fraction: Option<f64>,
    /// Analysis cell size, mm.
    pub cell_mm: f64,
    /// The solver's warnings, rendered (thin design space, filter
    /// under-resolved, ...; the setup's own face-id references aside:
    /// this crate builds them).
    pub warnings: Vec<String>,
    /// The faceting.
    pub facet: FacetStats,
    /// Pieces of the body (each joined to a support).
    pub pieces: usize,
    /// Material dropped as joined to no support, mm³.
    pub detached_volume_mm3: f64,
    /// Enclosed voids filled, mm³.
    pub cavity_volume_mm3: f64,
    /// Cells made solid to close edge and corner contacts (one-node
    /// hinges) into face connections.
    pub closed_cells: usize,
    /// Edge or corner contacts left open (their cells outside the part).
    pub open_contacts: usize,
}

/// The design body.
#[derive(Debug, Clone)]
pub struct GenerativeResult {
    /// The faceted body, in the design space's units, unnamed.
    pub brep: BRep,
    /// Per face of `brep` (in its face order): the design space face
    /// its triangles lie on, or `None` for the design's new surface.
    pub sources: Vec<Option<FaceHandle>>,
    /// What the design came to.
    pub report: GenerativeReport,
    /// The design's id: a digest of everything its solve reads (32 hex
    /// digits), the same for the same design in any process; it names
    /// the design's history ([`crate::design_history`]).
    pub design: String,
}

/// A design space meshed for the solver and the setup built on it.
pub struct Prepared {
    /// The analysis model, its faces' sources in the design space.
    pub named: NamedModel,
    /// The setup the solver reads.
    pub setup: StructuralSetup,
    /// The setup's faces resolved on the model.
    pub resolution: Resolution,
}

/// Meshes `space` (in millimetres) and builds `request`'s setup on it.
///
/// # Errors
///
/// [`GenerativeError`]: a bad request, a design space the adapter
/// cannot mesh, a face not on the design space.
pub fn prepare(space: &BRep, request: &GenerativeRequest) -> Result<Prepared, GenerativeError> {
    request.check()?;
    let mm = request.mm_per_unit;
    let in_mm = if mm == 1.0 {
        space.clone()
    } else {
        space
            .transform(&Transform3::scale_uniform(mm))
            .map_err(|e| GenerativeError::Model(e.to_string()))?
    };
    let named = model_from_breps_named(
        std::slice::from_ref(&in_mm),
        &TessellationOptions::display(),
        &FeatureTolerances::default(),
        Provenance::Popkern {
            document: "generative design space".into(),
            revision: 0,
        },
        &|_, _| None,
    )
    .map_err(|e| GenerativeError::Model(e.to_string()))?;
    let id_of = |h: FaceHandle| {
        named
            .face_sources
            .iter()
            .position(|&(_, f)| f == h)
            .map(|i| FaceId(i as u32))
    };
    let setup = request.setup(&id_of)?;
    let resolution = resolve_faces(&named.model, &setup.faces, &setup.tolerance)?;
    Ok(Prepared {
        named,
        setup,
        resolution,
    })
}

/// The design body of a prepared design space: its faces' sources as
/// analysis face ids (`None`: the design's own surface), and the
/// surface it was faceted from (in the body's units), from which the
/// same body facets again.
pub(crate) struct Solved {
    pub brep: BRep,
    pub sources: Vec<Option<u32>>,
    pub report: GenerativeReport,
    pub surface: crate::cache::StoredSurface,
    /// The fields the solve went through (PK-959).
    pub history: Option<crate::history::DesignHistory>,
}

/// Optimises `prepared`, then facets the design's body in the design
/// space's units (`mm_per_unit`, `tolerance`).
pub(crate) fn solve(
    prepared: &Prepared,
    mm_per_unit: f64,
    tolerance: Tolerance,
    observer: &mut dyn FnMut(Progress<'_>) -> ControlFlow<()>,
) -> Result<(Solved, Design), GenerativeError> {
    let model = &prepared.named.model;
    let faces = &prepared.resolution.faces;
    let mut recorder = crate::history::Recorder::default();
    let mut recording = recorder.observing(observer);
    let observer = &mut recording;
    // A design space whose facets sit an ulp off grid planes (round
    // coordinates a boolean left a hair off) reads short in the solver's
    // clipping and declines VolumeMismatch (PK-963). Recorded fallback:
    // the same solve on a grid shifted a thousandth of a cell. The first
    // grid stays the default: a part's faces on cell faces come out
    // exact there, within the shift on the shifted one.
    let (design, shifted) =
        match solve_generative_setup_observed(model, &prepared.setup, faces, observer) {
            Err(AnalysisError::DegenerateGeometry {
                detail: Degeneracy::VolumeMismatch,
                ..
            }) => {
                let mut setup = prepared.setup.clone();
                let filter = setup.generative.as_ref().map_or(0.0, |g| g.filter_mm);
                let h = setup.options.target_cell_size.unwrap_or(filter / 3.0);
                setup.options.grid_shift = GRID_SHIFT.map(|f| f * h);
                let (design, _) = solve_generative_setup_observed(model, &setup, faces, observer)?;
                (design, Some(GRID_SHIFT))
            }
            other => (other?.0, None),
        };
    drop(recording);
    let history = recorder.finish(&design, model, mm_per_unit);
    let closing = crate::surface::body_surface(&design);
    let surface = &closing.mesh;
    let body = design_model(model, surface, design.cell_size);
    let owners: Vec<u32> = body
        .solids
        .first()
        .map_or_else(Vec::new, |s| s.triangle_face.iter().map(|f| f.0).collect());
    // Supported and loaded faces, by analysis face id.
    let ids = |names: &[String]| -> Vec<u32> {
        prepared
            .resolution
            .faces
            .iter()
            .filter(|r| names.contains(&r.name))
            .flat_map(|r| r.faces.iter().map(|f| f.0))
            .collect()
    };
    let supported: Vec<u32> = prepared
        .setup
        .supports
        .iter()
        .flat_map(|s| ids(&s.faces))
        .collect();
    let loads: Vec<Vec<u32>> = prepared.setup.loads.iter().map(|l| ids(&l.faces)).collect();
    let pieces = crate::components::sort(
        &surface.vertices,
        &surface.triangles,
        &owners,
        &|o| supported.contains(&o),
        &loads,
    )?;
    let triangles: Vec<[u32; 3]> = pieces
        .triangles
        .iter()
        .map(|&t| surface.triangles[t])
        .collect();
    let owners: Vec<u32> = pieces.triangles.iter().map(|&t| owners[t]).collect();
    let vertices: Vec<[f64; 3]> = surface
        .vertices
        .iter()
        .map(|p| p.map(|x| x / mm_per_unit))
        .collect();
    let faceted = facet(
        &vertices,
        &triangles,
        &owners,
        &FacetOptions::for_tolerance(&tolerance),
        tolerance,
    )?;
    let faces = model.faces.len() as u32;
    let sources = faceted
        .owners
        .iter()
        .map(|&o| (o < faces).then_some(o))
        .collect();
    let space_volume: f64 = design.entries.iter().map(|e| e.volume).sum();
    let mm3 = mm_per_unit * mm_per_unit * mm_per_unit;
    let report = GenerativeReport {
        compliance_nmm: design.compliance,
        volume_fraction: design.volume_fraction,
        body_fraction: if space_volume > 0.0 {
            faceted.stats.mesh_volume * mm3 / space_volume
        } else {
            0.0
        },
        iterations: design.history.len(),
        converged: design.converged,
        grey_fraction: design.grey_fraction,
        member_thickness_mm: design.member_thickness,
        unsupported_fraction: design.unsupported_fraction,
        cell_mm: design.cell_size,
        warnings: design
            .report
            .warnings
            .iter()
            .filter(|w| !matches!(w, Warning::FaceIdReference { .. }))
            .map(|w| format!("{w:?}"))
            .chain(shifted.map(|s| {
                format!(
                    "AnalysisGridShifted {{ cells: {s:?}, reason: \"the design space read short on the first grid (VolumeMismatch, PK-963)\" }}"
                )
            }))
            .collect(),
        facet: faceted.stats,
        pieces: pieces.kept,
        detached_volume_mm3: pieces.detached_volume,
        cavity_volume_mm3: pieces.cavity_volume,
        closed_cells: closing.closed_cells,
        open_contacts: closing.open_contacts,
    };
    Ok((
        Solved {
            brep: faceted.brep,
            sources,
            report,
            surface: crate::cache::StoredSurface {
                vertices,
                triangles,
                owners,
                faces,
            },
            history,
        },
        design,
    ))
}

impl Prepared {
    /// Analysis face ids back to the design space's face handles.
    pub(crate) fn handles(&self, sources: &[Option<u32>]) -> Vec<Option<FaceHandle>> {
        sources
            .iter()
            .map(|s| s.and_then(|i| self.named.face_sources.get(i as usize).map(|&(_, h)| h)))
            .collect()
    }
}

/// The design `request` asks for on `space`, each update reported to
/// `observer` (which may stop the run). Uncached: see
/// [`crate::cache::generate_cached`].
///
/// # Errors
///
/// [`GenerativeError`]: a bad request, a design space the adapter
/// cannot mesh, a solver decline (typed, `Cancelled` on a stop), or a
/// design surface that does not make a valid body.
pub fn generate(
    space: &BRep,
    request: &GenerativeRequest,
    observer: &mut dyn FnMut(Progress<'_>) -> ControlFlow<()>,
) -> Result<(GenerativeResult, Design), GenerativeError> {
    let prepared = prepare(space, request)?;
    let id = crate::cache::design_id(crate::cache::key(&prepared, request, space));
    let (solved, design) = solve(&prepared, request.mm_per_unit, space.tolerance, observer)?;
    Ok((
        GenerativeResult {
            sources: prepared.handles(&solved.sources),
            brep: solved.brep,
            report: solved.report,
            design: id,
        },
        design,
    ))
}
