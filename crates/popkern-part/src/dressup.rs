//! Advanced dress-up: the fillet builder form (continuity requests,
//! caps, variable radii), blend reports, and slider-bound helpers —
//! the progressive-disclosure layer over [`Part::fillet`].

use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

use crate::error::PartError;
use crate::part::Part;
use crate::select::{EdgeQuery, edge_infos};

/// Requested blend-boundary continuity (re-exported kernel enum:
/// `G0` positional, `G1` tangent — the default, `G2` curvature).
pub use popkern_fillet::ContinuityLevel as Continuity;

/// End-cap style for open fillet ends (kernel enum).
pub use popkern_fillet::CapKind;

/// A fillet recipe: everything beyond "constant radius" that the
/// simple [`Part::fillet`] hides.
///
/// ```no_run
/// # use popkern_part::prelude::*;
/// # use popkern_part::dressup::Continuity;
/// # fn f(part: Part) -> Result<Part, PartError> {
/// part.fillet_with(
///     |f| f.radius(0.1).continuity(Continuity::G2),
///     |e| e.straight().direction(Dir::Z),
/// )
/// # }
/// ```
#[derive(Debug, Clone)]
pub struct FilletRecipe {
    radius: FilletRadius,
    options: FilletOptions,
}

impl FilletRecipe {
    pub(crate) fn new() -> FilletRecipe {
        FilletRecipe {
            radius: FilletRadius::Constant(1.0),
            options: FilletOptions::default(),
        }
    }

    /// Constant radius for every selected edge.
    pub fn radius(mut self, r: f64) -> FilletRecipe {
        self.radius = FilletRadius::Constant(r);
        self
    }

    /// Linearly varying radius from each edge's start to its end.
    pub fn variable_radius(mut self, start: f64, end: f64) -> FilletRecipe {
        self.radius = FilletRadius::Variable { start, end };
        self
    }

    /// Requested blend-boundary continuity (default `G1`; `G2`
    /// routes through the curvature-continuous constructions).
    pub fn continuity(mut self, level: Continuity) -> FilletRecipe {
        self.options.continuity = level;
        self
    }

    /// End-cap style for open fillet ends.
    pub fn cap(mut self, kind: CapKind) -> FilletRecipe {
        self.options.cap_kind = kind;
        self
    }

    /// PK-185 partial service for non-completable mid-face strip
    /// interference: instead of declining the whole operation, drop
    /// both members of every interfering pair and serve the
    /// remaining fillets. Each dropped pair is reported through a
    /// `FilletPairDropped` warning carrying the reason the pair
    /// could not be seam-completed. Default off (transactional
    /// decline).
    pub fn serve_non_interfering(mut self) -> FilletRecipe {
        self.options.interference_policy = popkern_fillet::InterferencePolicy::ServeNonInterfering;
        self
    }

    /// PK-189 conic cross-section: blend with a conic (rational
    /// quadratic) profile instead of the circular rolling-ball arc.
    /// `rho` is the shoulder fullness in (0, 1) — 0.5 is parabolic,
    /// smaller is flatter, larger is fuller; the circular value is
    /// `cos(θ/2)/(1+cos(θ/2))` for section opening angle θ, and
    /// supplying it reproduces the circular blend to roundoff.
    pub fn conic(mut self, rho: f64) -> FilletRecipe {
        self.options.profile = popkern_fillet::BlendProfile::Conic { rho };
        self
    }

    fn describe(&self) -> String {
        let radius = match self.radius {
            FilletRadius::Constant(r) => format!("radius {r}"),
            FilletRadius::Variable { start, end } => format!("radius {start}→{end}"),
            ref other => format!("{other:?}"),
        };
        format!("{radius}, {:?} requested", self.options.continuity)
    }
}

/// Outcome data from one blend operation on a part: the kernel's
/// per-edge warnings, corner patch reports, and the measured
/// continuity evidence, preserved instead of dropped.
#[derive(Debug, Clone)]
pub struct BlendReport {
    /// Number of edges blended.
    pub edges_blended: usize,
    /// Kernel warnings, stringified with their edge handles.
    pub warnings: Vec<String>,
    /// Corner patch count (vertex blends inserted).
    pub corner_patches: usize,
    /// The continuity level the operation was asked for. Chamfers
    /// record `G0` — their boundary is a designed crease.
    pub requested_continuity: Continuity,
    /// Measured continuity evidence over the operation's live seams
    /// (PK-188): per-seam G0/G1/G2 residuals and the achieved level
    /// classified against explicit gates, so a silent downgrade is
    /// visible next to `requested_continuity`. `None` for chamfers
    /// (nothing smooth was requested) — never `None` for a fillet.
    ///
    /// The evidence is a *record of the measurement at blend time*:
    /// its handles are provenance, not live references, and stay
    /// meaningful in the report even after later operations rebuild
    /// the body.
    pub continuity: Option<popkern_fillet::BlendContinuityEvidence>,
}

/// Seam sampling density for blend-time continuity evidence. The
/// apex gates use 100 on a single corner; 33 keeps whole-part
/// measurement cheap while still straddling every patch span.
pub(crate) const EVIDENCE_SAMPLES_PER_SEAM: usize = 33;

/// Measure the continuity evidence for the faces `post` gained over
/// `pre_faces` (the blend delta), at the commercial gates.
pub(crate) fn measure_blend_evidence(
    post: &popkern_core::topology::BRep,
    pre_faces: &std::collections::BTreeSet<popkern_core::topology::FaceHandle>,
    requested: Continuity,
) -> popkern_fillet::BlendContinuityEvidence {
    let new_faces: Vec<_> = post
        .faces
        .keys()
        .filter(|f| !pre_faces.contains(f))
        .collect();
    popkern_fillet::measure_blend_continuity_evidence(
        post,
        &new_faces,
        requested,
        &popkern_fillet::ContinuityGates::commercial(),
        EVIDENCE_SAMPLES_PER_SEAM,
    )
}

impl Part {
    /// Fillet with full control: a recipe (radius profile,
    /// continuity, caps) applied to every selected edge in one
    /// kernel operation.
    ///
    /// # Errors
    ///
    /// As [`Part::fillet`]: empty selection unless permitted, kernel
    /// blend diagnostics on failure.
    pub fn fillet_with(
        &self,
        recipe: impl FnOnce(FilletRecipe) -> FilletRecipe,
        select: impl FnOnce(EdgeQuery) -> EdgeQuery,
    ) -> Result<Part, PartError> {
        let recipe = recipe(FilletRecipe::new());
        let query = select(EdgeQuery::new());
        let infos = edge_infos(self.brep());
        let mut edges = query.resolve(&infos);
        let selection = query.describe();
        let seams: std::collections::HashSet<_> = infos
            .iter()
            .filter(|i| i.is_seam)
            .map(|i| i.handle)
            .collect();
        edges.retain(|e| !seams.contains(e));
        if edges.is_empty() {
            if query.empty_permitted() {
                return Ok(self.clone());
            }
            return Err(PartError::EmptySelection {
                selection,
                entity: "edges",
                available: infos.len(),
            });
        }
        let specs: Vec<FilletSpec> = edges
            .iter()
            .map(|&edge| FilletSpec {
                edge,
                radius: recipe.radius.clone(),
            })
            .collect();
        let pre_faces: std::collections::BTreeSet<_> = self.brep().faces.keys().collect();
        let result = fillet_edges(self.brep(), &specs, &recipe.options, self.tolerance(), None)
            .map_err(|source| PartError::Blend {
                selection: selection.clone(),
                source,
            })?;
        let requested = recipe.options.continuity;
        let evidence = measure_blend_evidence(&result.brep, &pre_faces, requested);
        let report = BlendReport {
            edges_blended: edges.len(),
            warnings: result
                .edge_warnings
                .iter()
                .map(|(eh, w)| format!("{eh:?}: {w:?}"))
                .collect(),
            corner_patches: result.corner_patches.len(),
            requested_continuity: requested,
            continuity: Some(evidence),
        };
        self.rebuild_after_blend(
            result.brep,
            format!(
                "{} on {} edges: {selection}",
                recipe.describe(),
                edges.len()
            ),
            report,
        )
    }

    /// A conservative upper bound for the fillet radius over the
    /// selected edges: half the shortest selected edge's chord (the
    /// slider bound per the UI-forward rules; the kernel may accept
    /// slightly more or decline slightly less depending on the
    /// surrounding geometry).
    pub fn max_fillet_radius(&self, select: impl FnOnce(EdgeQuery) -> EdgeQuery) -> Option<f64> {
        let query = select(EdgeQuery::new());
        let infos = edge_infos(self.brep());
        let edges = query.resolve(&infos);
        let min_chord = infos
            .iter()
            .filter(|i| edges.contains(&i.handle) && i.chord_length > 0.0)
            .map(|i| i.chord_length)
            .fold(f64::INFINITY, f64::min);
        min_chord.is_finite().then_some(min_chord / 2.0)
    }
}
