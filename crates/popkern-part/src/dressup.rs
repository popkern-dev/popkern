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

    fn describe(&self) -> String {
        let radius = match self.radius {
            FilletRadius::Constant(r) => format!("radius {r}"),
            FilletRadius::Variable { start, end } => format!("radius {start}→{end}"),
            ref other => format!("{other:?}"),
        };
        format!("{radius}, {:?} requested", self.options.continuity)
    }
}

/// Outcome data from the most recent blend on a part: the kernel's
/// per-edge warnings, corner patch reports, and repair diagnostics,
/// preserved instead of dropped.
#[derive(Debug, Clone)]
pub struct BlendReport {
    /// Number of edges blended.
    pub edges_blended: usize,
    /// Kernel warnings, stringified with their edge handles.
    pub warnings: Vec<String>,
    /// Corner patch count (vertex blends inserted).
    pub corner_patches: usize,
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
        let result = fillet_edges(self.brep(), &specs, &recipe.options, self.tolerance(), None)
            .map_err(|source| PartError::Blend {
                selection: selection.clone(),
                source,
            })?;
        let report = BlendReport {
            edges_blended: edges.len(),
            warnings: result
                .edge_warnings
                .iter()
                .map(|(eh, w)| format!("{eh:?}: {w:?}"))
                .collect(),
            corner_patches: result.corner_patches.len(),
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
