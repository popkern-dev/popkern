//! Describing the session's document for the page after a change: the
//! full envelope (meshes, report, parameters, assembly view), and the
//! placement-only envelope an assembly parameter that drives no part
//! earns (PK-725, design `six-axis-arm.md` §9.8): a joint slider moves
//! occurrences, so the page keeps its meshes and its report and
//! re-places the meshes. Re-sending the arm's report and its 9976
//! parameter rows cost a joint step 25 ms over its 5 ms budget.

use std::collections::HashMap;

use popkern_assembly::Document;
use popkern_core::time::monotonic_clock_available;
use popkern_rhai::{
    Outcome, document_feature_count, document_first_failure, document_report_value,
};

use crate::envelope::{AssemblyView, DefinitionView, Envelope};
use crate::mesh_buffers::MeshPurpose;
use crate::session::Session;
use crate::session_document::{
    DisplayCache, assembly_param_rows, document_params, instance_view, mesh_document,
};

/// What the page holds from the last served full describe: part meshes
/// by definition, and the envelope rows a placement-only change leaves
/// as they were. A mutable borrow of the document, a run or an opened
/// file clears it.
#[derive(Default)]
pub(crate) struct DisplayState {
    /// Display meshes of an assembly's parts (see `DisplayCache`).
    pub(crate) cache: DisplayCache,
    /// The last served assembly describe's part-dependent rows.
    placement_base: Option<PlacementBase>,
}

/// The part-dependent rows of a served assembly envelope.
struct PlacementBase {
    definitions: Vec<DefinitionView>,
    features: usize,
}

impl DisplayState {
    /// A state holding nothing.
    pub(crate) fn new() -> DisplayState {
        DisplayState {
            cache: HashMap::new(),
            placement_base: None,
        }
    }

    /// Forget every mesh and row (the document may now be another).
    pub(crate) fn clear(&mut self) {
        self.cache.clear();
        self.placement_base = None;
    }
}

impl Session {
    /// The full envelope of the document as it stands (part meshes from
    /// the display cache): for a host that reads bodies and the report
    /// after a placement-only `set_param`, which carries neither.
    pub fn current_envelope(&mut self) -> Envelope {
        self.describe(Outcome::Served, Vec::new(), None)
    }

    /// The full envelope of the document after a change ending in
    /// `outcome`.
    pub(crate) fn describe(
        &mut self,
        outcome: Outcome,
        printed: Vec<String>,
        elapsed_ms: Option<f64>,
    ) -> Envelope {
        let clock_available = monotonic_clock_available();
        let wall_budget = self.wall_budget();
        // Every describe follows a document change: step meshes keyed
        // by (body, step) may now hold other geometry.
        self.step_meshes.clear();
        self.display.placement_base = None;
        let Some(document) = self.document.as_mut() else {
            return Envelope::bare(outcome, clock_available, wall_budget);
        };
        let mut outcome = outcome;
        if let Some(decline) = document_first_failure(document).filter(|_| outcome.is_served()) {
            outcome = decline;
        }
        let (bodies, assembly, mesh_error) = mesh_document(
            document,
            MeshPurpose::Display,
            Some(&mut self.display.cache),
        );
        self.last_display_triangles = bodies.iter().map(|b| b.indices.len() / 3).sum();
        if let Some(error) = mesh_error
            && outcome.is_served()
        {
            outcome = Outcome::GateFailure {
                step: "mesh".into(),
                detail: error,
            };
        }
        let report = document_report_value(document);
        let features = document_feature_count(document);
        let params = document_params(document);
        let film_context = Some(popkern_rhai::film::FilmContext::from_document(document));
        if outcome.is_served()
            && let Some(view) = &assembly
        {
            self.display.placement_base = Some(PlacementBase {
                definitions: view.definitions.clone(),
                features,
            });
        }
        Envelope {
            outcome,
            printed,
            elapsed_ms,
            clock_available,
            wall_budget,
            report: Some(report),
            features,
            params,
            bodies,
            bodies_unchanged: false,
            assembly,
            project: None,
            film_context,
        }
    }

    /// Whether setting `name` can only move occurrences: an assembly
    /// parameter that drives no part parameter, with a served full
    /// describe on screen.
    pub(crate) fn moves_placements_only(&self, name: &str) -> bool {
        let Some(Document::Assembly(assembly)) = &self.document else {
            return false;
        };
        self.display.placement_base.is_some()
            && !name.contains(':')
            && assembly.params().iter().any(|p| p.name == name)
            && !assembly.drives().iter().any(|d| d.param == name)
    }

    /// The envelope after a placement-only change: no bodies and no
    /// report (the page keeps the ones it has, `bodies_unchanged`), the
    /// occurrences at their new placements, and only the assembly's own
    /// parameter rows (no part row changed; the arm has ten thousand). Falls back to the full describe when nothing served
    /// is on screen.
    pub(crate) fn describe_placements(&mut self) -> Envelope {
        let clock_available = monotonic_clock_available();
        let wall_budget = self.wall_budget();
        let (Some(Document::Assembly(assembly)), Some(base)) =
            (&self.document, &self.display.placement_base)
        else {
            return self.describe(Outcome::Served, Vec::new(), None);
        };
        Envelope {
            outcome: Outcome::Served,
            printed: Vec::new(),
            elapsed_ms: None,
            clock_available,
            wall_budget,
            report: None,
            features: base.features,
            params: assembly_param_rows(assembly),
            bodies: Vec::new(),
            bodies_unchanged: true,
            assembly: Some(AssemblyView {
                name: assembly.name().to_string(),
                definitions: base.definitions.clone(),
                instances: assembly.tree().iter().map(instance_view).collect(),
            }),
            project: None,
            film_context: None,
        }
    }
}
