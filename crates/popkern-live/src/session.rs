//! One popkern live session: the model the last successful script
//! authored, the host budgets, and the verbs the page speaks.

use std::time::Duration;

use serde::Deserialize;

use popkern_assembly::{Document, InstancePath};
use popkern_core::time::{Stopwatch, monotonic_clock_available};
use popkern_rhai::outcome::from_run_error;
use popkern_rhai::{
    HostOptions, Outcome, ScriptHost, document_feature_count, document_first_failure,
    document_report_value, from_assembly_error,
};
use popkern_script::ContinuityCertification;
use popkern_script::Model;

use crate::envelope::{BodyMesh, Envelope, WallBudget};
use crate::mesh_buffers::MeshPurpose;
use crate::session_document::{document_params, mesh_document};

/// The budgets as the page spells them (JSON).
#[derive(Debug, Clone, Copy, Default, Deserialize)]
pub struct RunOptions {
    /// Maximum Rhai operations (`null` = unlimited).
    #[serde(default)]
    pub max_operations: Option<u64>,
    /// Wall-clock budget in milliseconds (`null` = none).
    #[serde(default)]
    pub timeout_ms: Option<f64>,
    /// `"Inline"` (default) or `"Deferred"`: whether fillets certify
    /// their seams during the run or later through [`Session::certify`].
    #[serde(default)]
    pub continuity_certification: ContinuityCertification,
}

impl RunOptions {
    /// The host options these budgets ask for.
    pub fn host(&self) -> HostOptions {
        HostOptions {
            max_operations: self.max_operations,
            timeout: self
                .timeout_ms
                .filter(|ms| ms.is_finite() && *ms >= 0.0)
                .map(|ms| Duration::from_secs_f64(ms / 1000.0)),
            continuity_certification: self.continuity_certification,
        }
    }
}

/// A session over one document (a part or an assembly).
pub struct Session {
    pub(crate) options: HostOptions,
    pub(crate) document: Option<Document>,
    /// Display-mesh triangles of the last described model (the analysis
    /// mesh's budget estimate).
    last_display_triangles: usize,
    /// The project runner and store, once a project has been built.
    pub(crate) project: Option<crate::session_project::ProjectState>,
}

impl Session {
    /// A session under `options`.
    pub fn new(options: HostOptions) -> Session {
        Session {
            options,
            document: None,
            last_display_triangles: 0,
            project: None,
        }
    }

    /// The budgets in force.
    pub fn options(&self) -> HostOptions {
        self.options
    }

    /// Replace the budgets for every later run (the page's time-limit
    /// control). The retained model is kept.
    pub fn set_options(&mut self, options: HostOptions) {
        self.options = options;
        if let Some(project) = self.project.as_mut() {
            project.set_options(options);
        }
    }

    /// The model the last successful run left, when it authored a part
    /// document (a failed run keeps the previous document so exports
    /// and parameter edits still answer).
    pub fn model(&self) -> Option<&Model> {
        match &self.document {
            Some(Document::Part(model)) => Some(model),
            _ => None,
        }
    }

    /// Mutable access to the retained part document's model.
    pub fn model_mut(&mut self) -> Option<&mut Model> {
        match &mut self.document {
            Some(Document::Part(model)) => Some(model),
            _ => None,
        }
    }

    /// The retained document of either kind.
    pub fn document(&self) -> Option<&Document> {
        self.document.as_ref()
    }

    /// Mutable access to the retained document.
    pub fn document_mut(&mut self) -> Option<&mut Document> {
        self.document.as_mut()
    }

    fn wall_budget(&self) -> WallBudget {
        WallBudget::describe(self.options.timeout.is_some(), monotonic_clock_available())
    }

    pub(crate) fn bare(&self, outcome: Outcome) -> Envelope {
        Envelope::bare(outcome, monotonic_clock_available(), self.wall_budget())
    }

    /// Run `source` (named `name` in diagnostics). On success the
    /// session adopts the script's model; on failure the previous
    /// model is kept and the envelope carries no bodies.
    pub fn run(&mut self, source: &str, name: &str) -> Envelope {
        let mut host = ScriptHost::with_options(self.options);
        match host.run(source, name) {
            Err(error) => self.bare(from_run_error(&error)),
            Ok(mut outcome) => {
                let Some(document) = outcome.take_document() else {
                    return self.bare(Outcome::GateFailure {
                        step: "model".into(),
                        detail: "the script created no model (call model() or assembly() and \
                                 keep it in a variable)"
                            .into(),
                    });
                };
                self.document = Some(document);
                self.describe(
                    Outcome::Served,
                    outcome.printed,
                    outcome.elapsed.map(|d| d.as_secs_f64() * 1000.0),
                )
            }
        }
    }

    /// Open a STEP file's product structure as the session's document
    /// (PK-631, design §8.3): a view-only assembly of imported parts —
    /// placed, named, re-exportable and savable, with no recipe. The
    /// envelope's `printed` carries the import report's summary; on
    /// failure the previous document is kept.
    pub fn open_step(&mut self, bytes: &[u8], name: &str) -> Envelope {
        match popkern_assembly::Assembly::from_step_bytes(bytes, name) {
            Err(error) => self.bare(from_assembly_error(&error, "open_step")),
            Ok((assembly, report)) => {
                let mut printed = vec![format!(
                    "imported {name}: {} part(s), {} occurrence(s), {} solid(s)",
                    report.parts, report.instances, report.solids
                )];
                printed.extend(
                    report
                        .refused_solids
                        .iter()
                        .map(|r| format!("refused solid: {r}")),
                );
                printed.extend(
                    report
                        .refused_instances
                        .iter()
                        .map(|r| format!("left out: {r}")),
                );
                self.document = Some(Document::Assembly(Box::new(assembly)));
                self.describe(Outcome::Served, printed, None)
            }
        }
    }

    /// Set a declared parameter and replay.
    /// An assembly takes `part:name` and replays only that part (every
    /// instance of it updates).
    pub fn set_param(&mut self, name: &str, value: f64) -> Envelope {
        let Some(document) = self.document.as_mut() else {
            return self.bare(no_model());
        };
        // A part parameter replays its part; an assembly parameter
        // replays the parts it drives inside `set_param` already.
        let part = name.split_once(':').map(|(part, _)| part.to_string());
        let assembly_param = part.is_none() && matches!(document, Document::Assembly(_));
        let result = document.set_param(name, value).and_then(|_| {
            for (definition, model) in document.models_mut() {
                if assembly_param {
                    break;
                }
                if part.is_none() || definition == part {
                    model
                        .replay()
                        .map_err(|source| popkern_assembly::AssemblyError::Script {
                            definition: definition.unwrap_or_else(|| "<document>".into()),
                            source,
                        })?;
                }
            }
            Ok(())
        });
        match result {
            Err(error) => {
                let outcome = from_assembly_error(&error, &format!("set_param {name}"));
                self.describe(outcome, Vec::new(), None)
            }
            Ok(()) => self.describe(Outcome::Served, Vec::new(), None),
        }
    }

    /// Show or hide the occurrence at `path` (`stage/planet_1`) in an
    /// assembly. Visibility is document data (it is saved), but it
    /// changes no geometry, so nothing is re-meshed.
    pub fn set_visible(&mut self, path: &str, visible: bool) -> Result<(), Outcome> {
        let Some(Document::Assembly(assembly)) = self.document.as_mut() else {
            return Err(Outcome::GateFailure {
                step: "set_visible".into(),
                detail: "visibility applies to assembly instances; this document is one part"
                    .into(),
            });
        };
        assembly
            .resolve(&InstancePath::parse(path))
            .and_then(|at| assembly.set_visible(at, visible))
            .map_err(|error| from_assembly_error(&error, "set_visible"))
    }

    /// Certify deferred blend seams of the retained model for about
    /// `slice` (at least one seam; exactly one where no clock exists),
    /// upgrading its report in place. The step carries the report once
    /// nothing is left, so the page swaps it in whole.
    pub fn certify(&mut self, slice: Duration) -> CertifyStep {
        let Some(document) = self.document.as_mut() else {
            return CertifyStep::failed(no_model());
        };
        let clock = Stopwatch::start();
        let mut certified = 0;
        // Parts in definition order: each model's seams until it has
        // none left, within one slice.
        'parts: for (_, model) in document.models_mut() {
            loop {
                if certified > 0 && clock.elapsed().is_none_or(|spent| spent >= slice) {
                    break 'parts;
                }
                match model.certify_next_deferred_seam() {
                    Ok(Some(_)) => certified += 1,
                    Ok(None) => break,
                    Err(error) => {
                        return CertifyStep::failed(Outcome::GateFailure {
                            step: "certify".into(),
                            detail: error.to_string(),
                        });
                    }
                }
            }
        }
        let remaining = document
            .models_mut()
            .into_iter()
            .map(|(_, m)| m.deferred_seam_count())
            .sum();
        CertifyStep {
            certified,
            remaining,
            report: (remaining == 0).then(|| document_report_value(document)),
            failure: None,
        }
    }

    /// The document shown back as canonical Rhai.
    pub fn to_rhai(&mut self) -> Result<String, Outcome> {
        let Some(document) = self.document.as_mut() else {
            return Err(no_model());
        };
        popkern_rhai::document_to_rhai(document).map_err(|error| Outcome::GateFailure {
            step: "to_rhai".into(),
            detail: error.to_string(),
        })
    }

    /// The whole model as STEP bytes (every alive body its own solid)
    /// with `name` as the product name.
    /// An assembly exports its product structure.
    pub fn export_step(&mut self, name: &str) -> Result<Vec<u8>, Outcome> {
        let Some(document) = self.document.as_mut() else {
            return Err(no_model());
        };
        document
            .to_step_bytes(name)
            .map_err(|error| from_assembly_error(&error, "export_step"))
    }

    /// Interference and clearance between the assembly's part
    /// occurrences within `clearance` (design §8.1): the report, or a
    /// gate failure for a part document.
    pub fn interference(
        &mut self,
        clearance: f64,
    ) -> Result<popkern_assembly::InterferenceReport, Outcome> {
        let Some(Document::Assembly(assembly)) = self.document.as_mut() else {
            return Err(Outcome::GateFailure {
                step: "interference".into(),
                detail:
                    "interference is checked between assembly instances; this document is one part"
                        .into(),
            });
        };
        let options = popkern_assembly::InterferenceOptions {
            clearance: if clearance.is_finite() && clearance >= 0.0 {
                clearance
            } else {
                0.0
            },
            ..Default::default()
        };
        assembly
            .interference(&options)
            .map_err(|error| from_assembly_error(&error, "interference"))
    }

    /// The document as `.pkn` bytes.
    pub fn export_pkn(&self) -> Result<Vec<u8>, Outcome> {
        let Some(document) = self.document.as_ref() else {
            return Err(no_model());
        };
        document
            .to_pkn_bytes()
            .map_err(|error| from_assembly_error(&error, "export_pkn"))
    }

    /// The envelope for the retained model: report, parameters, and a
    /// mesh per alive body. A served outcome becomes a kernel decline
    /// when the last replay recorded a failed feature, and a gate
    /// failure when a body will not tessellate.
    /// The live bodies re-meshed at the zebra ANALYSIS bound
    /// ([`crate::mesh_buffers::analysis_options`], PK-577), in the same
    /// order as the last envelope's bodies.
    ///
    /// Fallback boundary: the analysis mesh has about four times the
    /// display mesh's triangles (halving the facet angle quadruples the
    /// density on curved faces). When that estimate exceeds
    /// [`crate::mesh_buffers::ANALYSIS_TRIANGLE_BUDGET`], the request
    /// declines typed and the zebra view keeps the display mesh. Mesh
    /// time grows about linearly with triangle count (PK-589), so the
    /// budget bounds both the buffers and the wait.
    pub fn analysis_meshes(&mut self) -> Result<Vec<BodyMesh>, Outcome> {
        let display_triangles: usize = self.last_display_triangles;
        let estimate = 4 * display_triangles;
        if estimate > crate::mesh_buffers::ANALYSIS_TRIANGLE_BUDGET {
            return Err(Outcome::GateFailure {
                step: "analysis_mesh".into(),
                detail: format!(
                    "the display mesh has {display_triangles} triangles; the 6° analysis mesh would be about {estimate}, over the {} budget: the zebra view keeps the display mesh",
                    crate::mesh_buffers::ANALYSIS_TRIANGLE_BUDGET
                ),
            });
        }
        let document = self.document.as_mut().ok_or_else(no_model)?;
        match mesh_document(document, MeshPurpose::ZebraAnalysis) {
            (bodies, _, None) => Ok(bodies),
            (_, _, Some(error)) => Err(Outcome::GateFailure {
                step: "mesh".into(),
                detail: error,
            }),
        }
    }

    pub(crate) fn describe(
        &mut self,
        outcome: Outcome,
        printed: Vec<String>,
        elapsed_ms: Option<f64>,
    ) -> Envelope {
        let clock_available = monotonic_clock_available();
        let wall_budget = self.wall_budget();
        let Some(document) = self.document.as_mut() else {
            return Envelope::bare(outcome, clock_available, wall_budget);
        };
        let mut outcome = outcome;
        if let Some(decline) = document_first_failure(document).filter(|_| outcome.is_served()) {
            outcome = decline;
        }
        let (bodies, assembly, mesh_error) = mesh_document(document, MeshPurpose::Display);
        self.last_display_triangles = bodies.iter().map(|b| b.indices.len() / 3).sum();
        if let Some(error) = mesh_error
            && outcome.is_served()
        {
            outcome = Outcome::GateFailure {
                step: "mesh".into(),
                detail: error,
            };
        }
        Envelope {
            outcome,
            printed,
            elapsed_ms,
            clock_available,
            wall_budget,
            report: Some(document_report_value(document)),
            features: document_feature_count(document),
            params: document_params(document),
            bodies,
            assembly,
            project: None,
        }
    }
}

/// One [`Session::certify`] slice.
#[derive(Debug, Clone, serde::Serialize)]
pub struct CertifyStep {
    /// Seams certified in this slice.
    pub certified: usize,
    /// Deferred seams still left.
    pub remaining: usize,
    /// The upgraded report, once `remaining` is zero.
    pub report: Option<serde_json::Value>,
    /// Why certification stopped, when it could not continue.
    pub failure: Option<Outcome>,
}

impl CertifyStep {
    fn failed(outcome: Outcome) -> CertifyStep {
        CertifyStep {
            certified: 0,
            remaining: 0,
            report: None,
            failure: Some(outcome),
        }
    }
}

fn no_model() -> Outcome {
    Outcome::GateFailure {
        step: "model".into(),
        detail: "no model in this session: run a script first".into(),
    }
}

/// The buffers of `bodies` by index (a convenience for hosts that hand
/// typed arrays over one at a time).
pub fn body_buffers(bodies: &[BodyMesh], index: usize) -> Option<&BodyMesh> {
    bodies.get(index)
}
