//! One playground session: the model the last successful script
//! authored, the host budgets, and the verbs the page speaks.

use std::time::Duration;

use serde::Deserialize;

use popkern_core::time::{Stopwatch, monotonic_clock_available};
use popkern_rhai::outcome::{first_failure, from_run_error, from_script_error, report_value};
use popkern_rhai::{HostOptions, Outcome, ScriptHost};
use popkern_script::ContinuityCertification;
use popkern_script::Model;

use crate::envelope::{BodyMesh, Envelope, ParamRow, WallBudget};
use crate::mesh_buffers::{MeshPurpose, body_mesh};

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

/// A session over one document.
pub struct Session {
    options: HostOptions,
    model: Option<Model>,
    /// Display-mesh triangles of the last described model (the analysis
    /// mesh's budget estimate).
    last_display_triangles: usize,
}

impl Session {
    /// A session under `options`.
    pub fn new(options: HostOptions) -> Session {
        Session {
            options,
            model: None,
            last_display_triangles: 0,
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
    }

    /// The model the last successful run left (a failed run keeps the
    /// previous one so exports and parameter edits still answer).
    pub fn model(&self) -> Option<&Model> {
        self.model.as_ref()
    }

    /// Mutable access to the retained model.
    pub fn model_mut(&mut self) -> Option<&mut Model> {
        self.model.as_mut()
    }

    fn wall_budget(&self) -> WallBudget {
        WallBudget::describe(self.options.timeout.is_some(), monotonic_clock_available())
    }

    fn bare(&self, outcome: Outcome) -> Envelope {
        Envelope::bare(outcome, monotonic_clock_available(), self.wall_budget())
    }

    /// Run `source` (named `name` in diagnostics). On success the
    /// session adopts the script's model; on failure the previous
    /// model is kept and the envelope carries no bodies.
    pub fn run(&mut self, source: &str, name: &str) -> Envelope {
        let mut host = ScriptHost::with_options(self.options);
        match host.run(source, name) {
            Err(error) => self.bare(from_run_error(&error)),
            Ok(outcome) => {
                let Some(model) = outcome.model else {
                    return self.bare(Outcome::GateFailure {
                        step: "model".into(),
                        detail: "the script created no model (call model() and keep it in a \
                                 variable)"
                            .into(),
                    });
                };
                self.model = Some(model);
                self.describe(
                    Outcome::Served,
                    outcome.printed,
                    outcome.elapsed.map(|d| d.as_secs_f64() * 1000.0),
                )
            }
        }
    }

    /// Set a declared parameter and replay.
    pub fn set_param(&mut self, name: &str, value: f64) -> Envelope {
        let Some(model) = self.model.as_mut() else {
            return self.bare(no_model());
        };
        let result = model
            .set_param(name, value)
            .and_then(|_| model.replay().map(|_| ()));
        match result {
            Err(error) => {
                let outcome = from_script_error(&error, &format!("set_param {name}"));
                self.describe(outcome, Vec::new(), None)
            }
            Ok(()) => self.describe(Outcome::Served, Vec::new(), None),
        }
    }

    /// Certify deferred blend seams of the retained model for about
    /// `slice` (at least one seam; exactly one where no clock exists),
    /// upgrading its report in place. The step carries the report once
    /// nothing is left, so the page swaps it in whole.
    pub fn certify(&mut self, slice: Duration) -> CertifyStep {
        let Some(model) = self.model.as_mut() else {
            return CertifyStep::failed(no_model());
        };
        let clock = Stopwatch::start();
        let mut certified = 0;
        loop {
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
            if clock.elapsed().is_none_or(|spent| spent >= slice) {
                break;
            }
        }
        let remaining = model.deferred_seam_count();
        CertifyStep {
            certified,
            remaining,
            report: (remaining == 0).then(|| report_value(model)),
            failure: None,
        }
    }

    /// The document shown back as canonical Rhai.
    pub fn to_rhai(&mut self) -> Result<String, Outcome> {
        let Some(model) = self.model.as_mut() else {
            return Err(no_model());
        };
        popkern_rhai::to_rhai(model).map_err(|error| Outcome::GateFailure {
            step: "to_rhai".into(),
            detail: error.to_string(),
        })
    }

    /// The active body as STEP bytes with `name` as the product name.
    pub fn export_step(&mut self, name: &str) -> Result<Vec<u8>, Outcome> {
        let Some(model) = self.model.as_mut() else {
            return Err(no_model());
        };
        model
            .to_step_bytes(name)
            .map_err(|error| from_script_error(&error, "export_step"))
    }

    /// The document as `.pkn` bytes.
    pub fn export_pkn(&self) -> Result<Vec<u8>, Outcome> {
        let Some(model) = self.model.as_ref() else {
            return Err(no_model());
        };
        model
            .to_pkn_bytes()
            .map_err(|error| from_script_error(&error, "export_pkn"))
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
        let model = self.model.as_mut().ok_or_else(no_model)?;
        match mesh_bodies(model, MeshPurpose::ZebraAnalysis) {
            (bodies, None) => Ok(bodies),
            (_, Some(error)) => Err(Outcome::GateFailure {
                step: "mesh".into(),
                detail: error,
            }),
        }
    }

    fn describe(
        &mut self,
        outcome: Outcome,
        printed: Vec<String>,
        elapsed_ms: Option<f64>,
    ) -> Envelope {
        let clock_available = monotonic_clock_available();
        let wall_budget = self.wall_budget();
        let Some(model) = self.model.as_mut() else {
            return Envelope::bare(outcome, clock_available, wall_budget);
        };
        let mut outcome = outcome;
        if let Some(decline) = first_failure(model).filter(|_| outcome.is_served()) {
            outcome = decline;
        }
        let (bodies, mesh_error) = mesh_bodies(model, MeshPurpose::Display);
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
            report: Some(report_value(model)),
            features: model.features().len(),
            params: params(model),
            bodies,
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

fn params(model: &Model) -> Vec<ParamRow> {
    model
        .params()
        .into_iter()
        .map(|p| ParamRow {
            name: p.name,
            value: p.value,
            unit: serde_json::to_value(p.unit).unwrap_or(serde_json::Value::Null),
            expression: p.expression,
            bounds: p.bounds.map(|(lo, hi)| [lo, hi]),
            read_only: p.read_only,
        })
        .collect()
}

/// The buffers of `bodies` by index (a convenience for hosts that hand
/// typed arrays over one at a time).
pub fn body_buffers(bodies: &[BodyMesh], index: usize) -> Option<&BodyMesh> {
    bodies.get(index)
}

/// Mesh every live body of `model` for `purpose`: the buffers, and the
/// first meshing error (a body that fails to mesh is skipped).
fn mesh_bodies(model: &mut Model, purpose: MeshPurpose) -> (Vec<BodyMesh>, Option<String>) {
    let mut bodies = Vec::new();
    let mut error = None;
    for body in model.alive_bodies() {
        let brep = match model.brep_of(body) {
            Ok(brep) => brep,
            Err(_) => continue,
        };
        let label = model.body_label(body);
        let label_face = |handle| {
            let name = brep.faces.get(handle)?.name?;
            let entity = brep.name_pool.get(name)?;
            model.spell_entity(entity).ok()
        };
        match body_mesh(&brep, body.id.0 as u64, label, label_face, purpose) {
            Ok(mesh) => bodies.push(mesh),
            Err(e) => {
                error.get_or_insert_with(|| e.to_string());
            }
        }
    }
    (bodies, error)
}
