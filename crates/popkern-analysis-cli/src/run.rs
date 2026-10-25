//! `popkern-analysis run part.step setup.toml`: one solve, a versioned
//! JSON report and a `.vtu` of the per-vertex fields.
//!
//! The report holds no timestamps or timings, so two runs of the same
//! inputs give identical bytes. A decline is a report too: `status` is
//! `"declined"` and `error` carries the typed kind and message.

use std::path::Path;

use popkern_analysis::assess::fillet::FilletStatus;
use popkern_analysis::assess::{AssessOptions, assess, label_by_growth, local_peak_growth};
use popkern_analysis::model::{AnalysisModel, FaceId};
use popkern_analysis::setup::SupportKind;
use popkern_analysis::setup::TemperatureSource;
use popkern_analysis::structural::from_setup::{problem_from_setup, solve_structural_with};
use popkern_analysis::structural::submodel::refine_raiser;
use popkern_analysis::thermal::from_setup::solve_thermal;
use popkern_analysis::{AnalysisError, SolveReport};
use popkern_analysis_brep::resolve_faces;
use serde_json::{Value, json};

use crate::png::write_field_png;
use crate::setup_file::parse_setup;
use crate::step_model::{model_from_step_file, sha256_hex};
use crate::thermal_file::parse_thermal_setup;
use crate::vtu::{VtuFields, write_vtu};

#[path = "run_generative.rs"]
mod generative;
#[path = "run_modal.rs"]
mod modal;
#[path = "run_thermal.rs"]
mod thermal;
use thermal::run_thermal;

/// Report schema identifier; bump on any incompatible change and keep
/// `docs/schema/report-v1.schema.json` in step.
pub const REPORT_SCHEMA: &str = "popkern-analysis-report/1";

/// The popkern commit this build came from: `POPKERN_COMMIT` at build
/// time when the release tooling sets it, else `"in-tree"` (the solver
/// builds in the popkern workspace since PK-948; the pin is retired).
pub const POPKERN_COMMIT: &str = match option_env!("POPKERN_COMMIT") {
    Some(commit) => commit,
    None => "in-tree",
};

/// How a run is carried out.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct RunOptions {
    /// Solve again at half the cell size around each solid's peak (a
    /// submodel, PH-047) and label peaks that grow as singularities
    /// (PH-038). On the C8 shoulder the submodel has 0.34x the first
    /// solve's dofs, where the whole part at h/2 had 6.8x.
    pub refinement_check: bool,
    /// Also render the main field (von Mises, or temperature) to this
    /// PNG path (PH-041).
    pub png: Option<std::path::PathBuf>,
    /// Also write the structural surface samples (every clipped piece's
    /// centroid and corners, with stress and displacement) to this CSV
    /// path, for point-by-point comparison with another solver.
    pub samples: Option<std::path::PathBuf>,
    /// Re-solve each fillet or bore whose peak is not resolved in a
    /// submodel around it (PH-063).
    pub refine: bool,
    /// A generative run's design surface (the iso-surface of its
    /// material at ½, closed) as binary STL at this path.
    pub design_mesh: Option<std::path::PathBuf>,
}

/// What one run produced.
#[derive(Debug, Clone, PartialEq)]
pub struct RunOutput {
    /// The JSON report.
    pub report: Value,
    /// The `.vtu` file, when the solve succeeded.
    pub vtu: Option<String>,
}

impl RunOutput {
    /// Whether the solve succeeded.
    pub fn solved(&self) -> bool {
        self.report["status"] == "solved"
    }
}

/// A typed failure: kind and message.
struct Failure {
    kind: &'static str,
    message: String,
}

impl Failure {
    fn new(kind: &'static str, e: impl std::fmt::Display) -> Failure {
        Failure {
            kind,
            message: e.to_string(),
        }
    }
}

impl From<AnalysisError> for Failure {
    fn from(e: AnalysisError) -> Failure {
        Failure::new(analysis_error_kind(&e), e)
    }
}

/// The variant name of an [`AnalysisError`], for the report.
pub fn analysis_error_kind(e: &AnalysisError) -> &'static str {
    e.kind()
}

/// The material's cited (or typical) linear expansion, 1/K.
fn material_expansion(m: &popkern_analysis::material::Material) -> Option<f64> {
    match m.thermal_expansion {
        popkern_analysis::material::Property::Cited { value, .. }
        | popkern_analysis::material::Property::Typical { value, .. } => Some(value),
        _ => None,
    }
}

fn to_value(v: &impl serde::Serialize) -> Result<Value, Failure> {
    serde_json::to_value(v).map_err(|e| Failure::new("Serialize", e))
}

fn model_summary(m: &AnalysisModel) -> Value {
    json!({
        "solids": m.solids.len(),
        "faces": m.faces.len(),
        "concave_edges": m.concave_edges.len(),
        "inferred_fillets": m.faces.iter().filter(|f| f.role.is_some()).count(),
    })
}

/// Runs the part at `step` with the setup at `setup`.
pub fn run(step: &Path, setup: &Path) -> RunOutput {
    run_with(step, setup, &RunOptions::default())
}

/// [`run`] with explicit run options.
pub fn run_with(step: &Path, setup: &Path, options: &RunOptions) -> RunOutput {
    let mut report = json!({
        "schema": REPORT_SCHEMA,
        "popkern_commit": POPKERN_COMMIT,
    });
    match run_inner(step, setup, options, &mut report) {
        Ok(vtu) => {
            report["status"] = json!("solved");
            RunOutput {
                report,
                vtu: Some(vtu),
            }
        }
        Err(f) => {
            report["status"] = json!("declined");
            report["error"] = json!({ "kind": f.kind, "message": f.message });
            RunOutput { report, vtu: None }
        }
    }
}

fn run_inner(
    step: &Path,
    setup: &Path,
    options: &RunOptions,
    report: &mut Value,
) -> Result<String, Failure> {
    let text = std::fs::read_to_string(setup).map_err(|e| Failure::new("SetupRead", e))?;
    let step_bytes = std::fs::read(step).map_err(|e| Failure::new("StepRead", e))?;
    report["input"] = json!({
        "step_sha256": sha256_hex(&step_bytes),
        "setup_sha256": sha256_hex(text.as_bytes()),
    });
    let kind: Value = toml::from_str(&text).map_err(|e| Failure::new("SetupSyntax", e))?;
    let thermal = kind.get("analysis").and_then(Value::as_str) == Some("thermal");
    report["analysis"] = json!(if thermal { "thermal" } else { "structural" });
    let model = model_from_step_file(step).map_err(|e| Failure::new("StepModel", e))?;
    report["model"] = model_summary(&model);
    let (vtu, field, name) = if thermal {
        if options.samples.is_some() {
            return Err(Failure::new(
                "Samples",
                "surface samples are written for structural runs only",
            ));
        }
        let (vtu, field) = run_thermal(&model, &text, options.refine, report)?;
        (vtu, field, "temperature_c")
    } else {
        let dir = setup.parent().unwrap_or_else(|| Path::new("."));
        run_structural(&model, &text, dir, options, report)?
    };
    if let Some(path) = &options.png {
        let (lo, hi) = write_field_png(&model, &field, path).map_err(|e| Failure::new("Png", e))?;
        report["result"]["png"] = json!({ "field": name, "min": lo, "max": hi });
    }
    Ok(vtu)
}

fn with_resolution(
    mut solve: SolveReport,
    warnings: &[popkern_analysis::report::Warning],
) -> SolveReport {
    let mut all = warnings.to_vec();
    all.append(&mut solve.warnings);
    solve.warnings = all;
    solve
}

/// Adds `face_names` (the setup names whose faces include the peak's
/// face) to every peak, so a report reads "peak on `hook`".
fn name_peaks(assessment: &mut Value, resolved: &[popkern_analysis::report::ResolvedFace]) {
    let Some(peaks) = assessment.get_mut("peaks").and_then(Value::as_array_mut) else {
        return;
    };
    for solid in peaks {
        for key in ["von_mises", "max_principal"] {
            let Some(peak) = solid.get_mut(key) else {
                continue;
            };
            let face = peak.get("face").and_then(Value::as_u64);
            let names: Vec<&str> = resolved
                .iter()
                .filter(|r| face.is_some_and(|f| r.faces.iter().any(|x| u64::from(x.0) == f)))
                .map(|r| r.name.as_str())
                .collect();
            peak["face_names"] = json!(names);
        }
    }
}

fn run_structural(
    model: &AnalysisModel,
    text: &str,
    setup_dir: &Path,
    run_options: &RunOptions,
    report: &mut Value,
) -> Result<(String, Vec<Vec<f64>>, &'static str), Failure> {
    let setup = parse_setup(text).map_err(|e| Failure::new("SetupFile", e))?;
    let resolution = resolve_faces(model, &setup.faces, &setup.tolerance)?;
    if setup.modal.is_some() {
        let (vtu, field) = modal::run_modal(model, &setup, &resolution, report)?;
        return Ok((vtu, field, "mode_shape"));
    }
    if setup.generative.is_some() {
        let (vtu, field) = generative::run_generative(
            model,
            &setup,
            &resolution,
            run_options.design_mesh.as_deref(),
            report,
        )?;
        return Ok((vtu, field, "density"));
    }
    // A thermoelastic setup naming a thermal setup: that solve first, on
    // the same part (a staggered solve).
    let temperature = match setup.thermal.as_ref().map(|t| &t.temperature) {
        Some(TemperatureSource::ThermalSetup(path)) => {
            let path = setup_dir.join(path);
            let text = std::fs::read_to_string(&path).map_err(|e| Failure::new("SetupRead", e))?;
            let thermal = parse_thermal_setup(&text).map_err(|e| Failure::new("ThermalFile", e))?;
            let faces = resolve_faces(model, &thermal.faces, &thermal.tolerance)?;
            let solved = solve_thermal(model, &thermal, &faces.faces)?;
            report["thermal_solve"] = json!({
                "setup": path.display().to_string(),
                "setup_sha256": sha256_hex(text.as_bytes()),
                "max_temperature_c": solved.max_temperature_c,
            });
            Some(solved.field)
        }
        _ => None,
    };
    let (result, material) =
        solve_structural_with(model, &setup, &resolution.faces, temperature.as_ref())?;
    let assess_options = AssessOptions {
        build_direction: setup.build_direction,
        ..AssessOptions::default()
    };
    let mut assessment =
        assess(model, &result, material, &assess_options).map_err(|e| Failure::new("Assess", e))?;
    let refinement = if run_options.refinement_check {
        // The same problem at half the cell size in a submodel around
        // each solid's peak (PH-047), driven by this solution; a decline
        // there is reported, not fatal.
        let h = result.report.cell_size.unwrap_or(0.0);
        // Clamps make edge singularities; rollers mostly sit on symmetry
        // planes and do not.
        let clamped: Vec<FaceId> = setup
            .supports
            .iter()
            .filter(|s| s.kind == SupportKind::Fixed)
            .flat_map(|s| &s.faces)
            .flat_map(|name| {
                resolution
                    .faces
                    .iter()
                    .filter(move |r| &r.name == name)
                    .flat_map(|r| r.faces.iter().copied())
            })
            .collect();
        let (problem, _) = problem_from_setup(model, &setup, &resolution.faces)?;
        match local_peak_growth(
            model,
            &problem,
            &result,
            &clamped,
            &setup.options,
            &assess_options,
        ) {
            Ok((growth, dof)) => {
                label_by_growth(&mut assessment.peaks, &growth, &assess_options);
                json!({
                    "status": "solved",
                    "cell_size": 0.5 * h,
                    "dof": dof,
                    "growth": to_value(&growth)?,
                })
            }
            Err(e) => json!({
                "status": "declined",
                "cell_size": 0.5 * h,
                "error": { "kind": analysis_error_kind(&e), "message": e.to_string() },
            }),
        }
    } else {
        Value::Null
    };
    let reactions: Vec<Value> = resolution
        .faces
        .iter()
        .map(|r| {
            let force: [f64; 3] = std::array::from_fn(|k| {
                r.faces
                    .iter()
                    .filter_map(|f| result.face_reaction.get(f.0 as usize))
                    .map(|v| v[k])
                    .sum()
            });
            json!({ "name": r.name, "force_n": force })
        })
        .collect();
    let submodels = if run_options.refine {
        let (problem, _) = problem_from_setup(model, &setup, &resolution.faces)?;
        let raisers: Vec<(FaceId, f64)> = assessment
            .fillets
            .iter()
            .chain(&assessment.bores)
            .filter(|c| !matches!(c.status, FilletStatus::Resolved))
            .map(|c| (c.face, c.radius))
            .collect();
        let entries: Vec<Value> = raisers
            .iter()
            .map(|&(face, radius)| {
                match refine_raiser(model, &problem, &result, face, radius, &setup.options) {
                    Ok(r) => to_value(&r),
                    Err(e) => Ok(json!({
                        "face": face,
                        "radius": radius,
                        "error": { "kind": analysis_error_kind(&e), "message": e.to_string() },
                    })),
                }
            })
            .collect::<Result<_, _>>()?;
        Value::Array(entries)
    } else {
        Value::Null
    };
    let samples = match &run_options.samples {
        Some(path) => {
            let (csv, count) =
                crate::samples::samples_csv(&result).map_err(|e| Failure::new("Samples", e))?;
            std::fs::write(path, csv).map_err(|e| Failure::new("Samples", e))?;
            json!({ "count": count })
        }
        None => Value::Null,
    };
    let solve = with_resolution(result.report.clone(), &resolution.warnings);
    report["solve"] = to_value(&solve)?;
    let direction_used = assessment.safety.build_direction.is_some();
    let mut assessment = to_value(&assessment)?;
    name_peaks(&mut assessment, &resolution.faces);
    report["result"] = json!({
        "material": material.name,
        "max_displacement_mm": result.max_displacement,
        "global_error": result.global_error,
        "face_reactions": reactions,
        "assessment": assessment,
        // A layered material checks its layers in build axes when a
        // direction is given; otherwise at the worst-case orientation.
        "build_direction_used": direction_used,
        "infill": setup.infill.as_ref().map(|i| {
            let (e, nu) = popkern_analysis::structural::infill::gyroid(
                i.density,
                material.youngs_modulus.value().unwrap_or(0.0),
                material.poisson_ratio.value().unwrap_or(0.0),
            );
            json!({
                "pattern": i.pattern,
                "density": i.density,
                "walls_mm": i.walls_mm,
                "line_width_mm": i.line_width_mm,
                "gyroid_cell_mm": popkern_analysis::structural::infill::Infill {
                    density: i.density,
                    walls_mm: i.walls_mm,
                    line_width_mm: i.line_width_mm,
                }
                .cell_mm(),
                "core_youngs_mpa": e,
                "core_poisson": nu,
            })
        }),
        "thermal": setup.thermal.as_ref().map(|t| json!({
            "reference_c": t.reference_c,
            "alpha_per_k": material_expansion(material),
            "uniform_c": match t.temperature {
                TemperatureSource::Uniform(c) => json!(c),
                TemperatureSource::ThermalSetup(_) => Value::Null,
            },
        })),
        "build_direction": setup.build_direction,
        "refinement": refinement,
        "samples": samples,
        "submodels": submodels,
    });
    let vtu = write_vtu(
        model,
        &VtuFields {
            displacement: Some(&result.vertex_displacement),
            stresses: Some(&result.stresses),
            temperature: None,
            heat_flux: None,
        },
    );
    let von_mises = result
        .stresses
        .iter()
        .map(|s| s.vertices.iter().map(|v| v.von_mises).collect())
        .collect();
    Ok((vtu, von_mises, "von_mises_mpa"))
}
