//! The `simulate` request (PK-948, docs/architecture/analysis-integration.md):
//! a structural or thermal setup, in the analysis setup file format,
//! solved on the session's part.
//!
//! The part's live bodies become one analysis model, tessellated at
//! popkern live's display bound. Each face carries its stable popkern
//! name, spelled as the display mesh's `face_labels` spell it, so a
//! setup that names the faces the panel picked resolves again after a
//! parametric edit. The reply carries the analysis surface (welded per
//! solid) with per-vertex fields, so the overlay draws fields on the
//! mesh they were sampled on.
//!
//! A solver or setup decline is a result, not a failure: `status` is
//! `"declined"` with the typed kind and message. A session without a
//! part document fails as every other request does.

use std::collections::BTreeMap;

use popkern_analysis::AnalysisError;
use popkern_analysis::assess::{AssessOptions, assess};
use popkern_analysis::material::PRESETS;
use popkern_analysis::model::{AnalysisModel, Provenance};
use popkern_analysis::recover::VertexStress;
use popkern_analysis::setup::TemperatureSource;
use popkern_analysis::structural::from_setup::solve_structural_with;
use popkern_analysis::thermal::from_setup::solve_thermal;
use popkern_analysis_brep::{FeatureTolerances, model_from_breps_named, resolve_faces};
use popkern_analysis_setup::{parse_setup, parse_thermal_setup};
use popkern_assembly::Document;
use popkern_core::geometry::Transform3;
use popkern_core::topology::BRep;
use popkern_history::unit::Unit;
use popkern_mesh::TessellationOptions;
use popkern_rhai::Outcome;
use serde_json::{Value, json};

use crate::session::{Session, no_model};

/// One solid's analysis surface and its per-vertex fields.
#[derive(Debug, Clone, PartialEq)]
pub struct FieldMesh {
    /// The live body this solid belongs to (`BodyMesh::id`).
    pub body: u64,
    /// Vertex positions in model units, xyz per vertex.
    pub positions: Vec<f32>,
    /// Triangle vertex indices.
    pub indices: Vec<u32>,
    /// Per-vertex scalar fields by name (`von_mises_mpa`,
    /// `max_principal_mpa`, `displacement_mm`, `zz_error_mpa`, or
    /// `temperature_c`).
    pub fields: BTreeMap<&'static str, Vec<f32>>,
    /// Per-vertex displacement, xyz in mm (structural only, else empty).
    pub displacement: Vec<f32>,
}

/// What a `simulate` request produced.
#[derive(Debug, Clone, PartialEq)]
pub struct SimulateResult {
    /// `status` (`"solved"` or `"declined"`), `analysis`, and on a
    /// solve the summary (`result`) and the solver's report (`solve`);
    /// on a decline `error.kind` and `error.message`.
    pub summary: Value,
    /// Field meshes, one per solid; empty on a decline.
    pub meshes: Vec<FieldMesh>,
}

/// A typed decline: the error kind and its message.
struct Declined {
    kind: &'static str,
    message: String,
}

impl Declined {
    fn new(kind: &'static str, e: impl std::fmt::Display) -> Declined {
        Declined {
            kind,
            message: e.to_string(),
        }
    }
}

impl From<AnalysisError> for Declined {
    fn from(e: AnalysisError) -> Declined {
        Declined::new(e.kind(), e)
    }
}

/// The material presets a setup may name, for the panel's dropdown:
/// `[{name, youngs_modulus_mpa, behaviour}]`.
pub fn materials() -> Value {
    PRESETS
        .iter()
        .map(|m| {
            json!({
                "name": m.name,
                "youngs_modulus_mpa": m.youngs_modulus.value(),
                "printed": m.strength.interlayer.value().is_some(),
            })
        })
        .collect()
}

impl Session {
    /// Solve `setup` (setup file text: structural, or `analysis =
    /// "thermal"`) on the session's part.
    pub fn simulate(&mut self, setup: &str) -> Result<SimulateResult, Outcome> {
        let document = self.document.as_mut().ok_or_else(no_model)?;
        let Document::Part(model) = document else {
            return Ok(declined(
                "structural",
                Declined::new(
                    "Unsupported",
                    "simulation runs on a part document; solve each part of an assembly on its own",
                ),
            ));
        };
        // The solver works in millimetres. A mm, cm or inch document
        // stores true SI lengths (metres); a default-unit document keeps
        // the millimetre-magnitude convention (PK-462, the rule the STEP
        // export's units follow).
        let mm_per_unit = match model.length_unit() {
            Unit::Millimeters | Unit::Centimeters | Unit::Inches => 1000.0,
            _ => 1.0,
        };
        let to_mm = Transform3::scale_uniform(mm_per_unit);
        let mut breps: Vec<BRep> = Vec::new();
        let mut body_ids = Vec::new();
        for body in model.alive_bodies() {
            let Ok(brep) = model.brep_of(body) else {
                continue;
            };
            let brep = if mm_per_unit == 1.0 {
                brep
            } else {
                match brep.transform(&to_mm) {
                    Ok(scaled) => scaled,
                    Err(e) => {
                        return Ok(declined(
                            if setup_is_thermal(setup) {
                                "thermal"
                            } else {
                                "structural"
                            },
                            Declined::new("Model", e),
                        ));
                    }
                }
            };
            breps.push(brep);
            body_ids.push(body.id.0 as u64);
        }
        let namer = |index: usize, handle| {
            let brep = breps.get(index)?;
            let name = brep.faces.get(handle)?.name?;
            let entity = brep.name_pool.get(name)?;
            model.spell_entity(entity).ok()
        };
        let thermal = setup_is_thermal(setup);
        let analysis = if thermal { "thermal" } else { "structural" };
        let named = match model_from_breps_named(
            &breps,
            &TessellationOptions::display(),
            &FeatureTolerances::default(),
            Provenance::Popkern {
                document: "popkern live".into(),
                revision: 0,
            },
            &namer,
        ) {
            Ok(named) => named,
            Err(e) => return Ok(declined(analysis, Declined::new("Model", e))),
        };
        let solids_body: Vec<u64> = named
            .model
            .solids
            .iter()
            .map(|s| {
                named
                    .model
                    .faces
                    .iter()
                    .find(|f| f.solid == s.id)
                    .and_then(|f| named.face_sources.get(f.id.0 as usize))
                    .and_then(|(brep, _)| body_ids.get(*brep))
                    .copied()
                    .unwrap_or(0)
            })
            .collect();
        let solved = if thermal {
            solve_thermal_setup(&named.model, setup)
        } else {
            solve_structural_setup(&named.model, setup)
        };
        Ok(match solved {
            Ok((result, solve, fields)) => SimulateResult {
                summary: json!({
                    "status": "solved",
                    "analysis": analysis,
                    "model": {
                        "solids": named.model.solids.len(),
                        "faces": named.model.faces.len(),
                        "named_faces": named.model.faces.iter().filter(|f| f.name.is_some()).count(),
                    },
                    "result": result,
                    "solve": solve,
                    // Mesh positions are in model units; fields in mm,
                    // MPa and degrees C.
                    "mm_per_unit": mm_per_unit,
                }),
                meshes: field_meshes(&named.model, &solids_body, fields, mm_per_unit),
            },
            Err(d) => declined(analysis, d),
        })
    }
}

fn declined(analysis: &str, d: Declined) -> SimulateResult {
    SimulateResult {
        summary: json!({
            "status": "declined",
            "analysis": analysis,
            "error": { "kind": d.kind, "message": d.message },
        }),
        meshes: Vec::new(),
    }
}

/// Whether the setup's top-level `analysis` key says thermal.
fn setup_is_thermal(setup: &str) -> bool {
    setup.lines().any(|line| {
        let code = line.split('#').next().unwrap_or("").trim();
        code.strip_prefix("analysis")
            .map(|rest| rest.trim_start())
            .and_then(|rest| rest.strip_prefix('='))
            .is_some_and(|value| value.trim().trim_matches('"') == "thermal")
    })
}

/// Per-solid per-vertex fields by name, and per-solid displacements.
type Fields = (Vec<BTreeMap<&'static str, Vec<f64>>>, Vec<Vec<[f64; 3]>>);

fn to_value(v: &impl serde::Serialize) -> Result<Value, Declined> {
    serde_json::to_value(v).map_err(|e| Declined::new("Serialize", e))
}

fn solve_structural_setup(
    model: &AnalysisModel,
    text: &str,
) -> Result<(Value, Value, Fields), Declined> {
    let setup = parse_setup(text).map_err(|e| Declined::new("SetupFile", e))?;
    if setup.modal.is_some() || setup.generative.is_some() {
        return Err(Declined::new(
            "Unsupported",
            "modal and generative runs are not in popkern live yet; run them with popkern-analysis",
        ));
    }
    if matches!(
        setup.thermal.as_ref().map(|t| &t.temperature),
        Some(TemperatureSource::ThermalSetup(_))
    ) {
        return Err(Declined::new(
            "Unsupported",
            "a thermoelastic setup naming a thermal setup file runs with popkern-analysis; use a uniform temperature here",
        ));
    }
    let resolution = resolve_faces(model, &setup.faces, &setup.tolerance)?;
    let (result, material) = solve_structural_with(model, &setup, &resolution.faces, None)?;
    let options = AssessOptions {
        build_direction: setup.build_direction,
        ..AssessOptions::default()
    };
    let assessment =
        assess(model, &result, material, &options).map_err(|e| Declined::new("Assess", e))?;
    let mut solve = result.report.clone();
    let mut warnings = resolution.warnings.clone();
    warnings.append(&mut solve.warnings);
    solve.warnings = warnings;
    let summary = json!({
        "material": material.name,
        "max_displacement_mm": result.max_displacement,
        "global_error": result.global_error,
        "assessment": with_peak_labels(model, to_value(&assessment)?),
        "resolved": to_value(&resolution.faces)?,
    });
    let fields = result
        .stresses
        .iter()
        .zip(&result.vertex_displacement)
        .map(|(stresses, displacement)| {
            let mut f = BTreeMap::new();
            let column = |pick: fn(&VertexStress) -> f64| {
                stresses.vertices.iter().map(pick).collect::<Vec<f64>>()
            };
            f.insert("von_mises_mpa", column(|v| v.von_mises));
            f.insert("max_principal_mpa", column(|v| v.max_principal));
            f.insert("zz_error_mpa", column(|v| v.error.unwrap_or(-1.0)));
            f.insert(
                "displacement_mm",
                displacement
                    .iter()
                    .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
                    .collect(),
            );
            f
        })
        .collect();
    Ok((
        summary,
        to_value(&solve)?,
        (fields, result.vertex_displacement.clone()),
    ))
}

/// Each peak's face, by the stable name the panel shows (`face_label`;
/// null for a face without one).
fn with_peak_labels(model: &AnalysisModel, mut assessment: Value) -> Value {
    let label = |face: Option<u64>| {
        face.and_then(|f| model.faces.get(usize::try_from(f).ok()?))
            .and_then(|f| f.name.clone())
    };
    if let Some(peaks) = assessment.get_mut("peaks").and_then(Value::as_array_mut) {
        for solid in peaks {
            for key in ["von_mises", "max_principal"] {
                if let Some(peak) = solid.get_mut(key) {
                    let name = label(peak.get("face").and_then(Value::as_u64));
                    peak["face_label"] = json!(name);
                }
            }
        }
    }
    assessment
}

fn solve_thermal_setup(
    model: &AnalysisModel,
    text: &str,
) -> Result<(Value, Value, Fields), Declined> {
    let setup = parse_thermal_setup(text).map_err(|e| Declined::new("ThermalFile", e))?;
    let resolution = resolve_faces(model, &setup.faces, &setup.tolerance)?;
    let result = solve_thermal(model, &setup, &resolution.faces)?;
    let mut solve = result.report.clone();
    let mut warnings = resolution.warnings.clone();
    warnings.append(&mut solve.warnings);
    solve.warnings = warnings;
    let summary = json!({
        "material": setup.material,
        "max_temperature_c": result.max_temperature_c,
        "max_heat_flux_w_per_m2": result.max_heat_flux_w_per_m2,
        "flux_in_w": result.heat.flux_in,
        "source_w": result.heat.source,
        "resolved": to_value(&resolution.faces)?,
    });
    let fields = result
        .vertex_temperature
        .iter()
        .map(|t| BTreeMap::from([("temperature_c", t.clone())]))
        .collect();
    Ok((summary, to_value(&solve)?, (fields, Vec::new())))
}

fn field_meshes(
    model: &AnalysisModel,
    bodies: &[u64],
    (fields, displacement): Fields,
    mm_per_unit: f64,
) -> Vec<FieldMesh> {
    model
        .solids
        .iter()
        .enumerate()
        .map(|(i, s)| FieldMesh {
            body: bodies.get(i).copied().unwrap_or(0),
            positions: s
                .vertices
                .iter()
                .flatten()
                .map(|&c| (c / mm_per_unit) as f32)
                .collect(),
            indices: s.triangles.iter().flatten().copied().collect(),
            fields: fields
                .get(i)
                .map(|f| {
                    f.iter()
                        .map(|(k, v)| (*k, v.iter().map(|&x| x as f32).collect()))
                        .collect()
                })
                .unwrap_or_default(),
            displacement: displacement
                .get(i)
                .map(|d| d.iter().flatten().map(|&c| c as f32).collect())
                .unwrap_or_default(),
        })
        .collect()
}
