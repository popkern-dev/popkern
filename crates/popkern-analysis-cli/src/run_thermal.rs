//! The thermal half of a run: solve, optional heat-flux submodels,
//! and the report's thermal result.

use popkern_analysis::AnalysisError;
use popkern_analysis::assess::AssessOptions;
use popkern_analysis::model::AnalysisModel;
use popkern_analysis::thermal::from_setup::{problem_from_setup as thermal_problem, solve_thermal};
use popkern_analysis::thermal::submodel::{refine_flux, unresolved_flux_raisers};
use popkern_analysis::thermal_setup::{ConvectionFaces, Emittance, HeatTransfer};
use popkern_analysis_brep::resolve_faces;
use serde_json::{Value, json};

use super::{Failure, analysis_error_kind, to_value, with_resolution};
use crate::thermal_file::{parse_thermal_setup, preset_name};
use crate::vtu::{VtuFields, write_vtu};

pub(super) fn run_thermal(
    model: &AnalysisModel,
    text: &str,
    refine: bool,
    report: &mut Value,
) -> Result<(String, Vec<Vec<f64>>), Failure> {
    let setup = parse_thermal_setup(text).map_err(|e| Failure::new("ThermalFile", e))?;
    let resolution = resolve_faces(model, &setup.faces, &setup.tolerance)?;
    if refine && setup.transient.is_some() {
        return Err(AnalysisError::Unsupported {
            what: "--refine on a transient (submodels are steady)",
        }
        .into());
    }
    let result = solve_thermal(model, &setup, &resolution.faces)?;
    // Unresolved bores and fillets carrying heat, each re-solved in a
    // submodel (judged as the structural bore check: 6 cells across the
    // radius, peak at least half the solid's).
    let submodels = if refine {
        let problem = thermal_problem(model, &setup, &resolution.faces)?;
        let assess = AssessOptions::default();
        let raisers = unresolved_flux_raisers(
            model,
            &result,
            assess.bore_resolved_cells,
            assess.stressed_fraction,
        )?;
        let entries: Vec<Value> = raisers
            .iter()
            .map(|&(face, radius)| {
                match refine_flux(model, &problem, &result, (face, radius), &setup.options) {
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
    let face_out: Vec<Value> = resolution
        .faces
        .iter()
        .map(|r| {
            let q: f64 = r
                .faces
                .iter()
                .filter_map(|f| result.heat.face_out.get(f.0 as usize))
                .sum();
            json!({ "name": r.name, "heat_out_w": q })
        })
        .collect();
    let solve = with_resolution(result.report.clone(), &resolution.warnings);
    report["solve"] = to_value(&solve)?;
    // Every convection condition with the coefficient actually used.
    let convection: Vec<Value> = setup
        .convection
        .iter()
        .map(|c| {
            let faces = match &c.faces {
                ConvectionFaces::Named(names) => json!(names),
                ConvectionFaces::AllOther => json!("all_other"),
            };
            let source = match c.h {
                HeatTransfer::Preset(p) => preset_name(p),
                HeatTransfer::Value(_) => "value",
            };
            let mut entry = json!({
                "faces": faces,
                "h_w_per_m2k": c.h.value(),
                "source": source,
                "ambient_c": c.ambient_c,
            });
            // Radiating entries only: their h is the convective part.
            match c.emittance {
                None => {}
                Some(Emittance::Value(v)) => entry["emittance"] = json!(v),
                Some(Emittance::Material) => entry["emittance"] = json!("material"),
            }
            entry
        })
        .collect();
    report["result"] = json!({
        "material": setup.material,
        "max_temperature_c": result.max_temperature_c,
        "max_heat_flux_w_per_m2": result.max_heat_flux_w_per_m2,
        "convection": convection,
        "flux_in_w": result.heat.flux_in,
        "source_w": result.heat.source,
        "face_heat": face_out,
        "submodels": submodels,
    });
    if !result.transient.is_empty() {
        // The fields and face heat above are the final state.
        report["result"]["transient"] = to_value(&result.transient)?;
    }
    let vtu = write_vtu(
        model,
        &VtuFields {
            displacement: None,
            stresses: None,
            temperature: Some(&result.vertex_temperature),
            heat_flux: Some(&result.vertex_heat_flux),
        },
    );
    Ok((vtu, result.vertex_temperature))
}
