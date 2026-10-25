//! A modal run: the lowest natural frequencies of a structural setup
//! with `[modal]`, its first mode shape in the VTU.

use popkern_analysis::model::AnalysisModel;
use popkern_analysis::setup::StructuralSetup;
use popkern_analysis::structural::from_setup::solve_modal_setup;
use popkern_analysis_brep::Resolution;
use serde_json::{Value, json};

use super::{Failure, to_value, with_resolution};
use crate::vtu::{VtuFields, write_vtu};

/// Solves the modes, fills `report`, and returns the VTU text and the
/// first mode's displacement magnitude per vertex (scaled to 1 at its
/// largest), per solid.
pub(super) fn run_modal(
    model: &AnalysisModel,
    setup: &StructuralSetup,
    resolution: &Resolution,
    report: &mut Value,
) -> Result<(String, Vec<Vec<f64>>), Failure> {
    let (result, material) = solve_modal_setup(model, setup, &resolution.faces)?;
    let solve = with_resolution(result.report.clone(), &resolution.warnings);
    report["solve"] = to_value(&solve)?;
    let free = setup.supports.is_empty();
    report["result"] = json!({
        "material": material.name,
        "mass_kg": result.mass_kg,
        "free": free,
        "modes": to_value(&result.modes)?,
    });
    // The first mode (the first elastic one for a free part), scaled so
    // its largest displacement is 1 mm.
    let pick = if free { 6 } else { 0 };
    let shape = result
        .shapes
        .get(pick)
        .or(result.shapes.first())
        .cloned()
        .unwrap_or_default();
    let norm = |d: &[f64; 3]| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    let largest = shape.iter().flatten().map(norm).fold(0.0, f64::max);
    let scale = if largest > 0.0 { 1.0 / largest } else { 0.0 };
    let scaled: Vec<Vec<[f64; 3]>> = shape
        .iter()
        .map(|s| s.iter().map(|d| d.map(|c| c * scale)).collect())
        .collect();
    let field: Vec<Vec<f64>> = scaled
        .iter()
        .map(|s| s.iter().map(norm).collect())
        .collect();
    let vtu = write_vtu(
        model,
        &VtuFields {
            displacement: Some(&scaled),
            stresses: None,
            temperature: None,
            heat_flux: None,
        },
    );
    Ok((vtu, field))
}
