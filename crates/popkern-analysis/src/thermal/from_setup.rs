//! A [`ThermalSetup`] with resolved faces → a [`ConductionProblem`].

use crate::error::AnalysisError;
use crate::material::{Property, preset};
use crate::model::{AnalysisModel, FaceId};
use crate::report::ResolvedFace;
use crate::thermal_setup::{ConvectionFaces, Emittance, ThermalSetup, ThermalSetupError};

use super::{ConductionProblem, FaceCondition, ThermalResult, faces_area, solve_conduction};

fn ids(resolved: &[ResolvedFace], names: &[String]) -> Result<Vec<FaceId>, AnalysisError> {
    let mut out = Vec::new();
    for name in names {
        let r = resolved
            .iter()
            .find(|r| &r.name == name)
            .ok_or_else(|| AnalysisError::StaleBoundaryCondition { name: name.clone() })?;
        out.extend(&r.faces);
    }
    out.sort_unstable();
    out.dedup();
    Ok(out)
}

/// The conduction problem of `setup` on `model` with faces `resolved`.
pub fn problem_from_setup(
    model: &AnalysisModel,
    setup: &ThermalSetup,
    resolved: &[ResolvedFace],
) -> Result<ConductionProblem, AnalysisError> {
    let conductivity = setup.validate()?;
    let material = preset(&setup.material).map_err(crate::setup::SetupError::from)?;
    let mut faces = vec![FaceCondition::Insulated; model.faces.len()];
    let mut set = |list: &[FaceId], c: FaceCondition| {
        for f in list {
            if let Some(slot) = faces.get_mut(f.0 as usize) {
                *slot = c;
            }
        }
    };
    for f in &setup.fixed {
        set(
            &ids(resolved, &f.faces)?,
            FaceCondition::Fixed {
                t0_c: f.temperature_c,
                gradient: [0.0; 3],
            },
        );
    }
    for f in &setup.flux {
        let list = ids(resolved, &f.faces)?;
        let area = faces_area(model, &list);
        set(
            &list,
            FaceCondition::Flux {
                w_per_mm2: f.power_w / area,
            },
        );
    }
    let mut all_other = None;
    for c in &setup.convection {
        let condition = match c.emittance {
            None => FaceCondition::Convection {
                h: c.h.value(),
                ambient_c: c.ambient_c,
            },
            Some(e) => FaceCondition::Radiating {
                h: c.h.value(),
                emittance: match e {
                    Emittance::Value(v) => v,
                    Emittance::Material => match material.emittance_default {
                        Property::Cited { value, .. } | Property::Typical { value, .. } => value,
                        _ => {
                            return Err(crate::setup::SetupError::from(
                                crate::material::MaterialError::MissingProperty {
                                    material: material.name,
                                    property: "emittance",
                                },
                            )
                            .into());
                        }
                    },
                },
                ambient_c: c.ambient_c,
            },
        };
        match &c.faces {
            ConvectionFaces::Named(names) => set(&ids(resolved, names)?, condition),
            ConvectionFaces::AllOther => all_other = Some(condition),
        }
    }
    if let Some(condition) = all_other {
        for slot in &mut faces {
            if *slot == FaceCondition::Insulated {
                *slot = condition;
            }
        }
    }
    let hdt_c = match material.heat_deflection_temp {
        Property::Cited { value, .. } | Property::Typical { value, .. } => Some(value),
        _ => None,
    };
    let mut source_w = 0.0;
    let mut body_sources = Vec::new();
    for s in &setup.sources {
        let Some(name) = &s.body else {
            source_w += s.power_w;
            continue;
        };
        let solid = model
            .solids
            .iter()
            .find(|m| m.name.as_deref() == Some(name.as_str()))
            .ok_or_else(|| ThermalSetupError::UnknownBody {
                name: name.clone(),
                known: model.solids.iter().filter_map(|m| m.name.clone()).collect(),
            })?;
        body_sources.push((solid.id, s.power_w));
    }
    Ok(ConductionProblem {
        conductivity,
        faces,
        source_w,
        body_sources,
        hdt_c,
    })
}

/// Solves `setup` on `model` with faces `resolved`, recording the
/// resolution in the report.
pub fn solve_thermal(
    model: &AnalysisModel,
    setup: &ThermalSetup,
    resolved: &[ResolvedFace],
) -> Result<ThermalResult, AnalysisError> {
    let problem = problem_from_setup(model, setup, resolved)?;
    let mut result = match setup.transient {
        None => solve_conduction(model, &problem, &setup.options)?,
        Some(t) => {
            let material = preset(&setup.material).map_err(crate::setup::SetupError::from)?;
            let value = |p: Property, property: &'static str| match p {
                Property::Cited { value, .. } | Property::Typical { value, .. } => Ok(value),
                _ => Err(crate::setup::SetupError::from(
                    crate::material::MaterialError::MissingProperty {
                        material: material.name,
                        property,
                    },
                )),
            };
            let transient = super::transient::TransientProblem {
                conduction: problem,
                density: value(material.density, "density")?,
                specific_heat: value(material.specific_heat, "specific heat")?,
                initial_c: t.initial_c,
                step_s: t.step_s,
                steps: t.steps(),
            };
            super::transient::solve_transient(model, &transient, &setup.options)?.last
        }
    };
    result.report.resolved_faces = resolved.to_vec();
    Ok(result)
}
