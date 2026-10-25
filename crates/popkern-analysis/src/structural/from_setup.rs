//! A [`StructuralSetup`] with resolved faces → an [`ElasticProblem`].
//!
//! A load's total force is spread uniformly over the area of its faces:
//! a traction F/A along the given direction, or a pressure F/A along the
//! inward normal. A roller constrains the face-normal displacement, so
//! its direction must be the normal of each (planar) face it is on;
//! anything else is declined, never silently re-interpreted.

use crate::error::AnalysisError;
use crate::material::{Material, Property, preset};
use crate::model::{AnalysisModel, FaceId, Surface};
use crate::report::ResolvedFace;
use crate::setup::{
    Axis, LoadDirection, SetupError, StructuralSetup, SupportKind, TemperatureSource,
};

use super::assemble::FaceMechanics;
use super::{ElasticProblem, StructuralResult};
use crate::element::vector_boundary::PieceSupport;

/// Largest angle (as 1 − |cos|) between a roller's direction and its
/// face's normal.
pub const ROLLER_NORMAL_TOL: f64 = 1.0e-6;

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

fn unit(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    v.map(|c| c / n)
}

fn axis(a: Axis) -> [f64; 3] {
    match a {
        Axis::PlusX => [1.0, 0.0, 0.0],
        Axis::MinusX => [-1.0, 0.0, 0.0],
        Axis::PlusY => [0.0, 1.0, 0.0],
        Axis::MinusY => [0.0, -1.0, 0.0],
        Axis::PlusZ => [0.0, 0.0, 1.0],
        Axis::MinusZ => [0.0, 0.0, -1.0],
    }
}

/// The elastic problem of `setup` on `model` with faces `resolved`, and
/// its material.
pub fn problem_from_setup(
    model: &AnalysisModel,
    setup: &StructuralSetup,
    resolved: &[ResolvedFace],
) -> Result<(ElasticProblem, &'static Material), AnalysisError> {
    setup.validate()?;
    let material = preset(&setup.material).map_err(SetupError::from)?;
    let youngs = material
        .youngs_modulus
        .value()
        .ok_or(SetupError::MissingMaterial)?;
    let poisson = material
        .poisson_ratio
        .value()
        .ok_or(SetupError::MissingMaterial)?;
    let mut faces = vec![FaceMechanics::Free; model.faces.len()];
    for (i, support) in setup.supports.iter().enumerate() {
        let list = ids(resolved, &support.faces)?;
        let mechanics = match support.kind {
            SupportKind::Fixed => PieceSupport::Fixed {
                u0: [0.0; 3],
                gradient: [[0.0; 3]; 3],
            },
            SupportKind::Roller { normal } => {
                let n = unit(normal);
                for f in &list {
                    let along = match model.faces.get(f.0 as usize).map(|f| &f.surface) {
                        Some(Surface::Plane { normal: m, .. }) => {
                            (n[0] * m[0] + n[1] * m[1] + n[2] * m[2]).abs()
                        }
                        _ => 0.0,
                    };
                    if 1.0 - along > ROLLER_NORMAL_TOL {
                        return Err(SetupError::RollerNotNormal {
                            used_by: format!("support {i}"),
                            face: f.0,
                        }
                        .into());
                    }
                }
                PieceSupport::Roller
            }
            SupportKind::SurfaceRoller => PieceSupport::Roller,
            SupportKind::DirectionalRoller { direction } => PieceSupport::Directional {
                direction: unit(direction),
            },
        };
        for f in list {
            if let Some(slot) = faces.get_mut(f.0 as usize) {
                *slot = FaceMechanics::Support(mechanics);
            }
        }
    }
    for load in &setup.loads {
        let list = ids(resolved, &load.faces)?;
        let area = super::super::thermal::faces_area(model, &list);
        let magnitude = load.total_force_n / area;
        let mechanics = match load.direction {
            LoadDirection::FaceNormal => FaceMechanics::Pressure(magnitude),
            LoadDirection::Axis(a) => FaceMechanics::Traction(axis(a).map(|c| c * magnitude)),
            LoadDirection::Vector(v) => FaceMechanics::Traction(unit(v).map(|c| c * magnitude)),
        };
        for f in list {
            if let Some(slot) = faces.get_mut(f.0 as usize) {
                *slot = mechanics;
            }
        }
    }
    Ok((
        ElasticProblem {
            youngs,
            poisson,
            faces,
        },
        material,
    ))
}

/// Solves `setup` on `model` with faces `resolved`, recording the
/// resolution in the report.
pub fn solve_structural(
    model: &AnalysisModel,
    setup: &StructuralSetup,
    resolved: &[ResolvedFace],
) -> Result<(StructuralResult, &'static Material), AnalysisError> {
    solve_structural_with(model, setup, resolved, None)
}

/// [`solve_structural`], with the temperature field a
/// `TemperatureSource::ThermalSetup` names, solved by the caller.
pub fn solve_structural_with(
    model: &AnalysisModel,
    setup: &StructuralSetup,
    resolved: &[ResolvedFace],
    field: Option<&crate::thermal::TemperatureField>,
) -> Result<(StructuralResult, &'static Material), AnalysisError> {
    let (problem, material) = problem_from_setup(model, setup, resolved)?;
    let infill = setup.infill.as_ref().map(|i| super::infill::Infill {
        density: i.density,
        walls_mm: i.walls_mm,
        line_width_mm: i.line_width_mm,
    });
    let uniform;
    let from_field;
    let load = match &setup.thermal {
        None => None,
        Some(spec) => {
            let alpha = match material.thermal_expansion {
                Property::Cited { value, .. } | Property::Typical { value, .. } => value,
                _ => {
                    return Err(crate::setup::SetupError::from(
                        crate::material::MaterialError::MissingProperty {
                            material: material.name,
                            property: "thermal expansion",
                        },
                    )
                    .into());
                }
            };
            let temperature: &dyn Fn([f64; 3]) -> Option<f64> = match (&spec.temperature, field) {
                (TemperatureSource::Uniform(c), _) => {
                    let c = *c;
                    uniform = move |_: [f64; 3]| Some(c);
                    &uniform
                }
                (TemperatureSource::ThermalSetup(_), Some(f)) => {
                    from_field = move |p: [f64; 3]| f.temperature_at(p).ok().flatten();
                    &from_field
                }
                (TemperatureSource::ThermalSetup(_), None) => {
                    return Err(AnalysisError::Unsupported {
                        what: "a thermal setup's field must be solved by the caller first",
                    });
                }
            };
            Some(super::thermoelastic::ThermalLoad {
                alpha,
                reference_c: spec.reference_c,
                temperature,
            })
        }
    };
    let mut result = super::solve_with(
        model,
        &problem,
        load.as_ref(),
        infill.as_ref(),
        &setup.options,
    )?;
    result.report.resolved_faces = resolved.to_vec();
    Ok((result, material))
}

/// A crisp design's final projection sharpness (Wang, Lazarov &
/// Sigmund 2011 continue to 128; 32 leaves under 6% grey on the
/// cantilevers measured).
const CRISP_BETA: f64 = 32.0;

/// The robust formulation's threshold offset: eroded at 0.7, dilated at
/// 0.3 (Wang, Lazarov & Sigmund 2011's examples). Measured on the crisp
/// 24 × 8 × 4 cantilever (r = 1.5 mm, 0.5 mm cells): thinnest members
/// 1.0 mm against 0.5, for 23% more compliance (PH-124, PH-128).
const ROBUST_ETA: f64 = 0.3;

/// CG's relative residual while a setup's design moves (PH-126): on the
/// crisp 12 × 4 × 2 cantilever 52% fewer CG iterations than at 1e-8 for
/// a compliance 0.21% lower; the final design is solved at the
/// analysis tolerance.
const UPDATE_TOLERANCE: f64 = 1e-3;

/// The generative design a setup asks for (`setup.generative`): the
/// stiffest part at its volume fraction under its supports and loads,
/// keeping the material near its named faces.
pub fn solve_generative_setup(
    model: &AnalysisModel,
    setup: &StructuralSetup,
    resolved: &[ResolvedFace],
) -> Result<(super::generative::Design, &'static Material), AnalysisError> {
    let spec = setup
        .generative
        .as_ref()
        .ok_or(AnalysisError::Unsupported {
            what: "a generative run without [generative]",
        })?;
    let (problem, material) = problem_from_setup(model, setup, resolved)?;
    let keep = ids(resolved, &spec.keep)?;
    let mut design =
        super::generative::GenerativeOptions::new(spec.volume_fraction, spec.filter_mm);
    design.max_iterations = spec.max_iterations;
    if let Some(k) = spec.keep_mm {
        design.keep_thickness = k;
    }
    // MMA throughout: it carries the projection and the print filter,
    // and without them it found stiffer designs than the optimality
    // criteria (PH-122).
    design.method = super::generative::UpdateMethod::Mma;
    if spec.crisp {
        design.projection_beta = CRISP_BETA;
        design.projection_kind = super::generative::projection::ProjectionKind::Threshold(0.5);
    }
    design.build_direction = spec.build_direction;
    design.update_tolerance = UPDATE_TOLERANCE;
    if spec.robust {
        design.robust = Some(ROBUST_ETA);
    }
    let options = generative_cells(model, spec.filter_mm, &setup.options)?;
    let mut result = super::generative::optimise(model, &problem, &design, &keep, &options)?;
    result.report.resolved_faces = resolved.to_vec();
    Ok((result, material))
}

/// A generative setup's options: without a cell size, cells of r/3
/// (`generative::FILTER_CELLS` per filter radius: coarser designs drift
/// under refinement, PH-124), declined as over the budget when the
/// part's box needs more than `max_active_cells` of them.
fn generative_cells(
    model: &AnalysisModel,
    filter_mm: f64,
    options: &crate::options::AnalysisOptions,
) -> Result<crate::options::AnalysisOptions, AnalysisError> {
    if options.target_cell_size.is_some() {
        return Ok(*options);
    }
    let h = filter_mm / super::generative::FILTER_CELLS;
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for v in model.solids.iter().flat_map(|s| &s.vertices) {
        for a in 0..3 {
            lo[a] = lo[a].min(v[a]);
            hi[a] = hi[a].max(v[a]);
        }
    }
    let chosen = crate::options::AnalysisOptions {
        target_cell_size: Some(h),
        ..*options
    };
    let grid = crate::grid::Grid::for_model(lo, hi, &chosen)
        .map_err(|e| super::internal("generative cells", e))?;
    let requested = grid.interior_cell_count();
    if requested > options.max_active_cells {
        return Err(AnalysisError::CellBudgetExceeded {
            requested,
            cap: options.max_active_cells,
        });
    }
    Ok(chosen)
}

/// The modes a modal setup asks for (`setup.modal`), with the
/// material's cited density; loads are ignored.
pub fn solve_modal_setup(
    model: &AnalysisModel,
    setup: &StructuralSetup,
    resolved: &[ResolvedFace],
) -> Result<(super::modal::ModalResult, &'static Material), AnalysisError> {
    let (problem, material) = problem_from_setup(model, setup, resolved)?;
    let density = match material.density {
        Property::Cited { value, .. } | Property::Typical { value, .. } => value,
        _ => {
            return Err(
                SetupError::from(crate::material::MaterialError::MissingProperty {
                    material: material.name,
                    property: "density",
                })
                .into(),
            );
        }
    };
    let count = setup.modal.unwrap_or(1);
    let infill = setup.infill.as_ref().map(|i| super::infill::Infill {
        density: i.density,
        walls_mm: i.walls_mm,
        line_width_mm: i.line_width_mm,
    });
    let mut result = super::modal::solve_modal_with(
        model,
        &problem,
        density,
        count,
        infill.as_ref(),
        &setup.options,
    )?;
    result.report.resolved_faces = resolved.to_vec();
    Ok((result, material))
}
