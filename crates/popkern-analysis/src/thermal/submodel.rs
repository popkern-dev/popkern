//! Thermal submodelling: a box around a bore or fillet re-solved at a
//! finer cell size for its heat-flux peak, driven by the coarse
//! temperature field.
//!
//! The solid is cut to the box (`model::cut_model`): its faces inside
//! keep their conditions, sources keep their density, and the caps on
//! the box's faces are held (Nitsche, per piece) to the coarse
//! temperature, affine at each piece from the field's value and
//! gradient there. As with structural submodels, the coarse field is
//! least accurate at the raiser and good a few radii away, so the box
//! should reach a few radii beyond the peak.

use crate::error::{AnalysisError, Degeneracy};
use crate::model::cut_model::{BoxFace, cut_model};
use crate::model::{AnalysisModel, FaceId, SolidId};
use crate::options::AnalysisOptions;
use crate::structural::submodel::{
    CAP_CLEARANCE_CELLS, SUBMODEL_CELL_BUDGET, SUBMODEL_CELLS, SUBMODEL_RADII, box_half, cut_budget,
};

use super::{ConductionProblem, FaceCondition, ThermalResult, solve_conduction_core};

/// A solved thermal submodel.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalSubmodel {
    /// The cut solid as its own model (its faces first the parent's
    /// that reach into the box, then the caps).
    pub model: AnalysisModel,
    /// The submodel's solution.
    pub result: ThermalResult,
    /// The submodel face of each parent face, by parent `FaceId`.
    pub face_map: Vec<Option<FaceId>>,
    /// Box faces that cut the solid, with their caps' submodel face.
    pub caps: Vec<(BoxFace, FaceId)>,
    /// The box cut to (each bound moved off mesh vertices).
    pub lo: [f64; 3],
    /// Upper corner of that box.
    pub hi: [f64; 3],
}

impl ThermalSubmodel {
    /// The largest heat flux magnitude (W/m²) and where, over the
    /// submodel's vertices on parent face `face`, at least `margin`
    /// from every cap.
    pub fn peak_flux_on(&self, face: FaceId, margin: f64) -> Option<(f64, [f64; 3])> {
        let sub = (*self.face_map.get(face.0 as usize)?)?;
        let mesh = self.model.solids.first()?;
        let flux = self.result.vertex_heat_flux.first()?;
        let clear = |p: [f64; 3]| {
            crate::model::cut_model::clear_of_caps(&self.caps, self.lo, self.hi, p, margin)
        };
        let mut best: Option<(f64, [f64; 3])> = None;
        for (t, f) in mesh.triangles.iter().zip(&mesh.triangle_face) {
            if *f != sub {
                continue;
            }
            for &v in t {
                let (Some(p), Some(q)) = (mesh.vertices.get(v as usize), flux.get(v as usize))
                else {
                    continue;
                };
                let m = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
                if clear(*p) && best.is_none_or(|b| m > b.0) {
                    best = Some((m, *p));
                }
            }
        }
        best
    }
}

/// Solves solid `solid` of `model` (loaded by `problem`, solved coarse
/// as `parent`) inside the box [`lo`, `hi`] with `options` (its cell
/// size the submodel's), its caps held to the coarse temperature.
pub fn solve_thermal_submodel(
    model: &AnalysisModel,
    problem: &ConductionProblem,
    parent: &ThermalResult,
    solid: SolidId,
    (lo, hi): ([f64; 3], [f64; 3]),
    options: &AnalysisOptions,
) -> Result<ThermalSubmodel, AnalysisError> {
    let degenerate = || AnalysisError::DegenerateGeometry {
        solid: Some(solid),
        detail: Degeneracy::BoxCut,
    };
    let mesh = model.solids.get(solid.0 as usize).ok_or_else(degenerate)?;
    let cut = cut_model(model, solid, lo, hi).map_err(|_| degenerate())?;
    cut_budget(&cut, options)?;
    let faces: Vec<FaceCondition> = cut
        .parent_face
        .iter()
        .map(|parent| match parent {
            Some(p) => problem
                .faces
                .get(p.0 as usize)
                .copied()
                .unwrap_or(FaceCondition::Insulated),
            None => FaceCondition::Prescribed,
        })
        .collect();
    // The solid's source density, W/mm³, carried into the cut volume.
    let total_volume: f64 = model.solids.iter().map(|s| s.enclosed_volume()).sum();
    let own: f64 = problem
        .body_sources
        .iter()
        .filter(|(id, _)| *id == solid)
        .map(|(_, w)| w)
        .sum();
    let density = problem.source_w / total_volume + own / mesh.enclosed_volume();
    let cut_volume = cut
        .model
        .solids
        .first()
        .map_or(0.0, |s| s.enclosed_volume());
    let sub_problem = ConductionProblem {
        conductivity: problem.conductivity,
        faces,
        source_w: 0.0,
        body_sources: if density == 0.0 {
            Vec::new()
        } else {
            vec![(SolidId(0), density * cut_volume)]
        },
        hdt_c: problem.hdt_c,
    };
    let result = solve_conduction_core(&cut.model, &sub_problem, options, Some(&parent.field))?;
    Ok(ThermalSubmodel {
        result,
        model: cut.model,
        face_map: cut.face_map,
        caps: cut.caps,
        lo: cut.lo,
        hi: cut.hi,
    })
}

/// A heat-flux raiser re-solved in a submodel.
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct RefinedFlux {
    /// The raiser's face (in the parent model).
    pub face: FaceId,
    /// Its radius, mm.
    pub radius: f64,
    /// The coarse solve's peak heat flux on the face, W/m².
    pub coarse_peak: f64,
    /// The submodel's peak heat flux on the face, W/m².
    pub peak_heat_flux: f64,
    /// Where the submodel's peak is.
    pub position: [f64; 3],
    /// The submodel's cell size, mm.
    pub cell_size: f64,
    /// Its active cells.
    pub cells_active: usize,
    /// Lower corner of the box.
    pub lo: [f64; 3],
    /// Upper corner of the box.
    pub hi: [f64; 3],
}

/// Re-solves the heat-flux raiser on `face` (radius `radius`) of the
/// coarse solution `parent` in a submodel, sized as the structural
/// `refine_raiser`: a box [`SUBMODEL_RADII`] radii (at least two coarse
/// cells) around the face's coarse peak, at [`SUBMODEL_CELLS`] cells per
/// radius, within [`SUBMODEL_CELL_BUDGET`].
pub fn refine_flux(
    model: &AnalysisModel,
    problem: &ConductionProblem,
    parent: &ThermalResult,
    (face, radius): (FaceId, f64),
    options: &AnalysisOptions,
) -> Result<RefinedFlux, AnalysisError> {
    let missing = |detail: String| AnalysisError::Internal {
        stage: "submodel",
        detail,
    };
    let info = model
        .faces
        .get(face.0 as usize)
        .ok_or_else(|| missing(format!("no face {}", face.0)))?;
    let solid = info.solid;
    let mesh = model
        .solids
        .get(solid.0 as usize)
        .ok_or_else(|| missing(format!("no solid {}", solid.0)))?;
    let flux = parent
        .vertex_heat_flux
        .get(solid.0 as usize)
        .map(|f| f.as_slice())
        .unwrap_or_default();
    let mut coarse: Option<(f64, [f64; 3])> = None;
    for (t, f) in mesh.triangles.iter().zip(&mesh.triangle_face) {
        if *f != face {
            continue;
        }
        for &v in t {
            if let (Some(p), Some(q)) = (mesh.vertices.get(v as usize), flux.get(v as usize)) {
                let m = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
                if coarse.is_none_or(|b| m > b.0) {
                    coarse = Some((m, *p));
                }
            }
        }
    }
    let (coarse_peak, centre) =
        coarse.ok_or_else(|| missing(format!("no vertex on face {}", face.0)))?;
    let coarse_h = parent.report.cell_size.unwrap_or(radius);
    let h = radius / SUBMODEL_CELLS;
    let half = box_half(SUBMODEL_RADII, radius, coarse_h);
    let sub_options = AnalysisOptions {
        target_cell_size: Some(h),
        max_active_cells: SUBMODEL_CELL_BUDGET,
        grid_shift: [0.0; 3],
        ..*options
    };
    let box_ = (centre.map(|c| c - half), centre.map(|c| c + half));
    let sub = solve_thermal_submodel(model, problem, parent, solid, box_, &sub_options)?;
    let (peak, position) = sub
        .peak_flux_on(face, CAP_CLEARANCE_CELLS * h)
        .ok_or_else(|| {
            missing(format!(
                "no submodel vertex on face {} clear of the caps",
                face.0
            ))
        })?;
    Ok(RefinedFlux {
        face,
        radius,
        coarse_peak,
        peak_heat_flux: peak,
        position,
        cell_size: h,
        cells_active: sub.result.report.cells_active,
        lo: sub.lo,
        hi: sub.hi,
    })
}

/// The concave fillets and bores of `model` (`assess::concave_raisers`)
/// with fewer than `resolved_cells` coarse cells across their radius and
/// a coarse peak heat flux at least `flux_fraction` of their solid's:
/// the raisers whose flux peak reads low (T3's bore at 2 cells per
/// radius: 0.77..0.87 of the wall value).
pub fn unresolved_flux_raisers(
    model: &AnalysisModel,
    parent: &ThermalResult,
    resolved_cells: f64,
    flux_fraction: f64,
) -> Result<Vec<(FaceId, f64)>, AnalysisError> {
    let h = parent.report.cell_size.unwrap_or(f64::INFINITY);
    let raisers = crate::assess::concave_raisers(model).map_err(|e| AnalysisError::Internal {
        stage: "submodel",
        detail: e.to_string(),
    })?;
    let magnitude = |q: &[f64; 3]| (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]).sqrt();
    let solid_peak: Vec<f64> = parent
        .vertex_heat_flux
        .iter()
        .map(|f| f.iter().map(magnitude).fold(0.0, f64::max))
        .collect();
    let mut out = Vec::new();
    for (face, radius) in raisers {
        if radius / h >= resolved_cells {
            continue;
        }
        let Some(solid) = model.faces.get(face.0 as usize).map(|f| f.solid.0 as usize) else {
            continue;
        };
        let (Some(mesh), Some(flux)) =
            (model.solids.get(solid), parent.vertex_heat_flux.get(solid))
        else {
            continue;
        };
        let peak = mesh
            .triangles
            .iter()
            .zip(&mesh.triangle_face)
            .filter(|(_, f)| **f == face)
            .flat_map(|(t, _)| t.iter())
            .filter_map(|&v| flux.get(v as usize))
            .map(magnitude)
            .fold(0.0, f64::max);
        if peak > 0.0 && peak >= flux_fraction * solid_peak.get(solid).copied().unwrap_or(0.0) {
            out.push((face, radius));
        }
    }
    Ok(out)
}

#[cfg(test)]
#[path = "submodel_tests.rs"]
mod tests;
