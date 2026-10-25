//! Submodelling (PH-063): a box around a stress raiser re-solved at a
//! finer cell size, driven by the coarse solution.
//!
//! The solid is cut to the box (`model::box_cut`): its faces inside
//! keep their loads and supports, and the caps on the box's faces are
//! held (Nitsche, per piece) to the coarse displacement field, affine at
//! each piece from the field's value and gradient there. Everything else
//! is the ordinary solve. The coarse field is least accurate at the
//! raiser itself and good a few radii away (St Venant), so the box
//! should reach a few raiser radii beyond the peak.

use crate::error::{AnalysisError, Degeneracy};
use crate::model::cut_model::{BoxFace, CutModel, clear_of_caps, cut_model};
use crate::model::{AnalysisModel, FaceId, SolidId};
use crate::options::AnalysisOptions;

use super::assemble::FaceMechanics;
use super::{ElasticProblem, StructuralResult};

/// A solved submodel.
#[derive(Debug, Clone, PartialEq)]
pub struct Submodel {
    /// The cut solid as its own model (one solid; its faces first the
    /// parent's that reach into the box, then the caps).
    pub model: AnalysisModel,
    /// The submodel's solution.
    pub result: StructuralResult,
    /// The submodel face of each parent face, by parent `FaceId`.
    pub face_map: Vec<Option<FaceId>>,
    /// Box faces that cut the solid (axis, upper), with their caps'
    /// submodel face.
    pub caps: Vec<(BoxFace, FaceId)>,
    /// The box cut to (each bound moved off mesh vertices).
    pub lo: [f64; 3],
    /// Upper corner of that box.
    pub hi: [f64; 3],
}

impl Submodel {
    /// The largest von Mises (and where) over the submodel's surface
    /// samples on parent face `face`, at least `margin` from every cap.
    pub fn peak_on(&self, face: FaceId, margin: f64) -> Option<(f64, [f64; 3])> {
        self.peak_by(face, margin, |s| s.von_mises)
    }

    /// [`Submodel::peak_on`] of the largest principal stress.
    pub fn principal_peak_on(&self, face: FaceId, margin: f64) -> Option<(f64, [f64; 3])> {
        self.peak_by(face, margin, |s| s.max_principal)
    }

    fn peak_by(
        &self,
        face: FaceId,
        margin: f64,
        value: impl Fn(&crate::recover::VertexStress) -> f64,
    ) -> Option<(f64, [f64; 3])> {
        let sub = (*self.face_map.get(face.0 as usize)?)?;
        let clear = |p: [f64; 3]| clear_of_caps(&self.caps, self.lo, self.hi, p, margin);
        self.result
            .stresses
            .first()?
            .surface
            .iter()
            .filter(|s| s.face == sub && clear(s.position))
            .map(|s| (value(&s.stress), s.position))
            .fold(None, |best: Option<(f64, [f64; 3])>, x| match best {
                Some(b) if b.0 >= x.0 => Some(b),
                _ => Some(x),
            })
    }
}

impl Submodel {
    /// The submodel's surface samples on the parent's faces (parent
    /// ids), at least `clearance` from every cap: (position, von Mises,
    /// face).
    pub fn parent_samples(&self, clearance: f64) -> Vec<([f64; 3], f64, FaceId)> {
        let parent_of: Vec<Option<FaceId>> = {
            let mut back = vec![None; self.model.faces.len()];
            for (parent, sub) in self.face_map.iter().enumerate() {
                if let Some(slot) = sub.and_then(|s| back.get_mut(s.0 as usize)) {
                    *slot = Some(FaceId(parent as u32));
                }
            }
            back
        };
        let clear = |p: [f64; 3]| clear_of_caps(&self.caps, self.lo, self.hi, p, clearance);
        self.result
            .stresses
            .first()
            .map(|s| s.surface.as_slice())
            .unwrap_or_default()
            .iter()
            .filter(|s| clear(s.position))
            .filter_map(|s| {
                let face = (*parent_of.get(s.face.0 as usize)?)?;
                Some((s.position, s.stress.von_mises, face))
            })
            .collect()
    }
}

/// How a submodel's caps are driven by the coarse field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Caps {
    /// Held to the coarse displacement.
    Displacement,
    /// Loaded by the coarse traction, with a weak spring against rigid
    /// motion: right where the coarse solution is too stiff (a singular
    /// corner).
    Traction,
    /// Each cap held to the coarse displacement when a supported face of
    /// the parent lies beyond it (it carries the reaction), loaded by the
    /// coarse traction otherwise (it carries the load path): statically
    /// consistent where the coarse solution is too stiff, without the
    /// self-equilibrium all-traction caps lack.
    Mixed,
}

/// The weak spring of traction caps, as a fraction of E over the box's
/// half-size.
pub const CAP_SPRING: f64 = 1.0e-3;

/// Solves solid `solid` of `model` (loaded by `problem`, solved coarse
/// as `parent`) inside the box [`lo`, `hi`] with `options` (its cell
/// size the submodel's), its caps driven as `caps`.
#[allow(clippy::too_many_arguments)]
pub fn solve_submodel(
    model: &AnalysisModel,
    problem: &ElasticProblem,
    parent: &StructuralResult,
    solid: SolidId,
    lo: [f64; 3],
    hi: [f64; 3],
    options: &AnalysisOptions,
    caps_mode: Caps,
) -> Result<Submodel, AnalysisError> {
    let degenerate = || AnalysisError::DegenerateGeometry {
        solid: Some(solid),
        detail: Degeneracy::BoxCut,
    };
    let mesh = model.solids.get(solid.0 as usize).ok_or_else(degenerate)?;
    let cut = cut_model(model, solid, lo, hi).map_err(|_| degenerate())?;
    cut_budget(&cut, options)?;
    let toward_support = |(axis, upper): BoxFace| {
        let at = if upper { cut.hi[axis] } else { cut.lo[axis] };
        let beyond = |v: &[f64; 3]| if upper { v[axis] > at } else { v[axis] < at };
        mesh.triangles
            .iter()
            .zip(&mesh.triangle_face)
            .any(|(t, f)| {
                matches!(
                    problem.faces.get(f.0 as usize),
                    Some(FaceMechanics::Support(_))
                ) && t
                    .iter()
                    .any(|&k| mesh.vertices.get(k as usize).is_some_and(beyond))
            })
    };
    let mechanics: Vec<FaceMechanics> = cut
        .parent_face
        .iter()
        .enumerate()
        .map(|(f, parent)| match parent {
            Some(p) => problem
                .faces
                .get(p.0 as usize)
                .copied()
                .unwrap_or(FaceMechanics::Free),
            None => {
                let side = cut.caps.iter().find(|c| c.1.0 as usize == f).map(|c| c.0);
                match caps_mode {
                    Caps::Displacement => FaceMechanics::Prescribed,
                    Caps::Traction => FaceMechanics::PrescribedTraction,
                    Caps::Mixed if side.is_some_and(toward_support) => FaceMechanics::Prescribed,
                    Caps::Mixed => FaceMechanics::PrescribedTraction,
                }
            }
        })
        .collect();
    let sub_problem = ElasticProblem {
        youngs: problem.youngs,
        poisson: problem.poisson,
        faces: mechanics,
    };
    let result = super::solve::solve_core(
        &cut.model,
        &sub_problem,
        options,
        crate::recover::Fallback::InCell,
        None,
        crate::moments::Splitting::Exact,
        &super::solve::Extras {
            prescribed: Some(super::assemble::Prescribed {
                field: &parent.field,
                spring: CAP_SPRING * problem.youngs
                    / (0..3)
                        .map(|a| 0.5 * (cut.hi[a] - cut.lo[a]))
                        .fold(f64::INFINITY, f64::min)
                        .max(f64::MIN_POSITIVE),
            }),
            ..super::solve::Extras::default()
        },
    )?;
    Ok(Submodel {
        result,
        model: cut.model,
        face_map: cut.face_map,
        caps: cut.caps,
        lo: cut.lo,
        hi: cut.hi,
    })
}

/// A raiser's submodel box reaches this many raiser radii from its
/// coarse peak (C6: 3 radii already within 1.2% of the whole plate at
/// the same cell size, 5 within 0.7%).
pub const SUBMODEL_RADII: f64 = 4.0;

/// A raiser's submodel resolves its radius with this many cells: past
/// the 6 a bore needs (PH-058) and the 4 a fillet needs.
pub const SUBMODEL_CELLS: f64 = 8.0;

/// Peaks are read this many submodel cells clear of the caps, where
/// the coarse field holds the solution.
pub const CAP_CLEARANCE_CELLS: f64 = 3.0;

/// Most cells a raiser's submodel may use: the interior cells of the
/// grid covering the cut part at the submodel's cell size (as a derived
/// budget counts them; a box of 8 radii at 8 cells per radius is at most
/// 64³ = 262,144). The box reaches two coarse cells when those are wider
/// than [`SUBMODEL_RADII`] radii, so the caps sit where the coarse
/// field is meaningful. A submodel past the budget is declined
/// (`cut_budget`). Its explicit cell size bypasses the solve's own
/// budget (PH-104): living_hinges' 0.4 mm raisers under 2.6 mm coarse
/// cells took 13 radii and 3.5M cells each, and their peaks moved by up
/// to a factor of 3 when the box was clipped to 6.25 radii, so no
/// clipped answer stands. Every other gallery submodel uses at most
/// 174k cells.
pub const SUBMODEL_CELL_BUDGET: usize = 1_000_000;

/// A submodel box's half-width: `radii` raiser radii, widened to two
/// coarse cells.
pub fn box_half(radii: f64, radius: f64, coarse_h: f64) -> f64 {
    (radii * radius).max(2.0 * coarse_h)
}

/// `CellBudgetExceeded` when the grid covering the cut part at the
/// cell size of `options` has more interior cells than its
/// `max_active_cells` (the coarse cells are too wide for the raiser:
/// refine the coarse solve first).
pub fn cut_budget(cut: &CutModel, options: &AnalysisOptions) -> Result<(), AnalysisError> {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for v in cut.model.solids.iter().flat_map(|s| &s.vertices) {
        for a in 0..3 {
            lo[a] = lo[a].min(v[a]);
            hi[a] = hi[a].max(v[a]);
        }
    }
    let grid =
        crate::grid::Grid::for_model(lo, hi, options).map_err(|e| AnalysisError::Internal {
            stage: "submodel",
            detail: e.to_string(),
        })?;
    let cells = grid.interior_cell_count();
    if cells > options.max_active_cells {
        return Err(AnalysisError::CellBudgetExceeded {
            requested: cells,
            cap: options.max_active_cells,
        });
    }
    Ok(())
}

/// A stress raiser re-solved in a submodel.
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct RefinedRaiser {
    /// The raiser's face (in the parent model).
    pub face: FaceId,
    /// Its radius, mm.
    pub radius: f64,
    /// The coarse solve's peak von Mises on the face, MPa.
    pub coarse_peak: f64,
    /// The submodel's peak von Mises on the face, MPa.
    pub peak_von_mises: f64,
    /// The submodel's largest principal stress on the face, MPa (what
    /// handbook K_t values give).
    pub peak_principal: f64,
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

/// Re-solves the raiser on `face` (radius `radius`) of the coarse
/// solution `parent` in a submodel: a box [`SUBMODEL_RADII`] radii
/// (at least two coarse cells) around the face's coarse peak, at
/// [`SUBMODEL_CELLS`] cells per radius, its caps held to the coarse
/// displacement (a resolved concentration barely changes the part's
/// stiffness: boss_plate 0.992 of CalculiX, against 0.981 with mixed
/// caps). `options` are the coarse solve's, with the budget
/// [`SUBMODEL_CELL_BUDGET`].
pub fn refine_raiser(
    model: &AnalysisModel,
    problem: &ElasticProblem,
    parent: &StructuralResult,
    face: FaceId,
    radius: f64,
    options: &AnalysisOptions,
) -> Result<RefinedRaiser, AnalysisError> {
    let policy = (SUBMODEL_RADII, SUBMODEL_CELLS);
    refine_raiser_with(model, problem, parent, (face, radius), policy, options)
}

/// [`refine_raiser`] for the raiser (face, radius) with the box reaching
/// `policy.0` radii and `policy.1` cells per radius.
pub fn refine_raiser_with(
    model: &AnalysisModel,
    problem: &ElasticProblem,
    parent: &StructuralResult,
    (face, radius): (FaceId, f64),
    policy: (f64, f64),
    options: &AnalysisOptions,
) -> Result<RefinedRaiser, AnalysisError> {
    let (radii, cells) = policy;
    let info = model
        .faces
        .get(face.0 as usize)
        .ok_or(AnalysisError::Internal {
            stage: "submodel",
            detail: format!("no face {}", face.0),
        })?;
    let solid = info.solid;
    let samples = parent
        .stresses
        .get(solid.0 as usize)
        .map(|s| s.surface.as_slice())
        .unwrap_or_default();
    let (coarse_peak, centre) = samples
        .iter()
        .filter(|s| s.face == face)
        .map(|s| (s.stress.von_mises, s.position))
        .fold((0.0, None), |best, x| {
            if x.0 > best.0 { (x.0, Some(x.1)) } else { best }
        });
    let centre = centre.ok_or(AnalysisError::Internal {
        stage: "submodel",
        detail: format!("no surface sample on face {}", face.0),
    })?;
    let coarse_h = parent.report.cell_size.unwrap_or(radius);
    let h = radius / cells;
    let half = box_half(radii, radius, coarse_h);
    let lo = centre.map(|c| c - half);
    let hi = centre.map(|c| c + half);
    let sub_options = AnalysisOptions {
        target_cell_size: Some(h),
        max_active_cells: SUBMODEL_CELL_BUDGET,
        grid_shift: [0.0; 3],
        ..*options
    };
    let sub = solve_submodel(
        model,
        problem,
        parent,
        solid,
        lo,
        hi,
        &sub_options,
        Caps::Displacement,
    )?;
    let (peak, position) =
        sub.peak_on(face, CAP_CLEARANCE_CELLS * h)
            .ok_or(AnalysisError::Internal {
                stage: "submodel",
                detail: format!("no submodel sample on face {} clear of the caps", face.0),
            })?;
    let peak_principal = sub
        .principal_peak_on(face, CAP_CLEARANCE_CELLS * h)
        .map_or(peak, |p| p.0);
    Ok(RefinedRaiser {
        face,
        radius,
        coarse_peak,
        peak_von_mises: peak,
        peak_principal,
        position,
        cell_size: h,
        cells_active: sub.result.report.cells_active,
        lo: sub.lo,
        hi: sub.hi,
    })
}

#[cfg(test)]
#[path = "submodel_tests.rs"]
mod tests;
