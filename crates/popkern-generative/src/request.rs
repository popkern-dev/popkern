//! What a generative design asks for, and the analysis setup it makes:
//! supports, loads and kept faces as design space face handles, lengths
//! in the design space's units.

use popkern_analysis::AnalysisError;
use popkern_analysis::face_ref::{Expect, FaceRef, Tolerance as SelectorTolerance};
use popkern_analysis::model::FaceId;
use popkern_analysis::options::AnalysisOptions;
use popkern_analysis::setup::{
    FaceDef, GenerativeSpec, LoadDirection, LoadSpec, StructuralSetup, SupportKind, SupportSpec,
};
use popkern_analysis::structural::generative::print::BuildDirection;
use popkern_core::topology::FaceHandle;

use crate::generate::GenerativeError;

/// A support on design space faces.
#[derive(Debug, Clone, PartialEq)]
pub struct SupportGroup {
    /// The faces.
    pub faces: Vec<FaceHandle>,
    /// What it holds.
    pub kind: SupportKind,
}

/// A total force spread over design space faces.
#[derive(Debug, Clone, PartialEq)]
pub struct LoadGroup {
    /// The faces.
    pub faces: Vec<FaceHandle>,
    /// Total force, N (positive; the sense is the direction's).
    pub force_n: f64,
    /// Direction the force acts in (non-zero).
    pub direction: [f64; 3],
}

/// What to design, on which faces. Lengths are in the design space's
/// units (`mm_per_unit` converts them).
#[derive(Debug, Clone, PartialEq)]
pub struct GenerativeRequest {
    /// Material preset name (`popkern_analysis::material::PRESETS`).
    pub material: String,
    /// Supports (at least one).
    pub supports: Vec<SupportGroup>,
    /// Loads (at least one).
    pub loads: Vec<LoadGroup>,
    /// Faces whose nearby material is kept.
    pub keep: Vec<FaceHandle>,
    /// Material to keep, as a fraction of the design space, in (0, 1].
    pub volume_fraction: f64,
    /// Density filter radius (the smallest member is about twice it).
    pub filter: f64,
    /// How deep kept material reaches (`None`: the filter radius).
    pub keep_depth: Option<f64>,
    /// Analysis cell size (`None`: a third of the filter radius).
    pub cell: Option<f64>,
    /// The build direction the design prints along without supports.
    pub build_direction: Option<BuildDirection>,
    /// Most updates.
    pub max_iterations: usize,
    /// Whether to thicken the thinnest members (robust formulation).
    pub robust: bool,
    /// Most active analysis cells (`None`: the solver's default).
    pub max_cells: Option<usize>,
    /// Millimetres per design space unit.
    pub mm_per_unit: f64,
}

impl GenerativeRequest {
    pub(crate) fn check(&self) -> Result<(), GenerativeError> {
        let bad = |what: &str| Err(GenerativeError::Request(what.to_string()));
        if !(self.mm_per_unit.is_finite() && self.mm_per_unit > 0.0) {
            return bad("the unit scale must be positive");
        }
        if self.supports.is_empty() || self.supports.iter().any(|s| s.faces.is_empty()) {
            return bad("every support needs a face, and the design at least one support");
        }
        if self.loads.is_empty() || self.loads.iter().any(|l| l.faces.is_empty()) {
            return bad("every load needs a face, and the design at least one load");
        }
        Ok(())
    }

    /// The setup the solver reads: one face definition per distinct
    /// face, by its analysis id.
    pub(crate) fn setup(
        &self,
        id_of: &dyn Fn(FaceHandle) -> Option<FaceId>,
    ) -> Result<StructuralSetup, GenerativeError> {
        let mut faces: Vec<FaceDef> = Vec::new();
        let mut name = |h: FaceHandle| -> Result<String, GenerativeError> {
            let id = id_of(h).ok_or_else(|| {
                GenerativeError::Request(format!("face {h:?} is not on the design space"))
            })?;
            let n = format!("f{:06}", id.0);
            if !faces.iter().any(|f| f.name == n) {
                faces.push(FaceDef {
                    name: n.clone(),
                    reference: FaceRef::Id(id),
                    expect: Expect::Exactly(1),
                });
            }
            Ok(n)
        };
        let mut names = |hs: &[FaceHandle]| -> Result<Vec<String>, GenerativeError> {
            hs.iter().map(|&h| name(h)).collect()
        };
        let supports = self
            .supports
            .iter()
            .map(|s| {
                Ok(SupportSpec {
                    faces: names(&s.faces)?,
                    kind: s.kind,
                })
            })
            .collect::<Result<Vec<_>, GenerativeError>>()?;
        let loads = self
            .loads
            .iter()
            .map(|l| {
                Ok(LoadSpec {
                    faces: names(&l.faces)?,
                    total_force_n: l.force_n,
                    direction: LoadDirection::Vector(l.direction),
                })
            })
            .collect::<Result<Vec<_>, GenerativeError>>()?;
        let keep = names(&self.keep)?;
        faces.sort_by(|a, b| a.name.cmp(&b.name));
        let mm = self.mm_per_unit;
        let mut options = AnalysisOptions {
            target_cell_size: self.cell.map(|c| c * mm),
            ..AnalysisOptions::default()
        };
        if let Some(cap) = self.max_cells {
            options.max_active_cells = cap;
        }
        let build = self.build_direction.map(|b| {
            let mut v = [0.0; 3];
            v[b.axis] = if b.negative { -1.0 } else { 1.0 };
            v
        });
        let setup = StructuralSetup {
            material: self.material.clone(),
            build_direction: build,
            tolerance: SelectorTolerance::default(),
            faces,
            supports,
            loads,
            thermal: None,
            modal: None,
            infill: None,
            generative: Some(GenerativeSpec {
                volume_fraction: self.volume_fraction,
                filter_mm: self.filter * mm,
                keep,
                keep_mm: self.keep_depth.map(|k| k * mm),
                max_iterations: self.max_iterations,
                crisp: true,
                robust: self.robust,
                build_direction: self.build_direction,
            }),
            options,
        };
        setup
            .validate()
            .map_err(|e| GenerativeError::Solver(AnalysisError::InvalidSetup(e)))?;
        Ok(setup)
    }
}
