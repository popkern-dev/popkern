//! A structural analysis setup: named face references, supports, loads,
//! material and options.
//!
//! The setup file format (TOML) is parsed by `popkern-analysis-cli`,
//! the only crate that touches files. This module holds the semantic
//! types and their validation, which do not depend on any format.

use crate::face_ref::{Expect, FaceRef, Tolerance};
use crate::material::{MaterialError, preset};
use crate::options::{AnalysisOptions, OptionError};

/// A complete structural setup.
#[derive(Debug, Clone, PartialEq)]
pub struct StructuralSetup {
    /// Material preset name; must resolve with `material::preset`.
    pub material: String,
    /// Build direction for layer-aware checks, when the user gives it.
    pub build_direction: Option<[f64; 3]>,
    /// Default selector tolerances.
    pub tolerance: Tolerance,
    /// Named face references, sorted by name.
    pub faces: Vec<FaceDef>,
    /// Supports.
    pub supports: Vec<SupportSpec>,
    /// Loads.
    pub loads: Vec<LoadSpec>,
    /// A thermal strain (`None`: isothermal).
    pub thermal: Option<ThermalSpec>,
    /// Natural frequencies instead of a static solve: how many modes
    /// (loads are ignored; without supports the part is free).
    pub modal: Option<usize>,
    /// A printed part's infill (`None`: solid).
    pub infill: Option<InfillSpec>,
    /// A generative design instead of a static solve: where to keep
    /// material for the stiffest part at a volume
    /// (`structural::generative`).
    pub generative: Option<GenerativeSpec>,
    /// Discretisation and solver options.
    pub options: AnalysisOptions,
}

/// A printed part's infill: solid walls around a homogenised core.
#[derive(Debug, Clone, PartialEq)]
pub struct InfillSpec {
    /// The slicer's pattern; only `"gyroid"` has homogenised properties.
    pub pattern: String,
    /// Relative density of the core, in (0, 1] (the measured or
    /// double-bead density, not the slicer's nominal one where they
    /// differ: Bean 2022).
    pub density: f64,
    /// Solid wall thickness, mm (perimeters × line width).
    pub walls_mm: f64,
    /// The gyroid's wall (bead) width, mm: with the density it sets the
    /// gyroid's unit-cell size (`structural::infill::Infill::cell_mm`).
    pub line_width_mm: f64,
}

/// A generative design: minimum compliance at a volume fraction.
#[derive(Debug, Clone, PartialEq)]
pub struct GenerativeSpec {
    /// Material to keep, as a fraction of the part's volume, in (0, 1].
    pub volume_fraction: f64,
    /// Density filter radius, mm (the smallest member is about twice
    /// it).
    pub filter_mm: f64,
    /// Named faces whose nearby material is kept (mounting faces).
    pub keep: Vec<String>,
    /// How deep the kept material reaches, mm (`None`: the filter
    /// radius).
    pub keep_mm: Option<f64>,
    /// Most updates.
    pub max_iterations: usize,
    /// Whether to push the design to solid and void (the threshold
    /// projection under MMA) rather than leave grey material.
    pub crisp: bool,
    /// Whether to thicken the design's thinnest members (the robust
    /// formulation; needs `crisp`).
    pub robust: bool,
    /// The build direction the design must print along without
    /// supports (`None`: no overhang control).
    pub build_direction: Option<crate::structural::generative::print::BuildDirection>,
}

/// Most modes a modal setup may ask for.
pub const MAX_MODES: usize = 50;

/// A thermal strain on a structural setup: the temperature over a
/// stress-free reference (`structural::thermoelastic`).
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalSpec {
    /// Stress-free temperature, °C.
    pub reference_c: f64,
    /// Where the temperature comes from.
    pub temperature: TemperatureSource,
}

/// The temperature of a thermoelastic solve.
#[derive(Debug, Clone, PartialEq)]
pub enum TemperatureSource {
    /// Uniform, °C.
    Uniform(f64),
    /// A thermal setup the caller solves on the same part first (its
    /// path as written in the setup; the library reads no files).
    ThermalSetup(String),
}

/// A named face reference that supports and loads use.
#[derive(Debug, Clone, PartialEq)]
pub struct FaceDef {
    /// Name used by supports and loads.
    pub name: String,
    /// How the faces are found.
    pub reference: FaceRef,
    /// How many faces it must match.
    pub expect: Expect,
}

/// A support on one or more named faces.
#[derive(Debug, Clone, PartialEq)]
pub struct SupportSpec {
    /// Names of [`FaceDef`] entries.
    pub faces: Vec<String>,
    /// What the support constrains.
    pub kind: SupportKind,
}

/// What a support constrains.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SupportKind {
    /// All displacement components.
    Fixed,
    /// Only the displacement along `normal`, which must be the normal
    /// of every (planar) face the support is on.
    Roller {
        /// Constrained direction (non-zero).
        normal: [f64; 3],
    },
    /// Only the displacement along each face's own normal: on a curved
    /// face it follows the surface (radial on a cylinder).
    SurfaceRoller,
    /// Only the displacement along `direction`, on any face; motion
    /// normal to it is free.
    DirectionalRoller {
        /// Constrained direction (non-zero).
        direction: [f64; 3],
    },
}

/// A total force spread over one or more named faces.
#[derive(Debug, Clone, PartialEq)]
pub struct LoadSpec {
    /// Names of [`FaceDef`] entries.
    pub faces: Vec<String>,
    /// Total force in newtons (positive; the sense is in `direction`).
    pub total_force_n: f64,
    /// Direction the force acts in.
    pub direction: LoadDirection,
}

/// The direction of a load.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LoadDirection {
    /// Along the inward face normal at each point (a pressure).
    FaceNormal,
    /// Along a signed coordinate axis.
    Axis(Axis),
    /// Along a given vector (non-zero; normalised by the solver).
    Vector([f64; 3]),
}

/// A signed coordinate axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Axis {
    /// +X.
    PlusX,
    /// −X.
    MinusX,
    /// +Y.
    PlusY,
    /// −Y.
    MinusY,
    /// +Z.
    PlusZ,
    /// −Z.
    MinusZ,
}

/// Why a setup is invalid.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum SetupError {
    /// The material name is empty.
    #[error("no material is given")]
    MissingMaterial,
    /// An infill pattern without homogenised properties.
    #[error("infill: homogenised properties are not available for `{0}` (only gyroid)")]
    UnsupportedInfill(String),
    /// An infill density outside (0, 1], or walls or line width not
    /// positive.
    #[error("infill: density must be in (0, 1], and walls_mm and line_width_mm positive")]
    BadInfill,
    /// Infill on a material that is not printed.
    #[error("infill: `{0}` is not a printed material")]
    InfillNotPrinted(String),
    /// A generative design's settings are out of range.
    #[error("generative: {0}")]
    BadGenerative(&'static str),
    /// A modal count outside 1..=[`MAX_MODES`].
    #[error("modal: {0} modes is outside 1..=50")]
    BadModeCount(usize),
    /// A thermal strain's temperatures are not valid.
    #[error("thermal: temperatures (reference {reference_c} °C) are not valid")]
    BadThermal {
        /// The reference given.
        reference_c: f64,
    },
    /// Two face definitions share a name.
    #[error("face name `{0}` is defined twice")]
    DuplicateFaceName(String),
    /// A support or load uses a face name with no definition.
    #[error("{used_by} uses face `{name}`, which is not defined")]
    UndefinedFace {
        /// Which entry uses it, e.g. `support 0` or `load 1`.
        used_by: String,
        /// The undefined name.
        name: String,
    },
    /// A support or load lists no faces.
    #[error("{0} lists no faces")]
    NoFaces(String),
    /// A direction vector is zero or not finite.
    #[error("{0} has a zero or non-finite direction")]
    ZeroDirection(String),
    /// A load's force is not a positive finite number.
    #[error("{used_by} has force {value} N; it must be positive and finite")]
    BadForce {
        /// Which load.
        used_by: String,
        /// The force given.
        value: f64,
    },
    /// A tolerance is not a positive finite number.
    #[error("{what} is {value}; it must be positive and finite")]
    NonPositiveTolerance {
        /// Which tolerance.
        what: String,
        /// The value given.
        value: f64,
    },
    /// A face reference's selector is malformed.
    #[error("face `{face}`: {what} is invalid ({value})")]
    BadSelector {
        /// The face definition.
        face: String,
        /// What is wrong.
        what: String,
        /// The offending value, when there is one.
        value: f64,
    },
    /// A face definition expects zero faces.
    #[error("face `{0}` expects zero faces")]
    ExpectsNothing(String),
    /// An option is out of range.
    #[error("option: {0}")]
    Option(#[from] OptionError),
    /// The material name does not resolve to a usable preset.
    #[error("material: {0}")]
    Material(#[from] MaterialError),
    /// A roller's direction is not the normal of a planar face it is on
    /// (only face-normal rollers are supported).
    #[error("{used_by}: roller direction is not the normal of planar face {face}")]
    RollerNotNormal {
        /// Which support.
        used_by: String,
        /// The face id.
        face: u32,
    },
}

fn nonzero_direction(v: &[f64; 3]) -> bool {
    v.iter().all(|c| c.is_finite()) && v.iter().any(|&c| c != 0.0)
}

fn positive_tolerance(what: &str, value: f64) -> Result<(), SetupError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(SetupError::NonPositiveTolerance {
            what: what.to_string(),
            value,
        })
    }
}

impl StructuralSetup {
    /// Checks the setup on its own, before any model is seen: the
    /// material resolves to a valid preset, names are defined once and
    /// used names exist, directions are non-zero,
    /// forces and tolerances are positive, and options are in range.
    pub fn validate(&self) -> Result<(), SetupError> {
        if self.material.trim().is_empty() {
            return Err(SetupError::MissingMaterial);
        }
        preset(&self.material)?.validate()?;
        positive_tolerance("angle tolerance", self.tolerance.angle_deg)?;
        positive_tolerance("distance tolerance", self.tolerance.distance)?;
        if let Some(d) = &self.build_direction
            && !nonzero_direction(d)
        {
            return Err(SetupError::ZeroDirection("build direction".to_string()));
        }
        for (i, def) in self.faces.iter().enumerate() {
            if self.faces[..i].iter().any(|d| d.name == def.name) {
                return Err(SetupError::DuplicateFaceName(def.name.clone()));
            }
            def.check()?;
        }
        for (i, support) in self.supports.iter().enumerate() {
            let used_by = format!("support {i}");
            self.check_names(&used_by, &support.faces)?;
            if let SupportKind::Roller { normal: v }
            | SupportKind::DirectionalRoller { direction: v } = &support.kind
                && !nonzero_direction(v)
            {
                return Err(SetupError::ZeroDirection(used_by));
            }
        }
        for (i, load) in self.loads.iter().enumerate() {
            let used_by = format!("load {i}");
            self.check_names(&used_by, &load.faces)?;
            if !(load.total_force_n.is_finite() && load.total_force_n > 0.0) {
                return Err(SetupError::BadForce {
                    used_by,
                    value: load.total_force_n,
                });
            }
            if let LoadDirection::Vector(v) = &load.direction
                && !nonzero_direction(v)
            {
                return Err(SetupError::ZeroDirection(used_by));
            }
        }
        if let Some(infill) = &self.infill {
            if infill.pattern != "gyroid" {
                return Err(SetupError::UnsupportedInfill(infill.pattern.clone()));
            }
            let positive = |x: f64| x.is_finite() && x > 0.0;
            if !(positive(infill.density)
                && infill.density <= 1.0
                && positive(infill.walls_mm)
                && positive(infill.line_width_mm))
            {
                return Err(SetupError::BadInfill);
            }
            if preset(&self.material)?
                .strength
                .interlayer
                .value()
                .is_none()
            {
                return Err(SetupError::InfillNotPrinted(self.material.clone()));
            }
        }
        if let Some(g) = &self.generative {
            let positive = |x: f64| x.is_finite() && x > 0.0;
            if !(positive(g.volume_fraction) && g.volume_fraction <= 1.0) {
                return Err(SetupError::BadGenerative(
                    "volume_fraction must be in (0, 1]",
                ));
            }
            if !positive(g.filter_mm) {
                return Err(SetupError::BadGenerative("filter_mm must be positive"));
            }
            if g.keep_mm.is_some_and(|k| !(k.is_finite() && k >= 0.0)) {
                return Err(SetupError::BadGenerative("keep_mm must be at least 0"));
            }
            if g.max_iterations == 0 {
                return Err(SetupError::BadGenerative(
                    "max_iterations must be at least 1",
                ));
            }
            if self.modal.is_some() || self.thermal.is_some() || self.infill.is_some() {
                return Err(SetupError::BadGenerative(
                    "a generative design takes no [modal], [thermal] or [infill]",
                ));
            }
            if !g.keep.is_empty() {
                self.check_names("generative keep", &g.keep)?;
            }
        }
        if let Some(count) = self.modal
            && !(1..=MAX_MODES).contains(&count)
        {
            return Err(SetupError::BadModeCount(count));
        }
        if let Some(t) = &self.thermal {
            let valid = |c: f64| c.is_finite() && c > crate::thermal_setup::ABSOLUTE_ZERO_C;
            let uniform = match t.temperature {
                TemperatureSource::Uniform(c) => c,
                TemperatureSource::ThermalSetup(_) => t.reference_c,
            };
            if !(valid(t.reference_c) && valid(uniform)) {
                return Err(SetupError::BadThermal {
                    reference_c: t.reference_c,
                });
            }
        }
        self.options.validate()?;
        Ok(())
    }

    fn check_names(&self, used_by: &str, names: &[String]) -> Result<(), SetupError> {
        if names.is_empty() {
            return Err(SetupError::NoFaces(used_by.to_string()));
        }
        match names
            .iter()
            .find(|n| !self.faces.iter().any(|d| &d.name == *n))
        {
            Some(name) => Err(SetupError::UndefinedFace {
                used_by: used_by.to_string(),
                name: name.clone(),
            }),
            None => Ok(()),
        }
    }
}

impl FaceDef {
    fn check(&self) -> Result<(), SetupError> {
        if matches!(self.expect, Expect::Exactly(0)) {
            return Err(SetupError::ExpectsNothing(self.name.clone()));
        }
        if let FaceRef::Select(selector) = &self.reference {
            selector.check().map_err(|d| SetupError::BadSelector {
                face: self.name.clone(),
                what: d.what.to_string(),
                value: d.value,
            })?;
        }
        Ok(())
    }
}
