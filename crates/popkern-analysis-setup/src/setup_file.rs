//! The setup file (TOML): its on-disk shape, parsing and writing.
//!
//! The file shape mirrors `docs/architecture/analysis-library-boundary.md`
//! ("Example `setup.toml`"). Every table rejects unknown keys. Parsing
//! converts into `popkern_analysis::StructuralSetup` and then runs its
//! format-independent validation, so a file that parses is a valid
//! setup.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

#[path = "setup_convert.rs"]
mod convert;
#[path = "setup_error.rs"]
mod error;
#[path = "setup_shared.rs"]
pub(crate) mod shared;

pub use convert::SetupFileError;
pub(crate) use convert::options as options_of;

use popkern_analysis::StructuralSetup;

/// The whole setup file.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SetupFile {
    /// `"structural"` when given; a thermal file is refused.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub analysis: Option<String>,
    /// Length unit; only `"mm"` is supported.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub units: Option<String>,
    /// Material preset name.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub material: Option<String>,
    /// Optional build direction.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub build_direction: Option<[f64; 3]>,
    /// Default selector tolerances.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub tolerance: Option<ToleranceFile>,
    /// Named face references (sorted by name).
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub faces: BTreeMap<String, FaceEntry>,
    /// `[[support]]` entries.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub support: Vec<SupportFile>,
    /// `[[load]]` entries.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub load: Vec<LoadFile>,
    /// `[thermal]`: a thermal strain.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub thermal: Option<ThermalLoadFile>,
    /// `[modal]`: natural frequencies instead of a static solve.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub modal: Option<ModalFile>,
    /// `[infill]`: solid walls around a homogenised core.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub infill: Option<InfillFile>,
    /// `[generative]`: the stiffest part at a volume fraction instead of
    /// a static solve.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub generative: Option<GenerativeFile>,
    /// Option overrides; unset fields keep their defaults.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub options: Option<OptionsFile>,
}

/// `[thermal]`: the stress-free reference and exactly one of a uniform
/// temperature or a thermal setup (a path relative to this file) solved
/// on the same part first.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ThermalLoadFile {
    /// Stress-free temperature, °C.
    pub reference_c: f64,
    /// Uniform temperature, °C.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub uniform_c: Option<f64>,
    /// A thermal setup file.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub setup: Option<String>,
}

/// `[infill]`.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct InfillFile {
    /// The slicer's pattern (`"gyroid"`).
    pub pattern: String,
    /// Relative density of the core, in (0, 1].
    pub density: f64,
    /// Solid wall thickness, mm.
    pub walls_mm: f64,
    /// The gyroid's bead width, mm (default 0.4, a 0.4 mm nozzle).
    #[serde(default = "default_line_width")]
    pub line_width_mm: f64,
}

/// A 0.4 mm nozzle's line, the common default.
fn default_line_width() -> f64 {
    0.4
}

/// `[generative]`.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct GenerativeFile {
    /// Material to keep, as a fraction of the part's volume, in (0, 1].
    pub volume_fraction: f64,
    /// Density filter radius, mm.
    pub filter_mm: f64,
    /// Named faces whose nearby material is kept.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub keep: Vec<String>,
    /// How deep kept material reaches, mm (default: `filter_mm`).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub keep_mm: Option<f64>,
    /// Most updates (default 300).
    #[serde(default = "default_max_iterations")]
    pub max_iterations: usize,
    /// Push the design to solid and void (default true).
    #[serde(default = "default_crisp")]
    pub crisp: bool,
    /// Print without supports along this axis: "+x", "-x", "+y", "-y",
    /// "+z" or "-z" (the direction the part grows from the bed).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub build_direction: Option<String>,
    /// Keep the design's thinnest members about twice as thick, at the
    /// cost of stiffness (the robust formulation; needs `crisp`;
    /// default false).
    #[serde(default, skip_serializing_if = "std::ops::Not::not")]
    pub robust: bool,
}

/// Enough for a crisp design's projection to reach its sharpest
/// (measured: 200–400 updates on a cantilever, 60–120 without).
fn default_max_iterations() -> usize {
    300
}

/// Crisp designs mesh to the volume asked for; grey ones do not.
fn default_crisp() -> bool {
    true
}

/// `[modal]`.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ModalFile {
    /// How many of the lowest modes.
    pub modes: usize,
}

/// `[tolerance]`: angle in degrees, distance in model units.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct ToleranceFile {
    /// Angle tolerance in degrees.
    pub angle_deg: f64,
    /// Distance tolerance.
    pub distance: f64,
}

/// One `[faces]` entry: exactly one of the reference keys, plus an
/// optional `expect`.
#[derive(Debug, Clone, PartialEq, Default, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct FaceEntry {
    /// Plane selector.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub plane: Option<PlaneFile>,
    /// Cylinder selector.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cylinder: Option<CylinderFile>,
    /// Sphere selector.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sphere: Option<SphereFile>,
    /// Nearest-face-to-point selector.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub nearest_to: Option<[f64; 3]>,
    /// Intersection of selectors.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub all_of: Option<Vec<SelectorFile>>,
    /// Union of selectors.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub any_of: Option<Vec<SelectorFile>>,
    /// Stable popkern face name.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
    /// Raw face id (fragile).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub id: Option<u32>,
    /// How many faces must match; default exactly one.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub expect: Option<ExpectFile>,
}

/// A selector nested inside `all_of` / `any_of`: exactly one key.
#[derive(Debug, Clone, PartialEq, Default, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SelectorFile {
    /// Plane selector.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub plane: Option<PlaneFile>,
    /// Cylinder selector.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cylinder: Option<CylinderFile>,
    /// Sphere selector.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sphere: Option<SphereFile>,
    /// Nearest-face-to-point selector.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub nearest_to: Option<[f64; 3]>,
    /// Intersection of selectors.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub all_of: Option<Vec<SelectorFile>>,
    /// Union of selectors.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub any_of: Option<Vec<SelectorFile>>,
}

/// `plane = { normal, through, angle_tol?, dist_tol? }`.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct PlaneFile {
    /// Outward normal.
    pub normal: [f64; 3],
    /// A point on the plane.
    pub through: [f64; 3],
    /// Angle tolerance in degrees.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub angle_tol: Option<f64>,
    /// Distance tolerance.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub dist_tol: Option<f64>,
}

/// `cylinder = { axis, radius, through?, angle_tol?, dist_tol? }`.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct CylinderFile {
    /// Axis direction.
    pub axis: [f64; 3],
    /// Radius.
    pub radius: f64,
    /// A point on the axis.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub through: Option<[f64; 3]>,
    /// Angle tolerance in degrees.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub angle_tol: Option<f64>,
    /// Distance tolerance.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub dist_tol: Option<f64>,
}

/// `sphere = { centre, radius, dist_tol? }`.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SphereFile {
    /// Centre.
    pub centre: [f64; 3],
    /// Radius.
    pub radius: f64,
    /// Distance tolerance.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub dist_tol: Option<f64>,
}

/// `expect = 2` (exactly) or `expect = { at_least = 2 }`.
#[derive(Debug, Clone, Copy, PartialEq, Deserialize, Serialize)]
#[serde(untagged)]
pub enum ExpectFile {
    /// Exactly this many faces.
    Exactly(u32),
    /// At least this many faces.
    AtLeast {
        /// Minimum count.
        at_least: u32,
    },
}

/// A `[[support]]` entry.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SupportFile {
    /// Face names from `[faces]`.
    pub faces: Vec<String>,
    /// `"fixed"` or `"roller"`.
    pub kind: SupportKindFile,
    /// A roller's direction, asserted to be the normal of every (planar)
    /// face it is on; refused on a fixed support.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub normal: Option<[f64; 3]>,
    /// A roller's constrained direction on any face (motion normal to
    /// it is free). A roller with neither key follows each face's own
    /// normal (radial on a cylinder).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub direction: Option<[f64; 3]>,
}

/// Support kind keyword.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum SupportKindFile {
    /// All displacement components fixed.
    Fixed,
    /// Displacement along one direction fixed.
    Roller,
}

/// A `[[load]]` entry.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct LoadFile {
    /// Face names from `[faces]`.
    pub faces: Vec<String>,
    /// Total force in newtons.
    pub force_n: f64,
    /// A vector, `"face_normal"`, or a signed axis such as `"-z"`.
    pub direction: DirectionFile,
}

/// Load direction: a vector or a keyword.
#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
#[serde(untagged)]
pub enum DirectionFile {
    /// A direction vector.
    Vector([f64; 3]),
    /// `"face_normal"` or `"+x"`, `"-x"`, …, `"-z"`.
    Keyword(String),
}

/// `[options]`: overrides of `AnalysisOptions` defaults.
#[derive(Debug, Clone, Copy, PartialEq, Default, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct OptionsFile {
    /// Grid cell size.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub target_cell_size: Option<f64>,
    /// Active cell budget.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_active_cells: Option<usize>,
    /// Sub-cell lattice depth.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub subcell_depth: Option<u8>,
    /// AgFEM threshold γ.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub aggregation_threshold: Option<f64>,
    /// Support penalty factor.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub support_penalty_factor: Option<f64>,
    /// Nitsche factor.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub nitsche_factor: Option<f64>,
    /// CG relative residual tolerance.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cg_rel_tol: Option<f64>,
    /// CG iteration budget.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cg_max_iters: Option<usize>,
    /// Grid origin shift, mm.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub grid_shift: Option<[f64; 3]>,
    /// Target coarsest-level cell count.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub coarsest_dofs: Option<usize>,
}

/// Parses and validates a setup file's text.
pub fn parse_setup(text: &str) -> Result<StructuralSetup, SetupFileError> {
    let file: SetupFile = toml::from_str(text).map_err(|e| SetupFileError::from_toml(&e))?;
    let setup = file.into_setup()?;
    setup.validate()?;
    Ok(setup)
}

/// Writes a setup as setup-file text that [`parse_setup`] reads back
/// to an equal setup.
pub fn write_setup(setup: &StructuralSetup) -> Result<String, SetupFileError> {
    toml::to_string(&SetupFile::from_setup(setup)).map_err(|e| SetupFileError::Write {
        detail: e.to_string(),
    })
}
