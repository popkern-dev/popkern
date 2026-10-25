//! Why a setup file was refused: one typed error per thing the file can
//! get wrong.

use popkern_analysis::SetupError;

/// Why a setup file was refused.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum SetupFileError {
    /// A key the schema does not define.
    #[error("unknown key: {detail}")]
    UnknownField {
        /// The parser's message, with location.
        detail: String,
    },
    /// A required key is absent.
    #[error("missing key `{field}`")]
    MissingField {
        /// The key.
        field: String,
    },
    /// The text is not valid TOML or a value has the wrong type.
    #[error("setup file: {detail}")]
    Syntax {
        /// The parser's message, with location.
        detail: String,
    },
    /// A unit other than millimetres.
    #[error("units `{0}` are not supported; use \"mm\"")]
    UnsupportedUnits(String),
    /// A face entry gives no reference, or more than one.
    #[error("face `{face}` must give exactly one reference key, found {found}")]
    FaceReferenceCount {
        /// The face name (or `face/all_of` for nested selectors).
        face: String,
        /// Number of reference keys given.
        found: usize,
    },
    /// A load direction keyword that is not recognised.
    #[error("load {load}: unknown direction `{value}`")]
    UnknownDirection {
        /// Load index.
        load: usize,
        /// The keyword given.
        value: String,
    },
    /// A build direction that is not a signed axis.
    #[error(
        "generative: unknown build_direction `{0}`; use \"+x\", \"-x\", \"+y\", \"-y\", \"+z\" or \"-z\""
    )]
    UnknownBuildDirection(String),
    /// A roller given both a normal and a direction.
    #[error("support {0}: a roller takes `normal` or `direction`, not both")]
    RollerNormalAndDirection(usize),
    /// `[thermal]` gives both or neither of `uniform_c` and `setup`.
    #[error("thermal: give exactly one of `uniform_c` and `setup`")]
    ThermalSourceCount,
    /// The file is for a different analysis.
    #[error("analysis `{found}` given to the {expected} reader")]
    WrongAnalysis {
        /// The analysis the reader handles.
        expected: &'static str,
        /// The `analysis` key in the file.
        found: String,
    },
    /// A fixed support with a normal or direction (which would be
    /// ignored).
    #[error("support {0}: a fixed support takes no normal or direction")]
    NormalOnFixed(usize),
    /// The setup parsed but is not valid.
    #[error(transparent)]
    Invalid(#[from] SetupError),
    /// Writing the setup back to TOML failed.
    #[error("cannot write setup: {detail}")]
    Write {
        /// The serialiser's message.
        detail: String,
    },
}

impl SetupFileError {
    /// Classifies a TOML parse error. serde's messages for unknown and
    /// missing keys are stable (`unknown field`, `missing field`).
    pub(crate) fn from_toml(e: &toml::de::Error) -> Self {
        let message = e.message();
        if message.starts_with("unknown field") {
            SetupFileError::UnknownField {
                detail: e.to_string(),
            }
        } else if let Some(rest) = message.strip_prefix("missing field `") {
            SetupFileError::MissingField {
                field: rest.split('`').next().unwrap_or(rest).to_string(),
            }
        } else {
            SetupFileError::Syntax {
                detail: e.to_string(),
            }
        }
    }
}
