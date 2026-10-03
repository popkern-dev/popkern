//! The one error type of the assembly layer.

use popkern_script::ScriptError;

/// Why an assembly operation was refused. Every variant carries the
/// names a user needs (instance, definition, parent) so a GUI can say
/// exactly what to change.
#[derive(Debug, thiserror::Error)]
pub enum AssemblyError {
    /// A sibling of the same parent already has this name.
    #[error("instance name '{name}' is already used under '{parent}'")]
    DuplicateInstanceName {
        /// The parent (`<root>` or a sub-assembly definition name).
        parent: String,
        /// The duplicated name.
        name: String,
    },
    /// Another definition already has this name.
    #[error("a part or sub-assembly named '{name}' already exists in this assembly")]
    DuplicateDefinitionName {
        /// The duplicated name.
        name: String,
    },
    /// A name is empty or contains a path separator.
    #[error("'{name}' is not a valid {what} name: {reason}")]
    InvalidName {
        /// The refused name.
        name: String,
        /// `instance`, `part` or `assembly`.
        what: &'static str,
        /// Why it was refused.
        reason: &'static str,
    },
    /// A definition id that this assembly does not hold.
    #[error("definition #{id} does not exist in this assembly")]
    UnknownDefinition {
        /// The id.
        id: u32,
    },
    /// No instance at this path.
    #[error("no instance at '{path}'")]
    UnknownInstance {
        /// The path asked for.
        path: String,
    },
    /// The placement is not a rigid motion.
    #[error(
        "the placement of '{instance}' is not a rigid motion (orthonormality defect {defect:.3e}{mirror}): \
         instances move and turn, they never scale, shear or mirror"
    )]
    NonRigidPlacement {
        /// The instance being placed (or `placement` when not yet named).
        instance: String,
        /// Largest deviation of `RᵀR` from the identity.
        defect: f64,
        /// `", mirrored"` when the determinant is negative.
        mirror: &'static str,
    },
    /// Placing this definition here would make it contain itself.
    #[error(
        "placing '{definition}' inside '{parent}' would make an assembly contain itself ({chain})"
    )]
    Cycle {
        /// The definition being placed.
        definition: String,
        /// The sub-assembly it was placed into.
        parent: String,
        /// The containment chain, `a > b > a`.
        chain: String,
    },
    /// A part definition with no alive body.
    #[error("part '{name}' has no body: its script recorded nothing that builds a solid")]
    EmptyDefinition {
        /// The definition name.
        name: String,
    },
    /// An operation that needs a part was given a sub-assembly, or the
    /// reverse.
    #[error("'{name}' is {actual}, not {expected}")]
    WrongDefinitionKind {
        /// The definition name.
        name: String,
        /// What the operation needs.
        expected: &'static str,
        /// What it got.
        actual: &'static str,
    },
    /// A part definition's recorder failed (replay, export, …).
    #[error("in part '{definition}': {source}")]
    Script {
        /// The definition whose model failed.
        definition: String,
        /// The recorder's error.
        #[source]
        source: ScriptError,
    },
    /// A kernel geometry operation failed (a placement could not be
    /// applied to a body).
    #[error("placing '{instance}': {detail}")]
    Geometry {
        /// The instance being placed.
        instance: String,
        /// The kernel's message.
        detail: String,
    },
    /// STEP writing or reading failed.
    #[error("STEP: {0}")]
    Step(String),
    /// The `.pkn` document is malformed or of an unsupported version.
    #[error("assembly document: {0}")]
    Document(String),
}
