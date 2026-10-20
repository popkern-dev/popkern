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
    /// An assembly parameter or drive was refused.
    #[error("assembly parameter '{name}': {detail}")]
    Param {
        /// The parameter (assembly or part) concerned.
        name: String,
        /// Why.
        detail: String,
    },
    /// A mate refers to something it cannot constrain.
    #[error("joint '{joint}': the {role} entity is {found}")]
    JointEntity {
        /// The joint.
        joint: String,
        /// Which argument: `axis` or `zero`.
        role: &'static str,
        /// What it resolved to.
        found: String,
    },
    #[error("mate {mate}: {detail}")]
    Mate {
        /// The mate (or entity) in words.
        mate: String,
        /// Why.
        detail: String,
    },
    /// The mates cannot all hold; nothing moved.
    #[error(
        "the mates conflict: {} cannot hold together (largest residual {:.3e}); nothing was moved",
        conflict_list(report),
        report.residual
    )]
    MateConflict {
        /// The solve's report (conflicting and redundant mates, free motions).
        report: Box<crate::mate_solve::MateReport>,
    },
    /// A released joint's solved angle is outside its limits; nothing
    /// moved (PK-808).
    #[error(
        "joint '{joint}' would turn to {angle:.3} degrees, outside its limits [{}, {}]; nothing was moved",
        limits[0],
        limits[1]
    )]
    JointLimit {
        /// The joint.
        joint: String,
        /// The angle the solve reached, degrees.
        angle: f64,
        /// Its limits, degrees.
        limits: [f64; 2],
    },
    /// STEP writing or reading failed.
    #[error("STEP: {0}")]
    Step(String),
    /// The `.pkn` document is malformed or of an unsupported version.
    #[error("assembly document: {0}")]
    Document(String),
}

fn conflict_list(report: &crate::mate_solve::MateReport) -> String {
    let names: Vec<String> = report
        .conflicting
        .iter()
        .filter_map(|i| report.mates.get(*i).cloned())
        .collect();
    if names.is_empty() {
        "the mates".to_string()
    } else {
        names.join("; ")
    }
}

impl AssemblyError {
    /// The structured error payload every skin raises: a part's
    /// recorder failure keeps its own class and fields (its message
    /// prefixed with the part, a `part` field added); every other
    /// variant is `kind: "AssemblyError"` with a snake-case `reason`
    /// and the variant's fields.
    pub fn payload(&self) -> serde_json::Value {
        use serde_json::json;
        let message = self.to_string();
        let mut value = match self {
            AssemblyError::Script { definition, source } => {
                let mut inner = source.payload();
                if let Some(map) = inner.as_object_mut() {
                    map.insert("part".into(), json!(definition));
                }
                inner
            }
            AssemblyError::DuplicateInstanceName { parent, name } => {
                json!({ "reason": "duplicate_instance_name", "parent": parent, "name": name })
            }
            AssemblyError::DuplicateDefinitionName { name } => {
                json!({ "reason": "duplicate_definition_name", "name": name })
            }
            AssemblyError::InvalidName { name, what, reason } => {
                json!({ "reason": "invalid_name", "name": name, "what": what, "detail": reason })
            }
            AssemblyError::UnknownDefinition { id } => {
                json!({ "reason": "unknown_definition", "id": id })
            }
            AssemblyError::UnknownInstance { path } => {
                json!({ "reason": "unknown_instance", "path": path })
            }
            AssemblyError::NonRigidPlacement {
                instance,
                defect,
                mirror,
            } => json!({
                "reason": "non_rigid_placement", "instance": instance,
                "defect": if defect.is_finite() { json!(defect) } else { json!(null) },
                "mirrored": !mirror.is_empty(),
            }),
            AssemblyError::Cycle {
                definition,
                parent,
                chain,
            } => json!({
                "reason": "cycle", "definition": definition, "parent": parent, "chain": chain,
            }),
            AssemblyError::EmptyDefinition { name } => {
                json!({ "reason": "empty_definition", "name": name })
            }
            AssemblyError::WrongDefinitionKind {
                name,
                expected,
                actual,
            } => json!({
                "reason": "wrong_definition_kind", "name": name,
                "expected": expected, "actual": actual,
            }),
            AssemblyError::Geometry { instance, detail } => {
                json!({ "reason": "geometry", "instance": instance, "detail": detail })
            }
            AssemblyError::Param { name, detail } => {
                json!({ "reason": "param", "name": name, "detail": detail })
            }
            AssemblyError::JointEntity { joint, role, found } => {
                json!({ "reason": "joint_entity", "joint": joint, "role": role, "found": found })
            }
            AssemblyError::Mate { mate, detail } => {
                json!({ "reason": "mate", "mate": mate, "detail": detail })
            }
            AssemblyError::MateConflict { report } => {
                json!({ "reason": "mate_conflict", "report": report })
            }
            AssemblyError::JointLimit {
                joint,
                angle,
                limits,
            } => {
                json!({ "reason": "joint_limit", "joint": joint, "angle": angle, "limits": limits })
            }
            AssemblyError::Step(detail) => json!({ "reason": "step", "detail": detail }),
            AssemblyError::Document(detail) => json!({ "reason": "document", "detail": detail }),
        };
        if let Some(map) = value.as_object_mut() {
            map.entry("kind").or_insert_with(|| json!("AssemblyError"));
            map.insert("message".into(), json!(message));
        }
        value
    }
}
