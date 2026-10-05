//! `.pkn` assembly documents (design §6.2) and the one opener for
//! either kind of document.
//!
//! A version-5 envelope carries `tree` (a part document, read by
//! `popkern-history`) or `assembly`:
//!
//! ```json
//! { "version": 5, "naming_scheme": 2,
//!   "assembly": { "name": "gearbox",
//!     "definitions": [ { "name": "planet", "part": <FeatureTree> },
//!                      { "name": "stage", "instances": [ … ] } ],
//!     "root": [ { "name": "planet_1", "definition": 0,
//!                 "placement": [[…]], "visible": true } ] } }
//! ```
//!
//! Each part definition's tree is exactly the `tree` payload of a part
//! document, read through the same loader (index rebuild, label
//! normalisation), so a part moved into an assembly keeps its meaning.

use popkern_core::naming::{NAMING_SCHEME_VERSION, OLDEST_NAMING_SCHEME};
use popkern_history::{CURRENT_VERSION, FeatureTree};
use popkern_script::Model;
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::assembly::{Assembly, Definition, DefinitionContent, DefinitionId, Instance, Parent};
use crate::error::AssemblyError;
use crate::family::FamilyVariant;
use crate::names::validate_name;

/// The first `.pkn` version that can hold an assembly.
pub const FIRST_ASSEMBLY_VERSION: u32 = 5;

/// An opened `.pkn`: one part document or one assembly.
#[derive(Clone)]
pub enum Document {
    /// A single feature tree.
    Part(Box<Model>),
    /// Part definitions and placed instances.
    Assembly(Box<Assembly>),
}

impl Document {
    /// Open either kind of document from `.pkn` bytes.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Document`] for malformed JSON, an unsupported
    /// version or naming scheme, or an assembly whose structure is
    /// invalid (unknown definition, duplicate name, cycle);
    /// [`AssemblyError::Script`] when a part document cannot be read.
    pub fn from_pkn_bytes(data: &[u8]) -> Result<Document, AssemblyError> {
        let envelope: EnvelopeIn =
            serde_json::from_slice(data).map_err(|e| AssemblyError::Document(e.to_string()))?;
        match envelope.assembly {
            None => Model::from_pkn_bytes(data)
                .map(|model| Document::Part(Box::new(model)))
                .map_err(|source| AssemblyError::Script {
                    definition: "<document>".into(),
                    source,
                }),
            Some(assembly) => {
                if !(FIRST_ASSEMBLY_VERSION..=CURRENT_VERSION).contains(&envelope.version) {
                    return Err(AssemblyError::Document(format!(
                        "assembly document version {} not supported (reads {FIRST_ASSEMBLY_VERSION}..={CURRENT_VERSION})",
                        envelope.version
                    )));
                }
                // Each part tree carries its own scheme (PK-677); the
                // envelope bounds them.
                if !(OLDEST_NAMING_SCHEME..=NAMING_SCHEME_VERSION).contains(&envelope.naming_scheme)
                {
                    return Err(AssemblyError::Document(format!(
                        "naming scheme {} not supported (reads {OLDEST_NAMING_SCHEME}..={NAMING_SCHEME_VERSION})",
                        envelope.naming_scheme
                    )));
                }
                assembly
                    .into_assembly()
                    .map(|a| Document::Assembly(Box::new(a)))
            }
        }
    }

    /// Open a `.pkn` file of either kind.
    ///
    /// # Errors
    ///
    /// As [`Document::from_pkn_bytes`], and [`AssemblyError::Document`]
    /// when the file cannot be read.
    pub fn open(path: impl AsRef<std::path::Path>) -> Result<Document, AssemblyError> {
        let path = path.as_ref();
        let data = std::fs::read(path)
            .map_err(|e| AssemblyError::Document(format!("{}: {e}", path.display())))?;
        Document::from_pkn_bytes(&data)
    }
}

impl Assembly {
    /// The assembly as `.pkn` bytes.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Document`] when serialization fails.
    pub fn to_pkn_bytes(&self) -> Result<Vec<u8>, AssemblyError> {
        let out = EnvelopeOut {
            version: CURRENT_VERSION,
            naming_scheme: NAMING_SCHEME_VERSION,
            assembly: AssemblyOut {
                name: &self.name,
                definitions: self
                    .definitions
                    .iter()
                    .map(|d| DefinitionOut {
                        name: &d.name,
                        family: d.family.as_ref(),
                        part: d.model().map(Model::tree),
                        instances: match &d.content {
                            DefinitionContent::Assembly(instances) => Some(instances),
                            DefinitionContent::Part(_) => None,
                        },
                    })
                    .collect(),
                root: &self.root,
                mates: &self.mates,
                ground: self.ground.as_deref(),
                params: &self.params,
                drives: &self.drives,
            },
        };
        serde_json::to_vec(&out).map_err(|e| AssemblyError::Document(e.to_string()))
    }

    /// Save to `path`.
    ///
    /// # Errors
    ///
    /// As [`Assembly::to_pkn_bytes`], and [`AssemblyError::Document`]
    /// when the file cannot be written.
    pub fn save(&self, path: impl AsRef<std::path::Path>) -> Result<(), AssemblyError> {
        let path = path.as_ref();
        std::fs::write(path, self.to_pkn_bytes()?)
            .map_err(|e| AssemblyError::Document(format!("{}: {e}", path.display())))
    }
}

#[derive(Serialize)]
struct EnvelopeOut<'a> {
    version: u32,
    naming_scheme: u32,
    assembly: AssemblyOut<'a>,
}

#[derive(Serialize)]
struct AssemblyOut<'a> {
    name: &'a str,
    definitions: Vec<DefinitionOut<'a>>,
    root: &'a [Instance],
    #[serde(skip_serializing_if = "<[crate::mates::Mate]>::is_empty")]
    mates: &'a [crate::mates::Mate],
    #[serde(skip_serializing_if = "Option::is_none")]
    ground: Option<&'a str>,
    #[serde(skip_serializing_if = "<[crate::params::AssemblyParam]>::is_empty")]
    params: &'a [crate::params::AssemblyParam],
    #[serde(skip_serializing_if = "<[crate::params::Drive]>::is_empty")]
    drives: &'a [crate::params::Drive],
}

#[derive(Serialize)]
struct DefinitionOut<'a> {
    name: &'a str,
    #[serde(skip_serializing_if = "Option::is_none")]
    family: Option<&'a FamilyVariant>,
    #[serde(skip_serializing_if = "Option::is_none")]
    part: Option<&'a FeatureTree>,
    #[serde(skip_serializing_if = "Option::is_none")]
    instances: Option<&'a Vec<Instance>>,
}

#[derive(Deserialize)]
struct EnvelopeIn {
    version: u32,
    #[serde(default = "first_naming_scheme")]
    naming_scheme: u32,
    #[serde(default)]
    assembly: Option<AssemblyIn>,
}

fn first_naming_scheme() -> u32 {
    1
}

#[derive(Deserialize)]
struct AssemblyIn {
    name: String,
    definitions: Vec<DefinitionIn>,
    root: Vec<Instance>,
    #[serde(default)]
    mates: Vec<crate::mates::Mate>,
    #[serde(default)]
    ground: Option<String>,
    #[serde(default)]
    params: Vec<crate::params::AssemblyParam>,
    #[serde(default)]
    drives: Vec<crate::params::Drive>,
}

#[derive(Deserialize)]
struct DefinitionIn {
    name: String,
    #[serde(default)]
    family: Option<FamilyVariant>,
    #[serde(default)]
    part: Option<Value>,
    #[serde(default)]
    instances: Option<Vec<Instance>>,
}

impl AssemblyIn {
    /// Rebuild through the public construction path, so a loaded
    /// document meets exactly the invariants a scripted one does
    /// (names, uniqueness, acyclicity, known definitions).
    fn into_assembly(self) -> Result<Assembly, AssemblyError> {
        let mut assembly = Assembly::new(&self.name);
        let mut pending: Vec<(DefinitionId, Vec<Instance>)> = Vec::new();
        for definition in self.definitions {
            let id = match (definition.part, definition.instances) {
                (Some(tree), None) => {
                    let bytes = serde_json::to_vec(&serde_json::json!({
                        "version": CURRENT_VERSION,
                        "naming_scheme": NAMING_SCHEME_VERSION,
                        "tree": tree,
                    }))
                    .map_err(|e| AssemblyError::Document(e.to_string()))?;
                    let model =
                        Model::from_pkn_bytes(&bytes).map_err(|source| AssemblyError::Script {
                            definition: definition.name.clone(),
                            source,
                        })?;
                    assembly.add_loaded_part(&definition.name, model, definition.family)?
                }
                (None, Some(instances)) => {
                    let id = assembly.add_subassembly(&definition.name)?;
                    pending.push((id, instances));
                    id
                }
                _ => {
                    return Err(AssemblyError::Document(format!(
                        "definition '{}' must hold exactly one of `part` or `instances`",
                        definition.name
                    )));
                }
            };
            debug_assert_eq!(id.0 as usize + 1, assembly.definitions.len());
        }
        for (parent, instances) in pending {
            for instance in instances {
                place_loaded(&mut assembly, Parent::Definition(parent), instance)?;
            }
        }
        for instance in self.root {
            place_loaded(&mut assembly, Parent::Root, instance)?;
        }
        if let Some(ground) = &self.ground {
            assembly.set_ground(ground)?;
        }
        // Stored placements are the solved ones: the mates are kept, not
        // re-solved on open.
        assembly.mates = self.mates;
        // Stored part parameters already carry the driven values.
        assembly.params = self.params;
        assembly.drives = self.drives;
        assembly.mark_history_unrecorded();
        Ok(assembly)
    }
}

fn place_loaded(
    assembly: &mut Assembly,
    parent: Parent,
    instance: Instance,
) -> Result<(), AssemblyError> {
    validate_name(&instance.name, "instance")?;
    let at = assembly.place(
        parent,
        instance.definition,
        Some(&instance.name),
        instance.placement,
    )?;
    assembly.set_visible(at, instance.visible)
}

impl Assembly {
    /// A loaded part keeps its stored family and is not re-checked for
    /// alive bodies (the document has not been replayed yet).
    fn add_loaded_part(
        &mut self,
        name: &str,
        model: Model,
        family: Option<FamilyVariant>,
    ) -> Result<DefinitionId, AssemblyError> {
        validate_name(name, "part")?;
        if self.definition_named(name).is_some() {
            return Err(AssemblyError::DuplicateDefinitionName {
                name: name.to_string(),
            });
        }
        self.definitions.push(Definition {
            name: name.to_string(),
            family,
            content: DefinitionContent::Part(Box::new(model)),
        });
        Ok(DefinitionId(self.definitions.len() as u32 - 1))
    }
}
