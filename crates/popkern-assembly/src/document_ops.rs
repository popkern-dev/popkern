//! What a host does with an authored document of either kind: save it,
//! export it, view it whole, and set parameters. Hosts (the CLI,
//! popkern live, the bake tool) hold a [`Document`] and never branch
//! on its kind for these.
//!
//! Parameters of an assembly are qualified by their part definition:
//! `planet:width` sets `width` in the `planet` definition (and so in
//! every instance of it). A part document takes the bare name.

use popkern_part::Part;
use popkern_script::Model;

use crate::assembly::{Assembly, DefinitionContent};
use crate::document::Document;
use crate::error::AssemblyError;

impl Document {
    /// The document as `.pkn` bytes.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Script`] / [`AssemblyError::Document`] when
    /// serialization fails.
    pub fn to_pkn_bytes(&self) -> Result<Vec<u8>, AssemblyError> {
        match self {
            Document::Part(model) => model
                .to_pkn_bytes()
                .map_err(|source| AssemblyError::Script {
                    definition: "<document>".into(),
                    source,
                }),
            Document::Assembly(assembly) => assembly.to_pkn_bytes(),
        }
    }

    /// Save to `path`.
    ///
    /// # Errors
    ///
    /// As [`Document::to_pkn_bytes`], and [`AssemblyError::Document`]
    /// when the file cannot be written.
    pub fn save(&self, path: impl AsRef<std::path::Path>) -> Result<(), AssemblyError> {
        let path = path.as_ref();
        std::fs::write(path, self.to_pkn_bytes()?)
            .map_err(|e| AssemblyError::Document(format!("{}: {e}", path.display())))
    }

    /// STEP Part 21 bytes: a part document as one product holding every
    /// body's solid, an assembly as its product structure.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Script`] / [`AssemblyError::Step`].
    pub fn to_step_bytes(&mut self, name: &str) -> Result<Vec<u8>, AssemblyError> {
        match self {
            Document::Part(model) => {
                model
                    .to_step_bytes(name)
                    .map_err(|source| AssemblyError::Script {
                        definition: "<document>".into(),
                        source,
                    })
            }
            Document::Assembly(assembly) => assembly.to_step_bytes(name),
        }
    }

    /// Write [`Document::to_step_bytes`] to `path` (the file stem names
    /// the product).
    ///
    /// # Errors
    ///
    /// As [`Document::to_step_bytes`], and [`AssemblyError::Step`] when
    /// the file cannot be written.
    pub fn to_step(&mut self, path: impl AsRef<std::path::Path>) -> Result<(), AssemblyError> {
        let path = path.as_ref();
        let name = path
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| "document".into());
        let bytes = self.to_step_bytes(&name)?;
        std::fs::write(path, bytes)
            .map_err(|e| AssemblyError::Step(format!("{}: {e}", path.display())))
    }

    /// The whole document as one part (bodies concatenated, no boolean
    /// semantics): every alive body of a part document, every visible
    /// part occurrence of an assembly at its world placement. The
    /// render view.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Script`] / [`AssemblyError::Geometry`]; an
    /// assembly with no visible part is [`AssemblyError::EmptyDefinition`].
    pub fn to_part_all(&mut self) -> Result<Part, AssemblyError> {
        match self {
            Document::Part(model) => model.to_part_all().map_err(|source| AssemblyError::Script {
                definition: "<document>".into(),
                source,
            }),
            Document::Assembly(assembly) => assembly.to_part_all(),
        }
    }

    /// Every recorded model with its part definition's name (`None`
    /// for a part document), in definition order.
    pub fn models_mut(&mut self) -> Vec<(Option<String>, &mut Model)> {
        match self {
            Document::Part(model) => vec![(None, model.as_mut())],
            Document::Assembly(assembly) => assembly.part_models_mut(),
        }
    }

    /// Set a parameter: the bare name on a part document,
    /// `<part>:<name>` on an assembly.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Script`] for an unknown parameter or a failed
    /// set; [`AssemblyError::UnknownInstance`]-free: an unknown part is
    /// [`AssemblyError::Document`] naming the parts.
    pub fn set_param(&mut self, qualified: &str, value: f64) -> Result<(), AssemblyError> {
        let script = |definition: &str| {
            let definition = definition.to_string();
            move |source| AssemblyError::Script { definition, source }
        };
        match self {
            Document::Part(model) => model
                .set_param(qualified, value)
                .map_err(script("<document>")),
            Document::Assembly(assembly) => {
                let Some((part, name)) = qualified.split_once(':') else {
                    // Unqualified: an assembly parameter (PK-638),
                    // which replays every part it drives.
                    assembly.set_assembly_param(qualified, value)?;
                    if !assembly.mates().is_empty() {
                        assembly.solve_mates()?;
                    }
                    return Ok(());
                };
                if let Some(param) = assembly.driven_by(part, name) {
                    return Err(AssemblyError::Param {
                        name: qualified.to_string(),
                        detail: format!(
                            "driven by the assembly parameter '{param}': set that instead"
                        ),
                    });
                }
                let Some(id) = assembly.definition_named(part) else {
                    return Err(AssemblyError::Document(format!(
                        "no part named '{part}' (parts: {})",
                        assembly.part_names().join(", ")
                    )));
                };
                assembly
                    .part_model_mut(id)?
                    .set_param(name, value)
                    .map_err(script(part))
            }
        }
    }
}

impl Assembly {
    /// Every part definition's model with its name, in definition order.
    pub fn part_models_mut(&mut self) -> Vec<(Option<String>, &mut Model)> {
        self.definitions
            .iter_mut()
            .filter_map(|d| match &mut d.content {
                DefinitionContent::Part(model) => Some((Some(d.name.clone()), model.as_mut())),
                DefinitionContent::Assembly(_) => None,
            })
            .collect()
    }

    /// The part definitions' names, in definition order.
    pub fn part_names(&self) -> Vec<String> {
        self.definitions
            .iter()
            .filter(|d| matches!(d.content, DefinitionContent::Part(_)))
            .map(|d| d.name.clone())
            .collect()
    }

    /// Every visible part occurrence at its world placement, as one part
    /// (bodies concatenated): the render view of the assembly.
    ///
    /// # Errors
    ///
    /// As [`Assembly::flatten`]; [`AssemblyError::EmptyDefinition`]
    /// when nothing is visible.
    pub fn to_part_all(&mut self) -> Result<Part, AssemblyError> {
        let tolerance = self
            .part_models_mut()
            .into_iter()
            .map(|(_, m)| m.tolerance())
            .max_by(|a, b| a.linear().total_cmp(&b.linear()))
            .unwrap_or_default();
        let mut whole: Option<Part> = None;
        for body in self.flatten()?.into_iter().filter(|b| b.visible) {
            let part = Part::from_brep(body.brep, tolerance, "popkern-assembly").map_err(|e| {
                AssemblyError::Geometry {
                    instance: body.path.to_string(),
                    detail: e.to_string(),
                }
            })?;
            whole = Some(match whole {
                None => part,
                Some(acc) => acc
                    .concatenate(&part)
                    .map_err(|e| AssemblyError::Geometry {
                        instance: body.path.to_string(),
                        detail: e.to_string(),
                    })?,
            });
        }
        whole.ok_or_else(|| AssemblyError::EmptyDefinition {
            name: self.name.clone(),
        })
    }
}
