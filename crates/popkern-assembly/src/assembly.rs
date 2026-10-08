//! The assembly document: part definitions built once, and named
//! rigid instances of them, nestable through sub-assembly
//! definitions (design `docs/architecture/assemblies.md` §5).

use popkern_script::Model;
use serde::{Deserialize, Serialize};

use crate::error::AssemblyError;
use crate::names::{InstancePath, validate_name};
use crate::placement::Placement;

/// Index of a definition in its assembly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct DefinitionId(pub u32);

/// What a definition holds.
#[derive(Clone)]
pub enum DefinitionContent {
    /// A part: one recorded document; every alive body is the part.
    Part(Box<Model>),
    /// A sub-assembly: instances of other definitions.
    Assembly(Vec<Instance>),
}

/// A part or sub-assembly that instances place.
#[derive(Clone)]
pub struct Definition {
    /// Unique among the assembly's definitions; the STEP product name.
    pub name: String,
    /// The family this definition is a variant of (stage 2), when any.
    pub family: Option<crate::family::FamilyVariant>,
    /// The geometry (part) or the instances (sub-assembly).
    pub content: DefinitionContent,
}

impl Definition {
    /// `"a part"` or `"a sub-assembly"`.
    pub fn kind(&self) -> &'static str {
        match self.content {
            DefinitionContent::Part(_) => "a part",
            DefinitionContent::Assembly(_) => "a sub-assembly",
        }
    }

    /// The part's model, when this is a part.
    pub fn model(&self) -> Option<&Model> {
        match &self.content {
            DefinitionContent::Part(model) => Some(model),
            DefinitionContent::Assembly(_) => None,
        }
    }
}

/// One named, placed occurrence of a definition.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Instance {
    /// Unique among its siblings.
    pub name: String,
    /// What it places.
    pub definition: DefinitionId,
    /// Where it sits in its parent.
    pub placement: Placement,
    /// Shown in viewers (exports write every instance).
    #[serde(default = "visible_default")]
    pub visible: bool,
}

fn visible_default() -> bool {
    true
}

/// Where an instance is placed: the assembly root or a sub-assembly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Parent {
    /// The top level.
    Root,
    /// Inside a sub-assembly definition.
    Definition(DefinitionId),
}

/// A handle to one instance (its parent and sibling index): stable
/// because instances are never reordered.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct InstanceRef {
    /// The parent the instance was placed into.
    pub parent: Parent,
    /// Its index among the parent's instances.
    pub index: usize,
}

/// One occurrence in the expanded tree (sub-assembly occurrences
/// included), as a viewer's tree panel lists them.
#[derive(Debug, Clone, PartialEq)]
pub struct TreeRow {
    /// Path from the root.
    pub path: InstancePath,
    /// The placed definition.
    pub definition: DefinitionId,
    /// World placement (all parents applied).
    pub world: Placement,
    /// Shown: this instance and every ancestor visible.
    pub visible: bool,
    /// A leaf (part) occurrence rather than a sub-assembly.
    pub is_part: bool,
}

/// A product structure: definitions and the top-level instances.
#[derive(Clone)]
pub struct Assembly {
    pub(crate) name: String,
    pub(crate) definitions: Vec<Definition>,
    pub(crate) root: Vec<Instance>,
    /// Mates, in the order they were added (stage 3).
    pub(crate) mates: Vec<crate::mates::Mate>,
    /// The root instance the mate solve holds fixed (`None`: the first).
    pub(crate) ground: Option<String>,
    /// Assembly parameters (PK-638).
    pub(crate) params: Vec<crate::params::AssemblyParam>,
    /// Assembly parameters driving part parameters.
    pub(crate) drives: Vec<crate::params::Drive>,
    /// The construction steps (PK-666).
    pub(crate) history: crate::history::AssemblyHistory,
    /// Joints over the mates (six-axis-arm §4.3).
    pub(crate) joints: Vec<crate::joints::Joint>,
}

impl Assembly {
    /// An empty assembly called `name`.
    pub fn new(name: &str) -> Assembly {
        Assembly {
            name: name.to_string(),
            definitions: Vec::new(),
            root: Vec::new(),
            mates: Vec::new(),
            ground: None,
            params: Vec::new(),
            drives: Vec::new(),
            history: crate::history::AssemblyHistory::default(),
            joints: Vec::new(),
        }
    }

    /// The assembly's name (the STEP root product).
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Every definition, in insertion order (index = [`DefinitionId`]).
    pub fn definitions(&self) -> &[Definition] {
        &self.definitions
    }

    /// One definition.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::UnknownDefinition`].
    pub fn definition(&self, id: DefinitionId) -> Result<&Definition, AssemblyError> {
        self.definitions
            .get(id.0 as usize)
            .ok_or(AssemblyError::UnknownDefinition { id: id.0 })
    }

    /// The definition called `name`.
    pub fn definition_named(&self, name: &str) -> Option<DefinitionId> {
        self.definitions
            .iter()
            .position(|d| d.name == name)
            .map(|i| DefinitionId(i as u32))
    }

    /// A part definition's model, for edits (parameters, replay).
    ///
    /// # Errors
    ///
    /// [`AssemblyError::UnknownDefinition`];
    /// [`AssemblyError::WrongDefinitionKind`] for a sub-assembly.
    pub fn part_model_mut(&mut self, id: DefinitionId) -> Result<&mut Model, AssemblyError> {
        let definition = self
            .definitions
            .get_mut(id.0 as usize)
            .ok_or(AssemblyError::UnknownDefinition { id: id.0 })?;
        match &mut definition.content {
            DefinitionContent::Part(model) => Ok(model),
            DefinitionContent::Assembly(_) => Err(AssemblyError::WrongDefinitionKind {
                name: definition.name.clone(),
                expected: "a part",
                actual: "a sub-assembly",
            }),
        }
    }

    /// Add a part definition built by `model`.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::InvalidName`], [`AssemblyError::DuplicateDefinitionName`],
    /// [`AssemblyError::EmptyDefinition`] when the model has no alive body.
    pub fn add_part(&mut self, name: &str, model: Model) -> Result<DefinitionId, AssemblyError> {
        self.add_part_with(name, model, None)
    }

    /// Add a part definition that is a variant of a family (the
    /// family's name and full parameter set are recorded with it).
    ///
    /// # Errors
    ///
    /// As [`Assembly::add_part`].
    pub fn add_family_variant(
        &mut self,
        name: &str,
        model: Model,
        variant: crate::family::FamilyVariant,
    ) -> Result<DefinitionId, AssemblyError> {
        self.add_part_with(name, model, Some(variant))
    }

    fn add_part_with(
        &mut self,
        name: &str,
        mut model: Model,
        family: Option<crate::family::FamilyVariant>,
    ) -> Result<DefinitionId, AssemblyError> {
        self.check_new_definition(name, "part")?;
        if model.alive_bodies().is_empty() {
            return Err(AssemblyError::EmptyDefinition {
                name: name.to_string(),
            });
        }
        Ok(self.push_definition(Definition {
            name: name.to_string(),
            family,
            content: DefinitionContent::Part(Box::new(model)),
        }))
    }

    /// Add an empty sub-assembly definition; place into it with
    /// [`Parent::Definition`].
    ///
    /// # Errors
    ///
    /// [`AssemblyError::InvalidName`], [`AssemblyError::DuplicateDefinitionName`].
    pub fn add_subassembly(&mut self, name: &str) -> Result<DefinitionId, AssemblyError> {
        self.check_new_definition(name, "assembly")?;
        Ok(self.push_definition(Definition {
            name: name.to_string(),
            family: None,
            content: DefinitionContent::Assembly(Vec::new()),
        }))
    }

    fn check_new_definition(&self, name: &str, what: &'static str) -> Result<(), AssemblyError> {
        validate_name(name, what)?;
        if self.definition_named(name).is_some() {
            return Err(AssemblyError::DuplicateDefinitionName {
                name: name.to_string(),
            });
        }
        Ok(())
    }

    fn push_definition(&mut self, definition: Definition) -> DefinitionId {
        self.definitions.push(definition);
        DefinitionId(self.definitions.len() as u32 - 1)
    }

    /// The instances directly under `parent`.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::UnknownDefinition`];
    /// [`AssemblyError::WrongDefinitionKind`] when `parent` is a part.
    pub fn instances(&self, parent: Parent) -> Result<&[Instance], AssemblyError> {
        match parent {
            Parent::Root => Ok(&self.root),
            Parent::Definition(id) => {
                let definition = self.definition(id)?;
                match &definition.content {
                    DefinitionContent::Assembly(instances) => Ok(instances),
                    DefinitionContent::Part(_) => Err(AssemblyError::WrongDefinitionKind {
                        name: definition.name.clone(),
                        expected: "a sub-assembly",
                        actual: "a part",
                    }),
                }
            }
        }
    }

    fn instances_mut(&mut self, parent: Parent) -> Result<&mut Vec<Instance>, AssemblyError> {
        match parent {
            Parent::Root => Ok(&mut self.root),
            Parent::Definition(id) => {
                let definition = self
                    .definitions
                    .get_mut(id.0 as usize)
                    .ok_or(AssemblyError::UnknownDefinition { id: id.0 })?;
                match &mut definition.content {
                    DefinitionContent::Assembly(instances) => Ok(instances),
                    DefinitionContent::Part(_) => Err(AssemblyError::WrongDefinitionKind {
                        name: definition.name.clone(),
                        expected: "a sub-assembly",
                        actual: "a part",
                    }),
                }
            }
        }
    }

    pub(crate) fn parent_name(&self, parent: Parent) -> String {
        match parent {
            Parent::Root => "<root>".to_string(),
            Parent::Definition(id) => self
                .definition(id)
                .map(|d| d.name.clone())
                .unwrap_or_else(|_| format!("#{}", id.0)),
        }
    }

    /// Place `definition` under `parent`. `name` defaults to
    /// `<definition>_<n>`, `n` counting that definition's siblings
    /// from 1 (skipping taken names).
    ///
    /// # Errors
    ///
    /// [`AssemblyError::UnknownDefinition`], [`AssemblyError::InvalidName`],
    /// [`AssemblyError::DuplicateInstanceName`], [`AssemblyError::Cycle`],
    /// [`AssemblyError::WrongDefinitionKind`] when `parent` is a part.
    pub fn place(
        &mut self,
        parent: Parent,
        definition: DefinitionId,
        name: Option<&str>,
        placement: Placement,
    ) -> Result<InstanceRef, AssemblyError> {
        let placed = self.definition(definition)?.name.clone();
        if let Parent::Definition(container) = parent {
            self.check_acyclic(container, definition)?;
        }
        let siblings = self.instances(parent)?;
        let name = match name {
            Some(name) => {
                validate_name(name, "instance")?;
                if siblings.iter().any(|i| i.name == name) {
                    return Err(AssemblyError::DuplicateInstanceName {
                        parent: self.parent_name(parent),
                        name: name.to_string(),
                    });
                }
                name.to_string()
            }
            None => {
                let mut n = siblings
                    .iter()
                    .filter(|i| i.definition == definition)
                    .count()
                    + 1;
                loop {
                    let candidate = format!("{placed}_{n}");
                    if !siblings.iter().any(|i| i.name == candidate) {
                        break candidate;
                    }
                    n += 1;
                }
            }
        };
        let instances = self.instances_mut(parent)?;
        instances.push(Instance {
            name,
            definition,
            placement,
            visible: true,
        });
        let at = InstanceRef {
            parent,
            index: instances.len() - 1,
        };
        self.record_place(at, placement);
        Ok(at)
    }

    /// Refuse placing `definition` into `container` when `definition`
    /// is `container` or (transitively) contains it.
    fn check_acyclic(
        &self,
        container: DefinitionId,
        definition: DefinitionId,
    ) -> Result<(), AssemblyError> {
        let mut chain = vec![definition];
        if self.reaches(definition, container, &mut chain) {
            let names: Vec<String> = std::iter::once(container)
                .chain(chain)
                .map(|id| self.parent_name(Parent::Definition(id)))
                .collect();
            return Err(AssemblyError::Cycle {
                definition: self.parent_name(Parent::Definition(definition)),
                parent: self.parent_name(Parent::Definition(container)),
                chain: names.join(" > "),
            });
        }
        Ok(())
    }

    fn reaches(
        &self,
        from: DefinitionId,
        target: DefinitionId,
        chain: &mut Vec<DefinitionId>,
    ) -> bool {
        if from == target {
            return true;
        }
        let Ok(Definition {
            content: DefinitionContent::Assembly(instances),
            ..
        }) = self.definition(from)
        else {
            return false;
        };
        for instance in instances {
            chain.push(instance.definition);
            if self.reaches(instance.definition, target, chain) {
                return true;
            }
            chain.pop();
        }
        false
    }

    /// One instance by handle.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::UnknownInstance`].
    pub fn instance(&self, at: InstanceRef) -> Result<&Instance, AssemblyError> {
        self.instances(at.parent)?
            .get(at.index)
            .ok_or_else(|| AssemblyError::UnknownInstance {
                path: format!("{}[{}]", self.parent_name(at.parent), at.index),
            })
    }

    /// Replace an instance's placement.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::UnknownInstance`].
    pub fn set_placement(
        &mut self,
        at: InstanceRef,
        placement: Placement,
    ) -> Result<(), AssemblyError> {
        self.instance_mut(at)?.placement = placement;
        self.record_move(at, placement);
        Ok(())
    }

    /// Show or hide an instance (and everything under it).
    ///
    /// # Errors
    ///
    /// [`AssemblyError::UnknownInstance`].
    pub fn set_visible(&mut self, at: InstanceRef, visible: bool) -> Result<(), AssemblyError> {
        self.instance_mut(at)?.visible = visible;
        Ok(())
    }

    fn instance_mut(&mut self, at: InstanceRef) -> Result<&mut Instance, AssemblyError> {
        let parent = self.parent_name(at.parent);
        self.instances_mut(at.parent)?
            .get_mut(at.index)
            .ok_or_else(|| AssemblyError::UnknownInstance {
                path: format!("{parent}[{}]", at.index),
            })
    }
}
