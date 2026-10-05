//! Editing a part definition in place (PK-678): a project migration
//! that continues an assembly (`let a = previous();`) changes one
//! part's feature tree — a feature inserted at the history playhead,
//! a parameter, a fillet — and every instance of it follows, because
//! geometry lives only in definitions.

use popkern_script::Model;

use crate::assembly::{Assembly, DefinitionContent, DefinitionId};
use crate::error::AssemblyError;

impl Assembly {
    /// Replace part `id`'s model with `model`, its edited version,
    /// keeping every instance. When the assembly has mates they are
    /// solved again against the edited geometry. Transactional: a
    /// model with no alive body, or mates that no longer hold, leave
    /// the assembly exactly as it was.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::UnknownDefinition`];
    /// [`AssemblyError::WrongDefinitionKind`] for a sub-assembly or a
    /// family variant (a variant is generated from its family's script:
    /// edit the family); [`AssemblyError::EmptyDefinition`];
    /// [`AssemblyError::MateConflict`] / [`AssemblyError::Mate`] when
    /// the mates cannot be solved on the edited part.
    pub fn replace_part_model(
        &mut self,
        id: DefinitionId,
        mut model: Model,
    ) -> Result<(), AssemblyError> {
        let definition = self.definition(id)?;
        let name = definition.name.clone();
        if definition.family.is_some() {
            return Err(AssemblyError::WrongDefinitionKind {
                name,
                expected: "a part",
                actual: "a family variant (edit the family it is generated from)",
            });
        }
        if matches!(definition.content, DefinitionContent::Assembly(_)) {
            return Err(AssemblyError::WrongDefinitionKind {
                name,
                expected: "a part",
                actual: "a sub-assembly",
            });
        }
        if model.alive_bodies().is_empty() {
            return Err(AssemblyError::EmptyDefinition { name });
        }
        let previous = std::mem::replace(self.part_model_mut(id)?, model);
        if !self.mates.is_empty()
            && let Err(error) = self.solve_mates()
        {
            *self.part_model_mut(id)? = previous;
            return Err(error);
        }
        Ok(())
    }
}
