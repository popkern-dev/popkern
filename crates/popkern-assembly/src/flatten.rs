//! Flattening: world-placed geometry of every part occurrence, on
//! demand (design §6.4). Nothing here is stored: a definition's bodies
//! are read once per call and transformed per occurrence.

use std::collections::HashMap;

use popkern_core::topology::BRep;
use popkern_part::Part;

use crate::assembly::{Assembly, DefinitionContent, DefinitionId};
use crate::error::AssemblyError;
use crate::names::InstancePath;
use crate::placement::Placement;

/// One alive body of one part occurrence, in world coordinates.
pub struct FlatBody {
    /// The occurrence.
    pub path: InstancePath,
    /// The part definition.
    pub definition: DefinitionId,
    /// Index of the body among the definition's alive bodies.
    pub body: usize,
    /// The occurrence's world placement.
    pub world: Placement,
    /// The body moved to its world placement.
    pub brep: BRep,
    /// Shown in viewers.
    pub visible: bool,
}

impl Assembly {
    /// Every part definition's alive bodies as validated parts, in
    /// body order, keyed by definition (replaying a stale definition
    /// first).
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Script`] when a definition's model fails.
    pub fn definition_parts(&mut self) -> Result<HashMap<DefinitionId, Vec<Part>>, AssemblyError> {
        let mut parts = HashMap::new();
        for (i, definition) in self.definitions.iter_mut().enumerate() {
            if let DefinitionContent::Part(model) = &mut definition.content {
                let built = model
                    .to_parts_all()
                    .map_err(|source| AssemblyError::Script {
                        definition: definition.name.clone(),
                        source,
                    })?;
                parts.insert(DefinitionId(i as u32), built);
            }
        }
        Ok(parts)
    }

    /// Every alive body of every part occurrence, moved to its world
    /// placement, depth-first in instance order.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Script`] when a definition fails to build;
    /// [`AssemblyError::Geometry`] when a body refuses its placement.
    pub fn flatten(&mut self) -> Result<Vec<FlatBody>, AssemblyError> {
        let parts = self.definition_parts()?;
        let mut flat = Vec::new();
        for row in self.leaves() {
            let Some(bodies) = parts.get(&row.definition) else {
                continue;
            };
            for (body, part) in bodies.iter().enumerate() {
                let brep = if row.world.is_identity() {
                    part.brep().clone()
                } else {
                    part.brep().transform(row.world.transform()).map_err(|e| {
                        AssemblyError::Geometry {
                            instance: row.path.to_string(),
                            detail: e.to_string(),
                        }
                    })?
                };
                flat.push(FlatBody {
                    path: row.path.clone(),
                    definition: row.definition,
                    body,
                    world: row.world,
                    brep,
                    visible: row.visible,
                });
            }
        }
        Ok(flat)
    }
}
