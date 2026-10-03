//! Walking the instance tree: path resolution, the expanded tree a
//! viewer lists, and the part (leaf) occurrences.

use crate::assembly::{Assembly, DefinitionContent, Instance, InstanceRef, Parent, TreeRow};
use crate::error::AssemblyError;
use crate::names::InstancePath;
use crate::placement::Placement;

impl Assembly {
    /// Resolve a path from the root to the instance it names.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::UnknownInstance`].
    pub fn resolve(&self, path: &InstancePath) -> Result<InstanceRef, AssemblyError> {
        let unknown = || AssemblyError::UnknownInstance {
            path: path.to_string(),
        };
        let mut parent = Parent::Root;
        let mut found = None;
        for (depth, name) in path.0.iter().enumerate() {
            if depth > 0 {
                let Some(at) = found else {
                    return Err(unknown());
                };
                parent = Parent::Definition(self.instance(at)?.definition);
            }
            let index = self
                .instances(parent)
                .map_err(|_| unknown())?
                .iter()
                .position(|i| &i.name == name)
                .ok_or_else(unknown)?;
            found = Some(InstanceRef { parent, index });
        }
        found.ok_or_else(unknown)
    }

    /// The expanded tree, depth-first in insertion order: every
    /// occurrence with its world placement and effective visibility.
    pub fn tree(&self) -> Vec<TreeRow> {
        let mut rows = Vec::new();
        self.expand(
            &self.root,
            &InstancePath(Vec::new()),
            &Placement::identity(),
            true,
            &mut rows,
        );
        rows
    }

    fn expand(
        &self,
        instances: &[Instance],
        prefix: &InstancePath,
        parent_world: &Placement,
        parent_visible: bool,
        rows: &mut Vec<TreeRow>,
    ) {
        for instance in instances {
            let Ok(definition) = self.definition(instance.definition) else {
                continue;
            };
            let path = prefix.child(&instance.name);
            let world = instance.placement.then(parent_world);
            let visible = parent_visible && instance.visible;
            rows.push(TreeRow {
                path: path.clone(),
                definition: instance.definition,
                world,
                visible,
                is_part: matches!(definition.content, DefinitionContent::Part(_)),
            });
            if let DefinitionContent::Assembly(children) = &definition.content {
                self.expand(children, &path, &world, visible, rows);
            }
        }
    }

    /// Every part (leaf) occurrence, depth-first.
    pub fn leaves(&self) -> Vec<TreeRow> {
        self.tree().into_iter().filter(|row| row.is_part).collect()
    }
}
