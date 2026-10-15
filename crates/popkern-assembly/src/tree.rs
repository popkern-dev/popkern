//! Walking the instance tree: path resolution, the expanded tree a
//! viewer lists, and the part (leaf) occurrences.

use std::collections::HashMap;

use crate::assembly::{Assembly, DefinitionContent, Instance, InstanceRef, Parent, TreeRow};
use crate::error::AssemblyError;
use crate::names::InstancePath;
use crate::placement::Placement;

/// Where an instance sits while walking (`None`: not placed yet).
type PlacementOf<'a> = dyn Fn(InstanceRef, &Instance) -> Option<Placement> + 'a;

/// The parent being expanded.
struct Frame {
    parent: Parent,
    prefix: InstancePath,
    world: Placement,
    visible: bool,
}

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
        self.expand_root(&|_, instance| Some(instance.placement))
    }

    /// The expanded tree after the first `count` construction steps
    /// ([`Assembly::history`]): only the instances placed by then, at
    /// their placements then. Visibility is the document's.
    pub fn tree_after(&self, count: usize) -> Vec<TreeRow> {
        let placed: HashMap<InstanceRef, Placement> =
            self.placements_after(count).into_iter().collect();
        self.expand_root(&|at, _| placed.get(&at).copied())
    }

    /// The construction view after the first `count` steps: the
    /// [`Assembly::tree_after`] rows, plus the contents placed so far
    /// in every sub-assembly occurrence that is not placed yet, drawn at
    /// that occurrence's final world placement. Sub-assemblies are
    /// usually filled before they are placed (the six-axis arm fills
    /// its joint modules for ~570 steps before placing a link), and
    /// `tree_after` alone shows nothing until then (PK-831). Rows keep
    /// their final occurrence paths, so a step's own occurrences
    /// (the timeline's `instances`) name rows on screen.
    pub fn construction_tree_after(&self, count: usize) -> Vec<TreeRow> {
        let placed: HashMap<InstanceRef, Placement> =
            self.placements_after(count).into_iter().collect();
        let placement_of = |at: InstanceRef, _: &Instance| placed.get(&at).copied();
        let mut rows = self.expand_root(&placement_of);
        let mut shown: std::collections::HashSet<String> =
            rows.iter().map(|row| row.path.to_string()).collect();
        // Final rows come depth-first, parents before children, so an
        // unplaced occurrence inside another unplaced one is reached
        // after its parent's placed contents are already shown.
        for occurrence in self.tree() {
            if occurrence.is_part || shown.contains(&occurrence.path.to_string()) {
                continue;
            }
            let Ok(definition) = self.definition(occurrence.definition) else {
                continue;
            };
            let DefinitionContent::Assembly(children) = &definition.content else {
                continue;
            };
            let frame = Frame {
                parent: Parent::Definition(occurrence.definition),
                prefix: occurrence.path.clone(),
                world: occurrence.world,
                visible: occurrence.visible,
            };
            let mut inner = Vec::new();
            self.expand(children, &frame, &placement_of, &mut inner);
            shown.insert(occurrence.path.to_string());
            for row in inner {
                if shown.insert(row.path.to_string()) {
                    rows.push(row);
                }
            }
        }
        rows
    }

    fn expand_root(&self, placement_of: &PlacementOf) -> Vec<TreeRow> {
        let mut rows = Vec::new();
        let top = Frame {
            parent: Parent::Root,
            prefix: InstancePath(Vec::new()),
            world: Placement::identity(),
            visible: true,
        };
        self.expand(&self.root, &top, placement_of, &mut rows);
        rows
    }

    fn expand(
        &self,
        instances: &[Instance],
        frame: &Frame,
        placement_of: &PlacementOf,
        rows: &mut Vec<TreeRow>,
    ) {
        for (index, instance) in instances.iter().enumerate() {
            let at = InstanceRef {
                parent: frame.parent,
                index,
            };
            let Some(placement) = placement_of(at, instance) else {
                continue;
            };
            let Ok(definition) = self.definition(instance.definition) else {
                continue;
            };
            let path = frame.prefix.child(&instance.name);
            let world = placement.then(&frame.world);
            let visible = frame.visible && instance.visible;
            rows.push(TreeRow {
                path: path.clone(),
                definition: instance.definition,
                world,
                visible,
                is_part: matches!(definition.content, DefinitionContent::Part(_)),
            });
            if let DefinitionContent::Assembly(children) = &definition.content {
                let inner = Frame {
                    parent: Parent::Definition(instance.definition),
                    prefix: path,
                    world,
                    visible,
                };
                self.expand(children, &inner, placement_of, rows);
            }
        }
    }

    /// Every part (leaf) occurrence, depth-first.
    pub fn leaves(&self) -> Vec<TreeRow> {
        self.tree().into_iter().filter(|row| row.is_part).collect()
    }
}
