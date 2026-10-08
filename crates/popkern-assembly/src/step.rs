//! STEP export of an assembly as a product structure (design §6.3):
//! a product per used definition, the root product named after the
//! document, and an occurrence per instance placed through a
//! context-dependent shape representation (`popkern-io`
//! `export_assembly`). A definition's solids are written once however
//! many instances place it.

use std::collections::HashMap;
use std::sync::Arc;

use popkern_io::{
    AssemblyChild, OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind,
    StepDocument, StepExportOptions, StepProductNode, UnitContext,
};

use crate::assembly::{Assembly, DefinitionContent, DefinitionId, Instance};
use crate::error::AssemblyError;

impl Assembly {
    /// The product structure as a [`StepDocument`] whose root product
    /// is called `name`.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Script`] when a definition fails to build.
    pub fn to_step_document(&mut self, name: &str) -> Result<StepDocument, AssemblyError> {
        let mut parts = self.definition_parts()?;
        // Products in first-use order (depth-first from the root), so
        // the file reads top-down and unused definitions are left out.
        let mut order: Vec<DefinitionId> = Vec::new();
        self.collect_used(&self.root, &mut order);
        let index: HashMap<DefinitionId, usize> = order
            .iter()
            .enumerate()
            .map(|(i, id)| (*id, i + 1))
            .collect();
        let children = |instances: &[Instance]| -> Vec<AssemblyChild> {
            instances
                .iter()
                .filter_map(|instance| {
                    Some(AssemblyChild {
                        product_index: *index.get(&instance.definition)?,
                        transform: *instance.placement.transform(),
                        instance_name: Some(instance.name.clone()),
                    })
                })
                .collect()
        };
        let mut products = vec![product(name, children(&self.root), Vec::new())];
        for id in &order {
            let definition = self.definition(*id)?;
            let node = match &definition.content {
                DefinitionContent::Part(model) => {
                    let solids = parts
                        .remove(id)
                        .unwrap_or_default()
                        .iter()
                        .map(|part| part.brep().clone())
                        .collect();
                    let mut node = product(&definition.name, Vec::new(), solids);
                    node.metadata.solid_color = model.appearance().map(|a| a.rgb_f32());
                    node
                }
                DefinitionContent::Assembly(instances) => {
                    product(&definition.name, children(instances), Vec::new())
                }
            };
            products.push(node);
        }
        Ok(StepDocument {
            products,
            roots: vec![0],
            units: UnitContext::default(),
            opaque_bag: OpaqueBag::default(),
        })
    }

    fn collect_used(&self, instances: &[Instance], order: &mut Vec<DefinitionId>) {
        for instance in instances {
            if order.contains(&instance.definition) {
                continue;
            }
            order.push(instance.definition);
            if let Ok(definition) = self.definition(instance.definition)
                && let DefinitionContent::Assembly(children) = &definition.content
            {
                self.collect_used(children, order);
            }
        }
    }

    /// The assembly as STEP Part 21 bytes, root product `name`.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Script`] when a definition fails to build;
    /// [`AssemblyError::Step`] when the writer declines a body.
    pub fn to_step_bytes(&mut self, name: &str) -> Result<Vec<u8>, AssemblyError> {
        let doc = self.to_step_document(name)?;
        popkern_io::export_step(&doc, &StepExportOptions::default())
            .map_err(|e| AssemblyError::Step(e.to_string()))
    }

    /// Write [`Assembly::to_step_bytes`] to `path`, the file stem naming
    /// the root product.
    ///
    /// # Errors
    ///
    /// As [`Assembly::to_step_bytes`], and [`AssemblyError::Step`] when
    /// the file cannot be written.
    pub fn to_step(&mut self, path: impl AsRef<std::path::Path>) -> Result<(), AssemblyError> {
        let path = path.as_ref();
        let name = path
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| self.name.clone());
        let bytes = self.to_step_bytes(&name)?;
        std::fs::write(path, bytes)
            .map_err(|e| AssemblyError::Step(format!("{}: {e}", path.display())))
    }
}

fn product(
    name: &str,
    children: Vec<AssemblyChild>,
    solids: Vec<popkern_core::topology::BRep>,
) -> StepProductNode {
    StepProductNode {
        name: name.to_string(),
        id: name.to_string(),
        children,
        shapes: if solids.is_empty() {
            Vec::new()
        } else {
            vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids,
                context: RepresentationContext::default(),
                non_brep_items: Vec::new(),
            })]
        },
        metadata: ProductMetadata::default(),
    }
}
