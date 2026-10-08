//! STEP assembly import (design §8.3, PK-631): a STEP product tree
//! becomes an assembly — a product with solids a part definition (each
//! solid a `RawBody` snapshot of its own model: re-exportable,
//! savable, inspectable, with no recipe), a product with children a
//! sub-assembly, the root's occurrences the top-level instances.
//!
//! Nothing is guessed: a product whose solids all fail to record, an
//! occurrence whose placement is not rigid, or one that would close a
//! cycle is left out and listed in the [`StepImportReport`] with the
//! reason; the rest of the tree loads.

use std::collections::HashMap;

use popkern_history::SnapshotOrigin;
use popkern_io::{StepDocument, StepImportOptions};
use popkern_script::Model;
use serde::Serialize;

use crate::assembly::{Assembly, DefinitionId, Parent};
use crate::error::AssemblyError;
use crate::placement::Placement;

/// What an import loaded and what it left out.
#[derive(Debug, Clone, Default, Serialize)]
pub struct StepImportReport {
    /// Products in the file.
    pub products: usize,
    /// Part definitions created.
    pub parts: usize,
    /// Sub-assembly definitions created.
    pub subassemblies: usize,
    /// Occurrences placed (top-level and nested).
    pub instances: usize,
    /// Solids recorded.
    pub solids: usize,
    /// Solids the recorder refused (invalid geometry), by product.
    pub refused_solids: Vec<String>,
    /// Occurrences left out, with the reason.
    pub refused_instances: Vec<String>,
    /// The STEP reader's warnings.
    pub reader_warnings: usize,
}

impl Assembly {
    /// Import a STEP file's product structure as an assembly named
    /// `name` (the file stem, typically).
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Step`] when the file cannot be read at all, or
    /// when it holds no product with geometry.
    pub fn from_step_bytes(
        bytes: &[u8],
        name: &str,
    ) -> Result<(Assembly, StepImportReport), AssemblyError> {
        let imported = popkern_io::import_step(bytes, &StepImportOptions::default())
            .map_err(|e| AssemblyError::Step(e.to_string()))?;
        let mut report = StepImportReport {
            products: imported.document.products.len(),
            reader_warnings: imported.diagnostics.warnings().count(),
            ..StepImportReport::default()
        };
        let assembly = build(&imported.document, name, &mut report)?;
        let mut assembly = assembly;
        assembly.mark_history_unrecorded();
        Ok((assembly, report))
    }
}

/// A definition name from a product: its name, else its id, made valid
/// (no path separators) and unique.
fn definition_name(
    taken: &mut HashMap<String, usize>,
    product: &popkern_io::StepProductNode,
    index: usize,
) -> String {
    let raw = [product.name.trim(), product.id.trim()]
        .into_iter()
        .find(|s| !s.is_empty())
        .map(str::to_string)
        .unwrap_or_else(|| format!("product_{index}"));
    let clean: String = raw
        .chars()
        .map(|c| if c == '/' || c == ':' { '_' } else { c })
        .collect();
    unique(taken, clean)
}

fn unique(taken: &mut HashMap<String, usize>, base: String) -> String {
    let count = taken.entry(base.clone()).or_insert(0);
    *count += 1;
    if *count == 1 {
        base
    } else {
        format!("{base}_{count}")
    }
}

fn build(
    doc: &StepDocument,
    name: &str,
    report: &mut StepImportReport,
) -> Result<Assembly, AssemblyError> {
    let mut assembly = Assembly::new(name);
    let mut taken = HashMap::new();
    // Definitions reachable from the roots, children first so every
    // sub-assembly's parts exist when it is filled.
    let mut order: Vec<usize> = Vec::new();
    let mut visiting = vec![false; doc.products.len()];
    for &root in &doc.roots {
        visit(doc, root, &mut order, &mut visiting);
    }
    // A single root assembly IS the document: its occurrences become the
    // top-level instances, and it gets no definition of its own.
    let single_root_assembly =
        doc.roots.len() == 1 && !doc.products[doc.roots[0]].children.is_empty();
    let mut ids: HashMap<usize, DefinitionId> = HashMap::new();
    // A product's own solids, when it also has children, ride in a part
    // placed at its origin inside the sub-assembly.
    let mut own_parts: HashMap<usize, DefinitionId> = HashMap::new();
    for &index in &order {
        let product = &doc.products[index];
        let solids: Vec<&popkern_core::topology::BRep> = product
            .shapes
            .iter()
            .flat_map(|s| s.solids.iter())
            .collect();
        let part = if solids.is_empty() {
            None
        } else {
            let label = definition_name(&mut taken, product, index);
            record_part(
                &mut assembly,
                &label,
                &solids,
                product.metadata.solid_color,
                report,
            )?
        };
        if product.children.is_empty() {
            if let Some(id) = part {
                ids.insert(index, id);
            }
        } else if single_root_assembly && index == doc.roots[0] {
            if let Some(id) = part {
                own_parts.insert(index, id);
            }
        } else {
            let label = definition_name(&mut taken, product, index);
            let sub = assembly.add_subassembly(&label)?;
            report.subassemblies += 1;
            ids.insert(index, sub);
            if let Some(id) = part {
                own_parts.insert(index, id);
            }
        }
    }
    let place = |assembly: &mut Assembly,
                 parent: Parent,
                 child: &popkern_io::AssemblyChild,
                 siblings: &mut HashMap<String, usize>,
                 report: &mut StepImportReport| {
        let Some(&definition) = ids.get(&child.product_index) else {
            report.refused_instances.push(format!(
                "occurrence of product {}: the product has no usable geometry",
                child.product_index
            ));
            return;
        };
        let def_name = assembly
            .definition(definition)
            .map(|d| d.name.clone())
            .unwrap_or_default();
        let wanted = child
            .instance_name
            .as_deref()
            .map(|n| {
                n.chars()
                    .map(|c| if c == '/' || c == ':' { '_' } else { c })
                    .collect::<String>()
            })
            .filter(|n| !n.trim().is_empty())
            .unwrap_or_else(|| def_name.clone());
        let instance = unique(siblings, wanted.trim().to_string());
        let placement = match Placement::from_matrix(child.transform, &instance) {
            Ok(p) => p,
            Err(e) => {
                report.refused_instances.push(format!("{instance}: {e}"));
                return;
            }
        };
        match assembly.place(parent, definition, Some(&instance), placement) {
            Ok(_) => report.instances += 1,
            Err(e) => report.refused_instances.push(format!("{instance}: {e}")),
        }
    };
    for &index in &order {
        let product = &doc.products[index];
        if product.children.is_empty() {
            continue;
        }
        let Some(&sub) = ids.get(&index) else {
            continue;
        };
        let mut siblings = HashMap::new();
        if let Some(&own) = own_parts.get(&index) {
            let label = unique(
                &mut siblings,
                assembly
                    .definition(own)
                    .map(|d| d.name.clone())
                    .unwrap_or_default(),
            );
            assembly.place(
                Parent::Definition(sub),
                own,
                Some(&label),
                Placement::identity(),
            )?;
            report.instances += 1;
        }
        for child in &product.children {
            place(
                &mut assembly,
                Parent::Definition(sub),
                child,
                &mut siblings,
                report,
            );
        }
    }
    // The top level: a single root assembly's occurrences ARE the
    // document's; otherwise every root is placed at the origin.
    let mut siblings = HashMap::new();
    if single_root_assembly {
        let root = &doc.products[doc.roots[0]];
        if !root.name.trim().is_empty() {
            assembly.name = root.name.trim().to_string();
        }
        if let Some(&own) = own_parts.get(&doc.roots[0]) {
            let label = unique(
                &mut siblings,
                assembly
                    .definition(own)
                    .map(|d| d.name.clone())
                    .unwrap_or_default(),
            );
            assembly.place(Parent::Root, own, Some(&label), Placement::identity())?;
            report.instances += 1;
        }
        for child in &root.children {
            place(&mut assembly, Parent::Root, child, &mut siblings, report);
        }
    } else {
        for &root in &doc.roots {
            let child = popkern_io::AssemblyChild {
                product_index: root,
                transform: popkern_core::geometry::Transform3::identity(),
                instance_name: None,
            };
            place(&mut assembly, Parent::Root, &child, &mut siblings, report);
        }
    }
    if report.parts == 0 {
        // Name what the file carries instead of solids (a wireframe or
        // surface-only representation, typically).
        let mut kinds: Vec<String> = doc
            .products
            .iter()
            .flat_map(|p| p.shapes.iter())
            .map(|s| match &s.kind {
                popkern_io::ShapeRepKind::AdvancedBrep => {
                    "ADVANCED_BREP_SHAPE_REPRESENTATION".to_string()
                }
                popkern_io::ShapeRepKind::FacetedBrep => {
                    "FACETED_BREP_SHAPE_REPRESENTATION".to_string()
                }
                popkern_io::ShapeRepKind::Other(name) => name.clone(),
            })
            .collect();
        kinds.sort();
        kinds.dedup();
        let carried = if kinds.is_empty() {
            "no shape representation".to_string()
        } else {
            kinds.join(", ")
        };
        return Err(AssemblyError::Step(format!(
            "the file holds no product with solid geometry ({} product(s), {} solid(s) refused; it carries: {carried})",
            report.products,
            report.refused_solids.len()
        )));
    }
    Ok(assembly)
}

fn visit(doc: &StepDocument, index: usize, order: &mut Vec<usize>, visiting: &mut [bool]) {
    if order.contains(&index) || visiting.get(index).copied().unwrap_or(true) {
        return;
    }
    visiting[index] = true;
    for child in &doc.products[index].children {
        visit(doc, child.product_index, order, visiting);
    }
    visiting[index] = false;
    order.push(index);
}

/// One part definition holding every solid that records.
fn record_part(
    assembly: &mut Assembly,
    label: &str,
    solids: &[&popkern_core::topology::BRep],
    colour: Option<[f32; 3]>,
    report: &mut StepImportReport,
) -> Result<Option<DefinitionId>, AssemblyError> {
    // The part is built at the tolerance its solids were imported at
    // (the file's uncertainty when looser than the default, PK-639):
    // validating an imported body at the kernel default refused gaps
    // the file's own precision allows (measured: a 2.9e-5 shell gap in
    // a weldment arm written to 2.3e-3).
    let tolerance = solids
        .iter()
        .map(|s| s.tolerance)
        .max_by(|a, b| a.linear().total_cmp(&b.linear()))
        .unwrap_or_default();
    let mut model = Model::with_tolerance(tolerance);
    let mut recorded = 0;
    for (k, solid) in solids.iter().enumerate() {
        let origin = SnapshotOrigin::Import {
            source: format!("STEP product '{label}' solid {k}"),
        };
        match model.raw_body((*solid).clone(), origin) {
            Ok(_) => recorded += 1,
            Err(e) => report
                .refused_solids
                .push(format!("{label} solid {k}: {e}")),
        }
    }
    if recorded == 0 {
        return Ok(None);
    }
    report.solids += recorded;
    report.parts += 1;
    // The file's solid colour becomes the part's appearance (finish at
    // its defaults: STEP carries no metalness or roughness).
    if let Some([r, g, b]) = colour
        && let Ok(appearance) = popkern_history::appearance::Appearance::new(
            [f64::from(r), f64::from(g), f64::from(b)],
            0.0,
            0.5,
        )
    {
        model.set_appearance(Some(appearance));
    }
    Ok(Some(assembly.add_part(label, model)?))
}
