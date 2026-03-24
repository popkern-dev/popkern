//! Product structure mapping: STEP product entities to [`StepProductNode`] assembly tree.
//!
//! Traces the STEP product chain:
//! `PRODUCT` → `PRODUCT_DEFINITION_FORMATION` → `PRODUCT_DEFINITION` →
//! `PRODUCT_DEFINITION_SHAPE` → `SHAPE_DEFINITION_REPRESENTATION` →
//! shape representation → BRep solids.
//!
//! Assembly relationships come from `NEXT_ASSEMBLY_USAGE_OCCURRENCE`.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use popkern_core::geometry::Transform3;
use popkern_core::topology::BRep;

use crate::step::diagnostics::DiagnosticLog;
use crate::step::document::{
    AssemblyChild, RepresentationContext, ShapeRep, ShapeRepKind, StepProductNode,
};
use crate::step::mapping::geometry::GeometryCache;
use crate::step::parser::graph::EntityGraph;
use crate::step::parser::types::EntityId;
use crate::step::schema::product::{
    NextAssemblyUsageOccurrence, ProductDefinition, ProductDefinitionFormation,
    ProductDefinitionShape, ShapeDefinitionRepresentation, StepProduct,
};

/// Map STEP product entities into a flat list of [`StepProductNode`]s and root indices.
///
/// `brep_solids` is the output of [`super::topology::map_brep_solids`], used to
/// attach geometry to the correct products.
pub fn map_products(
    graph: &EntityGraph,
    brep_solids: &[(EntityId, BRep)],
    geo: &mut GeometryCache<'_>,
    log: &mut DiagnosticLog,
) -> (Vec<StepProductNode>, Vec<usize>) {
    // Index: PRODUCT_DEFINITION entity_id → product_index
    let mut pd_to_index: HashMap<EntityId, usize> = HashMap::new();
    let mut products: Vec<StepProductNode> = Vec::new();

    // Build a map: shape_rep_entity_id → list of BRep solids referencing items in that rep
    let brep_by_id: HashMap<EntityId, &BRep> = brep_solids.iter().map(|(id, b)| (*id, b)).collect();

    // 1. Find all PRODUCT_DEFINITION entities and trace to PRODUCT
    let pd_ids: Vec<(EntityId, ProductDefinition)> = graph
        .iter()
        .filter_map(|(&id, e)| {
            if e.type_name == "PRODUCT_DEFINITION" {
                ProductDefinition::from_params(&e.params)
                    .ok()
                    .map(|pd| (id, pd))
            } else {
                None
            }
        })
        .collect();

    for (pd_id, pd) in &pd_ids {
        // Trace PD → PDF → PRODUCT
        let (name, product_id) = resolve_product_info(graph, pd, log);

        let node = StepProductNode {
            name,
            id: product_id,
            children: Vec::new(),
            shapes: Vec::new(),
            metadata: Default::default(),
        };
        let idx = products.len();
        products.push(node);
        pd_to_index.insert(*pd_id, idx);
    }

    // 2. Link shapes: PD_SHAPE → SDR → shape_rep → BRep solids
    let sdr_list: Vec<ShapeDefinitionRepresentation> = graph
        .iter()
        .filter_map(|(_, e)| {
            if e.type_name == "SHAPE_DEFINITION_REPRESENTATION" {
                ShapeDefinitionRepresentation::from_params(&e.params).ok()
            } else {
                None
            }
        })
        .collect();

    // Also build PD_SHAPE → PD index
    let pds_to_pd: HashMap<EntityId, EntityId> = graph
        .iter()
        .filter_map(|(&id, e)| {
            if e.type_name == "PRODUCT_DEFINITION_SHAPE" {
                ProductDefinitionShape::from_params(&e.params)
                    .ok()
                    .map(|pds| (id, pds.definition))
            } else {
                None
            }
        })
        .collect();

    for sdr in &sdr_list {
        // sdr.definition → PD_SHAPE → PD
        let pd_id = match pds_to_pd.get(&sdr.definition) {
            Some(&pd_id) => pd_id,
            None => continue,
        };
        let product_idx = match pd_to_index.get(&pd_id) {
            Some(&idx) => idx,
            None => continue,
        };

        // Find which rep this is and what kind
        let rep_id = sdr.used_representation;
        let rep_entity = match graph.get(rep_id) {
            Some(e) => e,
            None => continue,
        };

        let kind = match rep_entity.type_name.as_str() {
            "ADVANCED_BREP_SHAPE_REPRESENTATION" => ShapeRepKind::AdvancedBrep,
            "FACETED_BREP_SHAPE_REPRESENTATION" => ShapeRepKind::FacetedBrep,
            "SHAPE_REPRESENTATION" => ShapeRepKind::Other("SHAPE_REPRESENTATION".into()),
            other => ShapeRepKind::Other(other.into()),
        };

        // Collect BRep solids that are items of this representation
        let items =
            crate::step::schema::product::StepShapeRepresentation::from_params(&rep_entity.params);
        let solids: Vec<BRep> = match items {
            Ok(sr) => sr
                .items
                .iter()
                .filter_map(|item_id| brep_by_id.get(item_id).map(|b| (*b).clone()))
                .collect(),
            Err(_) => Vec::new(),
        };

        let shape = Arc::new(ShapeRep {
            kind,
            solids,
            context: RepresentationContext::default(),
            non_brep_items: Vec::new(),
        });
        products[product_idx].shapes.push(shape);
    }

    // 3. Assembly relationships via NEXT_ASSEMBLY_USAGE_OCCURRENCE
    let nauo_list: Vec<(EntityId, NextAssemblyUsageOccurrence)> = graph
        .iter()
        .filter_map(|(&id, e)| {
            if e.type_name == "NEXT_ASSEMBLY_USAGE_OCCURRENCE" {
                NextAssemblyUsageOccurrence::from_params(&e.params)
                    .ok()
                    .map(|n| (id, n))
            } else {
                None
            }
        })
        .collect();

    let mut child_pds: HashSet<EntityId> = HashSet::new();
    for (_nauo_id, nauo) in &nauo_list {
        let parent_idx = pd_to_index.get(&nauo.relating).copied();
        let child_idx = pd_to_index.get(&nauo.related).copied();

        if let (Some(pi), Some(ci)) = (parent_idx, child_idx) {
            child_pds.insert(nauo.related);

            // Try to resolve transform from REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION
            let transform = resolve_nauo_transform(graph, geo, nauo, log);

            products[pi].children.push(AssemblyChild {
                product_index: ci,
                transform,
                instance_name: nauo.reference_designator.clone(),
            });
        }
    }

    // 4. Identify roots: PDs not referenced as children
    let roots: Vec<usize> = pd_ids
        .iter()
        .filter_map(|(pd_id, _)| {
            if !child_pds.contains(pd_id) {
                pd_to_index.get(pd_id).copied()
            } else {
                None
            }
        })
        .collect();

    (products, roots)
}

/// Trace PRODUCT_DEFINITION → PDF → PRODUCT to get name and id.
fn resolve_product_info(
    graph: &EntityGraph,
    pd: &ProductDefinition,
    _log: &mut DiagnosticLog,
) -> (String, String) {
    // PD → PDF (formation)
    let pdf = graph
        .get(pd.formation)
        .and_then(|e| ProductDefinitionFormation::from_params(&e.params).ok());

    if let Some(pdf) = pdf {
        // PDF → PRODUCT (of_product)
        if let Some(product) = graph
            .get(pdf.of_product)
            .and_then(|e| StepProduct::from_params(&e.params).ok())
        {
            return (product.name, product.id);
        }
    }

    (pd.id.clone(), pd.id.clone())
}

/// Attempt to resolve the assembly transform for a NAUO.
///
/// This is a simplified version that returns identity when no transform is found.
fn resolve_nauo_transform(
    _graph: &EntityGraph,
    _geo: &mut GeometryCache<'_>,
    _nauo: &NextAssemblyUsageOccurrence,
    _log: &mut DiagnosticLog,
) -> Transform3 {
    // Full implementation would trace:
    // NAUO → CONTEXT_DEPENDENT_SHAPE_REPRESENTATION → REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION
    // → ITEM_DEFINED_TRANSFORMATION → two AXIS2_PLACEMENT_3D → compute relative transform
    // For now, return identity.
    Transform3::identity()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::step::parser::types::{EntityInstance, Parameter};

    fn make_entity(id: EntityId, type_name: &str, params: Vec<Parameter>) -> EntityInstance {
        EntityInstance {
            id,
            type_name: type_name.into(),
            params,
            complex_parts: vec![],
        }
    }

    fn p_str(s: &str) -> Parameter {
        Parameter::String(s.into())
    }
    fn p_ref(id: u64) -> Parameter {
        Parameter::EntityRef(id)
    }
    fn p_list(items: Vec<Parameter>) -> Parameter {
        Parameter::List(items)
    }

    #[test]
    fn map_single_product() {
        let entities = vec![
            // PRODUCT
            make_entity(
                1,
                "PRODUCT",
                vec![
                    p_str("P001"),
                    p_str("TestWidget"),
                    p_str("A widget"),
                    p_list(vec![p_ref(100)]),
                ],
            ),
            // PDF
            make_entity(
                2,
                "PRODUCT_DEFINITION_FORMATION",
                vec![p_str("v1"), p_str(""), p_ref(1)],
            ),
            // PD
            make_entity(
                3,
                "PRODUCT_DEFINITION",
                vec![p_str("pd1"), p_str(""), p_ref(2), p_ref(100)],
            ),
            // PD_SHAPE
            make_entity(
                4,
                "PRODUCT_DEFINITION_SHAPE",
                vec![p_str(""), p_str(""), p_ref(3)],
            ),
            // SDR
            make_entity(
                5,
                "SHAPE_DEFINITION_REPRESENTATION",
                vec![p_ref(4), p_ref(6)],
            ),
            // Shape rep (empty)
            make_entity(
                6,
                "ADVANCED_BREP_SHAPE_REPRESENTATION",
                vec![p_str("rep"), p_list(vec![]), p_ref(100)],
            ),
        ];

        let mut log = DiagnosticLog::new();
        let graph = EntityGraph::build(entities, &mut log);
        let mut geo = GeometryCache::new(&graph);

        let (products, roots) = map_products(&graph, &[], &mut geo, &mut log);
        assert_eq!(products.len(), 1);
        assert_eq!(products[0].name, "TestWidget");
        assert_eq!(products[0].id, "P001");
        assert_eq!(roots, vec![0]);
    }
}
