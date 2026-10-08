//! Presentation mapping: STEP styling entities to product metadata.
//!
//! Traces `STYLED_ITEM` → `PRESENTATION_STYLE_ASSIGNMENT` →
//! `SURFACE_STYLE_USAGE` → `SURFACE_SIDE_STYLE` → `SURFACE_STYLE_FILL_AREA` →
//! `FILL_AREA_STYLE` → `FILL_AREA_STYLE_COLOUR` → `COLOUR_RGB`.
//!
//! Also handles `PRESENTATION_LAYER_ASSIGNMENT` for layer names.

use crate::step::diagnostics::DiagnosticLog;
use crate::step::document::StepProductNode;
use crate::step::parser::graph::EntityGraph;
use crate::step::parser::types::EntityId;
use crate::step::schema::presentation::{
    ColourRgb, DraughtingPreDefinedColour, FillAreaStyle, FillAreaStyleColour,
    PresentationLayerAssignment, PresentationStyleAssignment, StyledItem, SurfaceSideStyle,
    SurfaceStyleFillArea, SurfaceStyleUsage,
};

/// Walk presentation entities and attach colors and layers to products.
///
/// This is a best-effort pass; unrecognized style chains are silently skipped.
pub fn map_presentation(
    graph: &EntityGraph,
    products: &mut [StepProductNode],
    log: &mut DiagnosticLog,
) {
    map_styled_items(graph, products, log);
    map_layer_assignments(graph, products, log);
}

/// Trace STYLED_ITEM → colour chain → assign solid_color on matching products.
fn map_styled_items(
    graph: &EntityGraph,
    products: &mut [StepProductNode],
    log: &mut DiagnosticLog,
) {
    for (_, entity) in graph.iter() {
        if entity.type_name != "STYLED_ITEM" {
            continue;
        }
        let si = match StyledItem::from_params(&entity.params) {
            Ok(s) => s,
            Err(_) => continue,
        };

        // Try to resolve color from the first style assignment
        let color = si
            .styles
            .iter()
            .find_map(|&style_id| resolve_color_from_psa(graph, style_id, log));

        if let Some(rgb) = color {
            apply_color_to_owner(products, si.item, rgb);
        }
    }
}

/// Trace presentation style assignment → colour RGB.
fn resolve_color_from_psa(
    graph: &EntityGraph,
    psa_id: EntityId,
    _log: &mut DiagnosticLog,
) -> Option<[f32; 3]> {
    let psa_entity = graph.get(psa_id)?;
    let psa = PresentationStyleAssignment::from_params(&psa_entity.params).ok()?;

    for &style_id in &psa.styles {
        let style_entity = graph.get(style_id)?;
        if style_entity.type_name != "SURFACE_STYLE_USAGE" {
            continue;
        }
        let ssu = SurfaceStyleUsage::from_params(&style_entity.params).ok()?;

        let sss_entity = graph.get(ssu.style)?;
        let sss = SurfaceSideStyle::from_params(&sss_entity.params).ok()?;

        for &elem_id in &sss.styles {
            let elem_entity = graph.get(elem_id)?;
            if elem_entity.type_name != "SURFACE_STYLE_FILL_AREA" {
                continue;
            }
            let ssfa = SurfaceStyleFillArea::from_params(&elem_entity.params).ok()?;

            let fas_entity = graph.get(ssfa.fill_area)?;
            let fas = FillAreaStyle::from_params(&fas_entity.params).ok()?;

            for &fs_id in &fas.fill_styles {
                let fs_entity = graph.get(fs_id)?;
                if fs_entity.type_name != "FILL_AREA_STYLE_COLOUR" {
                    continue;
                }
                let fasc = FillAreaStyleColour::from_params(&fs_entity.params).ok()?;

                return resolve_colour(graph, fasc.fill_colour);
            }
        }
    }
    None
}

/// Resolve a colour entity (COLOUR_RGB or DRAUGHTING_PRE_DEFINED_COLOUR).
fn resolve_colour(graph: &EntityGraph, colour_id: EntityId) -> Option<[f32; 3]> {
    let entity = graph.get(colour_id)?;
    match entity.type_name.as_str() {
        "COLOUR_RGB" => {
            let c = ColourRgb::from_params(&entity.params).ok()?;
            Some([c.red as f32, c.green as f32, c.blue as f32])
        }
        "DRAUGHTING_PRE_DEFINED_COLOUR" => {
            let c = DraughtingPreDefinedColour::from_params(&entity.params).ok()?;
            Some(predefined_colour_rgb(&c.name))
        }
        _ => None,
    }
}

/// Map pre-defined colour names to RGB values.
fn predefined_colour_rgb(name: &str) -> [f32; 3] {
    match name.to_lowercase().as_str() {
        "red" => [1.0, 0.0, 0.0],
        "green" => [0.0, 1.0, 0.0],
        "blue" => [0.0, 0.0, 1.0],
        "yellow" => [1.0, 1.0, 0.0],
        "cyan" => [0.0, 1.0, 1.0],
        "magenta" => [1.0, 0.0, 1.0],
        "white" => [1.0, 1.0, 1.0],
        "black" => [0.0, 0.0, 0.0],
        _ => [0.5, 0.5, 0.5],
    }
}

/// Colour the product that owns styled item `item` (a solid it took
/// on import). A styled item naming a face or another entity no
/// product took as a solid colours nothing: guessing would paint every
/// product of an assembly the first colour in the file.
fn apply_color_to_owner(products: &mut [StepProductNode], item: EntityId, color: [f32; 3]) {
    for product in products.iter_mut() {
        if product.metadata.solid_items.contains(&item) && product.metadata.solid_color.is_none() {
            product.metadata.solid_color = Some(color);
        }
    }
}

/// Map PRESENTATION_LAYER_ASSIGNMENT entities to product layer metadata.
fn map_layer_assignments(
    graph: &EntityGraph,
    products: &mut [StepProductNode],
    _log: &mut DiagnosticLog,
) {
    for (_, entity) in graph.iter() {
        if entity.type_name != "PRESENTATION_LAYER_ASSIGNMENT" {
            continue;
        }
        let pla = match PresentationLayerAssignment::from_params(&entity.params) {
            Ok(p) => p,
            Err(_) => continue,
        };

        // Apply layer name to all products that don't have one yet
        if !pla.assigned_items.is_empty() {
            for product in products.iter_mut() {
                if product.metadata.layer.is_none() {
                    product.metadata.layer = Some(pla.name.clone());
                }
            }
        }
    }
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
    fn p_real(v: f64) -> Parameter {
        Parameter::Real(v)
    }
    fn p_list(items: Vec<Parameter>) -> Parameter {
        Parameter::List(items)
    }

    #[test]
    fn resolve_colour_rgb() {
        let entities = vec![make_entity(
            1,
            "COLOUR_RGB",
            vec![p_str("red"), p_real(1.0), p_real(0.0), p_real(0.0)],
        )];
        let mut log = DiagnosticLog::new();
        let graph = EntityGraph::build(entities, &mut log);
        let c = super::resolve_colour(&graph, 1);
        assert_eq!(c, Some([1.0, 0.0, 0.0]));
    }

    #[test]
    fn resolve_predefined_colour() {
        let entities = vec![make_entity(
            1,
            "DRAUGHTING_PRE_DEFINED_COLOUR",
            vec![p_str("blue")],
        )];
        let mut log = DiagnosticLog::new();
        let graph = EntityGraph::build(entities, &mut log);
        let c = super::resolve_colour(&graph, 1);
        assert_eq!(c, Some([0.0, 0.0, 1.0]));
    }

    #[test]
    fn map_presentation_assigns_color() {
        let entities = vec![
            // Colour chain
            make_entity(
                1,
                "COLOUR_RGB",
                vec![p_str(""), p_real(0.8), p_real(0.2), p_real(0.1)],
            ),
            make_entity(2, "FILL_AREA_STYLE_COLOUR", vec![p_str(""), p_ref(1)]),
            make_entity(
                3,
                "FILL_AREA_STYLE",
                vec![p_str(""), p_list(vec![p_ref(2)])],
            ),
            make_entity(4, "SURFACE_STYLE_FILL_AREA", vec![p_ref(3)]),
            make_entity(
                5,
                "SURFACE_SIDE_STYLE",
                vec![p_str(""), p_list(vec![p_ref(4)])],
            ),
            make_entity(
                6,
                "SURFACE_STYLE_USAGE",
                vec![Parameter::Enum("BOTH".into()), p_ref(5)],
            ),
            make_entity(
                7,
                "PRESENTATION_STYLE_ASSIGNMENT",
                vec![p_list(vec![p_ref(6)])],
            ),
            make_entity(
                8,
                "STYLED_ITEM",
                vec![p_str(""), p_list(vec![p_ref(7)]), p_ref(100)],
            ),
        ];

        let mut log = DiagnosticLog::new();
        let graph = EntityGraph::build(entities, &mut log);

        // The styled item names solid #100, which the first product
        // took; the second product owns another solid and stays plain.
        let product = |name: &str, item: u64| StepProductNode {
            name: name.into(),
            id: name.into(),
            children: vec![],
            shapes: vec![],
            metadata: crate::step::document::ProductMetadata {
                solid_items: vec![item],
                ..Default::default()
            },
        };
        let mut products = vec![product("owner", 100), product("other", 200)];

        map_presentation(&graph, &mut products, &mut log);
        let c = products[0].metadata.solid_color.unwrap();
        assert!((c[0] - 0.8).abs() < 0.01);
        assert!(
            products[1].metadata.solid_color.is_none(),
            "a styled item colours only the product that owns its target"
        );
    }
}
