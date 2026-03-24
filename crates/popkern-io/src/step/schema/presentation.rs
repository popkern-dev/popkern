//! Typed presentation and styling entity definitions from the STEP schema.

use super::{MappingError, as_entity_ref, as_enum, as_real, as_ref_list, as_string, param_at};
use crate::step::parser::types::{EntityId, Parameter};

/// STYLED_ITEM — an item with applied presentation styles.
#[derive(Debug, Clone)]
pub struct StyledItem {
    /// Label name.
    pub name: String,
    /// References to `PRESENTATION_STYLE_ASSIGNMENT` entities.
    pub styles: Vec<EntityId>,
    /// Reference to the styled geometric item.
    pub item: EntityId,
}

impl StyledItem {
    const TYPE: &str = "STYLED_ITEM";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            styles: as_ref_list(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            item: as_entity_ref(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// PRESENTATION_STYLE_ASSIGNMENT — a collection of presentation styles.
#[derive(Debug, Clone)]
pub struct PresentationStyleAssignment {
    /// References to style entities.
    pub styles: Vec<EntityId>,
}

impl PresentationStyleAssignment {
    const TYPE: &str = "PRESENTATION_STYLE_ASSIGNMENT";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            styles: as_ref_list(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// SURFACE_STYLE_USAGE — binding of a surface style to a side of a surface.
#[derive(Debug, Clone)]
pub struct SurfaceStyleUsage {
    /// Which side the style applies to (e.g. `BOTH`).
    pub side: String,
    /// Reference to the `SURFACE_SIDE_STYLE`.
    pub style: EntityId,
}

impl SurfaceStyleUsage {
    const TYPE: &str = "SURFACE_STYLE_USAGE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            side: as_enum(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            style: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// SURFACE_SIDE_STYLE — styles applied to one side of a surface.
#[derive(Debug, Clone)]
pub struct SurfaceSideStyle {
    /// Label name.
    pub name: String,
    /// References to surface style element entities.
    pub styles: Vec<EntityId>,
}

impl SurfaceSideStyle {
    const TYPE: &str = "SURFACE_SIDE_STYLE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            styles: as_ref_list(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// SURFACE_STYLE_FILL_AREA — a surface style defined by a fill area.
#[derive(Debug, Clone)]
pub struct SurfaceStyleFillArea {
    /// Reference to the `FILL_AREA_STYLE`.
    pub fill_area: EntityId,
}

impl SurfaceStyleFillArea {
    const TYPE: &str = "SURFACE_STYLE_FILL_AREA";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            fill_area: as_entity_ref(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// FILL_AREA_STYLE — a style that fills a region.
#[derive(Debug, Clone)]
pub struct FillAreaStyle {
    /// Label name.
    pub name: String,
    /// References to fill style entities.
    pub fill_styles: Vec<EntityId>,
}

impl FillAreaStyle {
    const TYPE: &str = "FILL_AREA_STYLE";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            fill_styles: as_ref_list(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// FILL_AREA_STYLE_COLOUR — a fill-area style defined by a single colour.
#[derive(Debug, Clone)]
pub struct FillAreaStyleColour {
    /// Label name.
    pub name: String,
    /// Reference to the colour entity.
    pub fill_colour: EntityId,
}

impl FillAreaStyleColour {
    const TYPE: &str = "FILL_AREA_STYLE_COLOUR";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            fill_colour: as_entity_ref(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// COLOUR_RGB — a colour defined by red, green, blue components in \[0, 1\].
#[derive(Debug, Clone)]
pub struct ColourRgb {
    /// Label name.
    pub name: String,
    /// Red component.
    pub red: f64,
    /// Green component.
    pub green: f64,
    /// Blue component.
    pub blue: f64,
}

impl ColourRgb {
    const TYPE: &str = "COLOUR_RGB";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            red: as_real(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            green: as_real(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
            blue: as_real(param_at(params, 3, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// DRAUGHTING_PRE_DEFINED_COLOUR — a pre-defined colour referenced by name.
#[derive(Debug, Clone)]
pub struct DraughtingPreDefinedColour {
    /// Colour name (e.g. `"red"`, `"green"`).
    pub name: String,
}

impl DraughtingPreDefinedColour {
    const TYPE: &str = "DRAUGHTING_PRE_DEFINED_COLOUR";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

/// PRESENTATION_LAYER_ASSIGNMENT — assignment of items to a presentation layer.
#[derive(Debug, Clone)]
pub struct PresentationLayerAssignment {
    /// Layer name.
    pub name: String,
    /// Layer description.
    pub description: String,
    /// References to the assigned items.
    pub assigned_items: Vec<EntityId>,
}

impl PresentationLayerAssignment {
    const TYPE: &str = "PRESENTATION_LAYER_ASSIGNMENT";

    /// Parse from raw STEP parameters.
    pub fn from_params(params: &[Parameter]) -> Result<Self, MappingError> {
        Ok(Self {
            name: as_string(param_at(params, 0, Self::TYPE)?, Self::TYPE)?,
            description: as_string(param_at(params, 1, Self::TYPE)?, Self::TYPE)?,
            assigned_items: as_ref_list(param_at(params, 2, Self::TYPE)?, Self::TYPE)?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p_str(s: &str) -> Parameter {
        Parameter::String(s.into())
    }
    fn p_ref(id: u64) -> Parameter {
        Parameter::EntityRef(id)
    }
    fn p_real(v: f64) -> Parameter {
        Parameter::Real(v)
    }

    #[test]
    fn colour_rgb_parse() {
        let params = vec![p_str("red"), p_real(1.0), p_real(0.0), p_real(0.0)];
        let c = ColourRgb::from_params(&params).unwrap();
        assert_eq!(c.name, "red");
        assert!((c.red - 1.0).abs() < f64::EPSILON);
        assert!((c.green).abs() < f64::EPSILON);
    }

    #[test]
    fn styled_item_parse() {
        let params = vec![
            p_str("style1"),
            Parameter::List(vec![p_ref(10), p_ref(11)]),
            p_ref(20),
        ];
        let si = StyledItem::from_params(&params).unwrap();
        assert_eq!(si.name, "style1");
        assert_eq!(si.styles, vec![10, 11]);
        assert_eq!(si.item, 20);
    }

    #[test]
    fn presentation_layer_assignment_parse() {
        let params = vec![
            p_str("Layer0"),
            p_str("default layer"),
            Parameter::List(vec![p_ref(1), p_ref(2), p_ref(3)]),
        ];
        let pla = PresentationLayerAssignment::from_params(&params).unwrap();
        assert_eq!(pla.name, "Layer0");
        assert_eq!(pla.assigned_items, vec![1, 2, 3]);
    }
}
