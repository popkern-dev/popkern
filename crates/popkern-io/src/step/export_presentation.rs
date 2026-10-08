//! Solid colours on export (AP214 cc2 presentation): each solid of a
//! product with a `solid_color` gets a `STYLED_ITEM` through the
//! standard chain
//! `PRESENTATION_STYLE_ASSIGNMENT → SURFACE_STYLE_USAGE(.BOTH.) →
//! SURFACE_SIDE_STYLE → SURFACE_STYLE_FILL_AREA → FILL_AREA_STYLE →
//! FILL_AREA_STYLE_COLOUR → COLOUR_RGB`, and the styled items are
//! collected in one `MECHANICAL_DESIGN_GEOMETRIC_PRESENTATION_REPRESENTATION`
//! over the shape context — the form XCAF and the import side
//! (`mapping::presentation`) read. A document without colours writes
//! nothing here, byte-identical to before.

use crate::step::document::StepDocument;
use crate::step::export::{ExportContext, enum_, list, real, ref_, str_param};
use crate::step::parser::types::EntityId;

impl ExportContext {
    /// Style every coloured product's solid items. `shape_item_ids` is
    /// parallel to `doc.products`; `ctx_id` the geometric context.
    pub(crate) fn export_solid_colours(
        &mut self,
        doc: &StepDocument,
        shape_item_ids: &[Vec<EntityId>],
        ctx_id: EntityId,
    ) {
        let mut styled = Vec::new();
        for (product, items) in doc.products.iter().zip(shape_item_ids) {
            let Some([r, g, b]) = product.metadata.solid_color else {
                continue;
            };
            if items.is_empty() {
                continue;
            }
            let colour = self.emit(
                "COLOUR_RGB",
                vec![
                    str_param(""),
                    real(f64::from(r)),
                    real(f64::from(g)),
                    real(f64::from(b)),
                ],
            );
            let fill_colour =
                self.emit("FILL_AREA_STYLE_COLOUR", vec![str_param(""), ref_(colour)]);
            let fill = self.emit(
                "FILL_AREA_STYLE",
                vec![str_param(""), list(vec![ref_(fill_colour)])],
            );
            let area = self.emit("SURFACE_STYLE_FILL_AREA", vec![ref_(fill)]);
            let side = self.emit(
                "SURFACE_SIDE_STYLE",
                vec![str_param(""), list(vec![ref_(area)])],
            );
            let usage = self.emit("SURFACE_STYLE_USAGE", vec![enum_("BOTH"), ref_(side)]);
            let assignment = self.emit(
                "PRESENTATION_STYLE_ASSIGNMENT",
                vec![list(vec![ref_(usage)])],
            );
            for &item in items {
                styled.push(self.emit(
                    "STYLED_ITEM",
                    vec![str_param("color"), list(vec![ref_(assignment)]), ref_(item)],
                ));
            }
        }
        if styled.is_empty() {
            return;
        }
        self.emit(
            "MECHANICAL_DESIGN_GEOMETRIC_PRESENTATION_REPRESENTATION",
            vec![
                str_param(""),
                list(styled.into_iter().map(ref_).collect()),
                ref_(ctx_id),
            ],
        );
    }
}
