//! The file's geometric uncertainty (PK-639), split from `import.rs`
//! (500-line cap).

use crate::step::parser::graph::EntityGraph;
use crate::step::parser::types::Parameter;
use crate::step::schema::product::UncertaintyMeasureWithUnit;

/// The file's LENGTH uncertainty in file length units (PK-639): the
/// first `UNCERTAINTY_MEASURE_WITH_UNIT` that is a length (typed
/// `LENGTH_MEASURE`, or in a length unit) in the lowest-numbered
/// context assigning one, converted from its own unit. Angular
/// uncertainties (`PLANE_ANGLE_MEASURE`, `SOLID_ANGLE_MEASURE`) are
/// skipped. Before PK-639 the typed value never parsed, so every file
/// fell back to the default.
pub(crate) fn extract_step_uncertainty(graph: &EntityGraph, mm_per_length: f64) -> Option<f64> {
    let mut ids: Vec<u64> = graph.iter().map(|(id, _)| *id).collect();
    ids.sort_unstable();
    for id in ids {
        let Some(entity) = graph.get(id) else {
            continue;
        };
        let parts = std::iter::once((&entity.type_name, &entity.params))
            .chain(entity.complex_parts.iter().map(|(n, p)| (n, p)));
        for (part_type, part_params) in parts {
            if part_type != "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT" {
                continue;
            }
            // A complex instance carries only this partial's own
            // attribute (the list); a standalone one also carries the
            // inherited context attributes before it.
            let refs: Vec<u64> = part_params
                .iter()
                .find_map(|p| match p {
                    Parameter::List(items) => Some(
                        items
                            .iter()
                            .filter_map(|i| match i {
                                Parameter::EntityRef(r) => Some(*r),
                                _ => None,
                            })
                            .collect(),
                    ),
                    _ => None,
                })
                .unwrap_or_default();
            if let Some(u) = resolve_uncertainty_value(graph, &refs, mm_per_length) {
                return Some(u);
            }
        }
    }
    None
}

/// The first length uncertainty among `refs`, in file length units.
fn resolve_uncertainty_value(graph: &EntityGraph, refs: &[u64], mm_per_length: f64) -> Option<f64> {
    for &ref_id in refs {
        let Some(ent) = graph.get(ref_id) else {
            continue;
        };
        let Ok(umu) = UncertaintyMeasureWithUnit::from_params(&ent.params) else {
            continue;
        };
        let unit_mm = crate::step::mapping::units::mm_per_length_unit(graph, umu.unit);
        let is_length = umu
            .measure
            .as_deref()
            .map_or(unit_mm.is_some(), |m| m == "LENGTH_MEASURE");
        if !is_length || !(umu.value > 0.0 && umu.value.is_finite()) {
            continue;
        }
        let value = umu.value * unit_mm.unwrap_or(mm_per_length) / mm_per_length;
        return Some(value);
    }
    None
}
