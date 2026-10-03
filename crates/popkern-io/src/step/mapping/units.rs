//! The file's length and plane-angle units (ISO 10303-41 unit model).
//!
//! A geometric representation context lists its units in a
//! `GLOBAL_UNIT_ASSIGNED_CONTEXT`. Each unit is an `SI_UNIT` (with an
//! optional prefix) or a `CONVERSION_BASED_UNIT` whose factor is a
//! measure in another unit (`INCH` = 25.4 mm, `DEGREE` = 0.01745… rad).
//! The reader produces millimetres and radians, so every length in the
//! file is scaled by [`FileUnits::mm_per_length`] and every angle by
//! [`FileUnits::radians_per_angle`].
//!
//! Measured on the CAx-IF corpus (docs/reports/step-assembly-corpus.md):
//! 24 of 34 files declare angles in degrees, so their cone semi-angles
//! were read as radians and every cone face meshed open, and 9 declare
//! inch or centimetre lengths, which imported at the wrong scale.

use std::collections::HashSet;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::parser::graph::EntityGraph;
use crate::step::parser::types::{EntityId, EntityInstance, Parameter};

/// The units the file's geometry is written in, as conversion factors
/// to the kernel's millimetres and radians.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct FileUnits {
    /// Millimetres per file length unit (25.4 for inches).
    pub mm_per_length: f64,
    /// Radians per file plane-angle unit (π/180 for degrees).
    pub radians_per_angle: f64,
}

impl Default for FileUnits {
    fn default() -> Self {
        FileUnits {
            mm_per_length: 1.0,
            radians_per_angle: 1.0,
        }
    }
}

/// Which quantity a unit measures.
#[derive(Debug, Clone, Copy, PartialEq)]
enum UnitKind {
    Length,
    Angle,
    Other,
}

/// The units of the context the file's solids are written in: the
/// context of its first `ADVANCED_BREP_SHAPE_REPRESENTATION` (or any
/// shape representation), else the first context that assigns units.
/// Contexts that disagree are reported; an unresolvable unit keeps the
/// millimetre/radian default and is reported.
pub(crate) fn resolve_file_units(graph: &EntityGraph, log: &mut DiagnosticLog) -> FileUnits {
    let mut contexts: Vec<(EntityId, FileUnits)> = Vec::new();
    let mut ids: Vec<EntityId> = graph.iter().map(|(id, _)| *id).collect();
    ids.sort_unstable();
    for id in &ids {
        let Some(entity) = graph.get(*id) else {
            continue;
        };
        let Some(params) = part(entity, "GLOBAL_UNIT_ASSIGNED_CONTEXT") else {
            continue;
        };
        let Some(Parameter::List(units)) = params.first() else {
            continue;
        };
        let mut found = FileUnits::default();
        for unit in units {
            let Parameter::EntityRef(unit_id) = unit else {
                continue;
            };
            match unit_factor(graph, *unit_id, &mut HashSet::new()) {
                Some((UnitKind::Length, f)) => found.mm_per_length = f,
                Some((UnitKind::Angle, f)) => found.radians_per_angle = f,
                Some((UnitKind::Other, _)) => {}
                None => log.warn(
                    Phase::MapGeometry,
                    Some(*unit_id),
                    "unit not resolved: assuming millimetres and radians for it",
                ),
            }
        }
        contexts.push((*id, found));
    }
    if contexts.is_empty() {
        return FileUnits::default();
    }
    let solid_context = ids.iter().filter_map(|id| graph.get(*id)).find_map(|e| {
        (e.type_name == "ADVANCED_BREP_SHAPE_REPRESENTATION").then(|| match e.params.get(2) {
            Some(Parameter::EntityRef(c)) => Some(*c),
            _ => None,
        })?
    });
    let chosen = solid_context
        .and_then(|c| contexts.iter().find(|(id, _)| *id == c))
        .unwrap_or(&contexts[0])
        .1;
    if contexts.iter().any(|(_, u)| *u != chosen) {
        log.warn(
            Phase::MapGeometry,
            None,
            format!(
                "representation contexts declare different units; reading all geometry in the \
                 solids' context ({} mm per length unit, {} rad per angle unit)",
                chosen.mm_per_length, chosen.radians_per_angle
            ),
        );
    }
    if chosen != FileUnits::default() {
        log.warn(
            Phase::MapGeometry,
            None,
            format!(
                "file units: {} mm per length unit, {} rad per angle unit (converted on import)",
                chosen.mm_per_length, chosen.radians_per_angle
            ),
        );
    }
    chosen
}

/// Millimetres per unit for a LENGTH unit entity (`None` for any other
/// kind or an unresolvable unit).
pub(crate) fn mm_per_length_unit(graph: &EntityGraph, unit: EntityId) -> Option<f64> {
    match unit_factor(graph, unit, &mut HashSet::new())? {
        (UnitKind::Length, f) => Some(f),
        _ => None,
    }
}

/// The parameters of `name` in a simple or complex entity.
fn part<'e>(entity: &'e EntityInstance, name: &str) -> Option<&'e [Parameter]> {
    if entity.type_name == name {
        return Some(&entity.params);
    }
    entity
        .complex_parts
        .iter()
        .find(|(n, _)| n == name)
        .map(|(_, p)| p.as_slice())
}

fn has(entity: &EntityInstance, name: &str) -> bool {
    part(entity, name).is_some()
}

/// A unit's kind and its factor to millimetres or radians.
fn unit_factor(
    graph: &EntityGraph,
    id: EntityId,
    seen: &mut HashSet<EntityId>,
) -> Option<(UnitKind, f64)> {
    if !seen.insert(id) {
        return None;
    }
    let entity = graph.get(id)?;
    let kind = if has(entity, "LENGTH_UNIT") {
        UnitKind::Length
    } else if has(entity, "PLANE_ANGLE_UNIT") {
        UnitKind::Angle
    } else {
        UnitKind::Other
    };
    if let Some(si) = part(entity, "SI_UNIT") {
        // (dimensions, prefix, name) as a simple instance; (prefix,
        // name) in a complex one, where NAMED_UNIT holds the dimensions.
        let si = &si[si.len().saturating_sub(2)..];
        let prefix = match si.first() {
            Some(Parameter::Enum(p)) => si_prefix(p)?,
            _ => 1.0,
        };
        let base = match si.get(1) {
            Some(Parameter::Enum(n)) if n == "METRE" => 1000.0,
            Some(Parameter::Enum(n)) if n == "RADIAN" => 1.0,
            Some(Parameter::Enum(_)) => 1.0,
            _ => return None,
        };
        let kind = match (kind, si.get(1)) {
            (UnitKind::Other, Some(Parameter::Enum(n))) if n == "METRE" => UnitKind::Length,
            (UnitKind::Other, Some(Parameter::Enum(n))) if n == "RADIAN" => UnitKind::Angle,
            (k, _) => k,
        };
        return Some((kind, prefix * base));
    }
    if let Some(conv) = part(entity, "CONVERSION_BASED_UNIT") {
        let Some(Parameter::EntityRef(measure_id)) = conv.get(1) else {
            return None;
        };
        let measure = graph.get(*measure_id)?;
        let value = measure_value(measure.params.first()?)?;
        let Some(Parameter::EntityRef(base_id)) = measure.params.get(1) else {
            return None;
        };
        let (base_kind, base) = unit_factor(graph, *base_id, seen)?;
        let kind = if kind == UnitKind::Other {
            base_kind
        } else {
            kind
        };
        let factor = value * base;
        return (factor.is_finite() && factor > 0.0).then_some((kind, factor));
    }
    None
}

/// A measure's numeric value: plain, or typed (`LENGTH_MEASURE(25.4)`).
fn measure_value(param: &Parameter) -> Option<f64> {
    match param {
        Parameter::Real(v) => Some(*v),
        Parameter::Int(v) => Some(*v as f64),
        Parameter::Typed(_, inner) => measure_value(inner),
        Parameter::List(items) if items.len() == 1 => measure_value(&items[0]),
        _ => None,
    }
}

fn si_prefix(prefix: &str) -> Option<f64> {
    Some(match prefix {
        "EXA" => 1e18,
        "PETA" => 1e15,
        "TERA" => 1e12,
        "GIGA" => 1e9,
        "MEGA" => 1e6,
        "KILO" => 1e3,
        "HECTO" => 1e2,
        "DECA" => 1e1,
        "DECI" => 1e-1,
        "CENTI" => 1e-2,
        "MILLI" => 1e-3,
        "MICRO" => 1e-6,
        "NANO" => 1e-9,
        "PICO" => 1e-12,
        "FEMTO" => 1e-15,
        "ATTO" => 1e-18,
        _ => return None,
    })
}

/// Convert a mapped document from the file's length unit to
/// millimetres: every solid scaled by `mm_per_length` (radii, pcurve
/// lengths and tolerances with it) and every occurrence placement's
/// translation scaled the same way (`S ∘ T ∘ S⁻¹`), so shared
/// representations stay shared.
pub(crate) fn convert_lengths_to_mm(
    products: &mut [crate::step::document::StepProductNode],
    mm_per_length: f64,
    log: &mut DiagnosticLog,
) {
    use popkern_core::geometry::Transform3;
    use std::collections::HashMap;
    use std::sync::Arc;

    if mm_per_length == 1.0 {
        return;
    }
    let scale = Transform3::scale_uniform(mm_per_length);
    let unscale = Transform3::scale_uniform(1.0 / mm_per_length);
    let mut converted: HashMap<
        *const crate::step::document::ShapeRep,
        Arc<crate::step::document::ShapeRep>,
    > = HashMap::new();
    for product in products.iter_mut() {
        for shape in product.shapes.iter_mut() {
            let key = Arc::as_ptr(shape);
            if let Some(done) = converted.get(&key) {
                *shape = done.clone();
                continue;
            }
            let mut rep = (**shape).clone();
            for solid in rep.solids.iter_mut() {
                match scale_solid(solid, &scale, mm_per_length) {
                    Ok(scaled) => *solid = scaled,
                    Err(e) => log.error(
                        Phase::MapGeometry,
                        None,
                        format!("length unit conversion failed for a solid: {e}"),
                    ),
                }
            }
            let rep = Arc::new(rep);
            converted.insert(key, rep.clone());
            *shape = rep;
        }
        for child in product.children.iter_mut() {
            child.transform = scale.compose(&child.transform).compose(&unscale);
        }
    }
}

fn scale_solid(
    solid: &popkern_core::topology::BRep,
    scale: &popkern_core::geometry::Transform3,
    k: f64,
) -> Result<popkern_core::topology::BRep, popkern_core::error::CoreError> {
    let mut out = solid.transform(scale)?;
    out.tolerance = solid.tolerance.scale(k);
    out.max_imprecision = solid.max_imprecision * k;
    for (_, v) in out.vertices.iter_mut() {
        v.tolerance = v.tolerance.map(|t| t * k);
    }
    for (_, e) in out.edges.iter_mut() {
        e.tolerance = e.tolerance.map(|t| t * k);
    }
    for (_, f) in out.faces.iter_mut() {
        f.tolerance = f.tolerance.map(|t| t * k);
    }
    Ok(out)
}

#[cfg(test)]
#[path = "units_tests.rs"]
mod tests;
