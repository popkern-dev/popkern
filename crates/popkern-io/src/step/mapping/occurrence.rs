//! Assembly occurrence placement and naming (PK-629): what the writer
//! in `export_assembly` emits, read back.
//!
//! A `NEXT_ASSEMBLY_USAGE_OCCURRENCE` is placed through
//! `PRODUCT_DEFINITION_SHAPE(…, NAUO)` ←
//! `CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(rel, pds)`, where `rel` is a
//! `REPRESENTATION_RELATIONSHIP(…, rep_1, rep_2)` carrying a
//! `REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(op)`. The operator
//! is an `ITEM_DEFINED_TRANSFORMATION(…, item_1, item_2)` (two axis
//! placements; the transform maps `item_1`'s frame onto `item_2`'s) or
//! a `CARTESIAN_TRANSFORMATION_OPERATOR_3D`. By convention `rep_1` is
//! the child's representation; a file that lists the parent first is
//! detected against the child's own representations and inverted.

use std::collections::{HashMap, HashSet};

use popkern_core::geometry::{Point3, Transform3, Vec3};

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::mapping::geometry::GeometryCache;
use crate::step::parser::graph::EntityGraph;
use crate::step::parser::types::{EntityId, Parameter};
use crate::step::schema::product::NextAssemblyUsageOccurrence;

/// The placement wiring of every occurrence in a file, indexed once.
pub(crate) struct OccurrenceIndex {
    /// `PRODUCT_DEFINITION_SHAPE` ids by the entity they describe.
    shapes_of: HashMap<EntityId, Vec<EntityId>>,
    /// The relationship entity of each context-dependent shape
    /// representation, by its `PRODUCT_DEFINITION_SHAPE`.
    relation_of: HashMap<EntityId, EntityId>,
}

impl OccurrenceIndex {
    /// Index `graph`.
    pub(crate) fn build(graph: &EntityGraph) -> OccurrenceIndex {
        let mut shapes_of: HashMap<EntityId, Vec<EntityId>> = HashMap::new();
        let mut relation_of = HashMap::new();
        for (&id, entity) in graph.iter() {
            match entity.type_name.as_str() {
                "PRODUCT_DEFINITION_SHAPE" => {
                    if let Some(Parameter::EntityRef(of)) = entity.params.get(2) {
                        shapes_of.entry(*of).or_default().push(id);
                    }
                }
                "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION" => {
                    if let (Some(Parameter::EntityRef(rel)), Some(Parameter::EntityRef(pds))) =
                        (entity.params.first(), entity.params.get(1))
                    {
                        relation_of.insert(*pds, *rel);
                    }
                }
                _ => {}
            }
        }
        for shapes in shapes_of.values_mut() {
            shapes.sort_unstable();
        }
        OccurrenceIndex {
            shapes_of,
            relation_of,
        }
    }

    /// The transform placing occurrence `nauo_id` in its parent, or
    /// `None` when the file carries no placement for it (logged).
    /// `child_reps` are the child product's own representations.
    pub(crate) fn placement(
        &self,
        graph: &EntityGraph,
        geo: &mut GeometryCache<'_>,
        nauo_id: EntityId,
        child_reps: &HashSet<EntityId>,
        log: &mut DiagnosticLog,
    ) -> Option<Transform3> {
        let relation = self
            .shapes_of
            .get(&nauo_id)
            .into_iter()
            .flatten()
            .find_map(|pds| self.relation_of.get(pds).copied());
        let Some(relation) = relation else {
            log.warn(
                Phase::MapProduct,
                Some(nauo_id),
                "assembly occurrence has no context-dependent shape representation: placed at \
                 the parent origin",
            );
            return None;
        };
        let entity = graph.get(relation)?;
        let partials = std::iter::once((entity.type_name.as_str(), entity.params.as_slice()))
            .chain(
                entity
                    .complex_parts
                    .iter()
                    .map(|(t, p)| (t.as_str(), p.as_slice())),
            );
        let mut reps = None;
        let mut operator = None;
        for (type_name, params) in partials {
            match type_name {
                "REPRESENTATION_RELATIONSHIP" => {
                    if let (Some(Parameter::EntityRef(a)), Some(Parameter::EntityRef(b))) =
                        (params.get(2), params.get(3))
                    {
                        reps = Some((*a, *b));
                    }
                }
                "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION" => {
                    if let Some(Parameter::EntityRef(op)) = params.first() {
                        operator = Some(*op);
                    }
                }
                _ => {}
            }
        }
        let Some(operator) = operator else {
            log.warn(
                Phase::MapProduct,
                Some(relation),
                "occurrence relationship carries no transformation: placed at the parent origin",
            );
            return None;
        };
        let transform = self.operator_transform(graph, geo, operator, log)?;
        // rep_2 being the child's own representation means the file
        // relates parent → child: the operator then maps the parent
        // frame onto the child's, the inverse of a placement.
        let reversed =
            reps.is_some_and(|(a, b)| child_reps.contains(&b) && !child_reps.contains(&a));
        if reversed {
            return transform.inverse().ok();
        }
        Some(transform)
    }

    fn operator_transform(
        &self,
        graph: &EntityGraph,
        geo: &mut GeometryCache<'_>,
        operator: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<Transform3> {
        let entity = graph.get(operator)?;
        match entity.type_name.as_str() {
            "ITEM_DEFINED_TRANSFORMATION" => {
                let (Some(Parameter::EntityRef(a)), Some(Parameter::EntityRef(b))) =
                    (entity.params.get(2), entity.params.get(3))
                else {
                    return None;
                };
                let from = frame(geo.resolve_placement(*a, log)?)?;
                let to = frame(geo.resolve_placement(*b, log)?)?;
                Some(to.compose(&from.inverse().ok()?))
            }
            "CARTESIAN_TRANSFORMATION_OPERATOR_3D" => {
                // (name, description, axis1, axis2, local_origin, scale, axis3)
                let direction =
                    |index: usize, geo: &mut GeometryCache<'_>, log: &mut DiagnosticLog| {
                        match entity.params.get(index) {
                            Some(Parameter::EntityRef(id)) => geo.resolve_direction(*id, log),
                            _ => None,
                        }
                    };
                let x = direction(2, geo, log).unwrap_or(Vec3::new(1.0, 0.0, 0.0));
                let z = direction(6, geo, log).unwrap_or(Vec3::new(0.0, 0.0, 1.0));
                let origin = match entity.params.get(4) {
                    Some(Parameter::EntityRef(id)) => geo.resolve_point(*id, log)?,
                    _ => Point3::new(0.0, 0.0, 0.0),
                };
                Transform3::from_axis_placement(origin, z, x).ok()
            }
            other => {
                log.warn(
                    Phase::MapProduct,
                    Some(operator),
                    format!(
                        "unsupported occurrence transformation {other}: placed at the parent origin"
                    ),
                );
                None
            }
        }
    }
}

fn frame((origin, z, x): (Point3, Vec3, Vec3)) -> Option<Transform3> {
    Transform3::from_axis_placement(origin, z, x).ok()
}

/// The instance name an occurrence carries: its reference designator,
/// else its name (unless it is the generic "Next assembly
/// relationship" some writers put there), else its id.
pub(crate) fn occurrence_name(nauo: &NextAssemblyUsageOccurrence) -> Option<String> {
    let generic =
        |s: &str| s.trim().is_empty() || s.eq_ignore_ascii_case("next assembly relationship");
    nauo.reference_designator
        .as_deref()
        .filter(|s| !s.trim().is_empty())
        .or(Some(nauo.name.as_str()).filter(|s| !generic(s)))
        .or(Some(nauo.id.as_str()).filter(|s| !s.trim().is_empty()))
        .map(str::to_string)
}
