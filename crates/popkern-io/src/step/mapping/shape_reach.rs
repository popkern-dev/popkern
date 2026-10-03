//! The solids a product's shape representation reaches (PK-635).
//!
//! A product's `SHAPE_DEFINITION_REPRESENTATION` names one
//! representation. Writers put the solids in it directly, or in another
//! representation linked to it:
//!
//! - **`SHAPE_REPRESENTATION_RELATIONSHIP`** (plain, no transformation):
//!   the product's `SHAPE_REPRESENTATION` → the
//!   `ADVANCED_BREP_SHAPE_REPRESENTATION` holding the solids (UG, NX,
//!   SolidWorks, Creo style). A relationship that also carries a
//!   `REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION` is an assembly
//!   OCCURRENCE (mapping::occurrence), never followed here — that would
//!   pull a child's solids into its parent.
//! - **`MAPPED_ITEM`**: an item standing for another representation
//!   (`REPRESENTATION_MAP(origin, mapped)`) placed at a target frame;
//!   its solids are moved by `target ∘ origin⁻¹`.
//!
//! The walk never enters a representation another product names as
//! its own, so products sharing a geometry representation through a
//! link do not swallow each other.

use std::collections::{HashMap, HashSet};

use popkern_core::geometry::Transform3;
use popkern_core::topology::BRep;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::mapping::geometry::GeometryCache;
use crate::step::parser::graph::EntityGraph;
use crate::step::parser::types::{EntityId, Parameter};

/// Plain shape-representation links, both directions.
pub(crate) struct ShapeReach {
    links: HashMap<EntityId, Vec<EntityId>>,
}

impl ShapeReach {
    /// Index every plain `SHAPE_REPRESENTATION_RELATIONSHIP`.
    pub(crate) fn build(graph: &EntityGraph) -> ShapeReach {
        let mut links: HashMap<EntityId, Vec<EntityId>> = HashMap::new();
        for (_, entity) in graph.iter() {
            let plain_srr = entity.type_name == "SHAPE_REPRESENTATION_RELATIONSHIP"
                && !entity
                    .complex_parts
                    .iter()
                    .any(|(t, _)| t == "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION");
            if !plain_srr {
                continue;
            }
            if let (Some(Parameter::EntityRef(a)), Some(Parameter::EntityRef(b))) =
                (entity.params.get(2), entity.params.get(3))
            {
                links.entry(*a).or_default().push(*b);
                links.entry(*b).or_default().push(*a);
            }
        }
        for v in links.values_mut() {
            v.sort_unstable();
            v.dedup();
        }
        ShapeReach { links }
    }

    /// Every solid reachable from `start`, not entering `foreign`
    /// representations (other products' own).
    pub(crate) fn solids(
        &self,
        graph: &EntityGraph,
        geo: &mut GeometryCache<'_>,
        start: EntityId,
        foreign: &HashSet<EntityId>,
        brep_by_id: &HashMap<EntityId, Vec<&BRep>>,
        log: &mut DiagnosticLog,
    ) -> Vec<BRep> {
        let mut out = Vec::new();
        let mut seen = HashSet::new();
        self.collect(
            graph,
            geo,
            start,
            foreign,
            brep_by_id,
            &Transform3::identity(),
            &mut seen,
            &mut out,
            log,
        );
        out
    }

    #[allow(clippy::too_many_arguments)]
    fn collect(
        &self,
        graph: &EntityGraph,
        geo: &mut GeometryCache<'_>,
        rep: EntityId,
        foreign: &HashSet<EntityId>,
        brep_by_id: &HashMap<EntityId, Vec<&BRep>>,
        placed: &Transform3,
        seen: &mut HashSet<EntityId>,
        out: &mut Vec<BRep>,
        log: &mut DiagnosticLog,
    ) {
        if !seen.insert(rep) {
            return;
        }
        let Some(entity) = graph.get(rep) else { return };
        let items =
            crate::step::schema::product::StepShapeRepresentation::from_params(&entity.params)
                .map(|r| r.items)
                .unwrap_or_default();
        for item in items {
            if let Some(bodies) = brep_by_id.get(&item) {
                for body in bodies {
                    out.push(place(body, placed, item, log));
                }
                continue;
            }
            if let Some((mapped, transform)) = mapped_item(graph, geo, item, log) {
                if foreign.contains(&mapped) {
                    continue;
                }
                let inner = placed.compose(&transform);
                self.collect(
                    graph, geo, mapped, foreign, brep_by_id, &inner, seen, out, log,
                );
            }
        }
        for &next in self.links.get(&rep).into_iter().flatten() {
            if !foreign.contains(&next) {
                self.collect(
                    graph, geo, next, foreign, brep_by_id, placed, seen, out, log,
                );
            }
        }
    }
}

fn place(body: &BRep, placed: &Transform3, item: EntityId, log: &mut DiagnosticLog) -> BRep {
    if *placed == Transform3::identity() {
        return body.clone();
    }
    match body.transform(placed) {
        Ok(moved) => moved,
        Err(e) => {
            log.warn(
                Phase::MapProduct,
                Some(item),
                format!("mapped item placement not applied: {e}"),
            );
            body.clone()
        }
    }
}

/// `MAPPED_ITEM(name, REPRESENTATION_MAP(origin, mapped), target)`: the
/// mapped representation and the transform `target ∘ origin⁻¹`.
fn mapped_item(
    graph: &EntityGraph,
    geo: &mut GeometryCache<'_>,
    item: EntityId,
    log: &mut DiagnosticLog,
) -> Option<(EntityId, Transform3)> {
    let entity = graph.get(item)?;
    if entity.type_name != "MAPPED_ITEM" {
        return None;
    }
    let (Some(Parameter::EntityRef(map_id)), Some(Parameter::EntityRef(target))) =
        (entity.params.get(1), entity.params.get(2))
    else {
        return None;
    };
    let map = graph.get(*map_id)?;
    let (Some(Parameter::EntityRef(origin)), Some(Parameter::EntityRef(mapped))) =
        (map.params.first(), map.params.get(1))
    else {
        return None;
    };
    let frame = |id: EntityId, geo: &mut GeometryCache<'_>, log: &mut DiagnosticLog| {
        let (o, z, x) = geo.resolve_placement(id, log)?;
        Transform3::from_axis_placement(o, z, x).ok()
    };
    let from = frame(*origin, geo, log)?;
    let to = frame(*target, geo, log)?;
    Some((*mapped, to.compose(&from.inverse().ok()?)))
}
