//! `SHELL_BASED_SURFACE_MODEL` items: one sheet body per listed shell
//! (split from `mod.rs` under the file-size budget), a closed shell a
//! sheet around a void.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::RegionKind;
use popkern_core::topology::brep::BRep;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::mapping::geometry::GeometryCache;
use crate::step::parser::graph::EntityGraph;
use crate::step::parser::types::EntityId;
use crate::step::schema::topology::ShellBasedSurfaceModel;

/// Build a sheet body for every shell of every surface model.
pub(super) fn map_surface_models(
    graph: &EntityGraph,
    geo: &mut GeometryCache<'_>,
    log: &mut DiagnosticLog,
    tolerance: Tolerance,
) -> Vec<(EntityId, BRep)> {
    let mut results = Vec::new();
    let surface_model_ids: Vec<EntityId> = graph
        .iter()
        .filter(|(_, e)| e.type_name == "SHELL_BASED_SURFACE_MODEL")
        .map(|(&id, _)| id)
        .collect();

    for item_id in surface_model_ids {
        let Some(entity) = graph.get(item_id) else {
            continue;
        };
        let model = match ShellBasedSurfaceModel::from_params(&entity.params) {
            Ok(m) => m,
            Err(e) => {
                log.error(Phase::MapTopology, Some(item_id), format!("{e}"));
                continue;
            }
        };

        for shell_id in model.shells {
            match super::build_shell_body(graph, geo, shell_id, true, log, tolerance) {
                Some(mut brep) => {
                    // A surface model's closed shell is a closed sheet: what
                    // it encloses is a void, not material (PK-1114).
                    if graph
                        .get(shell_id)
                        .is_some_and(|e| e.type_name == "CLOSED_SHELL")
                    {
                        for region in brep.solid_regions() {
                            brep.regions[region].kind = RegionKind::Void;
                        }
                    }
                    results.push((item_id, brep))
                }
                None => log.error(
                    Phase::MapTopology,
                    Some(shell_id),
                    "failed to build sheet body from SHELL_BASED_SURFACE_MODEL",
                ),
            }
        }
    }
    results
}
