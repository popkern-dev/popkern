//! `SHELL_BASED_SURFACE_MODEL` items: one sheet body per listed shell
//! (split from `mod.rs` under the file-size budget).

use popkern_core::tolerance::Tolerance;
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
                Some(brep) => results.push((item_id, brep)),
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
