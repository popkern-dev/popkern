//! `BREP_WITH_VOIDS` items: a solid with cavities (PK-880).
//!
//! The outer `CLOSED_SHELL` and each void become shells of ONE popkern
//! solid, the cavity model the kernel serves. A void is an
//! `ORIENTED_CLOSED_SHELL` (or, leniently, a bare `CLOSED_SHELL`); with
//! orientation `.F.` (the schema's rule) its shell is oriented as the
//! void region's outside, so its faces are used reversed and face into
//! the void — out of the material, as every outer face does. These
//! items were ignored before: a hollow part imported as nothing.

use std::collections::HashSet;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::brep::BRep;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::mapping::geometry::GeometryCache;
use crate::step::parser::graph::EntityGraph;
use crate::step::parser::types::{EntityId, Parameter};
use crate::step::schema::topology::{BrepWithVoids, ClosedShell};

/// Build one body per `BREP_WITH_VOIDS` item.
pub(super) fn map_solids_with_voids(
    graph: &EntityGraph,
    geo: &mut GeometryCache<'_>,
    log: &mut DiagnosticLog,
    tolerance: Tolerance,
) -> Vec<(EntityId, BRep)> {
    let ids: Vec<EntityId> = graph
        .iter()
        .filter(|(_, e)| e.type_name == "BREP_WITH_VOIDS")
        .map(|(&id, _)| id)
        .collect();
    let mut results = Vec::new();
    for id in ids {
        match build_solid_with_voids(graph, geo, id, log, tolerance) {
            Some(brep) => results.push((id, brep)),
            None => log.error(
                Phase::MapTopology,
                Some(id),
                "failed to build solid from BREP_WITH_VOIDS",
            ),
        }
    }
    results
}

/// The faces a `CLOSED_SHELL` lists.
fn shell_faces(
    graph: &EntityGraph,
    shell: EntityId,
    log: &mut DiagnosticLog,
) -> Option<Vec<EntityId>> {
    let entity = graph.get(shell)?;
    if entity.type_name != "CLOSED_SHELL" {
        log.error(
            Phase::MapTopology,
            Some(shell),
            format!("expected CLOSED_SHELL, got {}", entity.type_name),
        );
        return None;
    }
    match ClosedShell::from_params(&entity.params) {
        Ok(s) => Some(s.faces),
        Err(e) => {
            log.error(Phase::MapTopology, Some(shell), format!("{e}"));
            None
        }
    }
}

/// A void reference's closed shell and whether it is used reversed.
fn void_shell(
    graph: &EntityGraph,
    void: EntityId,
    log: &mut DiagnosticLog,
) -> Option<(EntityId, bool)> {
    let entity = graph.get(void)?;
    match entity.type_name.as_str() {
        "ORIENTED_CLOSED_SHELL" => match (entity.params.get(2), entity.params.get(3)) {
            (Some(Parameter::EntityRef(inner)), orientation) => Some((
                *inner,
                matches!(orientation, Some(Parameter::Enum(v)) if v == "F"),
            )),
            _ => {
                log.error(
                    Phase::MapTopology,
                    Some(void),
                    "ORIENTED_CLOSED_SHELL without a closed shell element",
                );
                None
            }
        },
        "CLOSED_SHELL" => Some((void, false)),
        other => {
            log.error(
                Phase::MapTopology,
                Some(void),
                format!("expected a void shell, got {other}"),
            );
            None
        }
    }
}

fn build_solid_with_voids(
    graph: &EntityGraph,
    geo: &mut GeometryCache<'_>,
    id: EntityId,
    log: &mut DiagnosticLog,
    tolerance: Tolerance,
) -> Option<BRep> {
    let item = match BrepWithVoids::from_params(&graph.get(id)?.params) {
        Ok(b) => b,
        Err(e) => {
            log.error(Phase::MapTopology, Some(id), format!("{e}"));
            return None;
        }
    };
    let mut faces = shell_faces(graph, item.outer, log)?;
    let mut sizes = vec![faces.len()];
    let mut reversed = HashSet::new();
    for void in item.voids {
        let (shell, flip) = void_shell(graph, void, log)?;
        let listed = shell_faces(graph, shell, log)?;
        if flip {
            reversed.extend(listed.iter().copied());
        }
        sizes.push(listed.len());
        faces.extend(listed);
    }
    let layout = super::ShellLayout {
        sizes: &sizes,
        reversed_uses: &reversed,
    };
    super::build_brep_from_faces(graph, geo, &faces, layout, false, log, tolerance)
}
