//! STEP shell edge-incidence policy before BRep builder wiring.

use std::collections::HashMap;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::parser::types::EntityId;

use super::face_collector::FaceData;

pub(super) struct ShellIncidencePolicy {
    pub treat_as_open: bool,
    pub downgraded_closed_shell: bool,
}

pub(super) fn open_shell_incidence_policy() -> ShellIncidencePolicy {
    ShellIncidencePolicy {
        treat_as_open: true,
        downgraded_closed_shell: false,
    }
}

pub(super) fn closed_shell_incidence_policy(
    face_data: &FaceData,
    log: &mut DiagnosticLog,
) -> Option<ShellIncidencePolicy> {
    let mut use_counts: HashMap<EntityId, usize> = HashMap::new();
    for face in &face_data.faces {
        for bound in &face.bounds {
            for oriented_edge in &bound.oriented_edges {
                *use_counts.entry(oriented_edge.edge_curve_id).or_default() += 1;
            }
        }
    }

    let mut entries: Vec<_> = use_counts.into_iter().collect();
    entries.sort_by_key(|(edge_curve_id, _)| *edge_curve_id);

    let mut downgraded_closed_shell = false;
    let mut declined = false;
    for (edge_curve_id, uses) in entries {
        match uses {
            2 => {}
            1 => {
                downgraded_closed_shell = true;
                log.warn(
                    Phase::MapTopology,
                    Some(edge_curve_id),
                    "STEP CLOSED_SHELL edge incidence is single-sided; downgrading imported shell to open boundary topology",
                );
            }
            _ => {
                declined = true;
                log.error(
                    Phase::MapTopology,
                    Some(edge_curve_id),
                    format!(
                        "STEP CLOSED_SHELL edge incidence is non-manifold; edge has {uses} face uses, expected 2"
                    ),
                );
            }
        }
    }

    if declined {
        return None;
    }

    Some(ShellIncidencePolicy {
        treat_as_open: downgraded_closed_shell,
        downgraded_closed_shell,
    })
}
