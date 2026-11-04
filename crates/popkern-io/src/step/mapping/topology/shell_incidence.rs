//! STEP shell edge-incidence policy before BRep builder wiring.

use std::collections::BTreeMap;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::parser::types::EntityId;

use super::face_collector::FaceData;
use super::shell_pieces::{ShellPieces, touching_pieces};

pub(super) struct ShellIncidencePolicy {
    pub treat_as_open: bool,
    pub downgraded_closed_shell: bool,
    /// The manifold pieces of a closed shell that touches itself
    /// (`shell_pieces`); `None` for one piece.
    pub pieces: Option<ShellPieces>,
}

pub(super) fn open_shell_incidence_policy() -> ShellIncidencePolicy {
    ShellIncidencePolicy {
        treat_as_open: true,
        downgraded_closed_shell: false,
        pieces: None,
    }
}

/// The policy for a `CLOSED_SHELL`'s faces. A shell that touches itself
/// is split into manifold pieces when `allow_pieces` (a body of one
/// shell) and declines otherwise, as does one whose contact edges no
/// split pairs.
pub(super) fn closed_shell_incidence_policy(
    face_data: &FaceData,
    allow_pieces: bool,
    log: &mut DiagnosticLog,
) -> Option<ShellIncidencePolicy> {
    let mut use_counts: BTreeMap<EntityId, usize> = BTreeMap::new();
    for face in &face_data.faces {
        for bound in &face.bounds {
            for oriented_edge in &bound.oriented_edges {
                *use_counts.entry(oriented_edge.edge_curve_id).or_default() += 1;
            }
        }
    }

    let closed = use_counts.values().all(|&u| u >= 2);
    let pieces = match closed && allow_pieces {
        true => match touching_pieces(face_data, &use_counts) {
            Ok(pieces) => pieces,
            Err(unpaired) => {
                log.error(
                    Phase::MapTopology,
                    Some(unpaired.edge),
                    format!(
                        "STEP CLOSED_SHELL is non-manifold and no split into closed pieces pairs \
                         this edge's uses: {} of them lie in one piece (a surface touching itself \
                         along an edge)",
                        unpaired.uses_in_piece
                    ),
                );
                return None;
            }
        },
        false => None,
    };
    if let Some(pieces) = pieces {
        let contact_edges = use_counts.values().filter(|&&u| u > 2).count();
        log.warn(
            Phase::MapTopology,
            None,
            format!(
                "STEP CLOSED_SHELL is non-manifold: imported as {} manifold pieces meeting at \
                 {contact_edges} contact edge(s) and their shared vertices, each contact recorded",
                pieces.count()
            ),
        );
        return Some(ShellIncidencePolicy {
            treat_as_open: false,
            downgraded_closed_shell: false,
            pieces: Some(pieces),
        });
    }

    let mut downgraded_closed_shell = false;
    let mut declined = false;
    for (edge_curve_id, uses) in use_counts {
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
        pieces: None,
    })
}
