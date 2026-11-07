//! Importing a self-touching shell as one general body (PK-1108).
//!
//! A `CLOSED_SHELL` listing the faces of solids that touch (two blocks
//! sharing an edge or a corner, as OCC writes them) maps as manifold
//! pieces meeting at recorded contacts (`shell_pieces`). Under
//! [`StepNonManifold::General`] the pieces are glued at those contacts
//! (`BRep::glue_contacts`): one general body whose contact edge's ring
//! holds every piece's faces in radial order. A glue that declines, or a
//! glued body that does not validate, keeps the pieces and says why.

use popkern_core::topology::BRep;

use crate::step::diagnostics::{DiagnosticLog, Phase};

/// How a STEP shell that touches itself imports.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum StepNonManifold {
    /// Manifold pieces meeting at recorded contacts.
    #[default]
    Pieces,
    /// One general body: the pieces glued at their contacts.
    General,
}

/// Glue `brep`'s recorded contacts when it has any, logging the outcome
/// against STEP entity `id`.
pub(crate) fn glue_pieces(brep: &mut BRep, id: u64, log: &mut DiagnosticLog) {
    if brep.topology_contacts.is_empty() {
        return;
    }
    let mut glued = brep.clone();
    let outcome = match glued.glue_contacts() {
        Ok(report) => {
            // A shell imported face by face glues into sheets: the cells
            // they enclose become its solids.
            glued.fill_enclosed_cells();
            let issues = glued.validate_all();
            if issues.is_empty() {
                let kept = report.kept.len();
                *brep = glued;
                Ok(format!(
                    "self-touching shell imported as one general body: {} contact(s) glued, {kept} kept recorded",
                    report.glued
                ))
            } else {
                Err(format!("glued body invalid: {}", issues[0]))
            }
        }
        Err(e) => Err(e.to_string()),
    };
    match outcome {
        Ok(message) => log.warn(Phase::MapTopology, Some(id), message),
        Err(why) => log.warn(
            Phase::MapTopology,
            Some(id),
            format!("self-touching shell kept as pieces with recorded contacts: {why}"),
        ),
    }
}
