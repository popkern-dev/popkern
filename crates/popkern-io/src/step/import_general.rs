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

/// Glue `brep` into one general body where its pieces touch: its
/// recorded contacts (a shell split into pieces), then whatever
/// coincides without being shared (`BRep::glue_coincident`: a writer's
/// cavities touching along an edge each with their own copy of it), then
/// the cells a face-by-face shell encloses (`BRep::fill_enclosed_cells`).
/// Kept when it validates; logged against STEP entity `id` either way.
pub(crate) fn glue_pieces(brep: &mut BRep, id: u64, log: &mut DiagnosticLog) {
    let mut glued = brep.clone();
    let mut steps = Vec::new();
    let outcome = (|| -> Result<(), String> {
        if !glued.topology_contacts.is_empty() {
            let report = glued.glue_contacts().map_err(|e| e.to_string())?;
            steps.push(format!(
                "{} contact(s) glued, {} kept recorded",
                report.glued,
                report.kept.len()
            ));
        }
        let coincident = glued.glue_coincident().map_err(|e| e.to_string())?;
        if coincident.vertices + coincident.edges > 0 {
            steps.push(format!(
                "{} coincident vertex(es) and {} edge(s) glued",
                coincident.vertices, coincident.edges
            ));
        }
        if steps.is_empty() {
            return Ok(());
        }
        // A shell imported face by face glues into sheets: the cells
        // they enclose become its solids.
        glued.fill_enclosed_cells();
        let issues = glued.validate_all();
        if issues.is_empty() {
            Ok(())
        } else {
            Err(format!("glued body invalid: {}", issues[0]))
        }
    })();
    match outcome {
        Ok(()) if steps.is_empty() => {}
        Ok(()) => {
            *brep = glued;
            log.warn(
                Phase::MapTopology,
                Some(id),
                format!(
                    "touching pieces imported as one general body: {}",
                    steps.join("; ")
                ),
            );
        }
        Err(why) => log.warn(
            Phase::MapTopology,
            Some(id),
            format!("touching pieces kept as they are: {why}"),
        ),
    }
}
