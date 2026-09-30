//! The bridged setback border's retry (PK-611).
//!
//! The setback surgery carves the corners-at-crossings border first
//! (with reflex-face bridges at a saddle). When a setback arc crosses
//! its strip's trimmed boundary mid-arc, the crossing carve declines:
//! a thin strip's arc, run between the corners its wide neighbours fix,
//! overruns the strip. The surgery then retries on the bridged border
//! (`SetbackJunctionPolicy::BridgeEveryJunction`): each strip set back
//! to a station of its own, a chord bridge on every convex support
//! face, carved by the ring carve.
//!
//! Only that decline retries, so every corner the corners border
//! serves keeps its construction bit for bit. When the retry declines
//! too, both reasons are reported, the overrun first: the rung names
//! the corner's decline class from it.
//!
//! # The variational fit serves the bridged border
//!
//! The bridged border has 2n sides. The KP20 fan does not serve it:
//! - the KP18 guide weights stop at valence 10, so an apex of six or
//!   more strips has no fan;
//! - on the pentagon h0.4 variable and h0.8 wide apexes the valence-10
//!   fan passes its seam gates on sampled evidence, but its sectors
//!   swirl into a visible dent about the star (up to 49 of 110 gated
//!   seams stop certifying on a Jacobian enclosure through zero);
//! - its attempt costs up to 40 s before a fold decline.
//!
//! So the KP20 strategy declines there by name (a cap-stage decline,
//! which hands the corner to the fit on the same carve), and the
//! border is solved without the KP20 assembly.

use popkern_core::geometry::Point3;
use popkern_core::topology::BRep;

use super::SetbackTopologyOutcome;
use crate::construct::kp20_paper_border_extract::{
    SetbackJunctionPolicy, kp20_paper_border_from_fillet_setback_with,
};
use crate::construct::setback_cap_presolve::SetbackCapSolve;
use crate::construct::setback_decline::{
    ARC_BOUNDARY_CROSSINGS, SetbackCapStrategy, SetbackDecline, SetbackStage,
};
use crate::fillet_surface::FilletSurfaceResult;

/// The KP20 strategy's decline on the bridged border.
pub(crate) const BRIDGED_FAN_DECLINE: &str =
    "the bridged border's 2n-sided KP20 fan is not attempted (the fit serves it)";

/// Run `attempt` on the surgery's own border, then on the bridged
/// border after an overrun decline. `attempt` restores the model on
/// every decline.
pub(super) fn with_bridged_retry(
    brep: &mut BRep,
    fillet_surfaces: &[FilletSurfaceResult],
    vertex_position: Point3,
    strategy: SetbackCapStrategy,
    attempt: impl Fn(
        &mut BRep,
        Option<&SetbackCapSolve>,
    ) -> Result<SetbackTopologyOutcome, SetbackDecline>,
) -> Result<SetbackTopologyOutcome, SetbackDecline> {
    let overrun = match attempt(brep, None) {
        Err(d) if d.stage == SetbackStage::Carve && d.reason.contains(ARC_BOUNDARY_CROSSINGS) => d,
        result => return result,
    };
    if strategy == SetbackCapStrategy::Kp20 {
        return Err(both(
            overrun,
            SetbackDecline {
                stage: SetbackStage::Cap,
                reason: BRIDGED_FAN_DECLINE.to_string(),
                fit_certified: None,
            },
        ));
    }
    let solve = match kp20_paper_border_from_fillet_setback_with(
        fillet_surfaces,
        vertex_position,
        SetbackJunctionPolicy::BridgeEveryJunction,
    ) {
        Ok((border, extract)) => SetbackCapSolve {
            border,
            extract,
            assembly: Err(BRIDGED_FAN_DECLINE.to_string()),
        },
        Err(e) => return Err(both(overrun, SetbackDecline::from(e))),
    };
    attempt(brep, Some(&solve)).map_err(|bridged| both(overrun, bridged))
}

/// The overrun decline followed by the bridged border's, at the
/// bridged border's stage (a cap decline there is still a cap the
/// next construction may serve on the same carve).
fn both(overrun: SetbackDecline, bridged: SetbackDecline) -> SetbackDecline {
    SetbackDecline {
        stage: bridged.stage,
        reason: format!("{}; bridged border: {}", overrun.reason, bridged.reason),
        fit_certified: bridged.fit_certified,
    }
}
