//! What a public boolean entry point serves under each non-manifold
//! policy (PK-1115).
//!
//! The pipeline answers a non-manifold result as manifold pieces in one
//! body, touching without sharing: an intermediate. Declined by
//! default; under `General` glued into one general body; under
//! `SplitToManifold` separate bodies, which a single-body entry point
//! cannot return ([`BoolError::SplitIntoBodies`]) and its `_bodies`
//! form does. A result whose pieces touch nowhere is one body under
//! every policy, as the default serves it.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::BoolError;
use crate::options::NonManifoldPolicy;
use crate::postprocess::split_pieces::touch_contacts;

/// One body from the pipeline's answer `pieces` under `policy`.
///
/// # Errors
///
/// The pipeline's error; [`BoolError::GeneralBodyDeclined`];
/// [`BoolError::SplitIntoBodies`] when split pieces touch.
pub(crate) fn serve_one(
    pieces: Result<BRep, BoolError>,
    policy: NonManifoldPolicy,
    tolerance: &Tolerance,
) -> Result<BRep, BoolError> {
    let body = pieces?;
    match policy {
        NonManifoldPolicy::Decline => Ok(body),
        NonManifoldPolicy::General => crate::general::glue_result(body, tolerance),
        NonManifoldPolicy::SplitToManifold => {
            let contacts = touch_contacts(&body, tolerance);
            if contacts.is_empty() {
                Ok(body)
            } else {
                Err(BoolError::SplitIntoBodies {
                    bodies: body.solid_regions().len(),
                    contacts,
                })
            }
        }
    }
}

/// The bodies from the pipeline's answer `pieces` under `policy`: one
/// body, except split pieces that touch, each its own body.
///
/// # Errors
///
/// As [`serve_one`], never [`BoolError::SplitIntoBodies`].
pub(crate) fn serve_bodies(
    pieces: Result<BRep, BoolError>,
    policy: NonManifoldPolicy,
    tolerance: &Tolerance,
) -> Result<Vec<BRep>, BoolError> {
    match policy {
        NonManifoldPolicy::SplitToManifold => {
            let body = pieces?;
            if touch_contacts(&body, tolerance).is_empty() {
                return Ok(vec![body]);
            }
            body.piece_bodies()
                .map_err(|e| BoolError::DegenerateInput(e.to_string()))
        }
        _ => serve_one(pieces, policy, tolerance).map(|b| vec![b]),
    }
}
