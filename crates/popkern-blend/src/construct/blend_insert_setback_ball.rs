//! The setback rung's rolling-ball step (PK-864): where the rounds can
//! carry the ball, the corner is carved and the carved region takes the
//! sphere directly, instead of fitting and certifying a KP20 cap that
//! the rolling-ball replacement would throw away (12.9 s of a rounded
//! tetrahedron's 13 s). Any refusal leaves the model as it was and the
//! corner continues down the rung as before.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::BRep;

use super::{SetbackRung, commit};
use crate::construct::setback_decline::SetbackCapStrategy;
use crate::types::{ContinuityLevel, CornerMethod};

impl SetbackRung<'_> {
    /// The rolling-ball cap when `cap` asks for it, no user setback is
    /// named, and the corner's rounds carry a ball; else `None` (the
    /// model untouched).
    pub(super) fn ball_if_carried(
        &self,
        brep: &mut BRep,
        cap: crate::types::CornerCapPreference,
        user_setback: bool,
    ) -> Option<super::InsertVertexBlendOutput> {
        let carried = cap.construction == crate::types::CornerConstruction::RollingBall
            && !user_setback
            && crate::rolling_ball_corner::rounds_carry_a_ball(
                self.fillet_surfaces_pre.iter().map(|f| &f.surface),
                self.tolerance.linear().max(1e-9),
            );
        if carried { self.ball(brep) } else { None }
    }

    /// The rolling-ball cap, or `None` (the model restored) when the
    /// carve or the ball check refuses.
    pub(super) fn ball(&self, brep: &mut BRep) -> Option<super::InsertVertexBlendOutput> {
        let outcome = self.surgery(brep, SetbackCapStrategy::RollingBall).ok()?;
        let mut served = match commit(brep, &self.commit_snapshot, &outcome, |f, _| {
            self.info(
                f,
                CornerMethod::RollingBallSphere,
                ContinuityLevel::G1,
                None,
            )
        }) {
            Ok(served) => served,
            Err(_) => {
                *brep = self.commit_snapshot.clone();
                return None;
            }
        };
        for info in &mut served.0 {
            info.setback_cap = false;
            if let Some(SurfaceGeometry::Sphere { center, radius, .. }) = brep
                .faces
                .get(info.face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry)
            {
                info.cap_sphere = Some((*center, *radius));
            }
        }
        Some(served)
    }
}
