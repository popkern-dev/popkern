//! The cap commit on a carved setback region: the KP20 fan or the
//! variational fit (PK-574), by strategy.

use super::*;
use crate::construct::kp20_paper_border_extract::Kp20SetbackArc;
use crate::construct::setback_decline::{CapCommitError, SetbackCapStrategy};

/// Commit the cap `strategy` names over `blend_face`. Returns the cap
/// faces and, for the fit, its certified serve.
#[allow(clippy::too_many_arguments)]
pub(super) fn commit_cap(
    brep: &mut BRep,
    blend_face: FaceHandle,
    side_edges: &[EdgeHandle],
    fillet_surfaces: &[FilletSurfaceResult],
    arcs: &[Kp20SetbackArc],
    vertex_position: Point3,
    presolved: Option<&crate::construct::setback_cap_presolve::SetbackCapSolve>,
    strategy: SetbackCapStrategy,
    tolerance: &Tolerance,
) -> Result<(Vec<FaceHandle>, Option<setback_fit_fan::FitServe>), CapCommitError> {
    match strategy {
        SetbackCapStrategy::Kp20 => {
            if crate::dispatch_overrides::kp20_cap_decline_forced() {
                return Err("KP20 cap forced to decline (dispatch override)".into());
            }
            // The ranked assemblies (PK-613): the fairest first, then each
            // alternate on a fan postcondition failure, from the same carve.
            let alternates = presolved.map_or(&[][..], |solve| &solve.alternates[..]);
            let snapshot = (!alternates.is_empty()).then(|| brep.clone());
            let commit = |brep: &mut BRep,
                          solve: Option<
                &crate::construct::setback_cap_presolve::SetbackCapSolve,
            >| {
                super::setback_fan_commit::commit_cap_fan(
                    brep,
                    blend_face,
                    side_edges,
                    fillet_surfaces,
                    vertex_position,
                    solve,
                )
            };
            let first = match commit(brep, presolved) {
                Ok(faces) => return Ok((faces, None)),
                Err(e) => e,
            };
            if let (Some(solve), Some(snapshot)) = (presolved, snapshot) {
                for alternate in alternates {
                    *brep = snapshot.clone();
                    let retry = crate::construct::setback_cap_presolve::SetbackCapSolve {
                        border: solve.border.clone(),
                        extract: solve.extract.clone(),
                        assembly: Ok(alternate.clone()),
                        alternates: Vec::new(),
                    };
                    if let Ok(faces) = commit(brep, Some(&retry)) {
                        return Ok((faces, None));
                    }
                }
                *brep = snapshot;
            }
            Err(CapCommitError::from(first))
        }
        SetbackCapStrategy::RollingBall => {
            // The carved region must be the ball's triangle: rounds of one
            // radius around it, its boundary on the ball (the same check
            // the after-the-fact replacement makes on a served cap).
            let cap: std::collections::BTreeSet<FaceHandle> = [blend_face].into_iter().collect();
            let tol = tolerance.linear().max(1e-9);
            let (center, radius) = crate::rolling_ball_corner::rounds_ball(brep, blend_face, tol)
                .ok_or_else(|| {
                CapCommitError::from("the carved corner's rounds carry no ball".to_string())
            })?;
            // The carve's rims are coarse cubics: the exact cross-section
            // arcs replace them, then the ball check runs on the result.
            let touched = crate::rolling_ball_corner::exact_ball_rims(
                brep,
                blend_face,
                center,
                radius,
                (100.0 * tol).max(1e-5 * radius),
            )
            .map_err(|why| {
                CapCommitError::from(format!("the carved corner is not the rolling ball: {why}"))
            })?;
            let ball = crate::rolling_ball_corner::corner_ball(brep, &cap, tol).map_err(|why| {
                CapCommitError::from(format!("the carved corner is not the rolling ball: {why}"))
            })?;
            // The ball bulges where the rounds do: outward from its centre
            // when the round faces point away from their axes.
            let (probe, outward) = crate::rolling_ball_corner::round_sense(brep, blend_face)
                .ok_or_else(|| {
                    CapCommitError::from("the carved corner has no round neighbour".to_string())
                })?;
            crate::rolling_ball_corner::set_ball_surface(brep, blend_face, &ball, probe, outward);
            brep.recompute_pcurves_for_lenient(&touched);
            Ok((vec![blend_face], None))
        }
        SetbackCapStrategy::VariationalFit {
            mode,
            r_min,
            rival,
            held,
            defer,
            ..
        } => super::setback_fit_fan::commit_fit_fan(
            brep,
            blend_face,
            fillet_surfaces,
            arcs,
            mode,
            super::setback_fit_fan::FitStakes {
                r_min,
                rival,
                bar: held,
                defer,
            },
            tolerance,
        )
        .map(|(faces, serve)| (faces, Some(serve))),
    }
}
