//! Single-patch G1 cap for 4-sided holes (PK-224).
//!
//! A 4-sided hole does not need a sector complex: one rectangular
//! degree-5 Gregory patch over the whole ring has NO internal seams
//! (the fan-crease class cannot exist) and meets each surviving
//! neighbour tangentially through its boundary ribbon rows. The
//! Gregory inner pairs exist precisely to let two adjacent sides'
//! ribbons disagree at a shared corner.
//!
//! # Corner reality
//!
//! At a hole corner the two neighbour walls meet at an angle; no cap
//! can be G1 to both at the corner point. Each rim side's quintic is
//! therefore built with ZERO end legs (`[S, S, A, B, E, E]` — the
//! curve geometry is unchanged, the parameterization has a cusp at
//! the corners), which kills the adjacent side's corner-leg
//! contribution to this side's cross field: the cross-derivative
//! along each open side is then ribbon + tangential terms only —
//! G1 to the neighbour by construction away from the corner points.
//! Design: `docs/architecture/hole-cap-single-patch-g1.md`.

use crate::error::CoreError;
use crate::geometry::gregory::{GregoryPatch, RectangularGregoryDeg5};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::topology::BRep;
use crate::topology::coherence::CoherenceValidationScope;

use super::hole_fill::DetectedHole;
use super::hole_rim_ribbon::{
    RimRibbonError, fit_quintic_zero_leg, neighbour_surface_geometry, rim_ribbon_directions,
};
use super::n_sided_fill_g1::{G1HoleFillError, boundary_edges_with_frames};

/// Mid-band fractions used for the measured cap↔neighbour G1 gate.
/// Corner bands are excluded by the structural corner argument (see
/// module doc); the open middle of each side is where a G1 cap must
/// actually be tangent to its wall.
const G1_SAMPLE_FRACTIONS: [f64; 5] = [0.25, 0.375, 0.5, 0.625, 0.75];

/// Measured cap↔neighbour G1 bar, `|sin|` of the normal angle.
const G1_GATE: f64 = 5.0e-2;

/// Why the single-patch cap declined. Stringly summarized into the
/// heal report like the other rungs.
#[derive(Debug, thiserror::Error)]
pub enum SinglePatchCapError {
    /// Only 4-sided holes are served by this rung.
    #[error("single-patch cap serves 4-sided holes; this hole has {sides}")]
    UnsupportedSideCount {
        /// Sides of the detected hole.
        sides: usize,
    },
    /// The boundary frames could not be built.
    #[error("single-patch cap: {0}")]
    Frames(#[from] G1HoleFillError),
    /// A side's quintic fit missed the sampled rim curve.
    #[error(
        "single-patch cap: side {side} quintic fit residual {residual:.3e} exceeds {budget:.3e}"
    )]
    FitResidual {
        /// Side index in ring order.
        side: usize,
        /// Worst sample distance to the fitted Bezier.
        residual: f64,
        /// Linear tolerance budget scaled for fitting.
        budget: f64,
    },
    /// A side's ribbon direction could not be oriented.
    #[error("single-patch cap: side {side} has no usable ribbon direction")]
    RibbonDegenerate {
        /// Side index in ring order.
        side: usize,
    },
    /// The measured cap↔neighbour G1 exceeded the gate mid-side.
    #[error(
        "single-patch cap: side {side} measured cap-neighbour G1 {worst:.3e} exceeds {gate:.1e}"
    )]
    G1Exceeded {
        /// Side index in ring order.
        side: usize,
        /// Worst |sin| over the mid-band samples.
        worst: f64,
        /// The gate.
        gate: f64,
    },
    /// The topology or pcurve layer declined.
    #[error("single-patch cap: {0}")]
    Topology(#[from] CoreError),
    /// The committed cap failed the coherence gate.
    #[error(
        "single-patch cap: coherence gate failed — {missing_curves} curve(s), {missing_pcurves} pcurve(s) missing, worst L2 {l2_worst:.3e}"
    )]
    CoherenceRejected {
        /// L1 violations on the cap face.
        missing_curves: usize,
        /// L3 violations on the cap face.
        missing_pcurves: usize,
        /// Worst sampled boundary disagreement.
        l2_worst: f64,
    },
    /// Post-commit validation failed.
    #[error("single-patch cap: post-commit validation failed: {0}")]
    PostValidationFailed(#[from] crate::error::ValidationError),
}

impl From<RimRibbonError> for SinglePatchCapError {
    fn from(error: RimRibbonError) -> Self {
        match error {
            RimRibbonError::FitResidual {
                side,
                residual,
                budget,
            } => Self::FitResidual {
                side,
                residual,
                budget,
            },
            RimRibbonError::RibbonDegenerate { side } => Self::RibbonDegenerate { side },
        }
    }
}

/// A committed single-patch cap.
#[derive(Debug, Clone)]
pub struct SinglePatchCapCommit {
    /// The new cap face.
    pub face: crate::topology::FaceHandle,
    /// Worst L2 disagreement from the coherence gate.
    pub worst_l2_gap: f64,
    /// Worst measured cap↔neighbour G1 (`|sin|`) over the mid-band
    /// samples of all four sides.
    pub worst_mid_band_g1: f64,
}

/// Build and transactionally commit the single-patch G1 cap for a
/// 4-sided hole. Any decline leaves `brep` untouched.
pub fn commit_single_patch_hole_cap(
    brep: &mut BRep,
    hole: &DetectedHole,
    tolerance: f64,
) -> Result<SinglePatchCapCommit, SinglePatchCapError> {
    if hole.boundary_edges.len() != 4 {
        return Err(SinglePatchCapError::UnsupportedSideCount {
            sides: hole.boundary_edges.len(),
        });
    }
    let frames = boundary_edges_with_frames(brep, hole, tolerance)?;

    // Quintic side rows with zero end legs, chained CCW.
    let mut boundaries = [[Point3::new(0.0, 0.0, 0.0); 6]; 4];
    let mut chords = [0.0_f64; 4];
    for (side, frame) in frames.iter().enumerate() {
        boundaries[side] = fit_quintic_zero_leg(frame, side, tolerance)?;
        chords[side] = frame.start.distance_to(&frame.end);
    }

    // Per-side ribbon directions at the 1/3 and 2/3 stations,
    // oriented away from the neighbour face across the rim.
    let mut ribbon = [[Vec3::new(0.0, 0.0, 0.0); 2]; 4];
    for (side, frame) in frames.iter().enumerate() {
        ribbon[side] = rim_ribbon_directions(brep, hole, frame, side)?;
    }

    let patch = assemble_patch(&boundaries, &chords, &ribbon, tolerance)
        .map_err(|e| SinglePatchCapError::Topology(CoreError::PreconditionFailed(e.to_string())))?;

    // Measured G1 gate BEFORE committing topology: sample each side's
    // mid band and compare the cap normal against the neighbour's.
    let cap_geometry =
        SurfaceGeometry::Gregory(GregoryPatch::RectangularDeg5(Box::new(patch.clone())));
    let mut worst_mid_band_g1 = 0.0_f64;
    for side in 0..4 {
        let mut worst = 0.0_f64;
        for &f in &G1_SAMPLE_FRACTIONS {
            let (u, v) = side_uv(side, f);
            let cap_n = cap_geometry.normal_at(u, v);
            let Some(cap_n) = cap_n.normalized() else {
                continue;
            };
            let rim_point = cap_geometry.point_at(u, v);
            let Some(neighbour) = neighbour_surface_geometry(brep, hole, side) else {
                continue;
            };
            let Some((nu, nv)) = neighbour.closest_point(&rim_point) else {
                continue;
            };
            let Some(wall_n) = neighbour.normal_at(nu, nv).normalized() else {
                continue;
            };
            let sin = cap_n.cross(&wall_n).length();
            worst = worst.max(sin);
        }
        worst_mid_band_g1 = worst_mid_band_g1.max(worst);
        if worst > G1_GATE {
            return Err(SinglePatchCapError::G1Exceeded {
                side,
                worst,
                gate: G1_GATE,
            });
        }
    }

    // Transactional commit: single-face ring cap + surface + pcurves
    // + coherence + validation on a candidate clone.
    let mut candidate = brep.clone();
    let ring = candidate.cap_boundary_ring(&hole.boundary_edges)?;
    if let Some(face) = candidate.faces.get_mut(ring.face) {
        face.surface = Some(SurfaceOnFace {
            geometry: cap_geometry.clone(),
            orientation: true,
        });
    }
    candidate.recompute_face_pcurves(ring.face)?;
    // The orientation flag is derived state (PK-244); reconcile it
    // from the committed winding rather than trusting the stamp.
    let _ = candidate.reconcile_orientation_flags();

    let scope = CoherenceValidationScope::Faces([ring.face].into_iter().collect());
    let report = candidate.validate_coherence(&scope);
    let worst_l2_gap = report
        .l2_face_boundary_disagreements
        .iter()
        .map(|d| d.worst_gap)
        .fold(0.0_f64, f64::max);
    if !report.is_empty() {
        return Err(SinglePatchCapError::CoherenceRejected {
            missing_curves: report.l1_missing_edge_curves.len(),
            missing_pcurves: report.l3_missing_pcurves.len(),
            l2_worst: worst_l2_gap,
        });
    }
    candidate.validate()?;

    *brep = candidate;
    Ok(SinglePatchCapCommit {
        face: ring.face,
        worst_l2_gap,
        worst_mid_band_g1,
    })
}

/// Assemble the degree-5 Gregory patch: chained boundaries + ribbon
/// inner pairs per corner.
fn assemble_patch(
    boundaries: &[[Point3; 6]; 4],
    chords: &[f64; 4],
    ribbon: &[[Vec3; 2]; 4],
    tolerance: f64,
) -> Result<RectangularGregoryDeg5, crate::geometry::gregory::GregoryBuildError> {
    // Ribbon member for side s at its start (station 0) and end
    // (station 1): a point over the side's own row shifted by the
    // ribbon direction. The tangential part is irrelevant to G1 (it
    // lies in the neighbour plane); the ribbon part carries it.
    let member = |s: usize, station: usize| -> Point3 {
        let q = &boundaries[s];
        let base = if station == 0 {
            Point3::new(
                0.5 * (q[1].x + q[2].x),
                0.5 * (q[1].y + q[2].y),
                0.5 * (q[1].z + q[2].z),
            )
        } else {
            Point3::new(
                0.5 * (q[3].x + q[4].x),
                0.5 * (q[3].y + q[4].y),
                0.5 * (q[3].z + q[4].z),
            )
        };
        let h = chords[s] / 3.0;
        base + ribbon[s][station] * h
    };
    // Corner incident sides (CCW chaining b0..b3, corner k = start
    // of side k): side k starts at corner k, side (k+3)%4 ends at
    // corner k. Member-to-pair-slot assignment is pinned by the
    // calibration test `deg5_pair_member_convention` — on the u=0
    // edge the blend limits to the P_u member, so the side lying
    // along u=0 (boundaries[0]) owns the P_u slots of corners 0/1,
    // and so on around the ring.
    let inner_cps: [(Point3, Point3); 4] = [
        (member(0, 0), member(3, 1)),
        (member(0, 1), member(1, 0)),
        (member(2, 1), member(1, 1)),
        (member(2, 0), member(3, 0)),
    ];
    RectangularGregoryDeg5::from_boundary_beziers(*boundaries, inner_cps, tolerance.max(1e-9))
}

/// `(u, v)` on the patch for fraction `f` along ring side `s`,
/// following the `from_boundary_beziers` convention (side 0 = u=0
/// edge walking v, side 1 = v=1 edge walking u, side 2 = u=1 edge
/// walking v backwards, side 3 = v=0 edge walking u backwards).
fn side_uv(side: usize, f: f64) -> (f64, f64) {
    match side {
        0 => (0.0, f),
        1 => (f, 1.0),
        2 => (1.0, 1.0 - f),
        _ => (1.0 - f, 0.0),
    }
}

#[cfg(test)]
#[path = "hole_single_patch_tests.rs"]
mod tests;
