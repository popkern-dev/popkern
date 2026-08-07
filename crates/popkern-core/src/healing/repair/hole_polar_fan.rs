//! Polar-wedge G1 fan cap for N-sided holes (PK-248).
//!
//! The single-patch rung (PK-224) serves 4-sided holes with no
//! internal seams; general-N holes previously fell to the PnS bi-4
//! cap (which declines its own internal fan-seam gap at commit) and
//! then to the G0 ruled fan. This rung fills the gap for every
//! N != 4: N degree-5 Gregory polar wedges sharing exact quintic
//! spoke curves, tangent
//! to the surviving walls through the shared rim-ribbon machinery,
//! G1 to each other across the spokes by midpoint-offset member
//! pairing, with the pole enclosure closed by the Shirman-Séquin
//! solve (see [`super::hole_polar_fan_wedges`] for the construction
//! and per-seam mechanisms).
//!
//! # Serving contract
//!
//! Both G1 claims are MEASURED before any topology is touched: the
//! cap↔wall mid-band per rim side and the cross-spoke mid-band per
//! spoke must clear the commercial 5e-2 `|sin|` bar or the builder
//! declines with the worst sample. The returned [`HoleCapGeometry`]
//! commits through the transactional fan operator
//! ([`super::hole_commit::commit_hole_cap`], `ProjectableFan` chart),
//! whose coherence gate and validation are the final authority.
//! Corner points are excluded from the G1 claim for the same
//! structural reason as the single-patch rung: at a rim corner the
//! two walls disagree, and at the pole the parameterization
//! collapses.

use crate::geometry::curve::CurveGeometry;
use crate::geometry::gregory::GregoryPatch;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{CurveOnEdge, NurbsCurve, Point3};
use crate::topology::BRep;

use super::hole_commit_geometry::{HoleCapChart, HoleCapGeometry, HoleCapSectorGeometry};
use super::hole_fill::DetectedHole;
use super::hole_polar_fan_wedges::build_polar_wedges;
use super::hole_rim_ribbon::{RimRibbonError, neighbour_surface_geometry};
use super::n_sided_fill_g1::{G1HoleFillError, boundary_edges_with_frames};

/// Fewest sides this rung serves. 4-sided holes belong to the
/// seam-free single-patch rung; every other side count takes the fan
/// (N = 3 measured: open_prism3 serves with machine-zero spoke seams
/// and on-edge wall tangency, same as N >= 5).
const MIN_SIDES: usize = 3;

/// Mid-band fractions for both measured G1 gates (corners excluded
/// by the structural corner argument — see module doc).
const G1_SAMPLE_FRACTIONS: [f64; 5] = [0.25, 0.375, 0.5, 0.625, 0.75];

/// Measured G1 bar, `|sin|` of the normal angle.
const G1_GATE: f64 = 5.0e-2;

/// Why the polar-fan cap declined. Stringly summarized into the heal
/// report like the other rungs.
#[derive(Debug, thiserror::Error)]
pub enum PolarFanCapError {
    /// The hole's side count is outside this rung's service class
    /// (fewer than 3 sides, or exactly 4 — the seam-free
    /// single-patch rung owns 4-sided holes).
    #[error(
        "polar fan cap serves 3- and 5+-sided holes (4-sided holes take the single-patch rung); this hole has {sides}"
    )]
    UnsupportedSideCount {
        /// Sides of the detected hole.
        sides: usize,
        /// The rung's minimum side count.
        min: usize,
    },
    /// The boundary frames could not be built.
    #[error("polar fan cap: {0}")]
    Frames(#[from] G1HoleFillError),
    /// A rim side's quintic fit missed the sampled rim curve.
    #[error("polar fan cap: side {side} quintic fit residual {residual:.3e} exceeds {budget:.3e}")]
    RimFit {
        /// Side index in ring order.
        side: usize,
        /// Worst sample distance to the fitted Bezier.
        residual: f64,
        /// Linear tolerance budget scaled for fitting.
        budget: f64,
    },
    /// A rim side's ribbon direction could not be oriented.
    #[error("polar fan cap: side {side} has no usable ribbon direction")]
    RibbonDegenerate {
        /// Side index in ring order.
        side: usize,
    },
    /// The ring geometry cannot support the fan construction.
    #[error("polar fan cap: {detail}")]
    DegenerateRing {
        /// What was degenerate, with the measured value.
        detail: String,
    },
    /// The pole vertex-enclosure solve was unavailable.
    #[error("polar fan cap: the pole vertex-enclosure solve is unavailable for this ring")]
    EnclosureUnsolvable,
    /// A wedge's Gregory assembly failed.
    #[error("polar fan cap: wedge {wedge} assembly failed: {detail}")]
    Assembly {
        /// Wedge index in ring order.
        wedge: usize,
        /// The constructor's reason.
        detail: String,
    },
    /// The measured cap↔wall G1 exceeded the gate mid-side.
    #[error("polar fan cap: side {side} measured cap-wall G1 {worst:.3e} exceeds {gate:.1e}")]
    WallG1Exceeded {
        /// Side index in ring order.
        side: usize,
        /// Worst |sin| over the mid-band samples.
        worst: f64,
        /// The gate.
        gate: f64,
    },
    /// The measured wedge↔wedge G1 exceeded the gate along a spoke.
    #[error("polar fan cap: spoke {spoke} measured cross-spoke G1 {worst:.3e} exceeds {gate:.1e}")]
    SpokeG1Exceeded {
        /// Spoke index (spoke k ends at ring vertex k).
        spoke: usize,
        /// Worst |sin| over the mid-band samples.
        worst: f64,
        /// The gate.
        gate: f64,
    },
}

impl From<RimRibbonError> for PolarFanCapError {
    fn from(error: RimRibbonError) -> Self {
        match error {
            RimRibbonError::FitResidual {
                side,
                residual,
                budget,
            } => Self::RimFit {
                side,
                residual,
                budget,
            },
            RimRibbonError::RibbonDegenerate { side } => Self::RibbonDegenerate { side },
        }
    }
}

/// Build the polar-wedge G1 fan for a detected hole (any side count
/// except 4 — see [`PolarFanCapError::UnsupportedSideCount`]),
/// gating both measured G1 claims before returning. The result
/// commits through [`super::hole_commit::commit_hole_cap`]; nothing
/// here touches topology.
pub fn polar_fan_g1_cap(
    brep: &BRep,
    hole: &DetectedHole,
    tolerance: f64,
) -> Result<HoleCapGeometry, PolarFanCapError> {
    let sides = hole.boundary_edges.len();
    if sides < MIN_SIDES || sides == 4 {
        return Err(PolarFanCapError::UnsupportedSideCount {
            sides,
            min: MIN_SIDES,
        });
    }
    let frames = boundary_edges_with_frames(brep, hole, tolerance)?;
    let set = build_polar_wedges(brep, hole, &frames, tolerance)?;

    let surfaces: Vec<SurfaceGeometry> = set
        .wedges
        .iter()
        .map(|wedge| {
            SurfaceGeometry::Gregory(GregoryPatch::RectangularDeg5(Box::new(wedge.clone())))
        })
        .collect();

    // Gate 1 — cap↔wall mid-band G1 per rim side (the u = 0 edge of
    // wedge k, walked by v).
    for (side, surface) in surfaces.iter().enumerate() {
        let mut worst = 0.0_f64;
        for &f in &G1_SAMPLE_FRACTIONS {
            let Some(cap_n) = surface.normal_at(0.0, f).normalized() else {
                continue;
            };
            let rim_point = surface.point_at(0.0, f);
            let Some(neighbour) = neighbour_surface_geometry(brep, hole, side) else {
                continue;
            };
            let Some((nu, nv)) = neighbour.closest_point(&rim_point) else {
                continue;
            };
            let Some(wall_n) = neighbour.normal_at(nu, nv).normalized() else {
                continue;
            };
            worst = worst.max(cap_n.cross(&wall_n).length());
        }
        if worst > G1_GATE {
            return Err(PolarFanCapError::WallG1Exceeded {
                side,
                worst,
                gate: G1_GATE,
            });
        }
    }

    // Gate 2 — cross-spoke mid-band G1. Spoke k is wedge k's v = 0
    // edge and wedge k−1's v = 1 edge; along both, u = 1 − s where s
    // is the spoke fraction from the pole.
    let n = surfaces.len();
    for spoke in 0..n {
        let ccw = &surfaces[spoke];
        let cw = &surfaces[(spoke + n - 1) % n];
        let mut worst = 0.0_f64;
        for &s in &G1_SAMPLE_FRACTIONS {
            let u = 1.0 - s;
            let (Some(n_ccw), Some(n_cw)) = (
                ccw.normal_at(u, 0.0).normalized(),
                cw.normal_at(u, 1.0).normalized(),
            ) else {
                continue;
            };
            worst = worst.max(n_ccw.cross(&n_cw).length());
        }
        if worst > G1_GATE {
            return Err(PolarFanCapError::SpokeG1Exceeded {
                spoke,
                worst,
                gate: G1_GATE,
            });
        }
    }

    let sectors = surfaces
        .into_iter()
        .enumerate()
        .map(|(k, surface)| HoleCapSectorGeometry {
            surface,
            rim_start: frames[k].start,
            rim_end: frames[k].end,
            spoke_at_rim_start: quintic_spoke_curve(&set.spokes[k]),
            spoke_at_rim_end: quintic_spoke_curve(&set.spokes[(k + 1) % n]),
        })
        .collect();

    Ok(HoleCapGeometry {
        pole: set.pole,
        sectors,
        chart: HoleCapChart::ProjectableFan,
    })
}

/// Wrap a clamped quintic Bezier spoke as an edge curve on `[0, 1]`.
fn quintic_spoke_curve(cps: &[Point3; 6]) -> CurveOnEdge {
    CurveOnEdge {
        geometry: CurveGeometry::Nurbs(NurbsCurve {
            degree: 5,
            knots: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            control_points: cps.to_vec(),
            weights: vec![1.0; 6],
        }),
        t_start: 0.0,
        t_end: 1.0,
    }
}
