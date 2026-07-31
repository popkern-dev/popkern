//! Normalize hole-filler output into committable cap geometry.
//!
//! Both fillers produce a fan of N sector patches around a pole, but
//! in different containers: the G¹ PnS cap returns bi-4 NURBS sectors
//! whose spoke iso-curves are shared exactly, while the Coons fan
//! returns per-sector surfaces plus corner points and rebuilds its
//! spokes as centroid lines. [`HoleCapGeometry`] is the single shape
//! [`super::hole_commit::commit_hole_cap`] consumes, so both fillers
//! commit through the same operator (PK-199 acceptance 314).

use crate::geometry::curve::CurveGeometry;
use crate::geometry::{CurveOnEdge, NSidedFill, NurbsCurve, Point3, SurfaceGeometry};

use super::coons_patch::CoonsPatchResult;
use super::hole_commit::HoleCommitError;

/// Which chart family the sector surfaces follow, deciding how the
/// commit authors pcurves on the new faces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HoleCapChart {
    /// Fan-shaped patch sets whose sectors admit 2-D projection
    /// (Coons-fan wedges, PnS bi-4 attempts). Pcurves come from
    /// `recompute_face_pcurves` plus the pole-corner snap; the
    /// coherence gate then measures whether the sectors genuinely
    /// meet — the PnS bi-4 complex declines here with its real
    /// internal fan-seam gap (9.7e-2 on the open-box fill; see the
    /// PK-199 resolution notes for the falsified tilings).
    ProjectableFan,
    /// Ruled pole↔rim sectors from [`super::hole_ruled_fan`]: every
    /// boundary curve lies on its sector by construction, so all
    /// pcurves are authored as exact parameter lines — no projection.
    RuledFanAnalytic,
}

/// One committable sector of a hole cap.
#[derive(Debug, Clone)]
pub struct HoleCapSectorGeometry {
    /// The sector surface.
    pub surface: SurfaceGeometry,
    /// Rim start position (a hole ring vertex).
    pub rim_start: Point3,
    /// Rim end position (the next hole ring vertex).
    pub rim_end: Point3,
    /// Spoke curve from the pole to `rim_start`.
    pub spoke_at_rim_start: CurveOnEdge,
    /// Spoke curve from the pole to `rim_end`.
    pub spoke_at_rim_end: CurveOnEdge,
}

/// A hole cap normalized from either filler.
#[derive(Debug, Clone)]
pub struct HoleCapGeometry {
    /// The pole all sectors meet at (becomes the fan's pole vertex).
    pub pole: Point3,
    /// One sector per hole side, in filler order.
    pub sectors: Vec<HoleCapSectorGeometry>,
    /// Chart family of the sector surfaces.
    pub chart: HoleCapChart,
}

impl HoleCapGeometry {
    /// Normalize a G¹ PnS fill.
    ///
    /// Sector k of an [`NSidedFill`] is a bi-4 NURBS whose control
    /// net has the pole at `cp(0,0)`, the rim as the `v = 1` iso from
    /// `cp(0,4)` to `cp(4,4)`, and the spoke at its rim start as the
    /// `u = 0` iso. The sector *quads overlap*: `cp(4,0)` is the
    /// previous ring vertex, not the pole (measured on the open-box
    /// fill — the collapsed-wedge reading was falsified at 0.73 from
    /// the pole), so the committed face is a trimmed triangle of the
    /// patch and the spoke at a sector's rim END is the *next*
    /// sector's `u = 0` iso, shared in 3D by the cap's eq (2)
    /// machinery. Pole corners and spoke handoffs are verified within
    /// `tolerance`, measured.
    ///
    /// # Errors
    ///
    /// [`HoleCommitError::UnsupportedSurfaceForm`] when a sector is
    /// not a 5×5 NURBS; [`HoleCommitError::EndpointMismatch`] when a
    /// pole corner misses the pole or consecutive sectors do not hand
    /// their rim over within `tolerance`.
    pub fn from_pns_fill(fill: &NSidedFill, tolerance: f64) -> Result<Self, HoleCommitError> {
        let count = fill.surfaces.len();
        let mut rims = Vec::with_capacity(count);
        let mut start_spokes = Vec::with_capacity(count);
        for (index, surface) in fill.surfaces.iter().enumerate() {
            let nurbs = match surface {
                SurfaceGeometry::Nurbs(nurbs) if nurbs.num_u == 5 && nurbs.num_v == 5 => nurbs,
                other => {
                    return Err(HoleCommitError::UnsupportedSurfaceForm {
                        sector: index,
                        detail: format!(
                            "expected a 5x5 bi-4 NURBS sector, found {}",
                            surface_kind_name(other)
                        ),
                    });
                }
            };
            let cp = |u: usize, v: usize| nurbs.control_points[5 * u + v];
            let pole_distance = cp(0, 0).distance_to(&fill.pole);
            if pole_distance > tolerance {
                return Err(HoleCommitError::EndpointMismatch {
                    sector: index,
                    role: "pole corner (0,0)",
                    distance: pole_distance,
                    budget: tolerance,
                });
            }
            rims.push((cp(0, 4), cp(4, 4)));
            start_spokes.push(bezier_spoke((0..5).map(|v| cp(0, v)).collect()));
        }
        let mut sectors = Vec::with_capacity(count);
        for index in 0..count {
            let next = (index + 1) % count;
            let handoff = rims[index].1.distance_to(&rims[next].0);
            if handoff > tolerance {
                return Err(HoleCommitError::EndpointMismatch {
                    sector: index,
                    role: "rim handoff to next sector",
                    distance: handoff,
                    budget: tolerance,
                });
            }
            sectors.push(HoleCapSectorGeometry {
                surface: fill.surfaces[index].clone(),
                rim_start: rims[index].0,
                rim_end: rims[index].1,
                spoke_at_rim_start: start_spokes[index].clone(),
                spoke_at_rim_end: start_spokes[next].clone(),
            });
        }
        Ok(Self {
            pole: fill.pole,
            sectors,
            chart: HoleCapChart::ProjectableFan,
        })
    }

    /// Normalize a Coons fan.
    ///
    /// The fan's spokes are straight centroid lines that the filler
    /// does not return, so they are rebuilt here from the patch
    /// corners: `corners[0]` / `corners[1]` are the rim endpoints and
    /// the pole is the mean of the rim starts — exactly the centroid
    /// the fan radiated from. The known ~1% `center_near` offset of
    /// each patch's collapsed top edge is not hidden: it surfaces in
    /// the commit's coherence gate as measured L2 disagreement.
    ///
    /// # Errors
    ///
    /// [`HoleCommitError::SectorCountMismatch`] when the fan silently
    /// produced fewer patches than the hole has sides (individual
    /// patch failures are skipped by the filler, never committed
    /// here); [`HoleCommitError::UnsupportedSurfaceForm`] when a
    /// spoke line cannot be constructed.
    pub fn from_coons_patches(
        patches: &[CoonsPatchResult],
        expected_sides: usize,
        tolerance: f64,
    ) -> Result<Self, HoleCommitError> {
        if patches.len() != expected_sides {
            return Err(HoleCommitError::SectorCountMismatch {
                expected: expected_sides,
                got: patches.len(),
            });
        }
        let inverse_count = 1.0 / patches.len() as f64;
        let mut pole = Point3::new(0.0, 0.0, 0.0);
        for patch in patches {
            let start = patch.corners[0];
            pole = Point3::new(
                pole.x + start.x * inverse_count,
                pole.y + start.y * inverse_count,
                pole.z + start.z * inverse_count,
            );
        }
        let mut sectors = Vec::with_capacity(patches.len());
        for (index, patch) in patches.iter().enumerate() {
            let rim_start = patch.corners[0];
            let rim_end = patch.corners[1];
            let spoke = |rim: Point3| {
                CurveOnEdge::line_segment(pole, rim, tolerance).ok_or_else(|| {
                    HoleCommitError::UnsupportedSurfaceForm {
                        sector: index,
                        detail: format!("spoke from pole {pole:?} to rim {rim:?} is degenerate"),
                    }
                })
            };
            sectors.push(HoleCapSectorGeometry {
                surface: SurfaceGeometry::Nurbs(patch.surface.clone()),
                rim_start,
                rim_end,
                spoke_at_rim_start: spoke(rim_start)?,
                spoke_at_rim_end: spoke(rim_end)?,
            });
        }
        Ok(Self {
            pole,
            sectors,
            chart: HoleCapChart::ProjectableFan,
        })
    }
}

/// Wrap a clamped degree-4 Bézier spoke as an edge curve on `[0, 1]`.
fn bezier_spoke(control_points: Vec<Point3>) -> CurveOnEdge {
    let weights = vec![1.0; control_points.len()];
    CurveOnEdge {
        geometry: CurveGeometry::Nurbs(NurbsCurve {
            degree: 4,
            knots: vec![0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            control_points,
            weights,
        }),
        t_start: 0.0,
        t_end: 1.0,
    }
}

fn surface_kind_name(surface: &SurfaceGeometry) -> String {
    match surface {
        SurfaceGeometry::Nurbs(nurbs) => format!("a {}x{} NURBS", nurbs.num_u, nurbs.num_v),
        other => format!("{other:?}").chars().take(40).collect(),
    }
}
