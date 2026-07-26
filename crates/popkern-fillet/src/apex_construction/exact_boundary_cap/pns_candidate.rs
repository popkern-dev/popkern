//! A PnS bi-4 candidate for the exact-boundary cap route.
//!
//! # Why
//!
//! The exact-boundary route covers a corner with a fan of Coons
//! wedges, one per sector. Each wedge meets its supports well — that
//! is what acceptance gates — but adjacent wedges share only a
//! boundary curve, not a tangent plane. Measured on `n_gonal_pyramid`
//! mixed-radius apexes, the internal cap-cap seams run up to 89° of
//! normal misalignment with 47–60° *mid-seam*, while the route
//! reports `achieved_continuity = G1` (PK-197). The G1 claim comes
//! from gating cap-vs-support coplanarity; cap-vs-cap is never
//! measured.
//!
//! The PnS bi-4 cap is G¹ across its sector seams **by construction**
//! — that is exactly what the paper's eq (2) enforces — so it is the
//! natural replacement for the fan on this class.
//!
//! # How the boundary maps
//!
//! Each exact-boundary sector owns two edges: a fillet-end arc (with
//! the fillet strip beside it) and a face bridge (with the planar or
//! curved support beside it). The PnS cap takes one sector per
//! *pair*, spanning arc-then-bridge, which lines up one-to-one with
//! the `N` cap faces the route already commits:
//!
//! ```text
//!   FillBoundaryEdge.start  = arc start        (fillet support)
//!   FillBoundaryEdge.mid    = arc/bridge join  (fillet support)
//!   FillBoundaryEdge.end    = bridge end       (bridge support)
//! ```
//!
//! # Status
//!
//! This is a **candidate**, not the dispatched path. PK-192's rule is
//! that the measured numbers must win before the route changes. They
//! do — measured over the same corner rim, the shipped fan's worst
//! internal-seam normal misalignment is 89.9° (n=3) and 89.6° (n=4)
//! while this cap's is 0.0° — see
//! `pns_cap_internal_seams_beat_the_shipped_fan`.
//!
//! The substitution happens at the **solved-surfaces** stage, before
//! `validate_exact_boundary_candidate` — so the PnS surfaces go
//! through the route's own acceptance gates, topology validation and
//! transactional rollback exactly as the Coons wedges did. Nothing is
//! overwritten at commit: a dispatch that skipped acceptance would be
//! precisely the dishonest-green failure this work exists to remove.
//!
//! If the ring cannot be built, or the fill declines, the Coons
//! surfaces are left in place and the reason is recorded — the cap
//! still ships, just without the continuity upgrade.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{FillBoundaryEdge, NSidedFillError, Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use super::sector::ExactBoundaryCapSectorInput;

/// Why a PnS candidate could not be built for a corner.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub(crate) enum PnsCandidateError {
    /// A sector's boundary edge carried no curve.
    #[error("PnS candidate: sector {sector} edge {which} has no curve geometry")]
    MissingCurve {
        /// Sector index.
        sector: usize,
        /// Which edge ("fillet_end" or "face").
        which: &'static str,
    },
    /// A support face carried no surface, so no tangent plane exists.
    #[error("PnS candidate: sector {sector} support face for {which} has no surface")]
    MissingSupport {
        /// Sector index.
        sector: usize,
        /// Which edge's support.
        which: &'static str,
    },
    /// A boundary sample would not project onto its support.
    #[error("PnS candidate: sector {sector} sample on {which} does not project onto its support")]
    ProjectionFailed {
        /// Sector index.
        sector: usize,
        /// Which edge.
        which: &'static str,
    },
    /// The two edges of a sector do not meet.
    #[error(
        "PnS candidate: sector {sector}'s arc and bridge do not share an endpoint (gap {gap:.3e})"
    )]
    SectorNotJoined {
        /// Sector index.
        sector: usize,
        /// Measured gap.
        gap: f64,
    },
    /// The fill primitive declined.
    #[error("PnS candidate: {0}")]
    Fill(#[from] NSidedFillError),
}

/// Build the PnS boundary ring for a corner from its exact-boundary
/// sectors.
///
/// Returns one [`FillBoundaryEdge`] per sector, in sector order.
pub(crate) fn pns_ring_from_sectors(
    brep: &BRep,
    sectors: &[ExactBoundaryCapSectorInput],
    tolerance: f64,
) -> Result<Vec<FillBoundaryEdge>, PnsCandidateError> {
    let mut ring = Vec::with_capacity(sectors.len());
    for sector in sectors {
        ring.push(sector_to_boundary_edge(brep, sector, tolerance)?);
    }
    chain_ring(&mut ring, tolerance)?;
    Ok(ring)
}

/// Orient the sides so each one continues from the last.
///
/// Two orientations are involved and they are independent. Inside a
/// sector, the arc and the bridge have to meet — that is resolved when
/// the side is built. Around the ring, side `k`'s end has to meet side
/// `k+1`'s start, and nothing about sector order guarantees it: a
/// sector whose internal join is perfect can still run the wrong way
/// round the corner. Without this pass the ring read as open by
/// 2.4e-2 on a pyramid apex and the whole substitution declined.
fn chain_ring(ring: &mut [FillBoundaryEdge], tolerance: f64) -> Result<(), PnsCandidateError> {
    let n = ring.len();
    if n < 2 {
        return Ok(());
    }
    // Seed on whichever end of side 0 is closer to side 1.
    let touch =
        |p: Point3, e: &FillBoundaryEdge| p.distance_to(&e.start).min(p.distance_to(&e.end));
    if touch(ring[0].end, &ring[1]) > touch(ring[0].start, &ring[1]) {
        reverse_side(&mut ring[0]);
    }
    for k in 1..n {
        let prev_end = ring[k - 1].end;
        if prev_end.distance_to(&ring[k].end) < prev_end.distance_to(&ring[k].start) {
            reverse_side(&mut ring[k]);
        }
        let gap = prev_end.distance_to(&ring[k].start);
        if gap > tolerance {
            return Err(PnsCandidateError::SectorNotJoined { sector: k, gap });
        }
    }
    let closing = ring[n - 1].end.distance_to(&ring[0].start);
    if closing > tolerance {
        return Err(PnsCandidateError::SectorNotJoined {
            sector: 0,
            gap: closing,
        });
    }
    Ok(())
}

/// Reverse a side: swap its ends and negate its tangents.
fn reverse_side(side: &mut FillBoundaryEdge) {
    std::mem::swap(&mut side.start, &mut side.end);
    std::mem::swap(&mut side.normal_start, &mut side.normal_end);
    std::mem::swap(&mut side.tangent_start, &mut side.tangent_end);
    side.tangent_start = side.tangent_start * -1.0;
    side.tangent_mid = side.tangent_mid * -1.0;
    side.tangent_end = side.tangent_end * -1.0;
}

/// Sample one sector's arc-then-bridge span into a PnS boundary side.
fn sector_to_boundary_edge(
    brep: &BRep,
    sector: &ExactBoundaryCapSectorInput,
    tolerance: f64,
) -> Result<FillBoundaryEdge, PnsCandidateError> {
    let idx = sector.sector_index;
    let arc = edge_span(brep, sector.fillet_end_edge, idx, "fillet_end")?;
    let bridge = edge_span(brep, sector.face_edge, idx, "face")?;
    let arc_support = support_surface(brep, sector.fillet_support_face).ok_or(
        PnsCandidateError::MissingSupport {
            sector: idx,
            which: "fillet_end",
        },
    )?;
    let bridge_support = support_surface(brep, sector.bridge_support_face).ok_or(
        PnsCandidateError::MissingSupport {
            sector: idx,
            which: "face",
        },
    )?;

    // Orient the pair so the arc's end meets the bridge's start. Edge
    // direction is stored per edge and says nothing about loop order,
    // so either or both may run backwards.
    let (arc_start, join_a, arc_end_t) = (arc.start, arc.end, arc.t_end);
    let (b_start, b_end) = (bridge.start, bridge.end);
    let candidates = [
        (
            arc_start,
            join_a,
            b_end,
            arc.t_start,
            arc_end_t,
            bridge.t_end,
        ),
        (
            arc_start,
            join_a,
            b_start,
            arc.t_start,
            arc_end_t,
            bridge.t_start,
        ),
        (
            join_a,
            arc_start,
            b_end,
            arc.t_end,
            arc.t_start,
            bridge.t_end,
        ),
        (
            join_a,
            arc_start,
            b_start,
            arc.t_end,
            arc.t_start,
            bridge.t_start,
        ),
    ];
    let (start, mid, end, t_start, t_mid, t_end) = candidates
        .into_iter()
        .min_by(|a, b| {
            let ga = a.1.distance_to(&nearest_bridge_end(a.1, b_start, b_end));
            let gb = b.1.distance_to(&nearest_bridge_end(b.1, b_start, b_end));
            ga.partial_cmp(&gb).unwrap_or(std::cmp::Ordering::Equal)
        })
        .unwrap_or(candidates[0]);
    let join_gap = mid.distance_to(&nearest_bridge_end(mid, b_start, b_end));
    if join_gap > tolerance {
        return Err(PnsCandidateError::SectorNotJoined {
            sector: idx,
            gap: join_gap,
        });
    }

    let (n_start, tan_start) = frame(&arc, arc_support, t_start, idx, "fillet_end")?;
    let (n_mid, tan_mid) = frame(&arc, arc_support, t_mid, idx, "fillet_end")?;
    let (n_end, tan_end) = frame(&bridge, bridge_support, t_end, idx, "face")?;

    Ok(FillBoundaryEdge {
        start,
        mid,
        end,
        normal_start: n_start,
        normal_mid: n_mid,
        normal_end: n_end,
        tangent_start: tan_start,
        tangent_mid: tan_mid,
        tangent_end: tan_end,
    })
}

fn nearest_bridge_end(p: Point3, a: Point3, b: Point3) -> Point3 {
    if p.distance_to(&a) <= p.distance_to(&b) {
        a
    } else {
        b
    }
}

/// A boundary edge's endpoints and parameter range.
struct EdgeSpan {
    start: Point3,
    end: Point3,
    t_start: f64,
    t_end: f64,
    curve: popkern_core::geometry::curve::CurveOnEdge,
}

fn edge_span(
    brep: &BRep,
    edge: EdgeHandle,
    sector: usize,
    which: &'static str,
) -> Result<EdgeSpan, PnsCandidateError> {
    let curve = brep
        .edges
        .get(edge)
        .and_then(|e| e.curve.clone())
        .ok_or(PnsCandidateError::MissingCurve { sector, which })?;
    let eval = |t: f64| (&curve.geometry as &dyn CurveEval).point_at(t);
    Ok(EdgeSpan {
        start: eval(curve.t_start),
        end: eval(curve.t_end),
        t_start: curve.t_start,
        t_end: curve.t_end,
        curve,
    })
}

fn support_surface(brep: &BRep, face: FaceHandle) -> Option<&SurfaceGeometry> {
    Some(&brep.faces.get(face)?.surface.as_ref()?.geometry)
}

/// Surface normal and curve tangent at a parameter on a boundary edge.
fn frame(
    span: &EdgeSpan,
    support: &SurfaceGeometry,
    t: f64,
    sector: usize,
    which: &'static str,
) -> Result<(Vec3, Vec3), PnsCandidateError> {
    let eval = |t: f64| (&span.curve.geometry as &dyn CurveEval).point_at(t);
    let point = eval(t);
    let (u, v) = support
        .closest_point(&point)
        .ok_or(PnsCandidateError::ProjectionFailed { sector, which })?;
    let normal = support
        .normal_at(u, v)
        .normalized()
        .ok_or(PnsCandidateError::ProjectionFailed { sector, which })?;
    let dt = 1e-6 * (span.t_end - span.t_start).abs().max(1e-12);
    let lo = (t - dt).max(span.t_start.min(span.t_end));
    let hi = (t + dt).min(span.t_start.max(span.t_end));
    let tangent = (eval(hi) - eval(lo))
        .normalized()
        .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    Ok((normal, tangent))
}

/// Replace solved Coons sector surfaces with PnS bi-4 surfaces over
/// the same rim.
///
/// Returns `Ok(n)` with the number of surfaces substituted, or the
/// reason none were. On any error the caller must leave `solved`
/// untouched — this function only mutates it once the whole ring has
/// been built and filled successfully, so a partial substitution
/// (some sectors PnS, some Coons, with no continuity between them)
/// cannot occur.
pub(crate) fn substitute_pns_sector_surfaces(
    brep: &BRep,
    sectors: &[ExactBoundaryCapSectorInput],
    solved: &mut super::solve::ExactBoundarySolvedSurfaces,
    tolerance: f64,
) -> Result<usize, PnsCandidateError> {
    let ring = pns_ring_from_sectors(brep, sectors, tolerance)?;
    let fill = popkern_core::geometry::fill_n_sided_boundary(&ring, tolerance)?;
    if fill.surfaces.len() != sectors.len() {
        return Err(PnsCandidateError::Fill(
            NSidedFillError::UnsupportedSideCount {
                got: fill.surfaces.len(),
            },
        ));
    }

    // Map by sector_index rather than position: `solved.surfaces` is
    // in solver order, which is not guaranteed to be planner order.
    let mut substituted = 0usize;
    for (k, sector) in sectors.iter().enumerate() {
        let SurfaceGeometry::Nurbs(nurbs) = &fill.surfaces[k] else {
            continue;
        };
        if let Some(target) = solved
            .surfaces
            .iter_mut()
            .find(|s| s.sector_index == sector.sector_index)
        {
            target.surface = nurbs.clone();
            substituted += 1;
        }
    }
    Ok(substituted)
}

#[cfg(test)]
#[path = "pns_candidate_tests.rs"]
mod tests;

/// Whether the PnS cap replaces the Coons fan on exact-boundary
/// corners.
///
/// Default **on**: the substitution is measured to remove up to 89.9°
/// of internal-seam normal misalignment and goes through the route's
/// own acceptance gates, so a corner that would have shipped a folded
/// fan either ships a tangent cap or ships exactly what it did before.
/// `POPKERN_PNS_CAP=0` restores the fan for A/B measurement.
pub(crate) fn pns_cap_substitution_enabled() -> bool {
    std::env::var("POPKERN_PNS_CAP")
        .map(|v| v != "0")
        .unwrap_or(true)
}
