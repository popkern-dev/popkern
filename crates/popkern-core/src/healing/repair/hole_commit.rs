//! Commit a filled hole into topology (PK-199).
//!
//! Both hole fillers — the Coons fan and the G¹ PnS cap — return
//! patch geometry only. This module is the committing side: it
//! normalizes either filler's output into [`HoleCapGeometry`], runs
//! the Euler fan operator [`crate::topology::BRep::cap_boundary_ring_fan`]
//! on a candidate clone, attaches sector surfaces, spoke curves, and
//! recomputed pcurves, gates on coherence and full validation, and
//! only then swaps the candidate in. A fill that cannot commit leaves
//! the input BRep untouched.

use crate::error::{CoreError, ValidationError};
use crate::euler::BoundaryRingFanResult;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::{Point3, SurfaceOnFace, Vec3};
use crate::topology::BRep;
use crate::topology::coherence::CoherenceValidationScope;
use popkern_math::PortableMath;

use super::hole_commit_geometry::{HoleCapChart, HoleCapGeometry};
use super::hole_fill::DetectedHole;
use super::hole_ruled_fan::author_ruled_sector_pcurves;

/// Why a normalized cap could not be committed.
#[derive(Debug, thiserror::Error)]
pub enum HoleCommitError {
    /// The geometry has a different sector count than the hole has sides.
    #[error("hole commit: geometry has {got} sector(s) for a {expected}-sided hole")]
    SectorCountMismatch {
        /// Sides of the detected hole.
        expected: usize,
        /// Sectors supplied by the filler.
        got: usize,
    },
    /// A filler surface was not in the form the normalizer supports.
    #[error("hole commit: sector {sector} surface is not a committable form: {detail}")]
    UnsupportedSurfaceForm {
        /// Index of the offending sector.
        sector: usize,
        /// What was found instead.
        detail: String,
    },
    /// A sector corner or spoke endpoint missed its expected position.
    #[error(
        "hole commit: sector {sector} {role} sits {distance:.3e} from its expected position (budget {budget:.3e})"
    )]
    EndpointMismatch {
        /// Index of the offending sector.
        sector: usize,
        /// Which endpoint missed (e.g. "pole corner", "spoke end").
        role: &'static str,
        /// Measured distance.
        distance: f64,
        /// The budget it had to clear.
        budget: f64,
    },
    /// No committed sector rim matched a filler sector's endpoints.
    #[error(
        "hole commit: committed sector {sector} (rim {start:?} → {end:?}) matched {matches} filler sector(s); exactly one is required"
    )]
    SectorMatchFailed {
        /// Index of the committed sector in fan order.
        sector: usize,
        /// Rim start position.
        start: Point3,
        /// Rim end position.
        end: Point3,
        /// How many filler sectors matched.
        matches: usize,
    },
    /// The topology layer declined (ring validation, pcurve authoring).
    #[error("hole commit: {0}")]
    Topology(#[from] CoreError),
    /// The committed cap failed the coherence gate.
    #[error(
        "hole commit: coherence gate failed on the new faces — {missing_curves} edge curve(s) missing, {missing_pcurves} pcurve(s) missing, worst L2 gap {l2_worst:.3e}"
    )]
    CoherenceRejected {
        /// L1 violations among the new faces' edges.
        missing_curves: usize,
        /// L3 violations among the new faces' half-edges.
        missing_pcurves: usize,
        /// Worst sampled L2 disagreement, 0.0 when none reported.
        l2_worst: f64,
    },
    /// The capped body failed full fast validation.
    #[error("hole commit: post-commit validation failed: {0}")]
    PostValidationFailed(#[from] ValidationError),
}

/// A committed hole cap.
#[derive(Debug, Clone)]
pub struct HoleCapCommit {
    /// The fan topology that now closes the hole.
    pub fan: BoundaryRingFanResult,
    /// Worst L2 disagreement reported by the coherence gate
    /// (0.0 when every sample cleared the coherence tolerance).
    pub worst_l2_gap: f64,
}

/// Commit a normalized hole cap through the Euler fan operator.
///
/// Transactional: every mutation happens on a candidate clone; the
/// input BRep is replaced only after the coherence gate and
/// `validate()` both pass on the candidate.
///
/// `tolerance` is the linear budget used to match filler sectors to
/// committed rim vertices and to verify spoke endpoints. It is
/// explicit, never global state.
///
/// # Errors
///
/// See [`HoleCommitError`]; every decline names the sector, role, and
/// measured distance where applicable, and leaves the BRep untouched.
pub fn commit_hole_cap(
    brep: &mut BRep,
    hole: &DetectedHole,
    geometry: &HoleCapGeometry,
    tolerance: f64,
) -> Result<HoleCapCommit, HoleCommitError> {
    let sides = hole.boundary_edges.len();
    if geometry.sectors.len() != sides {
        return Err(HoleCommitError::SectorCountMismatch {
            expected: sides,
            got: geometry.sectors.len(),
        });
    }

    let mut candidate = brep.clone();
    let fan = candidate.cap_boundary_ring_fan(&hole.boundary_edges, geometry.pole)?;

    let rim_positions = fan_rim_positions(&candidate, &fan)?;
    let assignments = match_sectors(geometry, &rim_positions, tolerance)?;

    attach_sector_geometry(&mut candidate, &fan, geometry, &rim_positions, &assignments)?;
    attach_spoke_curves(
        &mut candidate,
        &fan,
        geometry,
        &rim_positions,
        &assignments,
        tolerance,
    )?;
    match geometry.chart {
        HoleCapChart::RuledFanAnalytic => {
            for &face in &fan.sector_faces {
                author_ruled_sector_pcurves(&mut candidate, face, fan.pole_vertex, tolerance)?;
            }
        }
        HoleCapChart::ProjectableFan => {
            for &face in &fan.sector_faces {
                candidate.recompute_face_pcurves(face)?;
            }
            snap_pole_pcurve_endpoints(&mut candidate, &fan, geometry.pole, tolerance)?;
        }
    }

    let scope = CoherenceValidationScope::Faces(fan.sector_faces.iter().copied().collect());
    let report = candidate.validate_coherence(&scope);
    let worst_l2_gap = report
        .l2_face_boundary_disagreements
        .iter()
        .map(|d| d.worst_gap)
        .fold(0.0_f64, f64::max);
    if !report.is_empty() {
        return Err(HoleCommitError::CoherenceRejected {
            missing_curves: report.l1_missing_edge_curves.len(),
            missing_pcurves: report.l3_missing_pcurves.len(),
            l2_worst: worst_l2_gap,
        });
    }
    candidate.validate()?;

    *brep = candidate;
    Ok(HoleCapCommit { fan, worst_l2_gap })
}

/// Ring vertex positions in fan order: `positions[k]` is `w_k`, the
/// origin of rim half-edge k.
fn fan_rim_positions(
    brep: &BRep,
    fan: &BoundaryRingFanResult,
) -> Result<Vec<Point3>, HoleCommitError> {
    fan.rim_half_edges
        .iter()
        .map(|&rim| {
            let origin = brep
                .half_edges
                .get(rim)
                .map(|half_edge| half_edge.origin)
                .ok_or_else(|| CoreError::InvalidHandle("fan rim half-edge is not live".into()))?;
            brep.vertices
                .get(origin)
                .map(|vertex| vertex.position)
                .ok_or_else(|| CoreError::InvalidHandle("fan rim vertex is not live".into()).into())
        })
        .collect()
}

/// Per committed sector: index of the matching filler sector and
/// whether it matched forward (`rim_start` at `w_k`).
fn match_sectors(
    geometry: &HoleCapGeometry,
    rim_positions: &[Point3],
    tolerance: f64,
) -> Result<Vec<(usize, bool)>, HoleCommitError> {
    let count = rim_positions.len();
    let mut used = vec![false; count];
    let mut assignments = Vec::with_capacity(count);
    for k in 0..count {
        let w_start = rim_positions[k];
        let w_end = rim_positions[(k + 1) % count];
        let mut matches = Vec::new();
        for (index, sector) in geometry.sectors.iter().enumerate() {
            if used[index] {
                continue;
            }
            let forward = sector.rim_start.distance_to(&w_start) <= tolerance
                && sector.rim_end.distance_to(&w_end) <= tolerance;
            let reversed = sector.rim_start.distance_to(&w_end) <= tolerance
                && sector.rim_end.distance_to(&w_start) <= tolerance;
            if forward || reversed {
                matches.push((index, forward));
            }
        }
        if matches.len() != 1 {
            return Err(HoleCommitError::SectorMatchFailed {
                sector: k,
                start: w_start,
                end: w_end,
                matches: matches.len(),
            });
        }
        used[matches[0].0] = true;
        assignments.push(matches[0]);
    }
    Ok(assignments)
}

fn attach_sector_geometry(
    brep: &mut BRep,
    fan: &BoundaryRingFanResult,
    geometry: &HoleCapGeometry,
    rim_positions: &[Point3],
    assignments: &[(usize, bool)],
) -> Result<(), HoleCommitError> {
    let count = rim_positions.len();
    for (k, &(sector_index, _)) in assignments.iter().enumerate() {
        let sector = &geometry.sectors[sector_index];
        // The sector loop runs rim (w_k → w_{k+1}), spoke to pole,
        // spoke back — its Newell normal is the face normal the
        // orientation flag is measured against.
        let newell = newell_normal(&[
            rim_positions[k],
            rim_positions[(k + 1) % count],
            geometry.pole,
        ]);
        let surface_normal = sector.surface.normal_at(0.5, 0.5);
        let orientation = newell.dot(&surface_normal) >= 0.0;
        let face = brep
            .faces
            .get_mut(fan.sector_faces[k])
            .ok_or_else(|| CoreError::InvalidHandle("fan sector face is not live".into()))?;
        face.surface = Some(SurfaceOnFace {
            geometry: sector.surface.clone(),
            orientation,
        });
    }
    Ok(())
}

fn attach_spoke_curves(
    brep: &mut BRep,
    fan: &BoundaryRingFanResult,
    geometry: &HoleCapGeometry,
    rim_positions: &[Point3],
    assignments: &[(usize, bool)],
    tolerance: f64,
) -> Result<(), HoleCommitError> {
    for (k, &(sector_index, forward)) in assignments.iter().enumerate() {
        let sector = &geometry.sectors[sector_index];
        // Spoke k connects w_k to the pole; sector k carries it as
        // the spoke at whichever rim endpoint sits at w_k.
        let spoke_curve = if forward {
            &sector.spoke_at_rim_start
        } else {
            &sector.spoke_at_rim_end
        };
        let eval = &spoke_curve.geometry as &dyn crate::geometry::curve::CurveEval;
        let start = eval.point_at(spoke_curve.t_start);
        let end = eval.point_at(spoke_curve.t_end);
        let pole_gap = start.distance_to(&geometry.pole);
        let rim_gap = end.distance_to(&rim_positions[k]);
        let budget = tolerance;
        if pole_gap > budget || rim_gap > budget {
            return Err(HoleCommitError::EndpointMismatch {
                sector: sector_index,
                role: "spoke end",
                distance: pole_gap.max(rim_gap),
                budget,
            });
        }
        let edge = brep
            .edges
            .get_mut(fan.spoke_edges[k])
            .ok_or_else(|| CoreError::InvalidHandle("fan spoke edge is not live".into()))?;
        edge.curve = Some(spoke_curve.clone());
    }
    Ok(())
}

/// Snap the pole end of every new pcurve to the sector surface's
/// exact pole preimage.
///
/// Fan sector patches can overlap (measured on the PnS cap: sector
/// quads reach past their committed triangle), which makes the
/// point-to-surface distance field multi-modal. Projection-derived
/// pcurves then land the pole endpoint in a wrong local minimum —
/// measured 9.85e-2 off on the open-box fill while the true corner
/// preimage evaluates to the pole exactly. Every sane fan filler puts
/// the pole at a corner of the sector's parameter domain, so the
/// preimage is recovered by evaluating the domain corners and
/// declining, measured, when none lands on the pole. The coherence
/// gate re-verifies the result afterwards; this snap is never the
/// acceptance authority.
fn snap_pole_pcurve_endpoints(
    brep: &mut BRep,
    fan: &BoundaryRingFanResult,
    pole: Point3,
    tolerance: f64,
) -> Result<(), HoleCommitError> {
    for (sector, &face) in fan.sector_faces.iter().enumerate() {
        let Some(surface) = brep
            .faces
            .get(face)
            .and_then(|face_data| face_data.surface.as_ref())
            .map(|on_face| on_face.geometry.clone())
        else {
            continue;
        };
        let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
        if !(u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite()) {
            continue;
        }
        let corners = [
            (u_min, v_min),
            (u_min, v_max),
            (u_max, v_min),
            (u_max, v_max),
        ];
        // A sector can have MORE than one pole-preimage corner — the
        // polar-wedge chart collapses a whole domain edge at the pole
        // (PK-248), so both of that edge's corners evaluate to it.
        // Each pcurve endpoint must snap to ITS OWN nearest preimage;
        // snapping every spoke to one global corner would drag one
        // spoke's pcurve across the collapsed edge.
        let pole_corners: Vec<(f64, f64)> = corners
            .iter()
            .copied()
            .filter(|&(u, v)| surface.point_at(u, v).distance_to(&pole) <= tolerance)
            .collect();
        if pole_corners.is_empty() {
            let nearest = corners
                .iter()
                .map(|&(u, v)| surface.point_at(u, v).distance_to(&pole))
                .fold(f64::INFINITY, f64::min);
            return Err(HoleCommitError::EndpointMismatch {
                sector,
                role: "pole preimage (no domain corner evaluates to the pole)",
                distance: nearest,
                budget: tolerance,
            });
        }

        let loop_ = brep.faces[face].outer_loop;
        let start = brep.loops[loop_].half_edge;
        let mut cursor = start;
        loop {
            let half_edge = &brep.half_edges[cursor];
            let next = half_edge.next;
            let touches_pole = half_edge.origin == fan.pole_vertex
                || brep.half_edges[next].origin == fan.pole_vertex;
            if let Some(nurbs) = touches_pole
                .then(|| brep.half_edges[cursor].pcurve.as_mut())
                .flatten()
                .and_then(|pcurve| pcurve.as_nurbs_mut())
                && let (Some(first), Some(last)) = (
                    nurbs.control_points.first().copied(),
                    nurbs.control_points.last().copied(),
                )
            {
                let first_gap = surface.point_at(first.x, first.y).distance_to(&pole);
                let last_gap = surface.point_at(last.x, last.y).distance_to(&pole);
                let snap_index = if first_gap <= last_gap {
                    0
                } else {
                    nurbs.control_points.len() - 1
                };
                let current = nurbs.control_points[snap_index];
                let nearest = pole_corners
                    .iter()
                    .copied()
                    .min_by(|a, b| {
                        let da = (a.0 - current.x).phypot(a.1 - current.y);
                        let db = (b.0 - current.x).phypot(b.1 - current.y);
                        da.total_cmp(&db)
                    })
                    .unwrap_or(pole_corners[0]);
                nurbs.control_points[snap_index] = Point3::new(nearest.0, nearest.1, 0.0);
            }
            cursor = next;
            if cursor == start {
                break;
            }
        }
    }
    Ok(())
}

/// Newell plane normal of a polygon; not normalized.
fn newell_normal(points: &[Point3]) -> Vec3 {
    let mut normal = Vec3::new(0.0, 0.0, 0.0);
    for (index, current) in points.iter().enumerate() {
        let next = &points[(index + 1) % points.len()];
        normal.x += (current.y - next.y) * (current.z + next.z);
        normal.y += (current.z - next.z) * (current.x + next.x);
        normal.z += (current.x - next.x) * (current.y + next.y);
    }
    normal
}

#[cfg(test)]
#[path = "hole_commit_tests.rs"]
mod tests;
