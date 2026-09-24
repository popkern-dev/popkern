//! PK-522 mitre-corner settlement: corners of three or more
//! equal-radius strips closed by trimming every adjacent strip pair
//! to its exact mutual intersection.
//!
//! **Problem.** A multi-fillet corner that gets no corner blend
//! (`VertexMode::Mitred`) keeps the face-trim residue. Each support
//! face carries a flap: a quad `c → v1 → X → v2` through the original
//! corner `c`, where `X` is the rail crossing of the two strips that
//! meet on that face and `v1`, `v2` are scaffold vertices on one strip
//! each. The flap edges lie on the support plane but not on the strips
//! (unit-cube corner, r = 0.05: 3.5e-2 off-surface), so the corner
//! resolution gate (PK-516) declines the fillet.
//!
//! **Mitre, by construction.** Equal-radius cylinder strips around a
//! corner have concurrent axes: each axis is the locus at distance r
//! from its two support planes, so all of them pass through the point
//! at distance r from every support plane. Two such cylinders meet in
//! bisector-plane ellipses (PK-152). A point on all N strips is the
//! mitre point P: each adjacent strip pair's seam runs from its rail
//! crossing X to P, and each strip ends in a V between its two
//! seams. For N = 3 a mitre point always exists (the bisector planes
//! share a line through the axes' common point, and distance to every
//! axis grows along it from zero). For N ≥ 4 it exists only for
//! symmetric corners such as a regular pyramid apex. Otherwise the
//! corner declines `NoCommonPoint`.
//!
//! **Surgery, Euler operators only.** For each flap, `kef` on its
//! `v1 → X` edge merges it into the strip beyond. Each flap is merged
//! in the same rotational sense, so each scaffold vertex is left with
//! two edges flanked by the same two strips. `kve` contracts it into
//! `c`. The corner is left as `c` joined to every rail crossing by
//! one edge between the right strip pair. `c` moves to P, and every
//! such edge is re-authored as its exact ellipse arc with exact
//! sinusoidal pcurves in both strip charts. Counts: V −N, E −2N,
//! F −N (χ unchanged).
//!
//! **Invariants:** strips keep their face identity, and support faces
//! and rails are untouched. Every seam geometry is authored and gated
//! BEFORE surgery. **Fallback:** a decline restores the corner's
//! snapshot and is recorded, typed, on the corner. The
//! corner-resolution gate then declines the fillet with that reason.
//! **Tolerance:** the caller's blend tolerance, local and explicit.

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveOnEdge, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::types::{MitreCornerDecline, MitreCornerDeclineReason as Reason, MitreCornerReport};

use super::blend_seam_corner::strip_cylinder;
use super::blend_seam_corner_geometry as geometry;
use super::blend_seam_corner_repair::{commit_edge_authoring, neighbour_chart_point};
use super::corner_resolution::{UnresolvedCorner, UnresolvedCornerCause};
use super::terminal_arc_trim::StripCylinder;

#[path = "blend_mitre_corner_point.rs"]
mod point;
#[path = "blend_mitre_corner_surgery.rs"]
mod surgery;

/// Mitre every requested-mitre corner that is still unresolved.
/// Corners that serve are removed from `corners`; declines stay, with
/// the reason recorded in the returned report.
pub(crate) fn settle_mitre_corners(
    brep: &mut BRep,
    corners: &mut Vec<UnresolvedCorner>,
    tolerance: &Tolerance,
) -> MitreCornerReport {
    let mut report = MitreCornerReport::default();
    corners.retain_mut(|corner| {
        if corner.cause != UnresolvedCornerCause::MitredRequested
            || brep.vertices.get(corner.vertex).is_none()
        {
            return true;
        }
        let snapshot = brep.clone();
        match mitre_corner(brep, corner.vertex, tolerance) {
            Ok(()) => {
                report.corners_built += 1;
                false
            }
            Err(reason) => {
                *brep = snapshot;
                corner.mitre_decline = Some(reason);
                report.declines.push(MitreCornerDecline {
                    corner: corner.vertex,
                    reason,
                });
                true
            }
        }
    });
    report
}

/// One seam authored before surgery.
struct StagedSeam {
    rail_crossing: VertexHandle,
    strips: (FaceHandle, FaceHandle),
    frame: geometry::SeamFrame,
}

/// The strips around `corner`: every non-planar face bordering the
/// corner's face star. Each must be a cylinder.
fn corner_strips(
    brep: &BRep,
    corner: VertexHandle,
) -> Result<Vec<(FaceHandle, StripCylinder)>, Reason> {
    let face_of = |he: HalfEdgeHandle| -> Option<FaceHandle> {
        Some(brep.loops.get(brep.half_edges.get(he)?.loop_)?.face)
    };
    let star: Vec<FaceHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, he)| he.origin == corner)
        .filter_map(|(h, _)| face_of(h))
        .collect();
    let mut strips: Vec<(FaceHandle, StripCylinder)> = Vec::new();
    for &face in &star {
        let outer = brep
            .faces
            .get(face)
            .ok_or(Reason::UnrecognizedResidue)?
            .outer_loop;
        let start = brep
            .loops
            .get(outer)
            .ok_or(Reason::UnrecognizedResidue)?
            .half_edge;
        let mut cursor = start;
        for _ in 0..=brep.half_edges.len() {
            let he = brep
                .half_edges
                .get(cursor)
                .ok_or(Reason::UnrecognizedResidue)?;
            let neighbour = face_of(he.twin).ok_or(Reason::UnrecognizedResidue)?;
            let planar = brep
                .faces
                .get(neighbour)
                .and_then(|f| f.surface.as_ref())
                .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Plane { .. }));
            if !star.contains(&neighbour) && !planar && strips.iter().all(|(f, _)| *f != neighbour)
            {
                let cylinder =
                    strip_cylinder(brep, neighbour).ok_or(Reason::NonCylindricalStrip)?;
                strips.push((neighbour, cylinder));
            }
            cursor = he.next;
            if cursor == start {
                break;
            }
        }
    }
    if strips.len() < 3 {
        return Err(Reason::UnrecognizedResidue);
    }
    Ok(strips)
}

fn mitre_corner(
    brep: &mut BRep,
    corner: VertexHandle,
    tolerance: &Tolerance,
) -> Result<(), Reason> {
    // Typed geometric verdicts first: they hold whatever shape the
    // residue took.
    let cylinders = corner_strips(brep, corner)?;
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let radius = cylinders[0].1.radius;
    if cylinders
        .iter()
        .any(|(_, c)| (c.radius - radius).abs() > pos_tol.max(radius * 1e-9))
    {
        return Err(Reason::MismatchedRadii);
    }
    let corner_pos = brep
        .vertices
        .get(corner)
        .ok_or(Reason::UnrecognizedResidue)?
        .position;
    let axes: Vec<&StripCylinder> = cylinders.iter().map(|(_, c)| c).collect();
    let mitre_point =
        point::solve_mitre_point(&axes, corner_pos, pos_tol).ok_or(Reason::NoCommonPoint)?;

    // Phase 1: reduce the residue to an N-gon of rail crossings.
    let strip_faces: Vec<FaceHandle> = cylinders.iter().map(|(f, _)| *f).collect();
    let polygon = surgery::reduce_to_residue_polygon(brep, corner, &strip_faces)?;

    // Author every seam frame before building the star.
    let cyl_of = |face: FaceHandle| -> Result<&StripCylinder, Reason> {
        cylinders
            .iter()
            .find(|(f, _)| *f == face)
            .map(|(_, c)| c)
            .ok_or(Reason::UnrecognizedResidue)
    };
    let n = polygon.crossings.len();
    let mut staged: Vec<StagedSeam> = Vec::with_capacity(n);
    for (k, &crossing) in polygon.crossings.iter().enumerate() {
        let strips = (polygon.sides[(k + n - 1) % n], polygon.sides[k]);
        let x = brep
            .vertices
            .get(crossing)
            .ok_or(Reason::UnrecognizedResidue)?
            .position;
        let frame = geometry::solve_equal_radius_seam_frame(
            cyl_of(strips.0)?,
            cyl_of(strips.1)?,
            x,
            mitre_point,
            tolerance,
        )
        .map_err(Reason::SeamUnresolved)?;
        geometry::author_seam_curve(&frame, x, mitre_point, corner_pos, tolerance)
            .map_err(Reason::SeamUnresolved)?;
        staged.push(StagedSeam {
            rail_crossing: crossing,
            strips,
            frame,
        });
    }

    // Phase 2: star the polygon from the mitre point.
    let mitre = surgery::star_from_mitre_point(brep, &polygon, mitre_point)?;
    for seam in &staged {
        author_seam_edge(brep, mitre, seam, corner_pos, tolerance)?;
    }
    brep.validate().map_err(|_| Reason::ValidationFailed)?;
    Ok(())
}

/// The half-edge `from → to`, when an edge joins them.
fn half_edge_between(brep: &BRep, from: VertexHandle, to: VertexHandle) -> Option<HalfEdgeHandle> {
    brep.half_edges.iter().find_map(|(h, he)| {
        let twin = brep.half_edges.get(he.twin)?;
        (he.origin == from && twin.origin == to).then_some(h)
    })
}

/// Re-author the surviving corner edge to `seam.rail_crossing` as the
/// exact seam arc, with a pcurve in each strip, and gate it.
fn author_seam_edge(
    brep: &mut BRep,
    corner: VertexHandle,
    seam: &StagedSeam,
    stale_corner: Point3,
    tolerance: &Tolerance,
) -> Result<(), Reason> {
    let he_out =
        half_edge_between(brep, corner, seam.rail_crossing).ok_or(Reason::EulerRejected)?;
    let he_in = brep
        .half_edges
        .get(he_out)
        .ok_or(Reason::EulerRejected)?
        .twin;
    let edge: EdgeHandle = brep
        .half_edges
        .get(he_out)
        .ok_or(Reason::EulerRejected)?
        .edge;
    let pos_of = |v: VertexHandle| brep.vertices.get(v).map(|v| v.position);
    let p = pos_of(corner).ok_or(Reason::EulerRejected)?;
    let x = pos_of(seam.rail_crossing).ok_or(Reason::EulerRejected)?;
    let primary = brep.edges.get(edge).ok_or(Reason::EulerRejected)?.half_edge;
    let (p_start, p_end) = if primary == he_out { (p, x) } else { (x, p) };
    let curve: CurveOnEdge =
        geometry::author_seam_curve(&seam.frame, p_start, p_end, stale_corner, tolerance)
            .map_err(Reason::SeamUnresolved)?;

    let mut pcurves: Vec<(FaceHandle, HalfEdgeHandle, Pcurve)> = Vec::new();
    for he in [he_out, he_in] {
        let h = brep.half_edges.get(he).ok_or(Reason::EulerRejected)?;
        let face = brep.loops.get(h.loop_).ok_or(Reason::EulerRejected)?.face;
        if face != seam.strips.0 && face != seam.strips.1 {
            return Err(Reason::ValidationFailed);
        }
        let surface = brep
            .faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| s.geometry.clone())
            .ok_or(Reason::EulerRejected)?;
        // The rail crossing's chart point comes from the rail edge
        // beside it; the mitre point is unwrapped onto that branch.
        let x_uv = if h.origin == seam.rail_crossing {
            neighbour_chart_point(brep, face, h.prev, true)
        } else {
            neighbour_chart_point(brep, face, h.next, false)
        }
        .ok_or(Reason::SeamUnresolved(
            crate::types::BlendSeamCornerDeclineReason::PcurveClosureFailed,
        ))?;
        let p_uv = point::chart_point_near(&surface, p, x_uv);
        let (origin_pos, dest_pos, origin_uv, dest_uv) = if h.origin == seam.rail_crossing {
            (x, p, Some(x_uv), p_uv)
        } else {
            (p, x, p_uv, Some(x_uv))
        };
        let pcurve = geometry::author_side_pcurve(
            &surface,
            &seam.frame,
            origin_pos,
            dest_pos,
            &geometry::ChartAnchor { uv: origin_uv },
            &geometry::ChartAnchor { uv: dest_uv },
            tolerance,
        )
        .map_err(Reason::SeamUnresolved)?;
        pcurves.push((face, he, pcurve));
    }
    commit_edge_authoring(brep, edge, curve, pcurves);
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    if brep.measure_edge_deviation(edge) > pos_tol {
        return Err(Reason::ValidationFailed);
    }
    Ok(())
}
