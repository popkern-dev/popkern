//! Mitre-corner settlement (PK-522, PK-528): corners of three or more
//! fillet strips closed by trimming every adjacent strip pair to its
//! mutual intersection, all seams meeting at the strips' common point.
//!
//! **Problem.** A multi-fillet corner that gets no corner blend
//! (`VertexMode::Mitred`) keeps the face-trim residue. That is a star
//! of support-plane flaps through the original corner, whose edges lie
//! on the support planes but not on the strips (unit-cube corner,
//! r = 0.05: 3.5e-2 off-surface). The corner-resolution gate (PK-516)
//! declines it.
//!
//! **Mitre point.** P is the point near the corner on every strip,
//! and every adjacent strip pair's seam runs from its rail crossing to
//! P. Equal-radius cylinder strips have concurrent axes: each axis
//! lies at distance r from its two support planes. P is then solved
//! against the axes (`point`), and each seam is an exact
//! bisector-plane ellipse (PK-152). Any other strip set (mismatched
//! radii, curvature-continuous or variable-radius NURBS) has P solved
//! against the surfaces themselves, with seams marched (`marched`).
//! For N = 3, P generically exists. For N ≥ 4 it exists only at
//! symmetric corners; otherwise the corner declines `NoCommonPoint`.
//!
//! **Surgery, Euler operators only** (`surgery`). The flaps merge and
//! reduce to an N-gon of rail crossings, and the N-gon is starred from
//! P. Whatever shape the residue took (quad flaps at orthogonal
//! corners, stub chains at oblique ones), it reduces the same way.
//!
//! **Invariants:** strips keep their face identity, and support faces
//! and rails are untouched. The mitre point and every seam geometry
//! are solved and gated before the star is built. **Fallback:** a
//! decline restores the corner's snapshot and is recorded, typed, on
//! the corner. The corner-resolution gate then declines the fillet
//! with that reason. **Tolerance:** the caller's blend tolerance,
//! local and explicit.

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

#[path = "blend_mitre_corner_marched.rs"]
mod marched;
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

/// A seam's geometry, solved before the star is built.
enum SeamGeometry {
    /// Equal-radius cylinders: the exact bisector-plane ellipse.
    Exact(geometry::SeamFrame),
    /// Any other strip pair: marched samples, crossing first.
    Marched(marched::MarchedSeam),
}

/// One seam solved before the star is built.
struct StagedSeam {
    rail_crossing: VertexHandle,
    strips: (FaceHandle, FaceHandle),
    geometry: SeamGeometry,
    /// The seam's declared band: the crossing's recorded band (see
    /// [`crossing_band`]), widened by a marched seam to cover a
    /// sub-precision crossing residual on fitted strips.
    band: f64,
}

/// The tolerance band the seam at `crossing` inherits: the largest of
/// `pos_tol` and the tolerances recorded on the two strips and on the
/// edges at the crossing. Fitted strips (G2, variable-radius) carry
/// their fit accuracy as a recorded band, and their rails meet the
/// strips only within it (measured 6.8e-5 on a G2 strip at an oblique
/// corner). A seam starting at the crossing cannot be tighter than the
/// geometry it starts from, so the band is inherited explicitly and
/// recorded on the seam edge.
fn crossing_band(
    brep: &BRep,
    crossing: VertexHandle,
    strips: (FaceHandle, FaceHandle),
    pos_tol: f64,
) -> f64 {
    let faces = [strips.0, strips.1]
        .into_iter()
        .filter_map(|f| brep.faces.get(f)?.tolerance);
    let edges = brep
        .half_edges
        .values()
        .filter(|he| he.origin == crossing)
        .filter_map(|he| brep.edges.get(he.edge)?.tolerance);
    faces.chain(edges).fold(pos_tol, f64::max)
}

/// The strips around `corner`: every non-planar face bordering the
/// corner's face star.
fn corner_strips(
    brep: &BRep,
    corner: VertexHandle,
) -> Result<Vec<(FaceHandle, SurfaceGeometry)>, Reason> {
    let face_of = |he: HalfEdgeHandle| -> Option<FaceHandle> {
        Some(brep.loops.get(brep.half_edges.get(he)?.loop_)?.face)
    };
    let star: Vec<FaceHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, he)| he.origin == corner)
        .filter_map(|(h, _)| face_of(h))
        .collect();
    let mut strips: Vec<(FaceHandle, SurfaceGeometry)> = Vec::new();
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
            let surface = brep
                .faces
                .get(neighbour)
                .and_then(|f| f.surface.as_ref())
                .map(|s| s.geometry.clone())
                .ok_or(Reason::UnrecognizedResidue)?;
            let planar = matches!(surface, SurfaceGeometry::Plane { .. });
            if !star.contains(&neighbour) && !planar && strips.iter().all(|(f, _)| *f != neighbour)
            {
                strips.push((neighbour, surface));
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

/// The equal-radius cylinders of `strips`, when every strip is one.
fn equal_radius_cylinders(
    brep: &BRep,
    strips: &[(FaceHandle, SurfaceGeometry)],
    pos_tol: f64,
) -> Option<Vec<StripCylinder>> {
    let cylinders: Vec<StripCylinder> = strips
        .iter()
        .map(|(f, _)| strip_cylinder(brep, *f))
        .collect::<Option<_>>()?;
    let radius = cylinders.first()?.radius;
    cylinders
        .iter()
        .all(|c| (c.radius - radius).abs() <= pos_tol.max(radius * 1e-9))
        .then_some(cylinders)
}

fn mitre_corner(
    brep: &mut BRep,
    corner: VertexHandle,
    tolerance: &Tolerance,
) -> Result<(), Reason> {
    let strips = corner_strips(brep, corner)?;
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let corner_pos = brep
        .vertices
        .get(corner)
        .ok_or(Reason::UnrecognizedResidue)?
        .position;
    // The exact path applies when the strips are equal-radius
    // cylinders with concurrent axes; every other strip set is solved
    // on the surfaces.
    let cylinders = equal_radius_cylinders(brep, &strips, pos_tol);
    let exact_point = cylinders.as_ref().and_then(|cyl| {
        let axes: Vec<&StripCylinder> = cyl.iter().collect();
        point::solve_mitre_point(&axes, corner_pos, pos_tol)
    });
    let mitre_point = match exact_point {
        Some(p) => p,
        None => {
            let reach = strips
                .iter()
                .filter_map(|(f, _)| strip_cylinder(brep, *f).map(|c| c.radius))
                .fold(0.0_f64, f64::max)
                .max(corner_reach(brep, corner));
            let surfaces: Vec<&SurfaceGeometry> = strips.iter().map(|(_, g)| g).collect();
            marched::solve_general_mitre_point(&surfaces, corner_pos, 4.0 * reach, pos_tol)
                .ok_or(Reason::NoCommonPoint)?
        }
    };

    // Phase 1: reduce the residue to an N-gon of rail crossings.
    let strip_faces: Vec<FaceHandle> = strips.iter().map(|(f, _)| *f).collect();
    let polygon = surgery::reduce_to_residue_polygon(brep, corner, &strip_faces)?;

    // Solve every seam before building the star.
    let index_of = |face: FaceHandle| -> Result<usize, Reason> {
        strips
            .iter()
            .position(|(f, _)| *f == face)
            .ok_or(Reason::UnrecognizedResidue)
    };
    let n = polygon.crossings.len();
    let mut staged: Vec<StagedSeam> = Vec::with_capacity(n);
    for (k, &crossing) in polygon.crossings.iter().enumerate() {
        let pair = (polygon.sides[(k + n - 1) % n], polygon.sides[k]);
        let (a, b) = (index_of(pair.0)?, index_of(pair.1)?);
        let x = brep
            .vertices
            .get(crossing)
            .ok_or(Reason::UnrecognizedResidue)?
            .position;
        let band = crossing_band(brep, crossing, pair, pos_tol);
        let seam_geometry = match (&cylinders, exact_point) {
            (Some(cyl), Some(_)) => {
                let frame = geometry::solve_equal_radius_seam_frame(
                    &cyl[a],
                    &cyl[b],
                    x,
                    mitre_point,
                    tolerance,
                )
                .map_err(Reason::SeamUnresolved)?;
                geometry::author_seam_curve(&frame, x, mitre_point, corner_pos, tolerance)
                    .map_err(Reason::SeamUnresolved)?;
                SeamGeometry::Exact(frame)
            }
            _ => SeamGeometry::Marched(
                marched::march_mitre_seam(
                    &strips[a].1,
                    &strips[b].1,
                    x,
                    mitre_point,
                    pos_tol,
                    band,
                )
                .map_err(Reason::SeamUnresolved)?,
            ),
        };
        let band = match &seam_geometry {
            SeamGeometry::Marched(marched_seam) => marched_seam.band,
            SeamGeometry::Exact(_) => band,
        };
        staged.push(StagedSeam {
            rail_crossing: crossing,
            strips: pair,
            geometry: seam_geometry,
            band,
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

/// Longest edge at the corner: the scale the mitre point must stay
/// within when no strip radius bounds it.
fn corner_reach(brep: &BRep, corner: VertexHandle) -> f64 {
    let Some(c) = brep.vertices.get(corner).map(|v| v.position) else {
        return 0.0;
    };
    brep.half_edges
        .values()
        .filter(|he| he.origin == corner)
        .filter_map(|he| brep.half_edges.get(he.twin))
        .filter_map(|tw| brep.vertices.get(tw.origin))
        .map(|v| (v.position - c).length())
        .fold(0.0, f64::max)
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
    let curve: CurveOnEdge = match &seam.geometry {
        SeamGeometry::Exact(frame) => {
            let (p_start, p_end) = if primary == he_out { (p, x) } else { (x, p) };
            geometry::author_seam_curve(frame, p_start, p_end, stale_corner, tolerance)
        }
        // Samples run crossing → mitre; `he_out` runs mitre → crossing.
        SeamGeometry::Marched(marched_seam) => marched::seam_curve(marched_seam, primary == he_out),
    }
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
        // beside it and fixes the chart branch of the whole seam.
        let from_crossing = h.origin == seam.rail_crossing;
        let x_uv = if from_crossing {
            neighbour_chart_point(brep, face, h.prev, true)
        } else {
            neighbour_chart_point(brep, face, h.next, false)
        };
        let pcurve = match &seam.geometry {
            SeamGeometry::Exact(frame) => {
                let x_uv = x_uv.ok_or(Reason::SeamUnresolved(
                    crate::types::BlendSeamCornerDeclineReason::PcurveClosureFailed,
                ))?;
                let p_uv = point::chart_point_near(&surface, p, x_uv);
                let (origin_pos, dest_pos, origin_uv, dest_uv) = if from_crossing {
                    (x, p, Some(x_uv), p_uv)
                } else {
                    (p, x, p_uv, Some(x_uv))
                };
                geometry::author_side_pcurve(
                    &surface,
                    frame,
                    origin_pos,
                    dest_pos,
                    &geometry::ChartAnchor { uv: origin_uv },
                    &geometry::ChartAnchor { uv: dest_uv },
                    tolerance,
                )
            }
            SeamGeometry::Marched(marched_seam) => {
                marched::seam_pcurve(&surface, marched_seam, !from_crossing, x_uv)
            }
        }
        .map_err(Reason::SeamUnresolved)?;
        pcurves.push((face, he, pcurve));
    }
    commit_edge_authoring(brep, edge, curve, pcurves);
    // The band is the seam edge's own declaration; vertex tolerances
    // stay put (the kernel's containment rule keeps a vertex's
    // tolerance at or under every incident edge's).
    if seam.band > (tolerance.linear() * 10.0).max(1e-9)
        && let Some(e) = brep.edges.get_mut(edge)
    {
        e.tolerance = Some(e.tolerance.unwrap_or(0.0).max(seam.band));
    }
    if brep.measure_edge_deviation(edge) > seam.band {
        return Err(Reason::ValidationFailed);
    }
    Ok(())
}
