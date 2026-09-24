//! A strip's end running into a FAN of faces (PK-530, PK-534).
//!
//! Two fixtures reach this:
//! - A single fillet on a pyramid's slant edge ends at the apex, where
//!   it meets not one planar cap but every face around the apex that
//!   is not one of its supports.
//! - A fillet on a cube's vertical edge, after its top edges are
//!   rounded (blend on blend), ends against the two top-fillet strips,
//!   which already meet along their PK-152 seam.
//!
//! The chain-end settlement needs a single cap, so it declined
//! (`NoCommonPlanarCap`), and the pipeline served the stale L through
//! the corner: 1.3e-2 off at the pyramid, (√2 − 1)·r = 4.1e-2 at the
//! cube.
//!
//! **Exact end.** Name the faces at the corner, from the in-stub's cap
//! to the out-stub's, `F_0 … F_m`, and the fan edges between them
//! `f_1 … f_m`. The strip's end is the chain of strip ∩ `F_j` arcs:
//! `foot_a → q_1` on `F_0`, then `q_j → q_(j+1)` on `F_j`, then
//! `q_m → foot_b` on `F_m`. Here `q_j` is where the strip crosses
//! `f_j`, nearest the corner. The arcs are built one of two ways:
//! - exact ellipses for a cylinder strip against planes;
//! - marched otherwise, for fitted strips, curved fan faces, and
//!   curved fan edges (`surface_pair_march`).
//!
//! **Topology.** For m = 1 the corner vertex just moves to q_1. For
//! m ≥ 2 the corner splits (`split_vertex`, the general MEV): each
//! split moves the fan sector past `f_j` to a new vertex at
//! `q_(j+1)`, and its connector edge lands between the strip and
//! `F_j`, which is exactly the new end arc. Each fan edge keeps its
//! curve, trimmed to start at its `q_j` (`edges`).
//!
//! **Fallback.** The whole end is transactional. It is solved and
//! gated, and any decline restores the body. **Tolerance:** the
//! caller's blend tolerance. A fitted fan edge keeps its own recorded
//! band.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::types::TerminalArcTrimDeclineReason as Reason;

use super::terminal_arc_oblique::{ObliqueSection, compute_oblique_stub_arc, unwrap_near};
use super::terminal_arc_trim::{CornerCandidate, StripCylinder};

#[path = "terminal_fan_end_edges.rs"]
mod edges;
#[path = "terminal_fan_end_fitted.rs"]
mod fitted;

/// Largest fan the settlement will split (m fan edges).
const MAX_FAN_EDGES: usize = 8;

type Plane = (Point3, Vec3, Vec3, Vec3);

fn face_of(brep: &BRep, he: HalfEdgeHandle) -> Option<FaceHandle> {
    Some(brep.loops.get(brep.half_edges.get(he)?.loop_)?.face)
}

fn plane_of(brep: &BRep, face: FaceHandle) -> Option<Plane> {
    match brep.faces.get(face)?.surface.as_ref().map(|s| &s.geometry) {
        Some(SurfaceGeometry::Plane {
            origin,
            u_axis,
            v_axis,
            normal,
        }) => Some((*origin, *u_axis, *v_axis, *normal)),
        _ => None,
    }
}

/// The strip whose end runs into the fan.
pub(super) enum FanStrip<'a> {
    /// An analytic cylinder: exact ellipse arcs.
    Cylinder(&'a StripCylinder),
    /// A fitted strip (G2 or variable-radius NURBS): arcs marched onto
    /// the strip and each fan plane (`terminal_fan_end_fitted`).
    Fitted(&'a SurfaceGeometry),
}

/// Settle a fan end transactionally.
pub(super) fn repair_fan_end(
    brep: &mut BRep,
    c: &CornerCandidate,
    strip: FanStrip<'_>,
    tolerance: &Tolerance,
) -> Result<(), Reason> {
    let snapshot = brep.clone();
    let outcome = settle(brep, c, &strip, tolerance);
    if outcome.is_err() {
        *brep = snapshot;
    }
    outcome
}

fn settle(
    brep: &mut BRep,
    c: &CornerCandidate,
    strip: &FanStrip<'_>,
    tolerance: &Tolerance,
) -> Result<(), Reason> {
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    // Fan from the out-stub: [he_out, f_m, …, f_1, twin(he_in)].
    let twin_in = brep
        .half_edges
        .get(c.he_in)
        .ok_or(Reason::NoCommonPlanarCap)?
        .twin;
    let mut fan: Vec<HalfEdgeHandle> = Vec::new();
    let mut cursor = brep
        .half_edges
        .get(
            brep.half_edges
                .get(c.he_out)
                .ok_or(Reason::NoCommonPlanarCap)?
                .twin,
        )
        .ok_or(Reason::NoCommonPlanarCap)?
        .next;
    while cursor != twin_in {
        if fan.len() > MAX_FAN_EDGES {
            return Err(Reason::NoCommonPlanarCap);
        }
        fan.push(cursor);
        let twin = brep
            .half_edges
            .get(cursor)
            .ok_or(Reason::NoCommonPlanarCap)?
            .twin;
        cursor = brep
            .half_edges
            .get(twin)
            .ok_or(Reason::NoCommonPlanarCap)?
            .next;
    }
    if fan.is_empty() {
        // One oblique cap (m = 0): only a fitted strip reaches here;
        // a cylinder's single oblique cap is settled exactly by
        // `terminal_arc_oblique`.
        return match strip {
            FanStrip::Fitted(geom) => fitted::settle_single_cap(brep, c, geom, tolerance),
            FanStrip::Cylinder(_) => Err(Reason::NoCommonPlanarCap),
        };
    }
    // fan = [f_m, …, f_1]; reorder to f_1 … f_m.
    fan.reverse();
    let m = fan.len();
    // Faces F_0 … F_m: F_0 is the in-stub's cap (twin(he_in)'s loop),
    // and F_j is the loop of f_j's own outgoing half-edge, the wedge
    // between f_j and f_(j+1) on the out-stub side.
    let cap_first = face_of(brep, twin_in).ok_or(Reason::NoCommonPlanarCap)?;
    let mut faces = vec![cap_first];
    for &f in &fan {
        faces.push(face_of(brep, f).ok_or(Reason::NoCommonPlanarCap)?);
    }
    let cap_last = face_of(
        brep,
        brep.half_edges
            .get(c.he_out)
            .ok_or(Reason::NoCommonPlanarCap)?
            .twin,
    )
    .ok_or(Reason::NoCommonPlanarCap)?;
    if faces[m] != cap_last {
        return Err(Reason::NoCommonPlanarCap);
    }
    let surfaces: Vec<SurfaceGeometry> = faces
        .iter()
        .map(|f| {
            brep.faces
                .get(*f)
                .and_then(|x| x.surface.as_ref())
                .map(|x| x.geometry.clone())
                .ok_or(Reason::NoCommonPlanarCap)
        })
        .collect::<Result<_, _>>()?;
    // The exact path: a cylinder strip against planar faces only.
    let planes: Option<Vec<Plane>> = faces.iter().map(|f| plane_of(brep, *f)).collect();
    let exact = match (strip, &planes) {
        (FanStrip::Cylinder(cyl), Some(planes)) => Some((*cyl, planes)),
        _ => None,
    };
    let strip_geometry = match strip {
        FanStrip::Cylinder(cyl) => SurfaceGeometry::Cylinder {
            origin: cyl.origin,
            axis: cyl.axis,
            radius: cyl.radius,
        },
        FanStrip::Fitted(geom) => (*geom).clone(),
    };

    // Crossings q_j of the strip with each fan edge, from the corner.
    let corner_pos = brep
        .vertices
        .get(c.corner)
        .ok_or(Reason::NoCommonPlanarCap)?
        .position;
    let mut crossings: Vec<Point3> = Vec::with_capacity(m);
    for &f in &fan {
        let q = edges::fan_edge_crossing(
            brep,
            f,
            corner_pos,
            exact.map(|e| e.0),
            &strip_geometry,
            pos_tol,
        )
        .ok_or(Reason::ArcBranchUnresolved)?;
        crossings.push(q);
    }

    // Stub anchors (the feet).
    let foot_a = brep
        .half_edges
        .get(c.he_in)
        .and_then(|h| brep.vertices.get(h.origin))
        .ok_or(Reason::ArcBranchUnresolved)?
        .position;
    let foot_b = brep
        .half_edges
        .get(c.he_out)
        .and_then(|h| brep.half_edges.get(h.twin))
        .and_then(|tw| brep.vertices.get(tw.origin))
        .ok_or(Reason::ArcBranchUnresolved)?
        .position;
    // Chain of end points foot_a, q_1, …, q_m, foot_b; arc j lies on F_j.
    let mut chain = vec![foot_a];
    chain.extend(&crossings);
    chain.push(foot_b);
    // Solve every arc before touching topology.
    let arcs = match exact {
        Some((cyl, planes)) => FanArcs::Exact(exact_arcs(cyl, planes, &chain, pos_tol)?),
        None => {
            let recorded = fitted::recorded_band(brep, c, pos_tol);
            FanArcs::Marched(fitted::march_arcs(
                &strip_geometry,
                &surfaces,
                &chain,
                pos_tol,
                recorded,
            )?)
        }
    };

    // Topology: split the corner so every q_j has its own vertex.
    // After splitting, `vertices[j]` sits at q_(j+1) (0-based).
    let mut vertices: Vec<VertexHandle> = vec![c.corner];
    for j in 1..m {
        // Move [he_out, f_m, …, f_(j+1)] to a new vertex at q_(j+1).
        let result = brep
            .split_vertex(c.he_out, fan[j], crossings[j])
            .map_err(|_| Reason::FanEndRejected)?;
        vertices.push(result.vertex);
    }
    if let Some(v) = brep.vertices.get_mut(c.corner) {
        v.position = crossings[0];
    }
    for (j, &v) in vertices.iter().enumerate() {
        if let Some(vx) = brep.vertices.get_mut(v) {
            vx.position = crossings[j];
        }
    }

    // Strip-side half-edges of the chain, in order: he_in, the
    // connectors, he_out (each from chain[j] to chain[j + 1]).
    let mut strip_half_edges = vec![c.he_in];
    let mut cursor = brep
        .half_edges
        .get(c.he_in)
        .ok_or(Reason::FanEndRejected)?
        .next;
    for _ in 0..m {
        strip_half_edges.push(cursor);
        cursor = brep
            .half_edges
            .get(cursor)
            .ok_or(Reason::FanEndRejected)?
            .next;
    }
    if strip_half_edges[m] != c.he_out {
        return Err(Reason::FanEndRejected);
    }
    let mut staged = Vec::with_capacity(m + 1);
    for j in 0..=m {
        let he = strip_half_edges[j];
        if face_of(
            brep,
            brep.half_edges.get(he).ok_or(Reason::FanEndRejected)?.twin,
        ) != Some(faces[j])
        {
            return Err(Reason::FanEndRejected);
        }
        let (stub, band) = match (&arcs, &planes) {
            (FanArcs::Exact((sections, angles)), Some(planes)) => {
                let (o, u, v, _) = planes[j];
                (
                    compute_oblique_stub_arc(
                        brep,
                        he,
                        &sections[j],
                        angles[j],
                        angles[j + 1],
                        (o, u, v),
                        tolerance,
                    )?,
                    None,
                )
            }
            (FanArcs::Marched(arcs), _) => {
                let staged = fitted::stage_marched_arc(
                    brep,
                    he,
                    &strip_geometry,
                    &surfaces[j],
                    &arcs[j],
                    j == 0,
                )?;
                (staged, Some(arcs[j].band))
            }
            _ => return Err(Reason::FanEndRejected),
        };
        staged.push((stub, band));
    }
    for (stub, band) in staged {
        commit_stub(brep, c.strip_face, stub, band, pos_tol)?;
    }
    // Fan edges: the same curves, shortened to their crossings.
    for (j, &f) in fan.iter().enumerate() {
        edges::reauthor_fan_edge(brep, f, crossings[j], (faces[j], faces[j + 1]), tolerance)?;
    }
    brep.validate().map_err(|_| Reason::FanEndRejected)?;
    Ok(())
}

/// Commit one staged end arc: its curve (declaring a fitted strip's
/// marched band), and its strip and cap pcurves.
fn commit_stub(
    brep: &mut BRep,
    strip_face: FaceHandle,
    stub: super::terminal_arc_author::StagedStubArc,
    band: Option<f64>,
    pos_tol: f64,
) -> Result<(), Reason> {
    if let Some(e) = brep.edges.get_mut(stub.edge) {
        e.curve = Some(stub.curve);
        if let Some(band) = band.filter(|b| *b > pos_tol) {
            e.tolerance = Some(e.tolerance.unwrap_or(0.0).max(band));
        }
    }
    if let Some(he) = brep.half_edges.get_mut(stub.strip_half_edge) {
        he.pcurve = Some(stub.strip_pcurve.clone());
    }
    brep.set_face_use_pcurve(strip_face, stub.strip_half_edge, stub.strip_pcurve);
    if let Some(tw) = brep.half_edges.get_mut(stub.cap_half_edge) {
        tw.pcurve = Some(stub.cap_pcurve.clone());
    }
    let cap = face_of(brep, stub.cap_half_edge).ok_or(Reason::FanEndRejected)?;
    brep.set_face_use_pcurve(cap, stub.cap_half_edge, stub.cap_pcurve);
    Ok(())
}

/// Every arc of the fan end, solved before the topology changes.
enum FanArcs {
    /// Cylinder strip: each face's exact section and the chart angles
    /// of the chain.
    Exact((Vec<ObliqueSection>, Vec<f64>)),
    /// Fitted strip: each face's marched intersection.
    Marched(Vec<super::surface_pair_march::MarchedIntersection>),
}

/// The exact sections and chart angles of a cylinder strip's fan end.
fn exact_arcs(
    cyl: &StripCylinder,
    planes: &[Plane],
    chain: &[Point3],
    pos_tol: f64,
) -> Result<(Vec<ObliqueSection>, Vec<f64>), Reason> {
    let m = planes.len() - 1;
    let foot_a = chain[0];
    let sections: Vec<ObliqueSection> = planes
        .iter()
        .map(|p| ObliqueSection::new(cyl, p.0, p.3).ok_or(Reason::AnchorsOffSurface))
        .collect::<Result<_, _>>()?;
    for (j, section) in sections.iter().enumerate() {
        if section.miss(chain[j]) > pos_tol || section.miss(chain[j + 1]) > pos_tol {
            return Err(Reason::AnchorsOffSurface);
        }
    }
    // Chart angles along the chain, on the strip chart's branch at
    // foot_a, monotone and within a half-turn in total.
    let strip = SurfaceGeometry::Cylinder {
        origin: cyl.origin,
        axis: cyl.axis,
        radius: cyl.radius,
    };
    let (chart_u, _) = popkern_core::geometry::uv_projection::project_point_to_uv(&strip, &foot_a)
        .ok_or(Reason::ArcBranchUnresolved)?;
    let mut angles = vec![unwrap_near(sections[0].angle_of(foot_a), chart_u)];
    for (j, p) in chain.iter().enumerate().skip(1) {
        let prev = angles[j - 1];
        angles.push(unwrap_near(sections[j.min(m)].angle_of(*p), prev));
    }
    let total = angles[m + 1] - angles[0];
    if total.abs() >= std::f64::consts::PI
        || angles.windows(2).any(|w| (w[1] - w[0]) * total <= 0.0)
    {
        return Err(Reason::ArcBranchUnresolved);
    }

    Ok((sections, angles))
}
