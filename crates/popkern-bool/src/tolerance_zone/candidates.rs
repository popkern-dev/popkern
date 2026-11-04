//! Candidate vertices: every point an interference produces, with
//! where it lies as well as where it is (plan §4.3a).
//!
//! Three sources feed the first stage: the operands' own vertices, the
//! edge-face hits, and the ends of the face-face records. Each
//! candidate carries its incidence signature (the vertices, edges,
//! faces and records it lies on, with local parameters) and a
//! tolerance radius, the sum of its sources' tolerances (Jackson
//! §3.2): an operand vertex its own tolerance; an edge-face hit the
//! edge's plus the face's; a record end the two faces' plus the
//! record's own fit deviation there, its 3D end's distance off either
//! surface (a record fit off its surfaces is that much less certain
//! where it ends).

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};
use slotmap::Key;

use super::region::Bounds;
use crate::intersect::{IntersectionRecord, IntersectionResult};
use crate::intersect_edge_face::EdgeFaceHit;

/// One entity a candidate lies on.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Incidence {
    /// An operand vertex (operand 0 = A, 1 = B).
    Vertex {
        /// The operand.
        operand: usize,
        /// The vertex.
        vertex: VertexHandle,
    },
    /// An operand edge at parameter `t`.
    Edge {
        /// The operand.
        operand: usize,
        /// The edge.
        edge: EdgeHandle,
        /// The edge curve's parameter.
        t: f64,
    },
    /// An operand face.
    Face {
        /// The operand.
        operand: usize,
        /// The face.
        face: FaceHandle,
        /// Where on the face, in its chart and period branch; `None`
        /// where the producer carries no chart for it (a trim residue's
        /// other face).
        uv: Option<(f64, f64)>,
    },
    /// A face-face record at curve parameter `s`.
    Record {
        /// The record's index in the intersection result.
        record: usize,
        /// The record curve's parameter.
        s: f64,
        /// Which end of the record the point is: `Some(false)` its
        /// start, `Some(true)` its end, `None` its interior (a record
        /// another candidate lies on, a T-junction to be).
        end: Option<bool>,
    },
}

/// Which interference produced a candidate, in dimension order: the
/// order is the elimination precedence (a group anchors at its lowest).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) enum SourceKind {
    /// An operand's own vertex.
    OperandVertex,
    /// An operand edge piercing the other operand's face.
    EdgeFaceHit,
    /// A face-face record's end.
    RecordEndpoint,
}

/// A point some interference produced, before any vertex is committed.
#[derive(Debug, Clone)]
pub(crate) struct CandidateVertex {
    /// Where it is.
    pub position: Point3,
    /// Its tolerance radius.
    pub radius: f64,
    /// Its reach along an edge it lies on: the radius over the sine of
    /// the angle at which that edge crosses the other operand's face
    /// (`incidence::add_crossing_zones`), so at least the radius. A
    /// shallow crossing is that uncertain along the edge: two
    /// computations of one crossing (a hit and a record end) land that
    /// far apart on it.
    pub zone: f64,
    /// What produced it.
    pub source: SourceKind,
    /// What it lies on.
    pub incidence: Vec<Incidence>,
    /// Canonical order key: source kind, operand, then the producing
    /// entities' handles and parameters; independent of the order the
    /// interferences were computed in.
    pub key: (SourceKind, usize, u64, u64, u64),
}

/// Where `record` lies at `s` on its face of operand B (`on_b`) or A,
/// by that face's own uv curve, when that uv reproduces the record's
/// point there within `radius`: a trim residue lies on its owner face
/// alone and carries no chart for the other, and a uv curve's end can
/// miss its 3D end (measured: 3e-3 to 2e-2 along cone and sphere
/// walls, the edge-face hit at the same point exact), which is no chart
/// evidence of where on the face the point lies.
fn chart_uv(
    (a, b): (&BRep, &BRep),
    record: &IntersectionRecord,
    on_b: bool,
    s: f64,
    (position, radius): (Point3, f64),
) -> Option<(f64, f64)> {
    let (curve, brep, face) = match (record.residue_owner_a, on_b) {
        (Some(owner_a), _) if owner_a == on_b => return None,
        (_, false) => (&record.curve.curve_uv1, a, record.face_a),
        (_, true) => (&record.curve.curve_uv2, b, record.face_b),
    };
    let p = curve.point_at(s);
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    (surface.point_at(p.x, p.y).distance_to(&position) <= radius).then_some((p.x, p.y))
}

fn handle_bits(k: impl Key) -> u64 {
    k.data().as_ffi()
}

/// Every candidate of a boolean between `a` and `b`, sorted by key:
/// the intersection result's edge-face hits and the ends of `records`,
/// the records the graph build inserts (residues stripped, closed
/// records rotated, multi-turn records cut at their band lattice), and
/// the operand vertices within reach of them.
///
/// An operand vertex interferes only where the operands meet, and
/// every meeting leaves a hit or a record there (a vertex touching the
/// other operand is its edges' hits on the face it touches), so a
/// vertex farther than `cap` plus the radii from every hit, record end
/// and record curve stays a vertex of its own operand alone: it joins
/// no group, and the table leaves it out. Without the filter every
/// boolean pays for its operands' whole vertex sets (measured: a 100-
/// knuckle piano hinge built 3.1 s against 1.8 s, growing as the
/// square of the knuckles as each union re-read the growing body).
pub(crate) fn collect_candidates(
    a: &BRep,
    b: &BRep,
    hits: &[EdgeFaceHit],
    records: &[IntersectionRecord],
    linear: f64,
    cap: f64,
) -> Vec<CandidateVertex> {
    let operands = [a, b];
    let mut out = Vec::new();
    collect_sites(a, b, hits, records, linear, &mut out);
    let site_radius = out.iter().map(|c| c.radius).fold(0.0_f64, f64::max);
    let vertex_tol = operands
        .iter()
        .flat_map(|o| o.vertices.values())
        .map(|v| v.tolerance.unwrap_or(linear))
        .fold(linear, f64::max);
    let reach = cap.max(site_radius) + vertex_tol + linear;
    let mut sites: Vec<Bounds> = out
        .iter()
        .filter_map(|c| Bounds::of(&[c.position], reach))
        .collect();
    sites.extend(
        records
            .iter()
            .filter_map(|r| Bounds::of(&r.curve.curve_3d.control_points, reach)),
    );
    for (operand, brep) in operands.iter().enumerate() {
        for (vh, v) in &brep.vertices {
            if !sites.iter().any(|s| s.holds(&v.position)) {
                continue;
            }
            out.push(CandidateVertex {
                position: v.position,
                radius: v.tolerance.unwrap_or(linear),
                zone: v.tolerance.unwrap_or(linear),
                source: SourceKind::OperandVertex,
                incidence: vec![Incidence::Vertex {
                    operand,
                    vertex: vh,
                }],
                key: (SourceKind::OperandVertex, operand, handle_bits(vh), 0, 0),
            });
        }
    }
    out.sort_by_key(|c| c.key);
    out
}

/// The edge-face hits of `ix` that are interferences: those on their
/// face's trim (or at the end of a record on that face). The
/// intersector reports an edge crossing the face's SURFACE, and only
/// a crossing on the face is one. The graph splits operand edges only
/// at these when the table decides vertices: a hit off its face split
/// an edge the table's paves and common blocks never knew about
/// (measured: a flush gusset's rim edge split where it crossed a
/// neighbouring face's plane, beside the record running along it,
/// and the duplicate pairing that needs both split alike failed).
pub(crate) fn interfering_hits(
    (a, b): (&BRep, &BRep),
    ix: &IntersectionResult,
    records: &[IntersectionRecord],
    linear: f64,
) -> Vec<EdgeFaceHit> {
    let operands = [a, b];
    ix.edge_face_hits
        .iter()
        .filter(|hit| {
            let face_op = usize::from(hit.edge_of_a);
            at_record_end(records, face_op, hit, hit_radius(operands, hit, linear))
                || !outside_face(operands[face_op], hit.face, &hit.point)
        })
        .cloned()
        .collect()
}

/// A hit's tolerance radius: its edge's plus its face's.
fn hit_radius(operands: [&BRep; 2], hit: &EdgeFaceHit, linear: f64) -> f64 {
    let (edge_op, face_op) = if hit.edge_of_a { (0, 1) } else { (1, 0) };
    let edge_tol = operands[edge_op]
        .edges
        .get(hit.edge)
        .and_then(|e| e.tolerance)
        .unwrap_or(linear);
    let face_tol = operands[face_op]
        .faces
        .get(hit.face)
        .and_then(|f| f.tolerance)
        .unwrap_or(linear);
    edge_tol + face_tol
}

/// The interference sites of a boolean: its interfering edge-face hits
/// (`interfering_hits`) and the ends of `records`, pushed onto `out`.
fn collect_sites(
    a: &BRep,
    b: &BRep,
    hits: &[EdgeFaceHit],
    records: &[IntersectionRecord],
    linear: f64,
    out: &mut Vec<CandidateVertex>,
) {
    let operands = [a, b];
    for hit in hits {
        let (edge_op, face_op) = if hit.edge_of_a { (0, 1) } else { (1, 0) };
        let radius = hit_radius(operands, hit, linear);
        out.push(CandidateVertex {
            position: hit.point,
            radius,
            zone: radius,
            source: SourceKind::EdgeFaceHit,
            incidence: vec![
                Incidence::Edge {
                    operand: edge_op,
                    edge: hit.edge,
                    t: edge_parameter(operands[edge_op], hit.edge, &hit.point)
                        .unwrap_or(hit.t_edge),
                },
                Incidence::Face {
                    operand: face_op,
                    face: hit.face,
                    uv: operands[face_op]
                        .faces
                        .get(hit.face)
                        .and_then(|f| f.surface.as_ref())
                        .filter(|sf| {
                            sf.geometry
                                .point_at(hit.uv.0, hit.uv.1)
                                .distance_to(&hit.point)
                                <= radius
                        })
                        .map(|_| hit.uv),
                },
            ],
            key: (
                SourceKind::EdgeFaceHit,
                edge_op,
                handle_bits(hit.edge),
                handle_bits(hit.face),
                hit.t_edge.to_bits(),
            ),
        });
    }
    for (ri, record) in records.iter().enumerate() {
        let (s0, s1) = record.curve.curve_3d.parameter_range();
        for (end, s) in [(false, s0), (true, s1)] {
            let position = record.curve.curve_3d.point_at(s);
            let face_tol = |brep: &BRep, f: FaceHandle| {
                brep.faces
                    .get(f)
                    .and_then(|f| f.tolerance)
                    .unwrap_or(linear)
            };
            let radius = face_tol(a, record.face_a)
                + face_tol(b, record.face_b)
                + record_fit_deviation(a, b, record, end, &position);
            out.push(CandidateVertex {
                position,
                radius,
                zone: radius,
                source: SourceKind::RecordEndpoint,
                incidence: vec![
                    Incidence::Record {
                        record: ri,
                        s,
                        end: Some(end),
                    },
                    Incidence::Face {
                        operand: 0,
                        face: record.face_a,
                        uv: chart_uv((a, b), record, false, s, (position, radius)),
                    },
                    Incidence::Face {
                        operand: 1,
                        face: record.face_b,
                        uv: chart_uv((a, b), record, true, s, (position, radius)),
                    },
                ],
                key: (
                    SourceKind::RecordEndpoint,
                    0,
                    handle_bits(record.face_a),
                    handle_bits(record.face_b),
                    position.x.to_bits() ^ position.y.to_bits().rotate_left(21) ^ u64::from(end),
                ),
            });
        }
    }
}

/// Whether `hit` lies within `reach` of an end of a record on its
/// face: a record is clipped to the trims of the faces it lies on (a
/// trim residue to its owner's only), so its ends are on the face and
/// a hit there needs no trim query. The query rebuilds the face's loop
/// polygon per point (measured: 2.2 ms a boolean for eight hits on a
/// disc face of 250 edges).
fn at_record_end(
    records: &[IntersectionRecord],
    face_op: usize,
    hit: &EdgeFaceHit,
    reach: f64,
) -> bool {
    records.iter().any(|r| {
        let (face, on) = if face_op == 0 {
            (r.face_a, r.residue_owner_a != Some(false))
        } else {
            (r.face_b, r.residue_owner_a != Some(true))
        };
        let curve = &r.curve.curve_3d;
        let (s0, s1) = curve.parameter_range();
        on && face == hit.face
            && [s0, s1]
                .iter()
                .any(|&s| curve.point_at(s).distance_to(&hit.point) <= reach)
    })
}

/// Whether `p` lies outside `face`'s trim by the face's exact trim
/// query. The intersector reports an edge crossing the face's SURFACE
/// (measured: a ring's rim circle crossing a web end's plane beyond the
/// web), and only a crossing on the face is an interference; a query
/// that cannot answer keeps the hit.
fn outside_face(brep: &BRep, face: FaceHandle, p: &Point3) -> bool {
    let Some(surface) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
        return false;
    };
    let Some((u, v)) = surface.geometry.closest_point(p) else {
        return false;
    };
    brep.face_uv_trim_containment(face, &surface.geometry, u, v)
        == Some(popkern_core::validation::FaceUvContainment::Outside)
}

/// The parameter of `p` on `edge`'s own curve. A hit's `t_edge` is the
/// intersector's parameter, which can run the other way along the edge
/// (measured: 0.8116 against the curve's 0.1884 at one point), so
/// every edge incidence is read on the edge curve itself.
pub(crate) fn edge_parameter(brep: &BRep, edge: EdgeHandle, p: &Point3) -> Option<f64> {
    let c = brep.edges.get(edge)?.curve.as_ref()?;
    Some(popkern_core::geometry::curve::closest_point_in_window(
        &c.geometry,
        p,
        c.t_start,
        c.t_end,
    ))
}

/// How far a record's 3D end sits off its two faces' surfaces: the
/// record's own fit deviation there. (The 3D curve carries the clipped
/// parameter window while the pcurves keep their own domains, so the
/// surfaces, not the pcurve ends, are the reference.)
fn record_fit_deviation(
    a: &BRep,
    b: &BRep,
    record: &IntersectionRecord,
    _end: bool,
    position: &Point3,
) -> f64 {
    let mut worst = 0.0_f64;
    for (brep, face) in [(a, record.face_a), (b, record.face_b)] {
        let Some(surface) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let Some((u, v)) = surface.geometry.closest_point(position) else {
            continue;
        };
        let d = surface.geometry.point_at(u, v).distance_to(position);
        if d.is_finite() {
            worst = worst.max(d);
        }
    }
    worst
}
