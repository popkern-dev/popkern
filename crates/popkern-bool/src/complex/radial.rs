//! Radial order of the fins around every complex edge (PK-973).
//!
//! The ordering itself is the core's
//! ([`popkern_core::topology::radial_order`]): each fin's face leaves
//! the edge in `d = folded_normal × fin_tangent`, the directions sort
//! by angle about the edge tangent, ties go down the tie ladder, and
//! the order must agree at every station. This module reads a complex
//! edge as the engine's source (the face's primary owner's surface,
//! the station's parameters seeded from the intersection graph's
//! pcurves) and records each fin's side-facing: whether its face's
//! Ahead side faces the NEXT wedge (counter-clockwise about `t`),
//! which holds iff the face's outward normal has a positive component
//! along `t × d`.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceOnFace};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::BRep;
use popkern_core::topology::radial_order::{
    RadialFins, RadialFrame, RadialOrderDecline, RadialTolerance, TieLadder, order_fins,
};

use super::ids::{EdgeId, FinId};
use super::model::CellComplex;

/// What the ladder measured for one tied fin.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct TieEvidence {
    /// The fin.
    pub fin: FinId,
    /// Its first-order angle about the tangent from the station's
    /// reference fin, in [0, 2π) (or just below zero for a tie
    /// wrapping past 2π).
    pub angle: f64,
    /// Its signed cross-section curvature, where the rung ran.
    pub curvature: Option<f64>,
    /// The signed angle of its offset probe point from `angle`, where
    /// the rung ran.
    pub probe: Option<f64>,
}

/// Why radial ordering declined.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum RadialDecline {
    /// Two or more fins leave the edge in the same direction and no
    /// rung of the tie ladder separates them.
    Ambiguous {
        /// The edge.
        edge: EdgeId,
        /// Fins that stayed tied.
        fins: Vec<FinId>,
        /// The station (fraction of the edge from its start).
        station: f64,
        /// What the ladder measured for every fin of the tie.
        evidence: Vec<TieEvidence>,
    },
    /// The cyclic order differs between stations along the edge: faces
    /// cross inside it, an intersection curve the intersection stage
    /// missed.
    Inconsistent {
        /// The edge.
        edge: EdgeId,
        /// The first station's order and the first that disagrees.
        stations: Vec<(f64, Vec<FinId>)>,
    },
    /// A face's normal or the edge tangent could not be evaluated.
    Unevaluable {
        /// The edge.
        edge: EdgeId,
        /// The fin whose face failed, if the edge itself did not.
        fin: Option<FinId>,
    },
}

impl RadialDecline {
    /// The engine's decline about `edge`, its fin indices into `fins`.
    fn from_engine(edge: EdgeId, fins: &[FinId], d: RadialOrderDecline) -> Self {
        let ids = |v: Vec<usize>| v.into_iter().map(|i| fins[i]).collect();
        match d {
            RadialOrderDecline::Ambiguous {
                fins: tied,
                station,
                evidence,
            } => Self::Ambiguous {
                edge,
                fins: ids(tied),
                station,
                evidence: evidence
                    .into_iter()
                    .map(|e| TieEvidence {
                        fin: fins[e.fin],
                        angle: e.angle,
                        curvature: e.curvature,
                        probe: e.probe,
                    })
                    .collect(),
            },
            RadialOrderDecline::Inconsistent { stations } => Self::Inconsistent {
                edge,
                stations: stations.into_iter().map(|(s, o)| (s, ids(o))).collect(),
            },
            RadialOrderDecline::Unevaluable { fin } => Self::Unevaluable {
                edge,
                fin: fin.map(|i| fins[i]),
            },
        }
    }
}

/// Per-fin result of the radial stage.
#[derive(Debug, Clone, Default)]
pub(crate) struct RadialSides {
    /// For each fin: whether its face's Ahead side faces the next
    /// (counter-clockwise) wedge about the edge tangent.
    pub ahead_faces_next: Vec<bool>,
}

/// Geometry the radial stage needs for a fin's face.
pub(crate) struct FaceGeometry<'a> {
    /// The operand BReps.
    pub operands: &'a [&'a BRep],
    /// The intersection graph the complex lifts, when there is one:
    /// its half-edges' pcurves seed each fin's projection onto its face.
    pub graph: Option<&'a crate::merge::graph::IntersectionGraph>,
}

/// One station: the point, the unit tangent along start → end.
pub(super) fn station(c: &CellComplex, edge: EdgeId, s: f64) -> Option<(Point3, Vec3)> {
    let e = c.edges.get(edge.index())?;
    let ps = c.vertices.get(e.start.index())?.position;
    let pe = c.vertices.get(e.end.index())?.position;
    match &e.curve {
        Some(curve) => {
            let (t0, t1) = (curve.t_start, curve.t_end);
            let a = curve.geometry.point_at(t0);
            // Which parameter end is the edge's start (closed edges:
            // both; keep the curve's own direction).
            let reversed = e.start != e.end && a.distance_to(&pe) < a.distance_to(&ps);
            let s = if reversed { 1.0 - s } else { s };
            let t = t0 + (t1 - t0) * s;
            let p = curve.geometry.point_at(t);
            let mut d = curve.geometry.derivative_at(t, 1);
            if (t1 - t0) < 0.0 {
                d = d * -1.0;
            }
            if reversed {
                d = d * -1.0;
            }
            Some((p, d.normalized()?))
        }
        None => {
            let d = (pe - ps).normalized()?;
            Some((ps + (pe - ps) * s, d))
        }
    }
}

/// The surface of `fin`'s face (its primary owner's).
fn fin_surface<'a>(
    c: &CellComplex,
    geo: &FaceGeometry<'a>,
    fin: FinId,
) -> Option<&'a SurfaceOnFace> {
    let f = c.fins.get(fin.index())?;
    let owner = c.faces.get(f.face.index())?.owners.first()?;
    geo.operands
        .get(owner.operand)?
        .faces
        .get(owner.face)?
        .surface
        .as_ref()
}

/// One complex edge as the engine's source.
struct ComplexEdgeFins<'c, 'g> {
    c: &'c CellComplex,
    geo: &'c FaceGeometry<'g>,
    edge: EdgeId,
    fins: &'c [FinId],
}

impl RadialFins for ComplexEdgeFins<'_, '_> {
    fn fin_count(&self) -> usize {
        self.fins.len()
    }

    fn station(&self, s: f64) -> Option<(Point3, Vec3)> {
        station(self.c, self.edge, s)
    }

    fn frame(&self, i: usize, p: Point3, t: Vec3) -> Option<RadialFrame<'_>> {
        let fin = self.fins[i];
        let f = self.c.fins.get(fin.index())?;
        let surf = fin_surface(self.c, self.geo, fin)?;
        let reach = self.c.edges.get(f.edge.index())?.tolerance.max(1e-9) * 1e3;
        let uv = seeded_uv(self.geo, f.source_half_edge, surf, p, reach)
            .or_else(|| surf.geometry.closest_point(&p))?;
        let n = surf.folded_normal(uv.0, uv.1)?;
        let tf = if f.forward { t } else { t * -1.0 };
        let d = n.cross(&tf);
        let d = (d - t * d.dot(&t)).normalized()?;
        Some(RadialFrame {
            d,
            surface: &surf.geometry,
            uv,
        })
    }

    fn clearance(&self, i: usize, p: Point3) -> f64 {
        super::radial_stations::face_clearance(self.c, self.fins[i], p)
    }

    fn chord(&self) -> f64 {
        let e = &self.c.edges[self.edge.index()];
        self.c.vertices[e.start.index()]
            .position
            .distance_to(&self.c.vertices[e.end.index()].position)
    }

    fn all_planar(&self) -> bool {
        self.fins
            .iter()
            .all(|&fin| super::radial_stations::on_plane(fin_surface(self.c, self.geo, fin)))
    }
}

/// Pcurve stations a fin's projection seed is chosen from.
const SEED_STATIONS: usize = 16;

/// The station's parameters on `surf` by single-start projection from
/// the fin's graph pcurve station whose image lies nearest `p`, kept
/// when its foot lies within `reach`. A multi-start projection onto a
/// multi-turn NURBS wall costs in proportion to its turns, once per fin
/// per station (PK-1039: 0.58 s / 2.12 s of a 40 / 80 thread cut, now
/// 0.05 s / 0.08 s). The pcurve's parameter need not track the edge
/// curve's (record uv fits drift: PK-1031), hence the nearest station
/// rather than the one at the station's fraction; and either half of
/// the graph edge may carry this face's pcurve (the other carries its
/// partner face's, whose uv means nothing here: seeds 1e6 away).
fn seeded_uv(
    geo: &FaceGeometry<'_>,
    half_edge: Option<usize>,
    surf: &SurfaceOnFace,
    p: Point3,
    reach: f64,
) -> Option<(f64, f64)> {
    let graph = geo.graph?;
    let he = graph.half_edges.get(half_edge?)?;
    let off = |q: (f64, f64)| surf.geometry.point_at(q.0, q.1).distance_to(&p);
    let (u0, v0) = [Some(he), graph.half_edges.get(he.twin)]
        .into_iter()
        .flatten()
        .filter_map(|h| h.pcurve.as_ref())
        .filter(|pc| {
            let (a, b) = pc.parameter_range();
            a.is_finite() && b.is_finite()
        })
        .flat_map(|pc| {
            let (a, b) = pc.parameter_range();
            (0..=SEED_STATIONS)
                .map(move |k| pc.point_at(a + (b - a) * k as f64 / SEED_STATIONS as f64))
        })
        .map(|q| (q.x, q.y))
        .min_by(|x, y| off(*x).total_cmp(&off(*y)))?;
    popkern_core::geometry::uv_projection::project_point_to_uv_near(&surf.geometry, &p, u0, v0)
        .filter(|&q| off(q) <= reach)
}

/// Whether `fin`'s face's Ahead side faces the next wedge about its
/// edge. With `d = n × tf` normal to `t`, `t × d` is `+n` when the fin
/// runs along `t` and `−n` against it, so the outward normal points
/// into the next wedge exactly when the face's folds outward and the
/// fin runs forward, or neither: topology alone, no evaluation.
fn ahead_faces_next(c: &CellComplex, fin: FinId) -> bool {
    let f = &c.fins[fin.index()];
    c.faces[f.face.index()].folded_is_outward == f.forward
}

/// Order the fins of every edge with three or more fins, and record
/// every fin's side-facing (two-fin edges need it too). `tie_angle` is
/// the angular tolerance for first-order ties; the tie-ladder probe
/// distance is relative to each edge's length, floored at a hundred
/// linear tolerances — never an absolute constant.
pub(crate) fn order_radially(
    c: &mut CellComplex,
    geo: &FaceGeometry<'_>,
    tie_angle: f64,
    linear: f64,
) -> Result<RadialSides, RadialDecline> {
    order_radially_with(c, geo, tie_angle, linear, TieLadder::default())
}

/// `order_radially` with the given ladder rungs.
pub(crate) fn order_radially_with(
    c: &mut CellComplex,
    geo: &FaceGeometry<'_>,
    tie_angle: f64,
    linear: f64,
    ladder: TieLadder,
) -> Result<RadialSides, RadialDecline> {
    let mut sides = RadialSides {
        ahead_faces_next: vec![false; c.fins.len()],
    };
    let tolerance = RadialTolerance {
        tie_angle,
        linear,
        ladder,
    };
    for ei in 0..c.edges.len() {
        let edge = EdgeId::from_index(ei);
        let fins = c.edges[ei].radial.clone();
        // One or two fins: any cyclic order is the order, and the sides
        // need no geometry. Evaluating every fin's face cost a
        // projection onto it per edge (PK-1039: on a multi-turn thread
        // cut, 0.58 s / 2.12 s at 205 / 365 edges, the complex's whole
        // cost).
        if fins.len() >= 3 {
            let src = ComplexEdgeFins {
                c,
                geo,
                edge,
                fins: &fins,
            };
            let order = order_fins(&src, tolerance)
                .map_err(|d| RadialDecline::from_engine(edge, &fins, d))?;
            c.edges[ei].radial = order.into_iter().map(|i| fins[i]).collect();
        }
        for &fin in &c.edges[ei].radial {
            sides.ahead_faces_next[fin.index()] = ahead_faces_next(c, fin);
        }
        if !fins.is_empty() {
            c.edges[ei].radial_ordered = true;
        }
    }
    Ok(sides)
}

#[cfg(test)]
#[path = "radial_tests.rs"]
mod tests;
