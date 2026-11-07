//! A body's edge as a radial-order source, and the ring it orders.
//!
//! The order is taken about the edge's own direction (its first fin
//! runs the stored curve forward, `t_start` → `t_end`), so the ring
//! [`BRep::order_ring`] links runs counter-clockwise about it seen from
//! the edge's end, as the wedge reading ([`BRep::edge_wedges`]) needs:
//! the wedge after a forward fin lies on its face's front.
//!
//! [`BRep::edge_wedges`]: crate::topology::BRep::edge_wedges

use super::{RadialFins, RadialFrame, RadialOrderDecline, RadialTolerance, TieLadder, order_fins};
use crate::geometry::curve::{CurveEval, CurveOnEdge, closest_point_in_window};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle, VertexHandle};

/// Why a body's edge could not be ordered.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum RadialOrderError {
    /// The edge is not in the body.
    #[error("edge {0:?} is not in the body")]
    UnknownEdge(EdgeHandle),
    /// A fin's face carries no surface to read its direction from.
    #[error("fin {fin:?} of edge {edge:?} lies on a face without a surface")]
    NoSurface {
        /// The edge.
        edge: EdgeHandle,
        /// The fin.
        fin: HalfEdgeHandle,
    },
    /// The engine declined; its fin indices are into `fins`.
    #[error("radial order of edge {edge:?} declined: {decline:?}")]
    Declined {
        /// The edge.
        edge: EdgeHandle,
        /// The edge's fins in ring order before ordering.
        fins: Vec<HalfEdgeHandle>,
        /// What the engine measured.
        decline: RadialOrderDecline,
    },
}

/// Pcurve stations a fin's projection seed is chosen from.
const SEED_STATIONS: usize = 16;

/// One edge of a body, read for ordering.
struct EdgeFins<'a> {
    brep: &'a BRep,
    edge: EdgeHandle,
    curve: Option<&'a CurveOnEdge>,
    ends: (Point3, Point3),
    end_vertices: (VertexHandle, VertexHandle),
    fins: Vec<HalfEdgeHandle>,
    surfaces: Vec<&'a SurfaceOnFace>,
    /// Whether each fin runs the edge's direction.
    forward: Vec<bool>,
    /// How far a projection may land from the station.
    reach: f64,
}

impl<'a> EdgeFins<'a> {
    /// `fins` about `edge`'s geometry: its own ring, or a ring a glue
    /// is about to make (fins of a coincident edge read against this
    /// edge's direction).
    fn new(
        brep: &'a BRep,
        edge: EdgeHandle,
        fins: Vec<HalfEdgeHandle>,
    ) -> Result<Self, RadialOrderError> {
        let e = brep
            .edges
            .get(edge)
            .ok_or(RadialOrderError::UnknownEdge(edge))?;
        let mut surfaces = Vec::with_capacity(fins.len());
        for &fin in &fins {
            let surface = brep
                .fin_face(fin)
                .and_then(|f| brep.faces.get(f))
                .and_then(|f| f.surface.as_ref())
                .ok_or(RadialOrderError::NoSurface { edge, fin })?;
            surfaces.push(surface);
        }
        let start_vertex = brep
            .half_edges
            .get(e.half_edge)
            .map(|h| h.origin)
            .ok_or(RadialOrderError::UnknownEdge(edge))?;
        let end_vertex = brep.fin_end(e.half_edge).unwrap_or(start_vertex);
        let start = brep
            .vertices
            .get(start_vertex)
            .map(|v| v.position)
            .ok_or(RadialOrderError::UnknownEdge(edge))?;
        let end = brep.vertices.get(end_vertex).map_or(start, |v| v.position);
        let mut src = Self {
            brep,
            edge,
            curve: e.curve.as_ref(),
            ends: (start, end),
            end_vertices: (start_vertex, end_vertex),
            fins,
            surfaces,
            forward: Vec::new(),
            reach: e.tolerance.unwrap_or(brep.tolerance.linear()).max(1e-9) * 1e3,
        };
        src.forward = (0..src.fins.len())
            .map(|i| src.fin_runs_forward(i))
            .collect::<Option<_>>()
            .ok_or_else(|| src.declined(RadialOrderDecline::Unevaluable { fin: None }))?;
        Ok(src)
    }

    fn declined(&self, decline: RadialOrderDecline) -> RadialOrderError {
        RadialOrderError::Declined {
            edge: self.edge,
            fins: self.fins.clone(),
            decline,
        }
    }

    /// Whether fin `i` runs along the edge: read off topology where it
    /// decides ([`BRep::fin_is_forward`] for the edge's own fins, the
    /// origin for another edge's fin on an open edge), otherwise (a
    /// closed edge of more than two fins) off the fin's pcurve image at
    /// mid-edge.
    fn fin_runs_forward(&self, i: usize) -> Option<bool> {
        let fin = self.brep.half_edges.get(self.fins[i])?;
        if fin.edge == self.edge {
            if let Some(forward) = self.brep.fin_is_forward(self.fins[i]) {
                return Some(forward);
            }
        } else if self.end_vertices.0 != self.end_vertices.1 {
            return match fin.origin {
                o if o == self.end_vertices.0 => Some(true),
                o if o == self.end_vertices.1 => Some(false),
                _ => None,
            };
        }
        let (p, t) = self.station(0.5)?;
        pcurve_sense(self.brep, self.fins[i], &self.surfaces[i].geometry, (p, t))
    }

    /// The station's parameters on fin `i`'s face: projected from the
    /// fin's pcurve station nearest `p`, kept within `reach`, else the
    /// surface's closest point.
    fn uv(&self, i: usize, p: Point3) -> Option<(f64, f64)> {
        let surface = &self.surfaces[i].geometry;
        let seeded = self
            .brep
            .half_edges
            .get(self.fins[i])
            .and_then(|h| h.pcurve.as_ref())
            .and_then(|pc| {
                let (a, b) = pc.parameter_range();
                let k = nearest_station(pc, surface, p)?;
                let q = pc.point_at(a + (b - a) * k as f64 / SEED_STATIONS as f64);
                crate::geometry::uv_projection::project_point_to_uv_near(surface, &p, q.x, q.y)
            })
            .filter(|&(u, v)| surface.point_at(u, v).distance_to(&p) <= self.reach);
        seeded.or_else(|| surface.closest_point(&p))
    }
}

/// Whether `fin`'s pcurve image runs along `t` where it passes `p`: its
/// uv tangent mapped through the surface's partials at the pcurve
/// station nearest `p`. `None` without a pcurve or with a tangent
/// normal to `t`.
fn pcurve_sense(
    brep: &BRep,
    fin: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    (p, t): (Point3, Vec3),
) -> Option<bool> {
    let pc = brep.half_edges.get(fin)?.pcurve.as_ref()?;
    let (a, b) = pc.parameter_range();
    let k = nearest_station(pc, surface, p)?;
    let s = a + (b - a) * k as f64 / SEED_STATIONS as f64;
    let duv = pc.tangent_at(s);
    let q = pc.point_at(s);
    let along = surface.partial_u(q.x, q.y) * duv.x + surface.partial_v(q.x, q.y) * duv.y;
    let dot = along.dot(&t);
    (dot.abs() > 0.0).then_some(dot > 0.0)
}

/// The pcurve station whose image lies nearest `p`.
fn nearest_station(
    pc: &crate::geometry::pcurve::Pcurve,
    surface: &SurfaceGeometry,
    p: Point3,
) -> Option<usize> {
    let (a, b) = pc.parameter_range();
    if !(a.is_finite() && b.is_finite()) {
        return None;
    }
    (0..=SEED_STATIONS).min_by(|&x, &y| {
        let off = |k: usize| {
            let q = pc.point_at(a + (b - a) * k as f64 / SEED_STATIONS as f64);
            surface.point_at(q.x, q.y).distance_to(&p)
        };
        off(x).total_cmp(&off(y))
    })
}

impl RadialFins for EdgeFins<'_> {
    fn fin_count(&self) -> usize {
        self.fins.len()
    }

    fn station(&self, s: f64) -> Option<(Point3, Vec3)> {
        match self.curve {
            Some(c) => {
                let t = c.t_start + (c.t_end - c.t_start) * s;
                let d = c.geometry.derivative_at(t, 1);
                let d = if c.t_end < c.t_start { d * -1.0 } else { d };
                Some((c.geometry.point_at(t), d.normalized()?))
            }
            None => {
                let (a, b) = self.ends;
                Some((a + (b - a) * s, (b - a).normalized()?))
            }
        }
    }

    fn frame(&self, fin: usize, p: Point3, t: Vec3) -> Option<RadialFrame<'_>> {
        let uv = self.uv(fin, p)?;
        let n = self.surfaces[fin].folded_normal(uv.0, uv.1)?;
        let tf = if self.forward[fin] { t } else { t * -1.0 };
        let d = n.cross(&tf);
        let d = (d - t * d.dot(&t)).normalized()?;
        Some(RadialFrame {
            d,
            surface: &self.surfaces[fin].geometry,
            uv,
        })
    }

    fn clearance(&self, fin: usize, p: Point3) -> f64 {
        let brep = self.brep;
        let Some(face) = brep
            .fin_face(self.fins[fin])
            .and_then(|f| brep.faces.get(f))
        else {
            return f64::INFINITY;
        };
        std::iter::once(face.outer_loop)
            .chain(face.inner_loops.iter().copied())
            .flat_map(|l| loop_fins(brep, l))
            .filter_map(|h| brep.half_edges.get(h).map(|h| h.edge))
            .filter(|&e| e != self.edge)
            .filter_map(|e| edge_distance(brep, e, p))
            .fold(f64::INFINITY, f64::min)
    }

    fn chord(&self) -> f64 {
        self.ends.0.distance_to(&self.ends.1)
    }

    fn all_planar(&self) -> bool {
        self.surfaces
            .iter()
            .all(|s| matches!(s.geometry, SurfaceGeometry::Plane { .. }))
    }
}

/// The fins of loop `l` in order (bounded by the fin count, so a
/// broken cycle ends).
fn loop_fins(brep: &BRep, l: crate::topology::LoopHandle) -> Vec<HalfEdgeHandle> {
    let Some(start) = brep.loops.get(l).map(|l| l.half_edge) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..brep.half_edges.len() {
        let Some(fin) = brep.half_edges.get(cursor) else {
            break;
        };
        out.push(cursor);
        cursor = fin.next;
        if cursor == start {
            break;
        }
    }
    out
}

/// The distance from `p` to edge `e`.
fn edge_distance(brep: &BRep, e: EdgeHandle, p: Point3) -> Option<f64> {
    let edge = brep.edges.get(e)?;
    if let Some(cv) = &edge.curve {
        let t = closest_point_in_window(&cv.geometry, &p, cv.t_start, cv.t_end);
        return Some(cv.geometry.point_at(t).distance_to(&p));
    }
    let (a, b) = brep.edge_vertices(e)?;
    let (a, b) = (
        brep.vertices.get(a)?.position,
        brep.vertices.get(b)?.position,
    );
    let ab = b - a;
    let k = ((p - a).dot(&ab) / ab.dot(&ab).max(1e-300)).clamp(0.0, 1.0);
    Some((a + ab * k).distance_to(&p))
}

impl BRep {
    /// The fins of `edge` in radial order about its direction, starting
    /// at its first fin: each fin's successor is the next face met
    /// turning counter-clockwise about the edge seen from its end.
    /// One or two fins are their own order. Ordered at the tolerance's
    /// angular tie and linear scale, agreeing at every station.
    ///
    /// # Errors
    ///
    /// [`RadialOrderError`]: an unknown edge, a fin on a face without a
    /// surface, or the engine's decline (an unresolved tie with its
    /// evidence, an order that changes along the edge, an unevaluable
    /// face).
    pub fn radial_order(
        &self,
        edge: EdgeHandle,
        tolerance: &Tolerance,
    ) -> Result<Vec<HalfEdgeHandle>, RadialOrderError> {
        self.radial_order_of(edge, self.radial_fins(edge), tolerance)
    }

    /// `fins` in radial order about `edge`'s direction, starting at the
    /// first: `edge`'s own ring, or the ring gluing a coincident edge's
    /// fins onto it would make. One or two fins are their own order.
    ///
    /// # Errors
    ///
    /// As [`BRep::radial_order`].
    pub fn radial_order_of(
        &self,
        edge: EdgeHandle,
        fins: Vec<HalfEdgeHandle>,
        tolerance: &Tolerance,
    ) -> Result<Vec<HalfEdgeHandle>, RadialOrderError> {
        if !self.edges.contains_key(edge) {
            return Err(RadialOrderError::UnknownEdge(edge));
        }
        if fins.len() < 3 {
            return Ok(fins);
        }
        Ok(self
            .radial_order_with_sense(edge, fins, tolerance)?
            .into_iter()
            .map(|(fin, _)| fin)
            .collect())
    }

    /// [`BRep::radial_order_of`] with each fin's sense: whether it runs
    /// along `edge`'s direction. Every fin's face must carry a surface,
    /// short rings too (the sense of another edge's fin on a closed
    /// edge is read off its pcurve image).
    ///
    /// # Errors
    ///
    /// As [`BRep::radial_order`].
    pub(crate) fn radial_order_with_sense(
        &self,
        edge: EdgeHandle,
        fins: Vec<HalfEdgeHandle>,
        tolerance: &Tolerance,
    ) -> Result<Vec<(HalfEdgeHandle, bool)>, RadialOrderError> {
        let src = EdgeFins::new(self, edge, fins)?;
        let order: Vec<usize> = if src.fins.len() < 3 {
            (0..src.fins.len()).collect()
        } else {
            order_fins(
                &src,
                RadialTolerance {
                    tie_angle: tolerance.angular(),
                    linear: tolerance.linear(),
                    ladder: TieLadder::default(),
                },
            )
            .map_err(|d| src.declined(d))?
        };
        Ok(order
            .into_iter()
            .map(|i| (src.fins[i], src.forward[i]))
            .collect())
    }

    /// Whether `fin` runs along its edge's direction (the way the edge's
    /// first fin runs the stored curve): read off topology where it
    /// decides ([`BRep::fin_is_forward`]), otherwise (a closed edge of
    /// more than two fins, whose ends are one vertex) off the fin's
    /// pcurve image at mid-edge against the curve's tangent. `None`
    /// when neither decides (no curve, pcurve or surface).
    pub fn fin_sense(&self, fin: HalfEdgeHandle) -> Option<bool> {
        self.fin_is_forward(fin)
            .or_else(|| self.fin_sense_by_pcurve(fin))
    }

    /// [`BRep::fin_sense`]'s geometric reading alone.
    pub(crate) fn fin_sense_by_pcurve(&self, fin: HalfEdgeHandle) -> Option<bool> {
        let edge = self.edges.get(self.half_edges.get(fin)?.edge)?;
        let c = edge.curve.as_ref()?;
        let t = 0.5 * (c.t_start + c.t_end);
        let d = c.geometry.derivative_at(t, 1);
        let d = if c.t_end < c.t_start { d * -1.0 } else { d };
        let surface = self
            .fin_face(fin)
            .and_then(|f| self.faces.get(f))
            .and_then(|f| f.surface.as_ref())?;
        pcurve_sense(
            self,
            fin,
            &surface.geometry,
            (c.geometry.point_at(t), d.normalized()?),
        )
    }

    /// Link `edge`'s ring in its radial order ([`BRep::radial_order`]).
    /// The body is unchanged when ordering declines.
    ///
    /// # Errors
    ///
    /// As [`BRep::radial_order`].
    pub fn order_ring(
        &mut self,
        edge: EdgeHandle,
        tolerance: &Tolerance,
    ) -> Result<(), RadialOrderError> {
        let order = self.radial_order(edge, tolerance)?;
        self.link_ring(&order);
        Ok(())
    }
}

#[cfg(test)]
#[path = "brep_order_tests.rs"]
mod tests;
