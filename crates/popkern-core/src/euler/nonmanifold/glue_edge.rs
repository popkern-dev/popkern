//! glue_edges / unglue_edge: merge two coincident edges into one ring,
//! and split a run of a ring off into an edge of its own.
//!
//! Gluing is how separate pieces become one general body (sewing two
//! sheets, two solids along a shared edge): the glued-away edge's fins
//! join the kept edge's ring in radial order about it, so each new
//! wedge lies between geometric neighbours. Both sides of every wedge
//! must lie in one region (two regions meeting in a wedge would mean
//! the faces cross, which gluing does not decide), and shells of that
//! region meeting in a wedge join. The ends must already be one vertex
//! each ([`BRep::glue_vertices`]), and the curves coincide within the
//! linear tolerance.

use std::collections::{BTreeMap, BTreeSet};

use super::EulerError;
use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, closest_point_in_window};
use crate::topology::{BRep, EdgeHandle, FaceSide, HalfEdgeHandle, ShellHandle};

/// What [`BRep::glue_edges`] changed.
#[derive(Debug, Clone, PartialEq)]
pub struct GlueEdgeResult {
    /// The largest distance from the glued-away curve to the kept one,
    /// now the kept edge's tolerance where it exceeds the body's.
    pub distance: f64,
    /// The kept edge's ring, in radial order.
    pub ring: Vec<HalfEdgeHandle>,
    /// Shells that joined another of their region and went.
    pub merged_shells: Vec<ShellHandle>,
}

/// Samples taken along a curve to measure coincidence.
const SAMPLES: usize = 9;

impl BRep {
    /// Glue edge `gone` onto `keep`: its fins join `keep`'s ring in
    /// radial order and `gone` goes. Counts: edges −1, shells −1 per
    /// shell joined.
    ///
    /// # Errors
    ///
    /// Unknown or equal edges, [`EulerError::NotAFaceEdge`],
    /// [`EulerError::EdgeEndsDiffer`] (glue the vertices first),
    /// [`EulerError::NotCoincident`] with the largest distance,
    /// [`EulerError::RadialOrder`] when the merged ring's order is
    /// undecided, [`EulerError::RegionConflict`] naming the wedge whose
    /// sides lie in two regions.
    pub fn glue_edges(
        &mut self,
        keep: EdgeHandle,
        gone: EdgeHandle,
    ) -> Result<GlueEdgeResult, EulerError> {
        for e in [keep, gone] {
            if !self.edges.contains_key(e) {
                return Err(EulerError::UnknownEdge(e));
            }
            if self.radial_fins(e).is_empty() {
                return Err(EulerError::NotAFaceEdge { edge: e });
            }
        }
        if keep == gone {
            return Err(EulerError::SameEdge { edge: keep });
        }
        let ends = |e| self.edge_vertices(e);
        let (k, g) = (ends(keep), ends(gone));
        let same_ends = match (k, g) {
            (Some((a, b)), Some((c, d))) => (a == c && b == d) || (a == d && b == c),
            _ => false,
        };
        if !same_ends {
            return Err(EulerError::EdgeEndsDiffer { keep, gone });
        }
        let distance = self.curve_distance(gone, keep);
        let tolerance = self
            .tolerance
            .linear()
            .max(self.edges[keep].tolerance.unwrap_or(0.0))
            .max(self.edges[gone].tolerance.unwrap_or(0.0));
        if distance > tolerance {
            return Err(EulerError::NotCoincident {
                distance,
                tolerance,
            });
        }
        let mut fins = self.radial_fins(keep);
        fins.extend(self.radial_fins(gone));
        let order = self
            .radial_order_with_sense(keep, fins, &self.tolerance)
            .map_err(|e| EulerError::RadialOrder(Box::new(e)))?;
        let joins = self.wedge_joins(keep, &order)?;
        // Every check passed: rewire.
        let ring: Vec<HalfEdgeHandle> = order.iter().map(|&(f, _)| f).collect();
        for &f in &ring {
            self.half_edges[f].edge = keep;
        }
        self.link_ring(&ring);
        self.edges.remove(gone);
        if distance > self.tolerance.linear() {
            let t = self.edges[keep].tolerance.unwrap_or(0.0).max(distance);
            self.edges[keep].tolerance = Some(t);
        }
        let merged_shells = self.join_shells(joins);
        Ok(GlueEdgeResult {
            distance,
            ring,
            merged_shells,
        })
    }

    /// Split the consecutive run `fins` of `edge`'s ring off onto a new
    /// edge with the same curve, the inverse of [`BRep::glue_edges`]:
    /// each ring closes up in its order, and shells the split leaves in
    /// pieces separate. Returns the new edge. Counts: edges +1, shells
    /// +1 per piece separated.
    ///
    /// # Errors
    ///
    /// [`EulerError::UnknownEdge`], [`EulerError::NotARingRun`] when the
    /// fins are not a proper consecutive run of the ring,
    /// [`EulerError::UnglueSplitsRegion`] when a bounded region's shell
    /// would fall into pieces (deciding which piece bounds it is
    /// geometry); the body is then unchanged.
    pub fn unglue_edge(
        &mut self,
        edge: EdgeHandle,
        fins: &[HalfEdgeHandle],
    ) -> Result<EdgeHandle, EulerError> {
        if !self.edges.contains_key(edge) {
            return Err(EulerError::UnknownEdge(edge));
        }
        let ring = self.radial_fins(edge);
        let (run, rest) = split_run(&ring, fins).ok_or(EulerError::NotARingRun { edge })?;
        let mut work = self.clone();
        let senses: BTreeMap<HalfEdgeHandle, bool> = ring
            .iter()
            .map(|&f| (f, work.fin_is_forward(f).unwrap_or(true)))
            .collect();
        let shells = work.shells_about(edge);
        let mut new = work.edges[edge].clone();
        new.half_edge = HalfEdgeHandle::default();
        let new = work.insert_edge_tracked(new);
        for &f in &run {
            work.half_edges[f].edge = new;
        }
        work.set_first_fin(new, &run, &senses);
        work.set_first_fin(edge, &rest, &senses);
        work.link_ring(&run);
        work.link_ring(&rest);
        for shell in shells {
            let region = work.shells[shell].region;
            let made = work.separate_shell_pieces(shell);
            if !made.is_empty() && region != work.outside {
                return Err(EulerError::UnglueSplitsRegion { region });
            }
        }
        *self = work;
        Ok(new)
    }

    /// Make `edge`'s first fin one of `fins` that runs its stored curve
    /// forward, reversing the curve's parameter window when none does.
    fn set_first_fin(
        &mut self,
        edge: EdgeHandle,
        fins: &[HalfEdgeHandle],
        senses: &BTreeMap<HalfEdgeHandle, bool>,
    ) {
        let forward = fins.iter().copied().find(|f| senses.get(f) == Some(&true));
        match forward {
            Some(f) => self.edges[edge].half_edge = f,
            None => {
                if let Some(&f) = fins.first() {
                    self.edges[edge].half_edge = f;
                    if let Some(c) = self.edges[edge].curve.as_mut() {
                        std::mem::swap(&mut c.t_start, &mut c.t_end);
                    }
                }
            }
        }
    }

    /// The shells holding a side of a face on `edge`.
    fn shells_about(&self, edge: EdgeHandle) -> BTreeSet<ShellHandle> {
        self.radial_fins(edge)
            .iter()
            .filter_map(|&f| self.fin_face(f))
            .filter_map(|f| self.face_sides(f))
            .flat_map(|s| [s.back, s.front])
            .collect()
    }

    /// The shell pairs the wedges of `order` about `edge` join, or the
    /// wedge whose sides lie in two regions.
    pub(super) fn wedge_joins(
        &self,
        edge: EdgeHandle,
        order: &[(HalfEdgeHandle, bool)],
    ) -> Result<Vec<(ShellHandle, ShellHandle)>, EulerError> {
        let n = order.len();
        let mut joins = Vec::new();
        for i in 0..n {
            let ((f, ff), (g, fg)) = (order[i], order[(i + 1) % n]);
            let from_side = if ff { FaceSide::Front } else { FaceSide::Back };
            let to_side = if fg { FaceSide::Back } else { FaceSide::Front };
            let side_shell = |fin, side| {
                self.fin_face(fin)
                    .and_then(|face| self.face_side_shell(face, side))
            };
            let (Some(a), Some(b)) = (side_shell(f, from_side), side_shell(g, to_side)) else {
                return Err(EulerError::NotAFaceEdge { edge });
            };
            let (ra, rb) = (self.shells[a].region, self.shells[b].region);
            if ra != rb {
                return Err(EulerError::RegionConflict {
                    edge,
                    fins: (f, g),
                    regions: (ra, rb),
                });
            }
            if a != b {
                joins.push((a, b));
            }
        }
        Ok(joins)
    }

    /// Join each pair's shells (of one region), keeping a region's
    /// outer shell where it is among them. Returns the shells that went.
    pub(super) fn join_shells(
        &mut self,
        joins: Vec<(ShellHandle, ShellHandle)>,
    ) -> Vec<ShellHandle> {
        let mut gone = Vec::new();
        let mut alias: BTreeMap<ShellHandle, ShellHandle> = BTreeMap::new();
        let resolve = |alias: &BTreeMap<ShellHandle, ShellHandle>, mut s: ShellHandle| {
            while let Some(&t) = alias.get(&s) {
                s = t;
            }
            s
        };
        for (a, b) in joins {
            let (a, b) = (resolve(&alias, a), resolve(&alias, b));
            if a == b {
                continue;
            }
            let region = self.shells[a].region;
            let outer = self.regions.get(region).map(|r| r.shell);
            let (target, other) = if outer == Some(b) { (b, a) } else { (a, b) };
            self.merge_shell_into(target, other);
            alias.insert(other, target);
            gone.push(other);
        }
        gone
    }

    /// The largest distance from `from`'s curve (sampled) to `to`'s.
    pub(super) fn curve_distance(&self, from: EdgeHandle, to: EdgeHandle) -> f64 {
        let points = |e: EdgeHandle| -> Vec<Point3> {
            let edge = &self.edges[e];
            match &edge.curve {
                Some(c) => (0..SAMPLES)
                    .map(|k| {
                        let s = k as f64 / (SAMPLES - 1) as f64;
                        c.geometry.point_at(c.t_start + (c.t_end - c.t_start) * s)
                    })
                    .collect(),
                None => self
                    .edge_vertices(e)
                    .map(|(a, b)| {
                        let (a, b) = (self.vertices[a].position, self.vertices[b].position);
                        (0..SAMPLES)
                            .map(|k| a + (b - a) * (k as f64 / (SAMPLES - 1) as f64))
                            .collect()
                    })
                    .unwrap_or_default(),
            }
        };
        let target = &self.edges[to];
        let distance_to = |p: Point3| -> f64 {
            match &target.curve {
                Some(c) => {
                    let t = closest_point_in_window(&c.geometry, &p, c.t_start, c.t_end);
                    c.geometry.point_at(t).distance_to(&p)
                }
                None => self.edge_vertices(to).map_or(f64::INFINITY, |(a, b)| {
                    let (a, b) = (self.vertices[a].position, self.vertices[b].position);
                    let ab = b - a;
                    let k = ((p - a).dot(&ab) / ab.dot(&ab).max(1e-300)).clamp(0.0, 1.0);
                    (a + ab * k).distance_to(&p)
                }),
            }
        };
        points(from)
            .into_iter()
            .map(distance_to)
            .fold(0.0, f64::max)
    }
}

/// `fins` as a proper consecutive run of the cyclic `ring`, with the
/// rest, both in ring order.
fn split_run(
    ring: &[HalfEdgeHandle],
    fins: &[HalfEdgeHandle],
) -> Option<(Vec<HalfEdgeHandle>, Vec<HalfEdgeHandle>)> {
    let n = ring.len();
    let set: BTreeSet<HalfEdgeHandle> = fins.iter().copied().collect();
    if set.is_empty()
        || set.len() >= n
        || set.len() != fins.len()
        || !set.iter().all(|f| ring.contains(f))
    {
        return None;
    }
    // Start the run where the ring enters it.
    let start =
        (0..n).find(|&i| set.contains(&ring[i]) && !set.contains(&ring[(i + n - 1) % n]))?;
    let run: Vec<HalfEdgeHandle> = (0..set.len()).map(|k| ring[(start + k) % n]).collect();
    if !run.iter().all(|f| set.contains(f)) {
        return None;
    }
    let rest = (set.len()..n).map(|k| ring[(start + k) % n]).collect();
    Some((run, rest))
}
