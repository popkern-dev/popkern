//! glue_vertices / unglue_vertex: merge two coincident vertices into
//! one, and split a vertex into one per use.
//!
//! A vertex is shared by every use reaching it: two solids touching at
//! a corner are one vertex with two fans, a wire standing on a face one
//! vertex with a fan and a wire use. Gluing moves every fin leaving the
//! glued-away vertex onto the kept one. Shells of face sides connect
//! through their edges, not their vertices, so they stay as they are;
//! wire-only shells of one region meeting at the vertex join (a wire's
//! shell is one connected set), into that region's face shell there
//! when there is one, and an acorn glued onto a used vertex is absorbed.
//! Ungluing gives each use but the first a vertex of its own at the
//! same place, and splits the wire shells that disconnects.

use std::collections::BTreeSet;

use super::EulerError;
use crate::topology::{BRep, RegionHandle, ShellHandle, Vertex, VertexHandle};

/// What [`BRep::glue_vertices`] changed.
#[derive(Debug, Clone, PartialEq)]
pub struct GlueVertexResult {
    /// The distance between the two vertices, now the kept vertex's
    /// tolerance where it exceeds the body's linear tolerance.
    pub distance: f64,
    /// Shells that went: absorbed acorns and wire-only shells merged
    /// into another shell of their region.
    pub removed_shells: Vec<ShellHandle>,
}

impl BRep {
    /// Glue `gone` onto `keep`: every fin leaving `gone` leaves `keep`
    /// and `gone` goes. The two must lie within the linear tolerance
    /// (the larger of the body's and the vertices' own). Counts:
    /// vertices −1, shells −1 per absorbed acorn or merged wire shell,
    /// acorns −1 per absorbed acorn.
    ///
    /// # Errors
    ///
    /// [`EulerError::UnknownVertex`], [`EulerError::SameVertex`],
    /// [`EulerError::NotCoincident`] with the measured distance,
    /// [`EulerError::VertexNotInRegion`] naming the other vertex when an
    /// acorn is glued onto a vertex none of whose uses lie in the
    /// acorn's region.
    pub fn glue_vertices(
        &mut self,
        keep: VertexHandle,
        gone: VertexHandle,
    ) -> Result<GlueVertexResult, EulerError> {
        for v in [keep, gone] {
            if !self.vertices.contains_key(v) {
                return Err(EulerError::UnknownVertex(v));
            }
        }
        if keep == gone {
            return Err(EulerError::SameVertex { vertex: keep });
        }
        let (a, b) = (&self.vertices[keep], &self.vertices[gone]);
        let distance = a.position.distance_to(&b.position);
        let tolerance = self
            .tolerance
            .linear()
            .max(a.tolerance.unwrap_or(0.0))
            .max(b.tolerance.unwrap_or(0.0));
        if distance > tolerance {
            return Err(EulerError::NotCoincident {
                distance,
                tolerance,
            });
        }
        let acorns = [self.acorn_shell(keep), self.acorn_shell(gone)];
        // An acorn must meet a use of its region at the other vertex.
        for (acorn, other) in [(acorns[0], gone), (acorns[1], keep)] {
            if let Some(s) = acorn {
                let region = self.shells[s].region;
                let other_acorn = self.acorn_shell(other);
                let meets = match other_acorn {
                    Some(o) => self.shells[o].region == region,
                    None => self
                        .shells_at(other)
                        .iter()
                        .any(|&t| self.shells[t].region == region),
                };
                if !meets {
                    return Err(EulerError::VertexNotInRegion {
                        vertex: other,
                        region,
                    });
                }
            }
        }
        let mut removed = Vec::new();
        match acorns {
            // Two acorns: the glued-away one's shell goes.
            [Some(_), Some(g)] => removed.push(g),
            [Some(k), None] => removed.push(k),
            [None, Some(g)] => removed.push(g),
            [None, None] => {}
        }
        for &s in &removed {
            self.shells.remove(s);
        }
        for fin in self.half_edges.values_mut() {
            if fin.origin == gone {
                fin.origin = keep;
            }
        }
        self.vertices.remove(gone);
        if distance > self.tolerance.linear() {
            let t = self.vertices[keep].tolerance.unwrap_or(0.0).max(distance);
            self.vertices[keep].tolerance = Some(t);
        }
        if !self.half_edges.contains_key(self.vertices[keep].outgoing) {
            let first = self
                .half_edges
                .iter()
                .find(|(_, h)| h.origin == keep)
                .map(|(k, _)| k)
                .unwrap_or_default();
            self.vertices[keep].outgoing = first;
        }
        removed.extend(self.join_wire_shells_at(keep));
        Ok(GlueVertexResult {
            distance,
            removed_shells: removed,
        })
    }

    /// Split `vertex` into one vertex per use: the first fan (or, with
    /// no fan, the first wire edge) keeps it, every other fan and wire
    /// end gets a new vertex at the same place, and wire shells the
    /// split disconnects split. Returns the new vertices, none when the
    /// vertex has one use. Counts: vertices +1 per further use, shells
    /// +1 per wire piece split off.
    ///
    /// # Errors
    ///
    /// [`EulerError::UnknownVertex`].
    pub fn unglue_vertex(&mut self, vertex: VertexHandle) -> Result<Vec<VertexHandle>, EulerError> {
        let Some(v) = self.vertices.get(vertex).cloned() else {
            return Err(EulerError::UnknownVertex(vertex));
        };
        let mut groups = self.fan_groups(vertex, &self.face_uses_leaving(vertex));
        groups.extend(
            self.half_edges
                .iter()
                .filter(|(_, h)| h.origin == vertex && self.is_wire_edge(h.edge))
                .map(|(k, _)| vec![k]),
        );
        if groups.len() < 2 {
            return Ok(Vec::new());
        }
        let wire_shells: BTreeSet<ShellHandle> = groups
            .iter()
            .flatten()
            .filter_map(|&h| self.wire_edge_shell(self.half_edges[h].edge))
            .collect();
        let mut made = Vec::new();
        for group in &groups[1..] {
            let fresh = self.vertices.insert(Vertex {
                outgoing: group[0],
                ..v.clone()
            });
            for &h in group {
                self.half_edges[h].origin = fresh;
            }
            made.push(fresh);
        }
        self.vertices[vertex].outgoing = groups[0][0];
        for shell in wire_shells {
            self.split_wire_shell(shell);
        }
        Ok(made)
    }

    /// The shells reaching `vertex`: its acorn's, its wire edges', and
    /// both sides of its faces'.
    fn shells_at(&self, vertex: VertexHandle) -> BTreeSet<ShellHandle> {
        let mut out = BTreeSet::new();
        out.extend(self.acorn_shell(vertex));
        for fin in self.half_edges.values().filter(|h| h.origin == vertex) {
            if let Some(s) = self.wire_edge_shell(fin.edge) {
                out.insert(s);
            } else if let Some(f) = self
                .loops
                .get(fin.loop_)
                .and_then(|l| self.faces.get(l.face))
            {
                out.insert(f.back);
                out.insert(f.front);
            }
        }
        out
    }

    /// Join the wire-only shells meeting at `vertex`, region by region,
    /// into that region's face shell there or the first of them.
    /// Returns the shells that went.
    fn join_wire_shells_at(&mut self, vertex: VertexHandle) -> Vec<ShellHandle> {
        let at = self.shells_at(vertex);
        let regions: BTreeSet<RegionHandle> = at.iter().map(|&s| self.shells[s].region).collect();
        let mut removed = Vec::new();
        for region in regions {
            let mine: Vec<ShellHandle> = at
                .iter()
                .copied()
                .filter(|&s| self.shells[s].region == region)
                .collect();
            let wire_only: Vec<ShellHandle> = mine
                .iter()
                .copied()
                .filter(|&s| !self.shell_has_faces(s) && !self.shells[s].wire_edges.is_empty())
                .collect();
            let target = mine
                .iter()
                .copied()
                .find(|&s| self.shell_has_faces(s))
                .or(wire_only.first().copied());
            let Some(target) = target else {
                continue;
            };
            for s in wire_only.into_iter().filter(|&s| s != target) {
                let moved = std::mem::take(&mut self.shells[s].wire_edges);
                self.shells[target].wire_edges.extend(moved);
                self.shells.remove(s);
                removed.push(s);
            }
        }
        removed
    }
}
