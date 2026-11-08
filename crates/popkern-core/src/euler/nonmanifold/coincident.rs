//! Gluing what coincides (PK-1072): vertices within the linear
//! tolerance of each other, then edges sharing both ends whose curves
//! coincide.
//!
//! A writer may put two pieces of one body in contact without sharing
//! the contact's entities (OpenCASCADE writes two cavities touching
//! along an edge with an edge and two vertices of its own each). In a
//! valid body nothing distinct coincides, so the glue does nothing
//! there; where pieces touch it makes the general body they are: one
//! vertex, one ring in radial order, the shells of a region meeting in
//! a new wedge joined.

use std::collections::{BTreeMap, HashMap};

use super::EulerError;
use crate::topology::{BRep, EdgeHandle, VertexHandle};

/// What [`BRep::glue_coincident`] glued.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct CoincidentGlue {
    /// Vertices glued onto another.
    pub vertices: usize,
    /// Edges glued onto another.
    pub edges: usize,
}

impl BRep {
    /// Glue every vertex within the linear tolerance of another onto
    /// it, then every two edges sharing both ends whose curves coincide
    /// (edges between the same vertices along different curves stay
    /// apart). All or nothing: on an error the body is unchanged.
    ///
    /// # Errors
    ///
    /// The glue operators' declines (a ring whose order is undecided, a
    /// wedge between two regions).
    pub fn glue_coincident(&mut self) -> Result<CoincidentGlue, EulerError> {
        let mut work = self.clone();
        let mut out = CoincidentGlue::default();
        let tolerance = work.tolerance.linear();
        // Vertices, bucketed on a grid of the tolerance: a coincident pair
        // lies in one cell or a neighbouring one.
        let cell = |p: crate::geometry::Point3| {
            let k = |x: f64| (x / tolerance.max(1e-12)).floor() as i64;
            (k(p.x), k(p.y), k(p.z))
        };
        let mut grid: HashMap<(i64, i64, i64), Vec<VertexHandle>> = HashMap::new();
        for (v, data) in &work.vertices {
            grid.entry(cell(data.position)).or_default().push(v);
        }
        let keys: Vec<VertexHandle> = work.vertices.keys().collect();
        let mut alias: BTreeMap<VertexHandle, VertexHandle> = BTreeMap::new();
        for v in keys {
            if alias.contains_key(&v) || !work.vertices.contains_key(v) {
                continue;
            }
            let p = work.vertices[v].position;
            let (cx, cy, cz) = cell(p);
            let near: Vec<VertexHandle> = (-1..=1)
                .flat_map(|dx| (-1..=1).flat_map(move |dy| (-1..=1).map(move |dz| (dx, dy, dz))))
                .filter_map(|(dx, dy, dz)| grid.get(&(cx + dx, cy + dy, cz + dz)))
                .flatten()
                .copied()
                .filter(|&w| w != v && !alias.contains_key(&w))
                .collect();
            for w in near {
                let Some(q) = work.vertices.get(w).map(|d| d.position) else {
                    continue;
                };
                if p.distance_to(&q) <= tolerance {
                    work.glue_vertices(v, w)?;
                    alias.insert(w, v);
                    out.vertices += 1;
                }
            }
        }
        // Edges, grouped by their end vertices.
        loop {
            let mut by_ends: HashMap<(VertexHandle, VertexHandle), Vec<EdgeHandle>> =
                HashMap::new();
            for e in work.edges.keys() {
                if work.radial_fins(e).is_empty() {
                    continue;
                }
                if let Some((a, b)) = work.edge_vertices(e) {
                    by_ends
                        .entry(if a < b { (a, b) } else { (b, a) })
                        .or_default()
                        .push(e);
                }
            }
            let mut glued = false;
            for edges in by_ends.values().filter(|g| g.len() > 1) {
                for (i, &keep) in edges.iter().enumerate() {
                    for &gone in &edges[i + 1..] {
                        match work.glue_edges(keep, gone) {
                            Ok(_) => {
                                out.edges += 1;
                                glued = true;
                            }
                            // Different curves between the same vertices.
                            Err(EulerError::NotCoincident { .. }) => continue,
                            Err(e) => return Err(e),
                        }
                        break;
                    }
                    if glued {
                        break;
                    }
                }
                if glued {
                    break;
                }
            }
            if !glued {
                break;
            }
        }
        *self = work;
        Ok(out)
    }
}

#[cfg(test)]
#[path = "coincident_tests.rs"]
mod tests;
