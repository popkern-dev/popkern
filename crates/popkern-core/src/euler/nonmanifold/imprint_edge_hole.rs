//! Closing a hole into a disc (split from `imprint_edge.rs` under the
//! 500-line cap): an imprinted edge between two vertices of one hole,
//! a closed contact curve's second half, splits the hole's loop in two,
//! and the half wound as the face's outer loop bounds the new disc face.

use crate::topology::{BRep, FaceHandle, VertexHandle};

impl BRep {
    /// After an edge closed a hole of `face` in two: the loop wound as
    /// the face's outer loop is (by its signed chart area) bounds the
    /// disc, the new face's outer loop; the other stays the hole. The
    /// split hands the new face either half, so the two swap faces when
    /// it took the hole's.
    pub(super) fn settle_closed_hole(
        &mut self,
        face: FaceHandle,
        new: FaceHandle,
    ) -> Result<(), &'static str> {
        let surface = self
            .faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| s.geometry.clone())
            .ok_or("the face has no surface")?;
        let outer = self.faces[face].outer_loop;
        let new_loop = self.faces[new].outer_loop;
        let sign = |b: &Self, l| {
            b.loop_signed_uv_area(l, &surface)
                .map(|a| a.signed.signum())
        };
        if sign(self, new_loop)? == sign(self, outer)? {
            return Ok(());
        }
        // The hole's half: the face's inner loop sharing the new edge's
        // vertices with the new face's loop.
        let new_vertices: Vec<VertexHandle> = self.hole_loop_vertices(new_loop);
        let Some(&hole) = self.faces[face].inner_loops.iter().find(|&&l| {
            self.hole_loop_vertices(l)
                .iter()
                .any(|v| new_vertices.contains(v))
        }) else {
            return Err("the closed hole's other half is not a hole of the face");
        };
        self.faces[new].outer_loop = hole;
        self.loops[hole].face = new;
        if let Some(slot) = self.faces[face]
            .inner_loops
            .iter_mut()
            .find(|l| **l == hole)
        {
            *slot = new_loop;
        }
        self.loops[new_loop].face = face;
        Ok(())
    }

    /// The vertices loop `l` leaves, in order.
    fn hole_loop_vertices(&self, l: crate::topology::LoopHandle) -> Vec<VertexHandle> {
        let Some(start) = self.loops.get(l).map(|d| d.half_edge) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        let mut h = start;
        for _ in 0..=self.half_edges.len() {
            let Some(d) = self.half_edges.get(h) else {
                break;
            };
            out.push(d.origin);
            h = d.next;
            if h == start {
                break;
            }
        }
        out
    }
}
