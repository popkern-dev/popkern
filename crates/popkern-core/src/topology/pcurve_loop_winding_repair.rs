//! Loop-winding postcondition for full pcurve re-authoring.
//!
//! Open-edge pcurve directions follow endpoint matching (the half-edge's
//! own origin/destination), and a CLOSED edge's follow the per-edge
//! canonical bit (its `half_edge` runs the curve forward, the twin
//! backward). Both are the half-edge's 3-D traversal, so re-authored
//! pcurves are honest by construction: two faces sharing an edge use it
//! oppositely in 3-D, and a left-handed chart (u×v against the material
//! normal) shows that as a clockwise loop — which is what its flag says.
//!
//! The winding contract (outer loop counter-clockwise in UV exactly when
//! the orientation flag is `true`) then has one free side left: the
//! flag a producer AUTHORED. A fully re-authored outer loop that
//! measures against its face's flag means the authored flag is wrong,
//! and the flag follows the loop. This used to be the other way round —
//! the closed members' pcurves were reversed until the loop agreed with
//! the flag — on the premise that the canonical bit could not express a
//! chart's direction. That premise was wrong, and the flip made the
//! wrong side true: a shelled cone's opening wall carried an authored
//! flag against its material, both of its rims were rewritten to run
//! against their edges, and the served body read winding-honest while
//! its pcurves contradicted its topology (PK-878; the outward shell of
//! the same cone exported to OpenCASCADE at 4.97 for 0.78).

use std::collections::HashSet;

use crate::geometry::surface::SurfaceGeometry;

use super::brep::BRep;
use super::{FaceHandle, HalfEdgeHandle, LoopHandle};

impl BRep {
    /// After a full pcurve recompute of `face_handle`, give the face the
    /// flag its re-authored outer loop winds with (see module docs).
    /// Loops whose winding is unmeasurable (period-wrapping,
    /// pole-touching, area-degenerate — the audit's own guards) leave
    /// the flag as authored.
    ///
    /// `freshly_authored` is the set of half-edges whose pcurves THIS
    /// recompute call wrote. Only loops that were re-authored IN FULL
    /// qualify: a partially re-authored loop (a boolean `fix_stale`
    /// cycle refreshing one stale member among trusted neighbours)
    /// carries directions consumers already rely on. The from-scratch
    /// case this repair exists for (assembly copies topology without
    /// pcurves) always authors whole loops.
    pub(crate) fn fix_reauthored_loop_flag(
        &mut self,
        face_handle: FaceHandle,
        freshly_authored: &HashSet<HalfEdgeHandle>,
    ) {
        let Some(face) = self.faces.get(face_handle) else {
            return;
        };
        let Some(sof) = &face.surface else {
            return;
        };
        let orientation = sof.orientation;
        let surface = sof.geometry.clone();
        // OUTER loops only. The interior-left rule implies inner loops
        // wind opposite, but the winding-honesty audit has never
        // measured inner loops and existing producers (boolean hole
        // rims) carry hole pcurves whose direction their consumers read
        // by containment, not winding. Inner-loop winding gets a
        // contract before it gets a repair.
        let loops: Vec<(LoopHandle, bool)> = vec![(face.outer_loop, true)];

        for (lh, is_outer) in loops {
            let members = self.loop_members(lh);
            if members.is_empty() || !members.iter().all(|heh| freshly_authored.contains(heh)) {
                continue; // partially trusted loop: never re-orient
            }
            let desired_positive = is_outer == orientation;
            if self.loop_winding_verdict(lh, &surface, desired_positive) != Some(false) {
                continue; // honest or unmeasurable
            }
            let closed: Vec<HalfEdgeHandle> = members
                .iter()
                .copied()
                .filter(|&heh| self.half_edge_on_closed_edge(heh))
                .collect();
            if closed.is_empty() {
                // A dishonest all-open loop is left to the serving
                // audit: no closed member, so not this re-authoring
                // path's measured population.
                continue;
            }

            // The honest pcurves stand; the authored flag follows them.
            if let Some(sof) = self
                .faces
                .get_mut(face_handle)
                .and_then(|f| f.surface.as_mut())
            {
                sof.orientation = !sof.orientation;
            }
        }
    }

    /// `Some(true)` = winds as desired, `Some(false)` = measurably
    /// dishonest, `None` = no verdict (degenerate/wrapping/pole).
    fn loop_winding_verdict(
        &self,
        lh: LoopHandle,
        surface: &SurfaceGeometry,
        desired_positive: bool,
    ) -> Option<bool> {
        let area = self.loop_signed_uv_area(lh, surface).ok()?;
        if area.signed.abs() <= area.perimeter.powi(2) * 1.0e-12 {
            return None;
        }
        Some((area.signed > 0.0) == desired_positive)
    }

    /// Half-edge handles of one loop, in traversal order.
    fn loop_members(&self, lh: LoopHandle) -> Vec<HalfEdgeHandle> {
        let Some(loop_data) = self.loops.get(lh) else {
            return Vec::new();
        };
        let start = loop_data.half_edge;
        if !self.half_edges.contains_key(start) {
            return Vec::new();
        }
        let mut out = vec![start];
        let mut cur = self.half_edges[start].next;
        while cur != start
            && self.half_edges.contains_key(cur)
            && out.len() <= self.half_edges.len()
        {
            out.push(cur);
            cur = self.half_edges[cur].next;
        }
        out
    }

    /// Whether the half-edge's owning edge carries a CLOSED 3D curve
    /// (coincident endpoints — the direction-ambiguous class).
    fn half_edge_on_closed_edge(&self, heh: HalfEdgeHandle) -> bool {
        use crate::geometry::curve::CurveEval;
        let Some(he) = self.half_edges.get(heh) else {
            return false;
        };
        let Some(curve) = self.edges.get(he.edge).and_then(|e| e.curve.as_ref()) else {
            return false;
        };
        if self.half_edges[heh].pcurve.is_none() {
            return false;
        }
        let start = curve.geometry.point_at(curve.t_start);
        let end = curve.geometry.point_at(curve.t_end);
        start.distance_to(&end) < 1e-9
    }
}
