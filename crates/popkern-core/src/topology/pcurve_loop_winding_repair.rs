//! Loop-winding repair postcondition for pcurve recompute.
//!
//! Open-edge pcurve directions are determined by endpoint matching (the
//! half-edge's own origin/destination disambiguate them), but a CLOSED
//! edge's two curve endpoints coincide, so authoring falls back to the
//! per-edge canonical bit — and a per-EDGE bit cannot express the
//! direction a face's CHART needs: a left-handed chart (u×v against the
//! folded normal — the cylinder's caps) honestly traverses a shared rim
//! the same 3D way as its right-handed neighbour, so no single
//! canonical assignment winds both faces honestly. Measured: every
//! re-authored hollow-body skin wound one lateral or cap loop against
//! its orientation flag (`WindingDishonest` on shell/offset assembly
//! output), while constructor-authored bodies stayed clean only because
//! their authored pcurves are never re-derived.
//!
//! The honest disambiguator is the winding contract itself (the
//! interior-left rule the winding-honesty audit measures): the outer
//! loop winds counter-clockwise in UV exactly when the stored
//! orientation flag is `true`. After chart coherence, each measurable
//! fully-re-authored loop that contradicts its face's flag has its
//! closed-edge members' pcurves reversed (a lossless, exact
//! locus-preserving flip) and coherence re-run so open members re-chain
//! across period branches. The repair is transactional: if the loop
//! still measures dishonest afterwards, every pcurve of the face is
//! restored — a wrong verdict must not half-rewrite a chart.

use std::collections::HashSet;

use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;

use super::brep::BRep;
use super::{FaceHandle, HalfEdgeHandle, LoopHandle};

impl BRep {
    /// Repair loop windings on `face_handle` after pcurve recompute by
    /// flipping closed-edge member pcurves (see module docs). Loops
    /// whose winding is unmeasurable (period-wrapping, pole-touching,
    /// area-degenerate — the audit's own guards) are left untouched.
    ///
    /// `freshly_authored` is the set of half-edges whose pcurves THIS
    /// recompute call wrote. Only loops that were re-authored IN FULL
    /// qualify: a partially re-authored loop (a boolean `fix_stale`
    /// cycle refreshing one stale member among trusted neighbours)
    /// carries directions consumers already rely on. The from-scratch
    /// case this repair exists for (assembly copies topology without
    /// pcurves) always authors whole loops.
    pub(crate) fn fix_loop_pcurve_winding(
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
                // A dishonest all-open loop means the endpoint-matched
                // authoring itself disagrees with the flag — that is a
                // flag/topology defect for the audit to surface, not a
                // direction ambiguity this repair owns.
                continue;
            }

            // Transaction: snapshot every pcurve of the face (coherence
            // re-chaining may translate members of any of its loops).
            let snapshot: Vec<(HalfEdgeHandle, Option<Pcurve>)> = self
                .face_half_edges(face_handle)
                .into_iter()
                .map(|heh| (heh, self.half_edges[heh].pcurve.clone()))
                .collect();

            for &heh in &closed {
                if let Some(pc) = self.half_edges[heh].pcurve.clone() {
                    self.half_edges[heh].pcurve = Some(pc.reversed());
                }
            }
            self.close_periodic_loop_winding(face_handle);

            if self.loop_winding_verdict(lh, &surface, desired_positive) != Some(true) {
                for (heh, pc) in snapshot {
                    self.half_edges[heh].pcurve = pc;
                }
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

    /// All half-edges of a face across its outer and inner loops.
    fn face_half_edges(&self, face_handle: FaceHandle) -> Vec<HalfEdgeHandle> {
        let Some(face) = self.faces.get(face_handle) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        out.extend(self.loop_members(face.outer_loop));
        for &lh in &face.inner_loops {
            out.extend(self.loop_members(lh));
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
