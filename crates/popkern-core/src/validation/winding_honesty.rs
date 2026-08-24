//! Winding-flag honesty audit — the adjacent-face consistency half of
//! check 7 that `validate_face_orientation` defers.
//!
//! The manifold contract every flag consumer rests on (edge convexity
//! classification, mesher triangle winding, STEP `same_sense`) is the
//! interior-left rule: walking a face's own half-edges with the
//! flag-folded surface normal up, the face interior lies on the LEFT.
//! Equivalently: the outer loop traversed in topological order winds
//! counter-clockwise in UV exactly when the stored orientation flag is
//! `true`. A face violating this is *flag-dishonest*: its mesh
//! normals/winding come out mirrored against its neighbours, its
//! measured volume contribution flips sign, and classifiers that fold
//! the flag mirror their into-face frames (the PK-238 wrong-branch
//! class). This module measures the invariant; it does not repair.

use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle};

/// Sample count per half-edge pcurve when polygonizing the loop.
const SAMPLES_PER_HALF_EDGE: usize = 8;

/// Per-face verdict from [`BRep::audit_winding_honesty`].
#[derive(Debug, Clone)]
pub struct FaceWindingHonesty {
    /// The audited face.
    pub face: FaceHandle,
    /// Signed UV area of the outer loop traversed in topological
    /// order (shoelace over the pcurve polygon, periodic-unwrapped).
    /// Positive = counter-clockwise in UV.
    pub signed_uv_area: f64,
    /// The stored orientation flag.
    pub orientation: bool,
    /// `Some(true)` when the flag agrees with the loop winding,
    /// `Some(false)` when it contradicts it, `None` when the loop
    /// could not be measured (see `reason`).
    pub honest: Option<bool>,
    /// Why the verdict is indeterminate; empty when determinate.
    pub reason: &'static str,
}

impl BRep {
    /// Measure winding-flag honesty for every face carrying surface
    /// geometry. Faces without a surface are skipped entirely;
    /// faces whose outer loop cannot be polygonized (missing
    /// pcurves, degenerate UV area) report an indeterminate verdict
    /// with the reason recorded. Uses each face's own pcurves at
    /// their stored tolerance — no global state.
    pub fn audit_winding_honesty(&self) -> Vec<FaceWindingHonesty> {
        let mut out = self.audit_winding_tier1();
        self.resolve_by_tangency_witness(&mut out);
        out
    }

    fn audit_winding_tier1(&self) -> Vec<FaceWindingHonesty> {
        let mut out = Vec::new();
        for (face, data) in &self.faces {
            let Some(sof) = &data.surface else {
                continue;
            };
            let verdict = match self.outer_loop_signed_uv_area(face, &sof.geometry) {
                Ok(area) => {
                    // Degeneracy floor: perimeter-relative, so a
                    // seam-spine loop (net-zero enclosed area at any
                    // chart scale — e.g. a full sphere wall) reads
                    // indeterminate instead of amplifying roundoff
                    // into a verdict. Honest loops sit near the
                    // isoperimetric scale (perimeter/4)^2, ~1e10x
                    // above this floor.
                    let degenerate = area.signed.abs() <= area.perimeter.powi(2) * 1.0e-12;
                    FaceWindingHonesty {
                        face,
                        signed_uv_area: area.signed,
                        orientation: sof.orientation,
                        honest: if degenerate {
                            None
                        } else {
                            Some((area.signed > 0.0) == sof.orientation)
                        },
                        reason: if degenerate { "degenerate-uv-area" } else { "" },
                    }
                }
                Err(reason) => FaceWindingHonesty {
                    face,
                    signed_uv_area: 0.0,
                    orientation: sof.orientation,
                    honest: None,
                    reason,
                },
            };
            out.push(verdict);
        }
        out
    }

    /// Second tier: resolve tier-1-indeterminate faces by tangency
    /// witness. Across a G1 trim the folded normals of one manifold
    /// shell's faces must agree, so a face whose winding cannot be
    /// measured (pole-crossing chart, seam-spine loop) inherits its
    /// verdict from the best high-|cos| neighbour that HAS a
    /// winding verdict: honest ⇔ (normals agree) == (neighbour
    /// honest). Transverse neighbours stay below the witness floor
    /// and cannot vote. Iterates to a fixpoint so chains of
    /// deferred faces resolve through their first anchored member.
    fn resolve_by_tangency_witness(&self, verdicts: &mut [FaceWindingHonesty]) {
        const WITNESS_FLOOR: f64 = 0.9;
        for _ in 0..verdicts.len() {
            let known: std::collections::BTreeMap<FaceHandle, bool> = verdicts
                .iter()
                .filter_map(|v| v.honest.map(|h| (v.face, h)))
                .collect();
            let mut changed = false;
            for v in verdicts.iter_mut() {
                if v.honest.is_some() {
                    continue;
                }
                let Some((dot, neighbour_honest)) = self.best_tangency_witness(v.face, &known)
                else {
                    continue;
                };
                if dot.abs() < WITNESS_FLOOR {
                    continue;
                }
                v.honest = Some((dot > 0.0) == neighbour_honest);
                v.reason = "tangency-witness";
                changed = true;
            }
            if !changed {
                break;
            }
        }
    }

    /// Best (max-|cos|) folded-normal agreement between `face` and
    /// any boundary neighbour with a known verdict, sampled at edge
    /// midpoints. Returns `(signed dot, neighbour's verdict)`.
    fn best_tangency_witness(
        &self,
        face: FaceHandle,
        known: &std::collections::BTreeMap<FaceHandle, bool>,
    ) -> Option<(f64, bool)> {
        let folded_at =
            |f: FaceHandle, p: &crate::geometry::Point3| -> Option<crate::geometry::Vec3> {
                let sof = self.faces.get(f)?.surface.as_ref()?;
                let eval = &sof.geometry as &dyn SurfaceEval;
                let (u, v) = eval.closest_point(p)?;
                let n = eval.normal_at(u, v).normalized()?;
                Some(n * if sof.orientation { 1.0 } else { -1.0 })
            };
        let outer = self.faces.get(face)?.outer_loop;
        let start = self.loops.get(outer)?.half_edge;
        let mut best: Option<(f64, bool)> = None;
        let mut cur = start;
        for _ in 0..self.half_edges.len() + 1 {
            let he = self.half_edges.get(cur)?;
            let twin = self.half_edges.get(he.twin)?;
            let neighbour = self.loops.get(twin.loop_)?.face;
            if neighbour != face
                && let Some(&neighbour_honest) = known.get(&neighbour)
            {
                // Sample ON the edge curve when one exists: a vertex
                // chord midpoint can sag off strongly curved trims
                // and closest_point then projects to the wrong sheet
                // (a small-tube torus chord lands past the inner
                // equator), inverting the witness normal.
                let mid = match self.edges.get(he.edge).and_then(|e| e.curve.as_ref()) {
                    Some(curve) => {
                        use crate::geometry::CurveEval;
                        curve.geometry.point_at((curve.t_start + curve.t_end) * 0.5)
                    }
                    None => {
                        let p0 = self.vertices.get(he.origin)?.position;
                        let p1 = self.vertices.get(twin.origin)?.position;
                        crate::geometry::Point3::new(
                            (p0.x + p1.x) * 0.5,
                            (p0.y + p1.y) * 0.5,
                            (p0.z + p1.z) * 0.5,
                        )
                    }
                };
                if let (Some(nf), Some(nn)) = (folded_at(face, &mid), folded_at(neighbour, &mid)) {
                    let dot = nf.dot(&nn);
                    if best.is_none_or(|(d, _)| dot.abs() > d.abs()) {
                        best = Some((dot, neighbour_honest));
                    }
                }
            }
            cur = he.next;
            if cur == start {
                break;
            }
        }
        best
    }

    /// Collecting validation form of the winding-honesty audit: one
    /// [`ValidationError::WindingDishonest`] per face whose stored
    /// flag contradicts its measured outer-loop winding (the armed
    /// adjacent-face-consistency half of check 7 — PK-245).
    /// Indeterminate verdicts (unmeasurable loops with no tangency
    /// witness) are not errors. Serving boundaries arm this via
    /// [`BRep::validate_all`]; mid-pipeline states may legitimately
    /// carry transient dishonesty that operator tails reconcile, so
    /// the fast [`BRep::validate`] path does not include it.
    pub fn validate_winding_honesty_all(&self) -> Vec<crate::error::ValidationError> {
        self.audit_winding_honesty()
            .into_iter()
            .filter(|v| v.honest == Some(false))
            .map(|v| crate::error::ValidationError::WindingDishonest {
                face: v.face,
                signed_uv_area: v.signed_uv_area,
                orientation: v.orientation,
            })
            .collect()
    }

    /// Faces whose stored orientation flag contradicts their outer
    /// loop's UV winding. Convenience filter over
    /// [`BRep::audit_winding_honesty`].
    pub fn dishonest_winding_faces(&self) -> Vec<FaceHandle> {
        self.audit_winding_honesty()
            .into_iter()
            .filter(|v| v.honest == Some(false))
            .map(|v| v.face)
            .collect()
    }

    /// Derive every measurable face's orientation flag from the
    /// audit verdict: the flag is redundant state (a function of
    /// loop winding and surface geometry), so at a serving boundary
    /// — after pcurves are finalized — the honest value is set
    /// directly rather than trusted from the author. Faces whose
    /// verdict is indeterminate (no confident winding or tangency
    /// evidence) keep their authored flag. Returns the number of
    /// flags flipped. Iterates because flipping a winding-dishonest
    /// face changes the tangency evidence available to its
    /// pole-degenerate neighbours.
    pub fn reconcile_orientation_flags(&mut self) -> usize {
        let mut flipped = 0;
        for _ in 0..4 {
            let dishonest = self.dishonest_winding_faces();
            if dishonest.is_empty() {
                break;
            }
            for face in dishonest {
                if let Some(sof) = self.faces.get_mut(face).and_then(|f| f.surface.as_mut()) {
                    sof.orientation = !sof.orientation;
                    flipped += 1;
                }
            }
        }
        flipped
    }

    fn outer_loop_signed_uv_area(
        &self,
        face: FaceHandle,
        surface: &SurfaceGeometry,
    ) -> Result<SignedUvArea, &'static str> {
        let outer = self.faces.get(face).ok_or("missing-face")?.outer_loop;
        self.loop_signed_uv_area(outer, surface)
    }

    /// Signed UV area + perimeter of one loop's pcurve polygon,
    /// traversed in topological order (periodic-unwrapped shoelace).
    /// Shared between the winding-honesty audit and the loop-winding
    /// repair postcondition of pcurve recompute.
    pub(crate) fn loop_signed_uv_area(
        &self,
        loop_handle: crate::topology::LoopHandle,
        surface: &SurfaceGeometry,
    ) -> Result<SignedUvArea, &'static str> {
        let start = self
            .loops
            .get(loop_handle)
            .ok_or("missing-outer-loop")?
            .half_edge;
        let (period_u, period_v) = surface.periodicity();
        let mut polygon: Vec<(f64, f64)> = Vec::new();
        let mut cur = start;
        for _ in 0..self.half_edges.len() + 1 {
            let samples = self.half_edge_uv_samples(cur, surface)?;
            for uv in samples {
                push_unwrapped(&mut polygon, uv, period_u, period_v);
            }
            cur = self.half_edges.get(cur).ok_or("missing-half-edge")?.next;
            if cur == start {
                break;
            }
        }
        if polygon.len() < 3 {
            return Err("loop-under-3-points");
        }
        // Chart-validity guards. A loop that encircles or touches a
        // parameterization pole has no meaningful UV winding: the
        // unwrapped polygon fails to close across a period (polar-cap
        // loops), or samples at the pole carry arbitrary u (measured:
        // a vertex-blend sphere cap whose corner vertex IS the chart
        // pole read a spurious −2.19π area). Such faces defer to the
        // tangency-witness tier.
        let (first, last) = (polygon[0], polygon[polygon.len() - 1]);
        let wraps =
            |gap: f64, period: Option<f64>| period.is_some_and(|p| p > 0.0 && gap.abs() > 0.25 * p);
        if wraps(first.0 - last.0, period_u) || wraps(first.1 - last.1, period_v) {
            return Err("pole-encircling-loop");
        }
        let mut normal_mags: Vec<f64> = polygon
            .iter()
            .map(|&(u, v)| surface.normal_at(u, v).length())
            .collect();
        normal_mags.sort_by(|a, b| a.total_cmp(b));
        let median = normal_mags[normal_mags.len() / 2];
        if median > 0.0 && normal_mags[0] < median * 1.0e-6 {
            return Err("pole-touching-loop");
        }
        let mut signed = 0.0;
        let mut perimeter = 0.0;
        for i in 0..polygon.len() {
            let (u0, v0) = polygon[i];
            let (u1, v1) = polygon[(i + 1) % polygon.len()];
            signed += (u0 * v1 - u1 * v0) * 0.5;
            perimeter += ((u1 - u0).powi(2) + (v1 - v0).powi(2)).sqrt();
        }
        if !signed.is_finite() {
            return Err("non-finite-uv-samples");
        }
        Ok(SignedUvArea { signed, perimeter })
    }

    /// UV samples along one half-edge's pcurve, oriented in the
    /// half-edge's own traversal direction. There is no repo-wide
    /// pcurve direction contract, so the direction is verified
    /// against the half-edge's vertex positions through the surface
    /// (the same witness `classify_convexity` uses).
    fn half_edge_uv_samples(
        &self,
        half_edge: HalfEdgeHandle,
        surface: &SurfaceGeometry,
    ) -> Result<Vec<(f64, f64)>, &'static str> {
        let he = self.half_edges.get(half_edge).ok_or("missing-half-edge")?;
        let pcurve = he.pcurve.as_ref().ok_or("missing-pcurve")?;
        let (t0, t1) = pcurve.parameter_range();
        if !(t0.is_finite() && t1.is_finite()) {
            return Err("unbounded-pcurve");
        }
        let origin = self
            .vertices
            .get(he.origin)
            .ok_or("missing-vertex")?
            .position;
        let dest_handle = self.half_edges.get(he.twin).ok_or("missing-twin")?.origin;
        let dest = self
            .vertices
            .get(dest_handle)
            .ok_or("missing-vertex")?
            .position;
        let uv_start = pcurve.point_at(t0);
        let q_start = surface.point_at(uv_start.x, uv_start.y);
        // Zero-length edges contribute nothing either way; sample
        // forward. Otherwise reverse when the pcurve's start maps
        // nearer the half-edge's destination than its origin.
        let reversed = q_start.distance_to(&dest) + 1.0e-12 < q_start.distance_to(&origin);
        let mut samples = Vec::with_capacity(SAMPLES_PER_HALF_EDGE);
        for k in 0..SAMPLES_PER_HALF_EDGE {
            let frac = k as f64 / (SAMPLES_PER_HALF_EDGE - 1) as f64;
            let frac = if reversed { 1.0 - frac } else { frac };
            let uv = pcurve.point_at(t0 + (t1 - t0) * frac);
            if !(uv.x.is_finite() && uv.y.is_finite()) {
                return Err("non-finite-uv-samples");
            }
            samples.push((uv.x, uv.y));
        }
        Ok(samples)
    }
}

/// Signed UV area and perimeter of a polygonized loop (crate-visible
/// for the pcurve loop-winding repair postcondition).
pub(crate) struct SignedUvArea {
    /// Shoelace area; positive = counter-clockwise in UV.
    pub(crate) signed: f64,
    /// Polygon perimeter, for the perimeter-relative degeneracy floor.
    pub(crate) perimeter: f64,
}

/// Append a UV sample, shifting by whole periods so consecutive
/// samples stay chart-continuous across periodic seams.
fn push_unwrapped(
    polygon: &mut Vec<(f64, f64)>,
    uv: (f64, f64),
    period_u: Option<f64>,
    period_v: Option<f64>,
) {
    let Some(&(prev_u, prev_v)) = polygon.last() else {
        polygon.push(uv);
        return;
    };
    let unwrap = |value: f64, prev: f64, period: Option<f64>| -> f64 {
        let Some(p) = period else {
            return value;
        };
        if p <= 0.0 {
            return value;
        }
        let shift = ((prev - value) / p).round();
        value + shift * p
    };
    polygon.push((
        unwrap(uv.0, prev_u, period_u),
        unwrap(uv.1, prev_v, period_v),
    ));
}

#[cfg(test)]
#[path = "winding_honesty_tests.rs"]
mod winding_honesty_tests;
