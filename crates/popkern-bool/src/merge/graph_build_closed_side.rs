//! Closed-edge side-fidelity witness for operand insertion (split
//! from `graph_build.rs` under the 500-line file budget). See
//! [`closed_he_runs_curve_forward`].

use std::collections::HashMap;

use popkern_core::topology::{BRep, FaceHandle as FHandle, HalfEdgeHandle};

/// Per-operand context for the closed-edge side witness: lazily
/// built face AABBs (for ray parity) and a memo of measured outward
/// normal signs per face.
pub(super) struct ClosedSideContext<'a> {
    brep: &'a BRep,
    face_aabbs: Option<
        Vec<(
            FHandle,
            popkern_core::geometry::intersection::surface_subdivision::Aabb,
        )>,
    >,
    outward_sign: HashMap<FHandle, Option<f64>>,
}

impl<'a> ClosedSideContext<'a> {
    pub(super) fn new(brep: &'a BRep) -> Self {
        Self {
            brep,
            face_aabbs: None,
            outward_sign: HashMap::new(),
        }
    }

    /// Sign to fold the face's geometric normal into the OUTWARD
    /// material normal, measured by ray parity (a two-sided probe at
    /// the face's boundary-centroid projection: exactly one side must
    /// read Inside). The stored orientation flag is deliberately not
    /// consulted — measured on the angled tap: the blind-hole floor
    /// and the thread tool's top cap both carried flag-folded normals
    /// pointing INTO their material, and any flag-trusting witness
    /// voted them backwards with full confidence.
    fn outward_sign(&mut self, face: FHandle) -> Option<f64> {
        if let Some(&v) = self.outward_sign.get(&face) {
            return v;
        }
        let v = self.measure_outward_sign(face);
        self.outward_sign.insert(face, v);
        v
    }

    fn measure_outward_sign(&mut self, face: FHandle) -> Option<f64> {
        use popkern_core::geometry::SurfaceEval as _;
        let brep = self.brep;
        let sof = brep.faces.get(face)?.surface.as_ref()?;
        let c = loop_boundary_centroid(brep, brep.faces.get(face)?.outer_loop)?;
        let (u, v) = sof.geometry.closest_point(&c)?;
        let q = sof.geometry.point_at(u, v);
        let n = sof.geometry.normal_at(u, v).normalized()?;
        if self.face_aabbs.is_none() {
            self.face_aabbs = Some(crate::intersect::compute_face_aabbs(brep));
        }
        let aabbs = self.face_aabbs.as_ref()?;
        let diag = aabbs
            .iter()
            .find(|(f, _)| *f == face)
            .map(|(_, a)| a.diagonal())
            .unwrap_or(0.0);
        let eps = (diag * 1.0e-3).max(brep.tolerance.linear() * 10.0);
        let classify = |p: popkern_core::geometry::Point3| {
            crate::classify::point_in_solid(&p, brep, aabbs, &brep.tolerance).ok()
        };
        use crate::error::PointClassification::{Inside, Outside};
        let plus = classify(q + n * eps)?;
        let minus = classify(q + n * (-eps))?;
        match (plus, minus) {
            (Outside, Inside) => Some(1.0),
            (Inside, Outside) => Some(-1.0),
            // Both sides void (the centroid projected into a hole of
            // the face) or both material or on-boundary: no witness.
            _ => None,
        }
    }
}

/// Whether a CLOSED edge's canonical half-edge traverses the edge's
/// stored curve in its increasing-t direction, by the material-left
/// rule: at the curve midpoint, a face traversing the rim has its
/// material on the LEFT of the traversal (outward normal × tangent
/// points toward the face interior). The outward normal is measured
/// by ray parity, never read from the orientation flag; the interior
/// is proxied loop-locally (toward the half-edge's own loop centroid
/// for an outer loop, away from it for an inner loop — a whole-face
/// centroid crosses the hole at an inner rim and votes backwards).
/// Manifold twins traverse a shared edge oppositely, so one decisive
/// vote suffices; when both faces vote they must agree or there is
/// no verdict. `None` keeps the historical canonical→forward
/// mapping. Pcurves are never consulted: a closed edge's traversal
/// cannot be read from its pcurve direction (legacy producers store
/// edge-aligned pcurves — measured: a pcurve chain witness swapped
/// the chamfered rod's tip rims and stranded its thread cut).
pub(super) fn closed_he_runs_curve_forward(
    ctx: &mut ClosedSideContext,
    he_h: HalfEdgeHandle,
) -> Option<bool> {
    use popkern_core::geometry::curve::CurveEval;
    let brep = ctx.brep;
    let he = brep.half_edges.get(he_h)?;
    let twin = brep.half_edges.get(he.twin)?;
    let curve = brep.edges.get(he.edge)?.curve.as_ref()?;
    let tm = 0.5 * (curve.t_start + curve.t_end);
    let p = curve.geometry.point_at(tm);
    let tangent = curve.geometry.derivative_at(tm, 1).normalized()?;
    // Signed material-left cosine for the half-edge's face under the
    // hypothesis that the face traverses the curve in +tangent
    // direction; `None` when unmeasurable.
    let mut vote = |heh: HalfEdgeHandle| -> Option<f64> {
        use popkern_core::geometry::SurfaceEval as _;
        let lh = ctx.brep.half_edges.get(heh)?.loop_;
        let face = ctx.brep.loops.get(lh)?.face;
        let fd = ctx.brep.faces.get(face)?;
        let sof = fd.surface.as_ref()?;
        let sign = ctx.outward_sign(face)?;
        let (u, v) = sof.geometry.closest_point(&p)?;
        let n = sof.geometry.normal_at(u, v).normalized()? * sign;
        let left = n.cross(&tangent);
        let c = loop_boundary_centroid(ctx.brep, lh)?;
        let inner = fd.inner_loops.contains(&lh);
        let d = if inner { p - c } else { c - p };
        let scale = left.length() * d.length();
        if scale <= 0.0 {
            return None;
        }
        let cosine = left.dot(&d) / scale;
        if popkern_core::trace::trace_enabled("closed-side") {
            eprintln!(
                "      [closed-side-face] {face:?} outward-sign {sign:+.0} inner {inner} cos {cosine:.3}"
            );
        }
        if cosine.abs() < 0.1 {
            None
        } else {
            Some(cosine)
        }
    };
    let v_own = vote(he_h);
    let v_twin = vote(he.twin);
    let _ = twin;
    // BOTH faces must vote, in opposition. A single vote is not
    // trusted: operands can carry non-manifold closed edges (the
    // chamfered rod's tip rim reaches graph-merge precisely because
    // cellular declined it as non-manifold), and there the pairing
    // argument that would let one face speak for its twin does not
    // hold (measured: a lone weak vote swapped that rim and
    // stranded the thread cut).
    match (v_own, v_twin) {
        (Some(a), Some(b)) if (a > 0.0) != (b > 0.0) => {
            if a > 0.0 {
                Some(true)
            } else {
                // The material-left truth says the canonical side
                // traverses BACKWARD. Re-registering is only
                // warranted when the stored pcurve LIES about that
                // (claims forward): a dishonest pcurve makes the
                // chart embedding fall back to the pair's
                // curve-window direction, where the side mapping is
                // the only direction authority (measured: the angled
                // tap's bore-wall rim). With an HONEST pcurve the
                // embedding already runs the loop's true direction
                // and the historical mapping is load-bearing
                // convention — swapping it broke the closed-rim
                // shell fixture's trim welds (measured: a rim edge
                // lost its opposing weld partner and the offset
                // declined).
                let own_len = loop_half_edge_count(brep, he_h);
                let twin_len = loop_half_edge_count(brep, he.twin);
                if popkern_core::trace::trace_enabled("closed-side") {
                    eprintln!(
                        "      [closed-side-loops] own {own_len} twin {twin_len} pc-claim {:?}",
                        pcurve_claims_forward(brep, he_h)
                    );
                }
                match pcurve_claims_forward(brep, he_h) {
                    Some(true) if own_len == 1 && twin_len > 1 => Some(false),
                    _ => None,
                }
            }
        }
        _ => None,
    }
}

/// Number of half-edges in the loop containing `he_h`.
fn loop_half_edge_count(brep: &BRep, he_h: HalfEdgeHandle) -> usize {
    let Some(he) = brep.half_edges.get(he_h) else {
        return 0;
    };
    let mut n = 1;
    let mut cur = he.next;
    while cur != he_h && n <= brep.half_edges.len() {
        let Some(x) = brep.half_edges.get(cur) else {
            break;
        };
        cur = x.next;
        n += 1;
    }
    n
}

/// The stored pcurve's own direction claim: `Some(true)` when
/// traversing it t_start→t_end moves ALONG the edge curve's
/// increasing-t direction (mapped through the half-edge's face
/// surface). `None` when unmeasurable.
fn pcurve_claims_forward(brep: &BRep, he_h: HalfEdgeHandle) -> Option<bool> {
    use popkern_core::geometry::SurfaceEval as _;
    use popkern_core::geometry::curve::CurveEval as _;
    let he = brep.half_edges.get(he_h)?;
    let curve = brep.edges.get(he.edge)?.curve.as_ref()?;
    let pc = he.pcurve.as_ref()?;
    let face = brep.loops.get(he.loop_)?.face;
    let surface = brep.faces.get(face)?.surface.as_ref()?;
    let (t0, t1) = pc.parameter_range();
    if !(t0.is_finite() && t1.is_finite()) {
        return None;
    }
    let tm = t0 + (t1 - t0) * 0.5;
    let te = t0 + (t1 - t0) * 0.501;
    let a = pc.point_at(tm);
    let b = pc.point_at(te);
    let p1 = surface.geometry.point_at(a.x, a.y);
    let p2 = surface.geometry.point_at(b.x, b.y);
    let step = p2 - p1;
    if step.length_squared() <= 0.0 {
        return None;
    }
    let (tc, dist) = super::super::graph_split::nearest_param_on_curve(curve, &p1);
    if dist > step.length() + 1.0e-3 {
        return None;
    }
    let tangent = curve.geometry.derivative_at(tc, 1);
    let dot = step.dot(&tangent);
    if dot == 0.0 { None } else { Some(dot > 0.0) }
}

/// Mean of 3D samples over one loop's boundary (edge-curve samples,
/// which are direction-agnostic; vertex positions where an edge has
/// no curve) — the loop-local interior proxy for the material-left
/// witness.
fn loop_boundary_centroid(
    brep: &BRep,
    lh: popkern_core::topology::LoopHandle,
) -> Option<popkern_core::geometry::Point3> {
    use popkern_core::geometry::curve::CurveEval;
    let (mut sx, mut sy, mut sz, mut n) = (0.0, 0.0, 0.0, 0usize);
    let mut add = |p: popkern_core::geometry::Point3| {
        sx += p.x;
        sy += p.y;
        sz += p.z;
        n += 1;
    };
    let start = brep.loops.get(lh)?.half_edge;
    let mut cur = start;
    for _ in 0..brep.half_edges.len() + 1 {
        let he = brep.half_edges.get(cur)?;
        match brep.edges.get(he.edge).and_then(|e| e.curve.as_ref()) {
            Some(c) => {
                for k in 0..4 {
                    let t = c.t_start + (c.t_end - c.t_start) * (k as f64 + 0.5) / 4.0;
                    add(c.geometry.point_at(t));
                }
            }
            None => add(brep.vertices.get(he.origin)?.position),
        }
        cur = he.next;
        if cur == start {
            break;
        }
    }
    if n == 0 {
        return None;
    }
    Some(popkern_core::geometry::Point3::new(
        sx / n as f64,
        sy / n as f64,
        sz / n as f64,
    ))
}
