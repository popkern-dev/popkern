//! Edge-curve and pcurve attachment for the cellular assembly —
//! split from `cellular/assemble.rs` under the structural line budget.

use popkern_core::{
    constructors::builder::BRepBuilder,
    geometry::Point3,
    geometry::curve::{CurveGeometry, CurveOnEdge},
    topology::EdgeHandle,
};

use crate::merge::graph::IntersectionGraph;

/// Attach curve geometry to a freshly created edge.
///
/// For shared (intersection) edges and self-loops: the graph's actual
/// intersection curve, aligned to the canonical edge direction (na→nb)
/// — both incident faces share the exact same [`CurveOnEdge`],
/// enforcing the shared-geometry principle (spec §4.5, §9.2) for
/// watertightness. For original/split edges: adopt the stored curve
/// when its endpoints coincide with this edge's nodes — an untouched
/// curved boundary edge synthesized as a line degraded arcs to their
/// chords (a cylinder wall's bottom arc and the wedge ellipse both
/// flattened; the chord then generated phantom edge-face hits and the
/// mesher's boundary-first pass bailed on the off-surface curve).
/// Fragments of a SPLIT curved edge whose stored curve spans the whole
/// original go through subrange + end-heal; otherwise a line from the
/// endpoints.
pub(super) fn attach_edge_curve(
    builder: &mut BRepBuilder,
    graph: &IntersectionGraph,
    sample_ghe: usize,
    edge_h: EdgeHandle,
    na: usize,
    nb: usize,
) {
    let is_shared = matches!(
        graph.half_edges[sample_ghe].source,
        crate::merge::graph::EdgeSource::Shared { .. }
    );
    if is_shared || na == nb {
        let (curve_he, curve_ref) = if graph.half_edges[sample_ghe].curve.is_some() {
            (sample_ghe, graph.half_edges[sample_ghe].curve.as_ref())
        } else {
            let twin = graph.half_edges[sample_ghe].twin;
            (twin, graph.half_edges[twin].curve.as_ref())
        };

        if let Some(curve) = curve_ref {
            let curve_origin = graph.half_edges[curve_he].origin;
            let mut aligned = curve.clone();
            if curve_origin != na && na != nb {
                // Curve goes nb→na but we need na→nb: swap parameters.
                std::mem::swap(&mut aligned.t_start, &mut aligned.t_end);
            }
            if popkern_core::trace::trace_enabled("boolean") {
                eprintln!(
                    "  [assemble] set curve on key {sample_ghe}: t [{:.3},{:.3}]",
                    aligned.t_start, aligned.t_end
                );
            }
            builder.set_edge_curve(edge_h, aligned);
        } else if popkern_core::trace::trace_enabled("boolean") {
            eprintln!(
                "  [assemble] shared/self-loop edge key {sample_ghe} ({na},{nb}): NO CURVE on either half"
            );
        }
    } else {
        let p_start = graph.nodes[na].position;
        let p_end = graph.nodes[nb].position;
        let stored = graph.half_edges[sample_ghe].curve.as_ref().or_else(|| {
            graph.half_edges[graph.half_edges[sample_ghe].twin]
                .curve
                .as_ref()
        });
        let adopted = stored.and_then(|c| {
            use popkern_core::geometry::curve::CurveEval;
            let c0 = c.geometry.point_at(c.t_start);
            let c1 = c.geometry.point_at(c.t_end);
            let tol = 1.0e-6;
            if c0.distance_to(&p_start) < tol && c1.distance_to(&p_end) < tol {
                Some(c.clone())
            } else if c0.distance_to(&p_end) < tol && c1.distance_to(&p_start) < tol {
                let mut rev = c.clone();
                std::mem::swap(&mut rev.t_start, &mut rev.t_end);
                Some(rev)
            } else {
                crate::merge::assemble::subrange_healed_curve(c, &p_start, &p_end)
            }
        });
        if let Some(curve) = adopted {
            builder.set_edge_curve(edge_h, curve);
        } else {
            let dir = p_end - p_start;
            let len = dir.length();
            if len > 1e-12 {
                let line = CurveOnEdge {
                    geometry: CurveGeometry::Line {
                        origin: p_start,
                        direction: dir * (1.0 / len),
                    },
                    t_start: 0.0,
                    t_end: len,
                };
                builder.set_edge_curve(edge_h, line);
            }
        }
    }
}

/// The stored pcurve for traversal direction `ghe` in `face_ref`'s
/// chart — the same per-face UV-space correction graph-merge's
/// assembly applies (a Shared pair stores face_a's chart on the
/// forward half and face_b's on the twin; a fragment walking the
/// "wrong" direction must borrow the twin's chart REVERSED).
///
/// The cellular assembly historically wired NO pcurves at all and
/// left every border to the projection recompute, which replaces
/// exact authored charts with fitted approximations (measured on the
/// spring end grind: the wall's degree-1 seam/rim border lines came
/// back as 23-control-point fits, and the NEXT boolean's UV clip
/// counted a phantom second crossing against the wiggle, collapsed
/// the 5.6e-4 corner-connector trim to a point, and the face split
/// failed the whole subtraction).
pub(super) fn pcurve_for_face(
    graph: &IntersectionGraph,
    ghe: usize,
    face_ref: &crate::merge::graph::FaceRef,
) -> Option<popkern_core::geometry::pcurve::Pcurve> {
    let twin = graph.half_edges[ghe].twin;
    if graph.half_edges[ghe].face.as_ref() == Some(face_ref) {
        return graph.half_edges[ghe].pcurve.clone();
    }
    if graph.half_edges[twin].face.as_ref() == Some(face_ref) {
        return graph.half_edges[twin].pcurve.as_ref().map(|p| p.reversed());
    }
    graph.half_edges[ghe].pcurve.clone()
}

/// Orient `pc` to the BRep half-edge pcurve convention: the edge's
/// canonical half-edge traces the edge curve FORWARD, its twin
/// backward. Judged by sampled comparison — the graph carries at
/// least three pcurve direction conventions (record-forward on BOTH
/// halves at Shared insertion, traversal-direction on walker
/// re-authored loops, per-side image cuts on split pieces), so no
/// ownership rule can decide direction; the geometry can (measured:
/// the ownership-rule install shipped both hole-rim pcurves of
/// cube−cylinder direction-reversed, 0.6 stale on the PK-287 gate).
/// A half-edge whose edge has no curve, or whose face surface is
/// unavailable, keeps the pcurve as handed.
pub(super) fn orient_pcurve_to_half_edge(
    brep: &popkern_core::topology::BRep,
    he: popkern_core::topology::HalfEdgeHandle,
    pc: popkern_core::geometry::pcurve::Pcurve,
    surface: Option<&popkern_core::geometry::surface::SurfaceGeometry>,
) -> popkern_core::geometry::pcurve::Pcurve {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;
    let Some(surface) = surface else { return pc };
    let Some(he_data) = brep.half_edges.get(he) else {
        return pc;
    };
    let Some(edge) = brep.edges.get(he_data.edge) else {
        return pc;
    };
    let Some(curve) = edge.curve.as_ref() else {
        return pc;
    };
    let forward = edge.half_edge == he;
    let (pt0, pt1) = pc.parameter_range();
    let mut d_fwd = 0.0;
    let mut d_rev = 0.0;
    for frac in [0.25, 0.5, 0.75] {
        let t = if forward {
            curve.t_start + (curve.t_end - curve.t_start) * frac
        } else {
            curve.t_end + (curve.t_start - curve.t_end) * frac
        };
        let p3d = curve.geometry.point_at(t);
        let uv_f = pc.point_at(pt0 + (pt1 - pt0) * frac);
        let uv_r = pc.point_at(pt1 + (pt0 - pt1) * frac);
        d_fwd += surface.point_at(uv_f.x, uv_f.y).distance_to(&p3d);
        d_rev += surface.point_at(uv_r.x, uv_r.y).distance_to(&p3d);
    }
    if d_rev < d_fwd { pc.reversed() } else { pc }
}

/// Midpoint-of-curve sample used by the span-unification registry.
pub(super) fn span_sample_points(
    graph: &IntersectionGraph,
    key: usize,
) -> ((usize, usize), Point3, Point3, Point3) {
    use popkern_core::geometry::curve::CurveEval;
    let na = graph.half_edges[key].origin;
    let nb = graph.half_edges[graph.half_edges[key].twin].origin;
    let (pa, pb) = (graph.nodes[na].position, graph.nodes[nb].position);
    let mid = graph.half_edges[key]
        .curve
        .as_ref()
        .or_else(|| graph.half_edges[graph.half_edges[key].twin].curve.as_ref())
        .map(|c| c.geometry.point_at(0.5 * (c.t_start + c.t_end)))
        .unwrap_or_else(|| {
            Point3::new(
                0.5 * (pa.x + pb.x),
                0.5 * (pa.y + pb.y),
                0.5 * (pa.z + pb.z),
            )
        });
    ((na, nb), pa, pb, mid)
}
