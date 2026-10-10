//! Edge-face intersection detection (diagnostic only).
//!
//! Split out of `intersect.rs` to keep the parent under the 500-line
//! structural cap. These hits are informational and not yet integrated
//! into the split pipeline.

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::{
    tolerance::Tolerance,
    topology::{BRep, FaceHandle},
};

use super::intersect::bounded_parameter_range;

/// A detected edge-face intersection point.
#[derive(Debug, Clone)]
#[allow(dead_code)] // Fields used for diagnostic inspection and future integration.
pub(crate) struct EdgeFaceHit {
    /// The edge handle from the source solid.
    pub edge: popkern_core::topology::EdgeHandle,
    /// The face handle from the target solid.
    pub face: FaceHandle,
    /// The 3D intersection point.
    pub point: popkern_core::geometry::Point3,
    /// Parameter t on the edge curve at the intersection.
    pub t_edge: f64,
    /// UV parameters (u, v) on the face surface at the intersection.
    pub uv: (f64, f64),
    /// True when `edge` belongs to solid A (and `face` to solid B);
    /// false for the reverse. Without this the operand-scoped handles
    /// are ambiguous to consumers holding both solids.
    pub edge_of_a: bool,
}

/// Detect edge-face intersections between two solids (diagnostic only).
///
/// For each edge of solid A, checks if it intersects any face of solid B,
/// and vice versa. Uses `intersect_surface_curve` with the edge's curve
/// geometry and the face's surface geometry.
///
/// Returns a list of [`EdgeFaceHit`] records. These are informational and
/// are not yet integrated into the split pipeline.
pub(crate) fn detect_edge_face_intersections(
    a: &BRep,
    b: &BRep,
    aabbs_b: &[(FaceHandle, Aabb)],
    aabbs_a: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> Vec<EdgeFaceHit> {
    let mut hits = Vec::new();

    // Edges of A vs faces of B.
    detect_edges_vs_faces(a, b, aabbs_b, tolerance, true, &mut hits);

    // Edges of B vs faces of A.
    detect_edges_vs_faces(b, a, aabbs_a, tolerance, false, &mut hits);

    if !hits.is_empty() && popkern_core::trace::trace_enabled("boolean") {
        eprintln!(
            "[intersect] detected {} edge-face intersection(s)",
            hits.len()
        );
    }

    hits
}

/// Check edges of `edge_brep` against faces of `face_brep`.
fn detect_edges_vs_faces(
    edge_brep: &BRep,
    face_brep: &BRep,
    face_aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
    edge_of_a: bool,
    hits: &mut Vec<EdgeFaceHit>,
) {
    use popkern_core::geometry::CurveEval;
    use popkern_core::geometry::intersection::SurfaceCurveIntersection;
    use popkern_core::geometry::intersection::surface_curve::intersect_surface_curve_dense;

    use popkern_core::geometry::curve::CurveGeometry;

    for (edge_h, edge) in &edge_brep.edges {
        // Get endpoint positions from vertices (always up-to-date, even if
        // the edge curve geometry is stale after vertex translation).
        let he = match edge_brep.half_edges.get(edge.half_edge) {
            Some(h) => h,
            None => continue,
        };
        let v_start = match edge_brep.vertices.get(he.origin) {
            Some(v) => v.position,
            None => continue,
        };
        let twin = match edge_brep.half_edges.get(he.twin) {
            Some(t) => t,
            None => continue,
        };
        let v_end = match edge_brep.vertices.get(twin.origin) {
            Some(v) => v.position,
            None => continue,
        };

        // CURVED edges use their real geometry: the vertex-chord
        // stand-in crosses faces the edge itself never touches (an
        // ellipse edge's chord hit two pocket walls at y=0; the graph
        // edge was split at points NOT on the ellipse and the arc
        // degraded to its chords, shredding the mesher's boundary).
        // STRAIGHT/missing-curve edges keep the vertex-position line,
        // which also covers stale line geometry.
        let stored = edge.curve.as_ref().filter(|c| {
            !matches!(c.geometry, CurveGeometry::Line { .. })
                && match &c.geometry {
                    // A degree-1 NURBS is straight ONLY as a 2-point
                    // segment: a multi-CP polyline is a curved edge
                    // spelled piecewise, and reducing it to its
                    // vertex chord puts the hit up to the full
                    // sagitta off the real edge (PK-295, measured: a
                    // torus strip's 17-point seam polyline chord-hit
                    // a cutter wall 0.073 from the true crossing, the
                    // split guard rightly rejected the off-curve
                    // point, and the pierce node never existed).
                    CurveGeometry::Nurbs(n) => n.degree > 1 || n.control_points.len() > 2,
                    _ => true,
                }
        });
        // The edge AABB must CONTAIN the edge — a vertex-only box
        // collapses to the SEAM POINT for a closed rim (both
        // vertices coincide), so the rim×face pair was pruned before
        // intersection ever ran (measured: a bore's bottom rim
        // crossing a breached wall produced no edge-face hits, the
        // rim never split at the wall, and the fragment walk
        // dead-ended at the mouth). Same conservative-extent rule as
        // the face AABBs.
        let mut extent_points = vec![v_start, v_end];
        if let Some(c) = stored {
            crate::intersect::conservative_curve_extent(c, &mut extent_points);
        }
        let edge_aabb = match Aabb::from_points(&extent_points) {
            Some(a) => a,
            None => continue,
        };
        let (edge_curve, t_range) = match stored {
            // The solver takes an ASCENDING window: a curve stored
            // with t_start > t_end (a reversed arc) errored "t_range
            // invalid" and every crossing of that edge was silently
            // dropped (measured, PK-750: a sketch rod's top half-rims).
            Some(c) => (
                c.geometry.clone(),
                (c.t_start.min(c.t_end), c.t_start.max(c.t_end)),
            ),
            None => {
                let direction = v_end - v_start;
                if direction.length_squared() < tolerance.numerical_floor().powi(2) {
                    continue;
                }
                (
                    CurveGeometry::Line {
                        origin: v_start,
                        direction,
                    },
                    (0.0, 1.0),
                )
            }
        };

        for (face_h, face_aabb) in face_aabbs {
            if !edge_aabb.overlaps(face_aabb, tolerance.linear()) {
                continue;
            }

            let face = &face_brep.faces[*face_h];
            let surface_on_face = match face.surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };

            let ((mut u_min, mut u_max), (v_min, v_max)) =
                bounded_parameter_range(face_brep, *face_h, surface_on_face);
            // CLOSED-WALL CAVEAT (same family as the trim-box caveat
            // below): a periodic wall bounded only by rims has 1–2
            // loop vertices, ALL on the seam, so the vertex-derived u
            // window collapses to ±padding around the seam azimuth
            // and prunes legitimate far-side roots (measured: the
            // cross-drilled boss — each seam-line × wall pair lost
            // the root opposite the seam, the pierce curves never
            // split at the missing nodes, and the graze filter
            // cascaded ALL record pieces off both walls). A closed
            // wall spans the full period: widen u to one full turn.
            let distinct_positions = {
                let ps = face_brep.face_outer_loop_positions(*face_h);
                let mut distinct: Vec<popkern_core::geometry::Point3> = Vec::new();
                for p in ps {
                    if distinct
                        .iter()
                        .all(|q| q.distance_to(&p) > tolerance.linear())
                    {
                        distinct.push(p);
                    }
                }
                distinct.len()
            };
            {
                use popkern_core::geometry::SurfaceEval as _;
                if let (Some(period), true) = (
                    surface_on_face.parameter_periodicity().0,
                    distinct_positions < 3,
                ) && u_max - u_min < period
                {
                    let mid = 0.5 * (u_min + u_max);
                    u_min = mid - 0.5 * period;
                    u_max = mid + 0.5 * period;
                }
            }

            // An edge lying IN the surface has no crossings, only a
            // coincidence (see `edge_lies_on_surface`).
            if crate::intersect_edge_face_coincident::edge_lies_on_surface(
                &edge_curve,
                t_range,
                surface_on_face,
                tolerance,
            ) {
                if popkern_core::trace::trace_enabled("edge-face-hits") {
                    eprintln!("  [hit-skip] edge {edge_h:?} lies on face {face_h:?}'s surface");
                }
                continue;
            }

            // COMPLETE crossing enumeration (PK-679): the tuned grid's
            // 48-per-axis cap gives a 111-turn thread wall cells ~2.3
            // turns wide whose 3x3-station boxes alias to a few
            // azimuths and do not contain their patch — the rod
            // cylinder's seam line found 56 of its 110 flank pierces,
            // the tip cone's seam pierce none, and the cone's record
            // crossed its seam without a node (ChartOverwrap 7.854).
            // The dense budget enumerates every crossing; the full
            // workspace suite (9622) is unchanged under it, so the
            // PK-461 note that extra hits re-wire the M5 seam-pinch
            // serve no longer holds.
            let result = intersect_surface_curve_dense(
                surface_on_face,
                (u_min, u_max),
                (v_min, v_max),
                &edge_curve,
                t_range,
                tolerance,
            );
            if popkern_core::trace::trace_enabled("edge-face-hits") {
                match &result {
                    Ok(r) => eprintln!(
                        "  [hit-pair] edge {edge_h:?} ({:.2},{:.2},{:.2})-({:.2},{:.2},{:.2}) x face {face_h:?} u[{u_min:.3},{u_max:.3}] v[{v_min:.3},{v_max:.3}] -> {} isects",
                        v_start.x,
                        v_start.y,
                        v_start.z,
                        v_end.x,
                        v_end.y,
                        v_end.z,
                        r.intersections.len()
                    ),
                    Err(e) => {
                        eprintln!("  [hit-pair] edge {edge_h:?} x face {face_h:?} -> error {e}")
                    }
                }
            }

            // Acceptance is trim-honest: the padded range above gives the
            // SOLVER headroom, but a hit landing in the padding is not on
            // the face (measured: a seat tool's wall plane crossed a deck
            // outline edge 4 mm outside the wall's bounds; the phantom hit
            // split the cellular graph edge there and the polluted loop
            // broke the NEXT boolean on the body). Bound acceptance by the
            // face's exact vertex-derived parameter box plus tolerance.
            // CLOSED-RIM CAVEAT: a face bounded by closed edges (a
            // cylinder cap's full-circle rim) has 1–2 loop vertices,
            // so its vertex-derived box collapses and would reject
            // every legitimate hit (measured: star-socket corner
            // hits on a cylinder cap all filtered, 427 open edges).
            // Fewer than 3 vertices → no verdict, keep the padded
            // acceptance.
            // Distinct-position guard (same reason as the widened u
            // window above): a closed wall's loop visits its 1–2 seam
            // vertices twice, so the RAW position count lies about
            // the trim information available.
            let face_positions = face_brep.face_outer_loop_positions(*face_h);
            let trim_box = if distinct_positions >= 3 {
                Some(
                    popkern_core::geometry::uv_projection::exact_bounded_parameter_range(
                        &face_positions,
                        surface_on_face,
                    ),
                )
            } else {
                None
            };
            let uv_margin = tolerance.linear();
            // A face that carries its own SEAM (an edge whose two
            // half-edges both bound it) spans the whole u period, but
            // its seam vertices project to u = 0 only, never 2π: the
            // vertex box stops at the last vertex before the seam and
            // rejects the strip beyond it. Measured (PK-746): a fin
            // straddling a housing's seam put a root corner at
            // u = 6.2427; a boss rod's seam edge pierced the housing
            // at u = 6.2785, the hit was rejected, its ring carried
            // no node at the rod's seam, the rod's chart walk fell
            // back to 3D-Newell, and the union served without the
            // boss wall. Such a face's u is not filtered.
            let u_period = {
                use popkern_core::geometry::SurfaceEval as _;
                surface_on_face.parameter_periodicity().0
            }
            .filter(|_| !face_has_own_seam(face_brep, *face_h));
            let full_turn = u_period.is_none() && {
                use popkern_core::geometry::SurfaceEval as _;
                surface_on_face.parameter_periodicity().0.is_some()
            };

            if let Ok(result) = result {
                for isect in &result.intersections {
                    if let SurfaceCurveIntersection::Point { t, u, v, .. } = *isect {
                        // The hit must lie STRICTLY INSIDE the edge's own
                        // parameter window: endpoint hits are vertex-face
                        // touches, and hits beyond the window live on the
                        // basis EXTENSION, not the edge (measured: the
                        // sphere seam arc's quadratic extension crossed a
                        // DISJOINT cube's face plane; the phantom node
                        // split the seam and the union carried a 5.0-gap
                        // open loop).
                        let eps = tolerance.linear();
                        let (t_lo, t_hi) = (t_range.0.min(t_range.1), t_range.0.max(t_range.1));
                        if t < t_lo + eps
                            || t > t_hi - eps
                            || crate::intersect_edge_face_coincident::touches_at_endpoint(
                                &edge_curve,
                                (t_lo, t_hi),
                                t,
                                surface_on_face,
                                tolerance,
                            )
                        {
                            if popkern_core::trace::trace_enabled("edge-face-hits") {
                                eprintln!(
                                    "  [hit-reject] endpoint t={t:.4} window [{t_lo:.4},{t_hi:.4}] edge {edge_h:?} face {face_h:?}"
                                );
                            }
                            continue;
                        }
                        // 3D honesty: the curve point must lie ON the
                        // face's surface, proven by EITHER evaluation
                        // route. The solver's claimed (u,v) can be coarse
                        // (measured 0.3 off on a plane-splitter's sheet
                        // face while the 3D hit itself was exact), so a
                        // projection check backs it up — but projection
                        // alone must not VETO: on a helical thread flank
                        // the cold projection basin converges on the
                        // ADJACENT COIL TURN, 0.875 away, while the
                        // solver uv evaluates exactly on the surface
                        // (measured: the tip-chamfer cone's seam edge
                        // crossing the groove flank — the rejected pierce
                        // node left the record chain seam-crossing and
                        // the chart walk declined the whole face,
                        // PK-463).
                        {
                            use popkern_core::geometry::SurfaceEval as _;
                            let p_curve = edge_curve.point_at(t);
                            let on_by = tolerance.linear() * 10.0;
                            let d_solver = surface_on_face.point_at(u, v).distance_to(&p_curve);
                            let off_surface = d_solver > on_by
                                && popkern_core::geometry::uv_projection::project_point_to_uv(
                                    surface_on_face,
                                    &p_curve,
                                )
                                .map(|(pu, pv)| {
                                    surface_on_face.point_at(pu, pv).distance_to(&p_curve) > on_by
                                })
                                .unwrap_or(false);
                            if off_surface {
                                if popkern_core::trace::trace_enabled("edge-face-hits") {
                                    eprintln!(
                                        "  [hit-reject] off-surface t={t:.3} uv=({u:.3},{v:.3}) d_solver={d_solver:.6}"
                                    );
                                }
                                continue;
                            }
                        }
                        // The box is in the principal chart; on a
                        // u-periodic surface the hit's u is tested
                        // modulo the period (measured, PK-681: the
                        // mouth-chamfer cone's seam edge pierced the
                        // nut's crest band at u = −0.0 against a box of
                        // [1.84, 2π] — the same azimuth as its upper
                        // edge — and the rejected pierce left the crest
                        // record seam-crossing without a node; the
                        // cone's chart walk declined ChartOverwrap).
                        let u_in_box = |lo: f64, hi: f64| {
                            let inside = |uu: f64| uu >= lo - uv_margin && uu <= hi + uv_margin;
                            match u_period {
                                Some(p) => {
                                    let k = ((lo - u) / p).ceil();
                                    inside(u + k * p) || inside(u + (k - 1.0) * p)
                                }
                                None => inside(u),
                            }
                        };
                        if let Some(((eu_min, eu_max), (ev_min, ev_max))) = trim_box
                            && ((!full_turn && !u_in_box(eu_min, eu_max))
                                || v < ev_min - uv_margin
                                || v > ev_max + uv_margin)
                        {
                            if popkern_core::trace::trace_enabled("edge-face-hits") {
                                eprintln!(
                                    "  [hit-reject] trim-box uv=({u:.4},{v:.4}) box u[{eu_min:.4},{eu_max:.4}] v[{ev_min:.4},{ev_max:.4}] edge {edge_h:?} face {face_h:?}"
                                );
                            }
                            continue;
                        }
                        // TRIPLE-POINT REFINEMENT (PK-480): the hit
                        // evaluated on the FITTED edge at the 1D
                        // solve's parameter slides along a
                        // near-tangent edge by fit-slop/sin(crossing
                        // angle) (see `clip_triple_junction`); the
                        // records meeting this junction end on the
                        // refined three-surface point, and the edge
                        // split must land on the SAME point or the
                        // duplicate node dead-ends the partner chart.
                        let raw_point = edge_curve.point_at(t);
                        let point = crate::clip_triple_junction::refine_to_triple_point(
                            edge_brep,
                            edge_h,
                            surface_on_face,
                            raw_point,
                            tolerance,
                        );
                        let (t_edge, uv) = if point == raw_point {
                            (t, (u, v))
                        } else {
                            // Keep parameter and chart data consistent
                            // with the refined position.
                            let t_ref = (&edge_curve as &dyn CurveEval)
                                .closest_point(&point)
                                .map(|tt| tt.clamp(t_lo, t_hi))
                                .unwrap_or(t);
                            let uv_ref =
                                popkern_core::geometry::uv_projection::project_point_to_uv_near(
                                    surface_on_face,
                                    &point,
                                    u,
                                    v,
                                )
                                .unwrap_or((u, v));
                            (t_ref, uv_ref)
                        };
                        hits.push(EdgeFaceHit {
                            edge: edge_h,
                            face: *face_h,
                            point,
                            t_edge,
                            uv,
                            edge_of_a,
                        });
                    }
                }
            }
        }
    }
}

/// Whether `face` carries its own seam: an edge both of whose
/// half-edges lie on loops of this face.
fn face_has_own_seam(brep: &BRep, face: FaceHandle) -> bool {
    let Some(f) = brep.faces.get(face) else {
        return false;
    };
    std::iter::once(f.outer_loop)
        .chain(f.inner_loops.iter().copied())
        .any(|lh| {
            let Some(start) = brep.loops.get(lh).map(|l| l.half_edge) else {
                return false;
            };
            let mut h = start;
            for _ in 0..=brep.half_edges.len() {
                let Some(he) = brep.half_edges.get(h) else {
                    return false;
                };
                let twin_on_face = brep
                    .half_edges
                    .get(he.twin)
                    .and_then(|tw| brep.loops.get(tw.loop_))
                    .is_some_and(|l| l.face == face);
                if twin_on_face && he.twin != h {
                    return true;
                }
                h = he.next;
                if h == start {
                    break;
                }
            }
            false
        })
}
