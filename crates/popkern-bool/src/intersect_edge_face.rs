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
    detect_edges_vs_faces(a, b, aabbs_b, tolerance, &mut hits);

    // Edges of B vs faces of A.
    detect_edges_vs_faces(b, a, aabbs_a, tolerance, &mut hits);

    if !hits.is_empty() {
        #[cfg(debug_assertions)]
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
    hits: &mut Vec<EdgeFaceHit>,
) {
    use popkern_core::geometry::CurveEval;
    use popkern_core::geometry::intersection::{SurfaceCurveIntersection, intersect_surface_curve};

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
                    CurveGeometry::Nurbs(n) => n.degree > 1,
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
            Some(c) => (c.geometry.clone(), (c.t_start, c.t_end)),
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

            let result = intersect_surface_curve(
                surface_on_face,
                (u_min, u_max),
                (v_min, v_max),
                &edge_curve,
                t_range,
                tolerance,
            );
            if std::env::var_os("PK_HIT_TRACE").is_some()
                && let Ok(r) = &result
                && !r.intersections.is_empty()
            {
                eprintln!(
                    "  [hit-pair] face {face_h:?} u[{u_min:.3},{u_max:.3}] v[{v_min:.3},{v_max:.3}] -> {} isects",
                    r.intersections.len()
                );
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
                        if t < t_lo + eps || t > t_hi - eps {
                            continue;
                        }
                        // 3D honesty: the curve point must lie ON the
                        // face's surface. Judged by PROJECTION, not by the
                        // solver's claimed (u,v) — the solver uv is coarse
                        // (measured 0.3 off on a plane-splitter's sheet
                        // face while the 3D hit itself was exact).
                        {
                            use popkern_core::geometry::SurfaceEval as _;
                            let p_curve = edge_curve.point_at(t);
                            let off_surface =
                                popkern_core::geometry::uv_projection::project_point_to_uv(
                                    surface_on_face,
                                    &p_curve,
                                )
                                .map(|(pu, pv)| {
                                    surface_on_face.point_at(pu, pv).distance_to(&p_curve)
                                        > tolerance.linear() * 10.0
                                })
                                .unwrap_or(false);
                            if off_surface {
                                if std::env::var_os("PK_HIT_TRACE").is_some() {
                                    eprintln!(
                                        "  [hit-reject] off-surface t={t:.3} uv=({u:.3},{v:.3})"
                                    );
                                }
                                continue;
                            }
                        }
                        if let Some(((eu_min, eu_max), (ev_min, ev_max))) = trim_box
                            && (u < eu_min - uv_margin
                                || u > eu_max + uv_margin
                                || v < ev_min - uv_margin
                                || v > ev_max + uv_margin)
                        {
                            continue;
                        }
                        hits.push(EdgeFaceHit {
                            edge: edge_h,
                            face: *face_h,
                            point: edge_curve.point_at(t),
                            t_edge: t,
                            uv: (u, v),
                        });
                    }
                }
            }
        }
    }
}
