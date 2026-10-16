//! Post-fillet edge reconciliation.
//!
//! After topology operations and surface assignment, edges between
//! curved faces (cylinder, sphere, NURBS) and adjacent faces can
//! end up with 3-D curves that don't lie on both adjacent surfaces.
//! `multi_chord` and corner-blend `mef` calls write straight-line
//! curves that were correct *at the time* (the adjacent face was
//! still planar) but became inconsistent once the face surface was
//! later replaced by a cylindrical/curved fit.
//!
//! Symptom in OCCT-strict slicers (Bambu, FreeCAD): cylindrical
//! fillet ribbons show up with sharp boundaries because the
//! slicer detects the edge geometry doesn't lie on the cylinder
//! surface and refuses to render the face.
//!
//! This pass walks every edge, computes the analytic surface ×
//! surface intersection of its two adjacent faces, and replaces
//! the edge's 3-D curve with the intersection curve trimmed to
//! the existing vertex endpoints — when the trim succeeds and
//! the new curve hits both endpoints within tolerance.
//!
//! ## Scope
//!
//! Currently handles **cylinder × plane** only. Other combinations
//! (cylinder × cylinder, cylinder × NURBS, sphere × plane, etc.)
//! are out of scope for this first pass; those edges keep their
//! existing curves. Pcurves are also untouched here — they remain
//! whatever was assigned at construction time. A more complete
//! reconciliation pass would update both pcurves and possibly
//! topology (kev internal edges); that's tracked as future work.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::geometry::edge_triple::{EdgeAnchors, EdgeTripleOptions, compute_edge_geometry};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle};
use popkern_math::PortableMath;

#[path = "edge_reconcile_fit.rs"]
mod edge_reconcile_fit;
#[path = "edge_reconcile_snap.rs"]
mod edge_reconcile_snap;
use edge_reconcile_fit::try_curved_fit;
pub(crate) use edge_reconcile_snap::snap_vertex_and_adjacent_chords;

/// Walk every edge in `brep`. For each edge bordered on one side
/// by a cylindrical face and on the other by a planar face,
/// replace the edge's 3-D curve with the analytic cylinder × plane
/// intersection curve, trimmed to the existing vertex endpoints.
/// Edges where the intersection round-trip doesn't match the
/// vertex positions within tolerance are left alone.
pub(crate) fn reconcile_edge_curves(brep: &mut BRep) {
    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    let tolerance = brep.tolerance;
    for edge in edge_handles {
        try_reconcile_cyl_plane_edge(brep, edge, &tolerance);
    }
}

fn try_reconcile_cyl_plane_edge(brep: &mut BRep, edge: EdgeHandle, tolerance: &Tolerance) -> bool {
    let Some(edge_data) = brep.edges.get(edge) else {
        return false;
    };
    let he = edge_data.half_edge;
    let Some(he_data) = brep.half_edges.get(he) else {
        return false;
    };
    let twin = he_data.twin;
    let Some(twin_data) = brep.half_edges.get(twin) else {
        return false;
    };

    let Some(face_1) = brep
        .loops
        .get(he_data.loop_)
        .map(|loop_data| loop_data.face)
    else {
        return false;
    };
    let Some(face_2) = brep
        .loops
        .get(twin_data.loop_)
        .map(|loop_data| loop_data.face)
    else {
        return false;
    };
    if face_1 == face_2 {
        // Self-loop (seam/slit) edge: no pair intersection exists,
        // but the curve's honest home is ON the one surface — refit
        // by direct projection when it deviates (measured: a strip-
        // internal slit at the boss top seam crossing sat 2.5e-3
        // off its own strip, masked behind the first-error report).
        return super::seam_self_fit::try_same_face_seam_fit(brep, edge);
    }

    let Some(sof_1) = brep.faces.get(face_1).and_then(|f| f.surface.as_ref()) else {
        return false;
    };
    let Some(sof_2) = brep.faces.get(face_2).and_then(|f| f.surface.as_ref()) else {
        return false;
    };

    // Reconciliation is for STALE curves (chord lines written
    // against a face whose surface was later replaced). An edge
    // whose existing curve already lies on BOTH surfaces is
    // authoritative — rewriting it re-trims through the analytic
    // intersection, and on periodic intersections (cylinder ×
    // plane ellipses) the endpoint pair does not distinguish the
    // short arc from the complement: PK-097 saw runoff edges
    // replaced by the full ellipse the long way around.
    if let Some(curve) = edge_data.curve.as_ref() {
        // Strict: only curves EXACT on both surfaces qualify (the
        // runoff bisection lands ~1e-9; straight chord lines sag off
        // the cylinder by far more even on short edges and must
        // still be rewritten, pcurves included).
        let gate = 1.0e-8;
        let _ = tolerance;
        let on_surface = |sfc: &SurfaceGeometry, p: &popkern_core::geometry::Point3| -> bool {
            match sfc {
                SurfaceGeometry::Plane { origin, normal, .. } => {
                    let n = normal.length();
                    n > 1e-300 && ((*p - *origin).dot(normal) / n).abs() < gate
                }
                SurfaceGeometry::Cylinder {
                    origin,
                    axis,
                    radius,
                } => {
                    let d = *p - *origin;
                    let along = d.dot(axis);
                    let radial = d - *axis * along;
                    (radial.length() - radius).abs() < gate
                }
                SurfaceGeometry::Cone {
                    apex,
                    axis,
                    half_angle,
                } => {
                    let n = axis.length();
                    if n < 1e-300 {
                        return false;
                    }
                    let an = *axis * (1.0 / n);
                    let d = *p - *apex;
                    let station = d.dot(&an);
                    let radial = (d - an * station).length();
                    (radial - station.abs() * half_angle.ptan().abs()).abs() < gate
                }
                _ => false,
            }
        };
        let mut all_on = true;
        for i in 0..=4 {
            let t = curve.t_start + (curve.t_end - curve.t_start) * i as f64 / 4.0;
            let p = curve.geometry.point_at(t);
            if !(on_surface(&sof_1.geometry, &p) && on_surface(&sof_2.geometry, &p)) {
                all_on = false;
                break;
            }
        }
        if all_on {
            return false;
        }
    }

    // Whitelist of analytic surface-pair combos that
    // `intersect_surfaces` resolves without falling back to the
    // marching solver. Plane × plane is included: the "line is
    // already correct" assumption held for fillets, but CHAMFERS
    // move plane×plane edges (wall × bevel borders keep the
    // original corner line — measured 5.66e-2 = the setback
    // distance off both supports); the stale-only gate above keeps
    // every already-exact line untouched. NURBS-involving combos go
    // through a separate alternating-projection fit
    // (`approximate_curved_intersection`) — `intersect_surfaces`
    // would hit the marching path which can loop indefinitely on
    // small fillet face inputs.
    let analytic_combo = matches!(
        (&sof_1.geometry, &sof_2.geometry),
        (SurfaceGeometry::Plane { .. }, SurfaceGeometry::Plane { .. })
            | (
                SurfaceGeometry::Cylinder { .. },
                SurfaceGeometry::Plane { .. }
            )
            | (
                SurfaceGeometry::Plane { .. },
                SurfaceGeometry::Cylinder { .. }
            )
            | (
                SurfaceGeometry::Sphere { .. },
                SurfaceGeometry::Plane { .. }
            )
            | (
                SurfaceGeometry::Plane { .. },
                SurfaceGeometry::Sphere { .. }
            )
            | (
                SurfaceGeometry::Sphere { .. },
                SurfaceGeometry::Cylinder { .. }
            )
            | (
                SurfaceGeometry::Cylinder { .. },
                SurfaceGeometry::Sphere { .. }
            )
            // CONE × plane: circular-chamfer strips are cones; their
            // stale trim lines (the perp-anchor placeholder across an
            // interrupted rim's gap) must reconcile to the conic like
            // any other stale chord (PK-171: the notched flange's gap
            // chord shipped 0.53 off the cone because Cone pairs
            // silently skipped here). The edge-triple trig dispatcher
            // resolves plane×cone analytically (ellipse) with the
            // marching fallback for parabola/hyperbola sections.
            | (SurfaceGeometry::Cone { .. }, SurfaceGeometry::Plane { .. })
            | (SurfaceGeometry::Plane { .. }, SurfaceGeometry::Cone { .. })
    );
    let nurbs_combo = matches!(
        (&sof_1.geometry, &sof_2.geometry),
        (SurfaceGeometry::Cylinder { .. }, SurfaceGeometry::Nurbs(_))
            | (SurfaceGeometry::Nurbs(_), SurfaceGeometry::Cylinder { .. })
            | (SurfaceGeometry::Sphere { .. }, SurfaceGeometry::Nurbs(_))
            | (SurfaceGeometry::Nurbs(_), SurfaceGeometry::Sphere { .. })
            | (SurfaceGeometry::Plane { .. }, SurfaceGeometry::Nurbs(_))
            | (SurfaceGeometry::Nurbs(_), SurfaceGeometry::Plane { .. })
            // Strip x strip boundaries (the cross-seam join of two
            // NURBS strip halves): a chord between settled endpoints
            // must refit onto both strips like any other stale chord
            // (measured: the boss top seam's cross-strip Line sat at
            // a 2.77e-3 mid-chord sagitta forever).
            | (SurfaceGeometry::Nurbs(_), SurfaceGeometry::Nurbs(_))
    );
    if !analytic_combo && !nurbs_combo {
        return false;
    }

    // Stage 2b.57: NURBS combos route directly to `try_curved_fit`
    // (chord-projection alternating-projection). Going through
    // `compute_edge_geometry` works (with 2b.54's timeout +
    // fallback) but adds ~4 s to the popkern-blend test suite —
    // marching is genuinely slower than chord projection on small
    // NURBS fillet face inputs, even after 2b.56's AABB tightening.
    // The unified-API benefit doesn't justify the overhead for
    // this consumer; the marching catch-all stays available in
    // `compute_edge_geometry` for callers that prefer the
    // EdgeTriple shape over per-input optimisation.
    if nurbs_combo {
        return try_curved_fit(brep, edge, he, twin);
    }

    let v1 = he_data.origin;
    let v2 = twin_data.origin;
    let p1 = brep.vertices[v1].position;
    let p2 = brep.vertices[v2].position;

    // Futility pre-gate (PK-274): the acceptance below requires the
    // reconciled curve's trim endpoints to land within `snap_tol` of
    // the existing vertices, and every point of a surface×surface
    // intersection lies ON both surfaces — so a vertex farther than
    // `snap_tol` from either surface makes acceptance PROVABLY
    // impossible and the intersection solve is skipped. This is what
    // the analytic dispatchers already do (AnchorOffSurface, µs);
    // the marching catch-all has no such check and was measured
    // spending 0.75–24.7 s PER EDGE concluding NoIntersection on
    // vertex-blend-fallback corner seams whose anchors sit 4.1e-2 off both
    // surfaces (56.1 s of a 75.7 s matrix test, zero rewrites).
    // Closed-form distances only (Plane/Cylinder/Sphere); pairs with
    // a Cone keep the full dispatch — the cone distance is slant-
    // corrected and PK-171's notched-flange reconcile depends on
    // this path staying permissive.
    let snap_tol = 5.0e-3;
    let exact_distance =
        |sfc: &SurfaceGeometry, p: &popkern_core::geometry::Point3| -> Option<f64> {
            match sfc {
                SurfaceGeometry::Plane { origin, normal, .. } => {
                    let n = normal.length();
                    (n > 1e-300).then(|| ((*p - *origin).dot(normal) / n).abs())
                }
                SurfaceGeometry::Cylinder {
                    origin,
                    axis,
                    radius,
                } => {
                    let d = *p - *origin;
                    let along = d.dot(axis);
                    Some(((d - *axis * along).length() - radius).abs())
                }
                SurfaceGeometry::Sphere { center, radius, .. } => {
                    Some(((*p - *center).length() - radius).abs())
                }
                _ => None,
            }
        };
    for p in [&p1, &p2] {
        for sfc in [&sof_1.geometry, &sof_2.geometry] {
            if exact_distance(sfc, p).is_some_and(|d| d > snap_tol) {
                if popkern_core::trace::trace_enabled("edge-reconcile") {
                    eprintln!(
                        "  [reconcile] edge {edge:?}: skipped, endpoint provably off-surface \
                         (> snap gate {snap_tol:.1e})"
                    );
                }
                return false;
            }
        }
    }

    // Analytic combos go through the EdgeTriple factory (stage
    // 2b.51). The trig dispatcher resolves cyl/sphere/plane pairs
    // in microseconds.
    let anchors = EdgeAnchors::Open { start: p1, end: p2 };
    let opts = EdgeTripleOptions::default();
    let triple = match compute_edge_geometry(
        &sof_1.geometry,
        &sof_2.geometry,
        anchors,
        None,
        tolerance,
        &opts,
    ) {
        Ok(t) => t,
        Err(e) => {
            if popkern_core::trace::trace_enabled("edge-reconcile") {
                let kind = |g: &SurfaceGeometry| match g {
                    SurfaceGeometry::Plane { .. } => "Plane",
                    SurfaceGeometry::Cylinder { .. } => "Cyl",
                    SurfaceGeometry::Sphere { .. } => "Sphere",
                    SurfaceGeometry::Torus { .. } => "Torus",
                    _ => "?",
                };
                eprintln!(
                    "  [reconcile] edge {edge:?}: triple ERR {e:?} ({} x {})",
                    kind(&sof_1.geometry),
                    kind(&sof_2.geometry)
                );
            }
            return false;
        }
    };

    // Defensive G5 verification: EdgeTriple guarantees endpoint
    // accuracy by construction, but we double-check at the strict
    // BRep tolerance (anything looser opens a validate() OpenLoop).
    //
    // Vertices within the SNAP gate of the analytic curve are MOVED
    // onto it (with adjacent chord-line edges rebuilt): perp-split
    // anchors land ~1e-3 off the true intersection, and refusing
    // the re-trim left every corner chord at 2.5e-1/3.3e-2
    // off-surface (the flagship chain's dominant tessellation
    // opens). Anything beyond the snap gate is a genuine mismatch.
    let strict_tol = tolerance.linear().max(1e-9);
    let trim_start = triple.curve_3d.point_at(triple.t_start);
    let trim_end = triple.curve_3d.point_at(triple.t_end);
    let gap = trim_start.distance_to(&p1).max(trim_end.distance_to(&p2));
    if popkern_core::trace::trace_enabled("edge-reconcile") {
        eprintln!("  [reconcile] edge {edge:?}: analytic gap {gap:.3e}");
    }
    if gap > snap_tol {
        return false;
    }
    if gap > strict_tol {
        snap_vertex_and_adjacent_chords(brep, v1, trim_start, edge);
        snap_vertex_and_adjacent_chords(brep, v2, trim_end, edge);
    }

    // ARC-SIDE ADJUDICATION on periodic intersections: the endpoint
    // pair does not distinguish the short arc from its complement
    // (PK-097), and the trig dispatcher's pick can cross a gap the
    // topology routes around (measured: the notched flange's gap
    // chord reconciled to the SHORT r24 arc THROUGH the notch void —
    // on-surface, check-9 clean, but the strip ring then wound a
    // full period and tessellation cracked the whole rim). Choose
    // the side whose midpoint lies nearer the two faces' material
    // (their boundary-sample centroid).
    let (mut pick_t0, mut pick_t1) = (triple.t_start, triple.t_end);
    {
        use popkern_core::geometry::CurveEval as _;
        // Angular geometries evaluate at ANY angle regardless of the
        // stored trim, so the complement is always constructible with
        // a 2π period; other geometries have no complement.
        let angular = matches!(
            triple.curve_3d,
            popkern_core::geometry::CurveGeometry::Circle { .. }
                | popkern_core::geometry::CurveGeometry::Arc { .. }
                | popkern_core::geometry::CurveGeometry::Ellipse { .. }
        );
        let period = std::f64::consts::TAU;
        let span = triple.t_end - triple.t_start;
        if popkern_core::trace::trace_enabled("edge-reconcile") {
            eprintln!(
                "  [reconcile] edge {edge:?}: angular={angular} span={span:.3} kind={:?}",
                std::mem::discriminant(&triple.curve_3d)
            );
        }
        if angular && span < period * (1.0 - 1e-6) {
            // PRIOR-CURVE side vote: both candidate midpoints lie ON
            // the intersection circle, so face containment is
            // degenerate for them (measured: every vote came back
            // (false, false)). The edge's PRE-REWRITE curve carries
            // the topological side — the interrupted-rim split
            // authors the long arc, and a genuinely stale chord's
            // midpoint sits nearest the short span (preserving the
            // historical pick for the straight cases).
            let prior_mid = brep
                .edges
                .get(edge)
                .and_then(|e| e.curve.as_ref())
                .map(|c| c.geometry.point_at(0.5 * (c.t_start + c.t_end)));
            if let Some(pm) = prior_mid {
                let mid_short = triple
                    .curve_3d
                    .point_at(0.5 * (triple.t_start + triple.t_end));
                let comp_mid_t = 0.5 * (triple.t_start + triple.t_end) + 0.5 * period;
                let mid_comp = triple.curve_3d.point_at(comp_mid_t);
                if popkern_core::trace::trace_enabled("edge-reconcile") {
                    eprintln!(
                        "  [reconcile] edge {edge:?}: prior-mid vote d_short={:.3} d_comp={:.3}",
                        mid_short.distance_to(&pm),
                        mid_comp.distance_to(&pm)
                    );
                }
                if mid_comp.distance_to(&pm) < mid_short.distance_to(&pm) {
                    // Complement wins: sweep the other way round.
                    pick_t0 = triple.t_end;
                    pick_t1 = triple.t_start + period;
                    if popkern_core::trace::trace_enabled("edge-reconcile") {
                        eprintln!(
                            "  [reconcile] edge {edge:?}: complement arc chosen (prior side)"
                        );
                    }
                }
            }
        }
    }
    let side_flipped = (pick_t0, pick_t1) != (triple.t_start, triple.t_end);
    // Apply: 3D curve on edge, pcurve_a on the half-edge whose face
    // is face_1, pcurve_b on the twin (face_2's half-edge).
    let geometry = triple
        .curve_3d
        .canonical_line_in_range((pick_t0, pick_t1), tolerance.linear().max(1.0e-8))
        .unwrap_or_else(|| triple.curve_3d.clone());
    let new_curve = CurveOnEdge {
        geometry,
        t_start: pick_t0,
        t_end: pick_t1,
    };
    if let Some(edge_mut) = brep.edges.get_mut(edge) {
        edge_mut.curve = Some(new_curve);
    } else {
        return false;
    }
    if side_flipped {
        // The triple's pcurves describe the SHORT side; drop them and
        // let the lenient pcurve recompute derive the complement's.
        if let Some(he_mut) = brep.half_edges.get_mut(he) {
            he_mut.pcurve = None;
        }
        if let Some(twin_mut) = brep.half_edges.get_mut(twin) {
            twin_mut.pcurve = None;
        }
        brep.recompute_all_pcurves_lenient();
        return true;
    }
    // The triple's pcurves run the CURVE's way; a half-edge starting at
    // the curve's end takes its pcurve reversed (a chamfer's boss-wall
    // quarter arcs served both twins forward, 0.113 off, PK-819). A
    // closed edge's one vertex cannot tell: kept as authored.
    for (half_edge, pcurve) in [(he, &triple.pcurve_a), (twin, &triple.pcurve_b)] {
        let backward = half_edge_runs_backward(brep, half_edge, (pick_t0, pick_t1));
        if let Some(h) = brep.half_edges.get_mut(half_edge) {
            h.pcurve = Some(if backward {
                pcurve.reversed()
            } else {
                pcurve.clone()
            });
        }
    }
    true
}

/// Whether `half_edge` (on an open edge) starts at its edge curve's
/// `range.1` end: its origin vertex is clearly nearer that end.
fn half_edge_runs_backward(brep: &BRep, half_edge: HalfEdgeHandle, range: (f64, f64)) -> bool {
    let Some(h) = brep.half_edges.get(half_edge) else {
        return false;
    };
    let (Some(origin), Some(curve)) = (
        brep.vertices.get(h.origin).map(|v| v.position),
        brep.edges.get(h.edge).and_then(|e| e.curve.as_ref()),
    ) else {
        return false;
    };
    let (a, b) = (
        curve.geometry.point_at(range.0),
        curve.geometry.point_at(range.1),
    );
    let chord = a.distance_to(&b);
    let (da, db) = (origin.distance_to(&a), origin.distance_to(&b));
    chord > brep.tolerance.linear() && db + 1.0e-3 * chord < da
}
