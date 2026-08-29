//! Curved-support extensions for vertex-blend corner topology.
//!
//! The planar corner machinery identifies which fillets border a
//! face by trim-on-plane membership and finds mef anchors by exact
//! radius distance on the OUTER loop. Both assumptions fail on
//! curved supports (PK-141 boss-on-wall: the fillet chain lives on
//! the wall's HOLE loop, membership needs surface projection, and
//! tangent-offset anchors land ~3% of r off the chord distance).
//! This module supplies the surface-general pieces; the planar fast
//! paths in `blend_topo` are unchanged.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{Point3, SurfaceEval};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use popkern_core::tolerance::Tolerance;

use crate::construct::EdgeGeometry;

/// Half-edges of every loop of `face` (outer first, then holes).
/// The corner machinery must see hole loops: on a pierced wall the
/// fillet trim anchors live on the hole boundary, and an outer-only
/// walk reports "0 trim-endpoint anchors" for every corner there.
pub(crate) fn face_all_loop_half_edges(brep: &BRep, face: FaceHandle) -> Vec<HalfEdgeHandle> {
    let Some(face_data) = brep.faces.get(face) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    let mut loops = vec![face_data.outer_loop];
    loops.extend(face_data.inner_loops.iter().copied());
    for lh in loops {
        let Some(loop_data) = brep.loops.get(lh) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cursor = start;
        for _ in 0..brep.half_edges.len() + 1 {
            if !brep.half_edges.contains_key(cursor) {
                break;
            }
            out.push(cursor);
            cursor = brep.half_edges[cursor].next;
            if cursor == start {
                break;
            }
        }
    }
    out
}

/// Mean distance of 8 interior samples of `trim` to a non-planar
/// surface, via closest-point projection. `None` when the surface
/// kind offers no projection (membership is then unverifiable, and
/// the caller must not claim the trim for this face).
pub(crate) fn trim_mean_surface_distance(
    geometry: &SurfaceGeometry,
    trim: &popkern_core::geometry::NurbsCurve,
) -> Option<f64> {
    let (t0, t1) = trim.parameter_range();
    let n_samples = 8usize;
    let mut mean = 0.0_f64;
    for i in 0..n_samples {
        let t = t0 + (t1 - t0) * (i as f64 + 0.5) / (n_samples as f64);
        let p = trim.point_at(t);
        let (u, v) = project_point_to_uv(geometry, &p)?;
        mean += p.distance_to(&geometry.point_at(u, v));
    }
    Some(mean / n_samples as f64)
}

/// Trim-endpoint anchor tier for curved supports: for each fillet
/// whose trim lies on `face`, take the trim ENDPOINT nearest the
/// corner vertex (within `3r` reach) and match it to the nearest
/// face-loop vertex inside a curvature-tolerant band (`5%` of r —
/// measured tangent-offset placement error on the r=1 wall is 3.4%
/// of r, while the nearest wrong candidate, the chord-crossing
/// point, sits at `(√2−1)·r ≈ 41%`). Radius-exact tier 1 stays the
/// primary; this tier only runs when it found fewer than 2 anchors.
pub(crate) fn trim_endpoint_anchors(
    brep: &BRep,
    face: FaceHandle,
    vertex: VertexHandle,
    edge_geometries: &[EdgeGeometry],
    on_face_radii: &[f64],
) -> Vec<VertexHandle> {
    let Some(v_pos) = brep.vertices.get(vertex).map(|v| v.position) else {
        return Vec::new();
    };
    let Some(sof) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
        return Vec::new();
    };
    let trace = crate::direct_apex_trace::enabled();
    let face_owner_surface = format!("{:?}", sof.geometry);
    if trace {
        eprintln!(
            "[direct-apex-trace]     t3 face {face:?} surface {:?} vertex {vertex:?} at ({:.4},{:.4},{:.4})",
            std::mem::discriminant(&sof.geometry),
            v_pos.x,
            v_pos.y,
            v_pos.z
        );
    }
    let loop_hes = face_all_loop_half_edges(brep, face);
    let mut anchors: Vec<VertexHandle> = Vec::new();
    for (ei, eg) in edge_geometries.iter().enumerate() {
        // Same self-exclusion as `collect_radii_for_face`: on a
        // planar-support fillet a strip face never anchors on its
        // own trim endpoints (that corner belongs to
        // blend_seam_corner).
        if planar_supports(eg) && format!("{:?}", eg.surface) == face_owner_surface {
            continue;
        }
        for (ti, trim) in [&eg.trim_3d_1, &eg.trim_3d_2].into_iter().enumerate() {
            let Some(mean) = trim_mean_surface_distance(&sof.geometry, trim) else {
                continue;
            };
            let r = eg.effective_radius.max(1e-12);
            if mean > r * 0.1 {
                if trace {
                    eprintln!(
                        "[direct-apex-trace]     t3 eg{ei} trim{ti}: membership declined mean {mean:.3e} vs {:.3e}",
                        r * 0.1
                    );
                }
                continue;
            }
            let (t0, t1) = trim.parameter_range();
            let ends = [trim.point_at(t0), trim.point_at(t1)];
            let Some(end) = ends
                .iter()
                .filter(|p| v_pos.distance_to(p) <= r * 3.0)
                .min_by(|a, b| v_pos.distance_to(a).total_cmp(&v_pos.distance_to(b)))
            else {
                if trace {
                    eprintln!(
                        "[direct-apex-trace]     t3 eg{ei} trim{ti}: on-face but no endpoint within 3r (ends {:.3e} / {:.3e})",
                        v_pos.distance_to(&ends[0]),
                        v_pos.distance_to(&ends[1])
                    );
                }
                continue;
            };
            let hit = nearest_loop_vertex(brep, &loop_hes, vertex, *end, r * 0.05);
            if trace {
                eprintln!(
                    "[direct-apex-trace]     t3 eg{ei} trim{ti}: end d(V) {:.3e} -> loop match {hit:?}",
                    v_pos.distance_to(end)
                );
            }
            if let Some(vh) = hit
                && !anchors.contains(&vh)
            {
                anchors.push(vh);
            }
        }
    }
    let _ = on_face_radii;
    anchors
}

fn nearest_loop_vertex(
    brep: &BRep,
    loop_hes: &[HalfEdgeHandle],
    exclude: VertexHandle,
    target: Point3,
    band: f64,
) -> Option<VertexHandle> {
    let mut best: Option<(VertexHandle, f64)> = None;
    for &he in loop_hes {
        let vh = brep.half_edges.get(he)?.origin;
        if vh == exclude {
            continue;
        }
        let pos = brep.vertices.get(vh)?.position;
        let d = pos.distance_to(&target);
        if d <= band && best.is_none_or(|(_, bd)| d < bd) {
            best = Some((vh, d));
        }
    }
    best.map(|(vh, _)| vh)
}

/// Radii of fillets whose trim curves lie on `face`'s plane — i.e.
/// the fillets actually bordering this face at `vertex`.
///
/// A trim curve is *on* `face` iff all its sample points lie on
/// `face`'s plane (mean plane-distance < `tol_plane`). This is far
/// more robust than checking a single closest-to-V point whose
/// sampling error bleeds into the face normal direction. Mean-of-8
/// averages out sampling noise, while a trim curve that truly does
/// NOT lie on this plane shows a mean distance ≥ the dihedral
/// angle's offset (much larger than sampling noise).
pub(crate) fn collect_radii_for_face(
    brep: &BRep,
    face: FaceHandle,
    vertex: VertexHandle,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
) -> Vec<f64> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let Some(face_data) = brep.faces.get(face) else {
        return Vec::new();
    };
    let Some(sof) = face_data.surface.as_ref() else {
        return Vec::new();
    };
    // Planes keep the exact plane-distance fast path; every other
    // surface kind measures membership through closest-point
    // projection (curved supports were structurally excluded here,
    // which mitred EVERY corner on the pierced-wall chain).
    let plane = match &sof.geometry {
        SurfaceGeometry::Plane { origin, normal, .. } => Some((*origin, *normal)),
        _ => None,
    };
    // Curved membership excludes SELF-matches: a fillet's trims lie
    // on its own strip by construction, so a strip face claiming its
    // OWN fillet's trims enters the corner machinery as if it were a
    // support (measured: 59 welded opens on both planar swept-corner
    // hardening gates — the changed vertex-blend topology blocked
    // the downstream blend_seam_corner repair). A strip face may
    // still claim a PARTNER fillet's trims: on region-path bodies
    // the junction faces are the strips themselves (the boss chain's
    // one completing junction builds through exactly that claim).
    let face_owner_surface = format!("{:?}", sof.geometry);

    let tol_plane = (tolerance.linear() * 100.0).max(1e-4);

    let mut radii = Vec::with_capacity(edge_geometries.len());
    for eg in edge_geometries {
        if plane.is_none()
            && planar_supports(eg)
            && format!("{:?}", eg.surface) == face_owner_surface
        {
            continue;
        }
        for trim in [&eg.trim_3d_1, &eg.trim_3d_2] {
            let mean_dist = if let Some((plane_origin, plane_normal)) = plane {
                let (t0, t1) = trim.parameter_range();
                let n_samples = 8usize;
                let mut mean = 0.0_f64;
                for i in 0..n_samples {
                    let t = t0 + (t1 - t0) * (i as f64 + 0.5) / (n_samples as f64);
                    let p = trim.point_at(t);
                    let delta = p - plane_origin;
                    mean += delta.dot(&plane_normal).abs();
                }
                mean / n_samples as f64
            } else {
                match super::blend_topo_curved::trim_mean_surface_distance(&sof.geometry, trim) {
                    Some(mean) => mean,
                    None => continue,
                }
            };

            if mean_dist < tol_plane {
                // §5.3 item #8 — for variable-radius fillets the
                // trim-endpoint distance from V equals the per-vertex
                // radius, not the edge-wide average. Select the
                // endpoint radius that matches `vertex`; fall back to
                // `effective_radius` when V is interior to the chain.
                let r = if vertex == eg.analysis.vertex_start {
                    eg.radius_at_start
                } else if vertex == eg.analysis.vertex_end {
                    eg.radius_at_end
                } else {
                    eg.effective_radius
                };
                radii.push(r);
                break; // At most one of eg's two trims lies on `face`.
            }
        }
    }
    radii
}

/// Both of a fillet's support faces are planes. Self-claims (a strip
/// face anchoring on its own fillet's trims) are excluded exactly for
/// this class: planar-support corners belong to `blend_seam_corner`,
/// while curved-support junctions have only the strips at the corner
/// vertex and must anchor through them.
pub(crate) fn planar_supports(eg: &EdgeGeometry) -> bool {
    use popkern_core::geometry::surface::SurfaceGeometry as SG;
    matches!(eg.analysis.surface_1, SG::Plane { .. })
        && matches!(eg.analysis.surface_2, SG::Plane { .. })
}
