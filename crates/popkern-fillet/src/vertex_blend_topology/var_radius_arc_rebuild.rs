//! Phase 6c-D Stage 3a — analytic cap-tangent arc edge rebuild.
//!
//! After Stage 2's per-vertex snap places cap-loop endpoints onto
//! `sphere ∩ cylinder` great circles, the EDGE CURVES between
//! those endpoints still trace the OLD rolling-ball geometry (the
//! NURBS curves authored during fillet construction, before any
//! sphere-fit substitution). The cap-vs-cylinder G0 metric samples
//! these edge curves via pcurve evaluation on each side, so even
//! perfectly-snapped endpoints leave a residual G0 from the
//! curve interior drifting off the analytic intersection.
//!
//! This module rewrites each cap-tangent arc edge with an analytic
//! `CurveGeometry::Ellipse` arc on the `sphere ∩ cylinder` great
//! circle, with closed-form Linear pcurves on the cylinder side
//! and a closest-point-projected NurbsCurve pcurve on the cap
//! side. Mirrors `var_radius::set_cut_edge_geometry` (which fires
//! on the Sphere-cap rung path) but works for any cap surface
//! kind including Gregory.
//!
//! Activation: invoked from `reroute_neighbour_cylinder_edges`'s
//! Stage 2 fallback after the vertex snap pass succeeds. Cube
//! fixtures route through `apply_path_b_plan`'s
//! `set_cut_edge_geometry` instead, so are unaffected.

use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{CurveGeometry, Point3, Vec3, fitting};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

/// Phase 6c-D Stage 3a entry point. For each cap face, walk its
/// outer loop and rewrite every cap-tangent arc edge (twin face
/// is the cap's OWN cylinder) with an analytic Ellipse arc on
/// the `sphere ∩ cylinder` great circle. Uses the same
/// re-derivation logic as Stage 2 to identify each cap's own
/// sphere/cylinder pair.
///
/// Returns the number of edges rewritten.
pub(crate) fn rebuild_cap_tangent_arc_curves(
    brep: &mut BRep,
    cap_faces: &[FaceHandle],
    per_fillet_caps: &[(Point3, f64)],
) -> usize {
    if cap_faces.len() != per_fillet_caps.len() {
        return 0;
    }
    let mut rewritten = 0usize;

    // Phase A — collect rewrites per edge. Don't mutate during
    // iteration to keep cap-loop walking deterministic.
    type Job = (
        HalfEdgeHandle, // cap_side_he
        Point3,         // sphere_centre
        f64,            // sphere_radius
        Point3,         // cyl_origin
        Vec3,           // cyl_axis
        f64,            // cyl_radius
    );
    let mut jobs: Vec<Job> = Vec::new();
    let mut seen_edges: std::collections::BTreeSet<popkern_core::topology::EdgeHandle> =
        std::collections::BTreeSet::new();

    for &cap_face in cap_faces {
        let Some((sphere_centre, sphere_radius)) =
            re_derive_cap_sphere(brep, cap_face, per_fillet_caps)
        else {
            continue;
        };
        let Some(face_data) = brep.faces.get(cap_face) else {
            continue;
        };
        let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let hes = collect_loop_half_edges(brep, start);
        for &he_h in &hes {
            let Some(he) = brep.half_edges.get(he_h) else {
                continue;
            };
            let edge = he.edge;
            if seen_edges.contains(&edge) {
                continue;
            }
            let twin_face = brep
                .half_edges
                .get(he.twin)
                .and_then(|t| brep.loops.get(t.loop_))
                .map(|l| l.face);
            let Some(twin_face) = twin_face else { continue };
            let Some(SurfaceGeometry::Cylinder {
                origin: cyl_origin,
                axis: cyl_axis,
                radius: cyl_radius,
            }) = brep
                .faces
                .get(twin_face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| s.geometry.clone())
            else {
                continue;
            };
            // Only OWN cylinder: sphere centre must lie on its axis.
            let Some(axis_unit) = cyl_axis.normalized() else {
                continue;
            };
            let to_centre = Vec3::new(
                sphere_centre.x - cyl_origin.x,
                sphere_centre.y - cyl_origin.y,
                sphere_centre.z - cyl_origin.z,
            );
            let along = to_centre.dot(&axis_unit);
            let proj = Point3::new(
                cyl_origin.x + axis_unit.x * along,
                cyl_origin.y + axis_unit.y * along,
                cyl_origin.z + axis_unit.z * along,
            );
            if sphere_centre.distance_to(&proj) > sphere_radius * 1e-3 {
                continue;
            }
            seen_edges.insert(edge);
            jobs.push((
                he_h,
                sphere_centre,
                sphere_radius,
                cyl_origin,
                cyl_axis,
                cyl_radius,
            ));
        }
    }

    // Phase B — apply rewrites.
    for (he_h, sphere_centre, sphere_radius, cyl_origin, cyl_axis, cyl_radius) in jobs {
        let Some(he) = brep.half_edges.get(he_h) else {
            continue;
        };
        let edge_h = he.edge;
        let twin_h = he.twin;
        let next_h = he.next;
        let v_origin = he.origin;
        let Some(next_he) = brep.half_edges.get(next_h) else {
            continue;
        };
        let v_dest = next_he.origin;

        // Both endpoints should now be on sphere ∩ cylinder thanks
        // to Stage 2's snap. Read their post-snap positions and
        // orient the arc to match the EDGE's primary half-edge
        // direction so `edge.curve.point_at(t_start)` matches the
        // primary's origin vertex (validate's edge invariant).
        let primary_he = brep.edges.get(edge_h).map(|e| e.half_edge);
        let we_are_primary = primary_he == Some(he_h);
        let p_origin = brep.vertices[v_origin].position;
        let p_dest = brep.vertices[v_dest].position;
        let (p_start, p_end) = if we_are_primary {
            (p_origin, p_dest)
        } else {
            (p_dest, p_origin)
        };

        let Some(arc) = build_analytic_arc(
            p_start,
            p_end,
            sphere_centre,
            sphere_radius,
            cyl_origin,
            cyl_axis,
            cyl_radius,
        ) else {
            continue;
        };

        // Cap-side surface (Gregory or Sphere) for pcurve refit.
        let cap_face = brep.loops.get(he.loop_).map(|l| l.face);
        let cap_surf = cap_face.and_then(|f| {
            brep.faces
                .get(f)
                .and_then(|fd| fd.surface.as_ref())
                .map(|s| s.geometry.clone())
        });

        // Cylinder UV: (theta, height). The Linear pcurve goes
        // straight across at fixed v=v_cyl, parameterised by theta.
        let axis_unit = match cyl_axis.normalized() {
            Some(u) => u,
            None => continue,
        };
        let v_cyl = (Vec3::new(
            p_start.x - cyl_origin.x,
            p_start.y - cyl_origin.y,
            p_start.z - cyl_origin.z,
        ))
        .dot(&axis_unit);
        let pcurve_cyl = Pcurve::linear(
            (arc.theta_a, v_cyl),
            (arc.theta_b, v_cyl),
            arc.theta_a,
            arc.theta_b,
        );

        // Cap-side pcurve: refit by closest-point projection from
        // analytic arc samples. For Sphere caps a Linear pcurve is
        // exact (great circle is a meridian / equator), but for
        // Gregory caps we refit a degree-3 NURBS pcurve.
        let pcurve_cap = cap_surf
            .as_ref()
            .and_then(|surf| refit_arc_pcurve_on_surface(&arc, surf, arc.theta_a, arc.theta_b));

        // Write the analytic 3-D curve to the edge.
        if let Some(edge_data) = brep.edges.get_mut(edge_h) {
            edge_data.curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Ellipse {
                    center: arc.circle_centre,
                    u_basis: arc.u_basis,
                    v_basis: arc.v_basis,
                    start_angle: arc.theta_a,
                    end_angle: arc.theta_b,
                },
                t_start: arc.theta_a,
                t_end: arc.theta_b,
            });
            rewritten += 1;
        }

        // Wire the pcurves. We're iterating in the cap's outer loop,
        // so `he_h` is always the cap-side half-edge and `twin_h` is
        // the cylinder-side. Primary/twin is independent — that's
        // about edge curve orientation (handled above by p_start /
        // p_end swapping), not about which side is which surface.
        if let Some(pc) = pcurve_cap {
            brep.half_edges[he_h].pcurve = Some(pc.into());
        }
        brep.half_edges[twin_h].pcurve = Some(pcurve_cyl);
    }

    rewritten
}

struct AnalyticArc {
    circle_centre: Point3,
    u_basis: Vec3,
    v_basis: Vec3,
    theta_a: f64,
    theta_b: f64,
}

fn build_analytic_arc(
    p_origin: Point3,
    p_dest: Point3,
    sphere_centre: Point3,
    sphere_radius: f64,
    cyl_origin: Point3,
    cyl_axis: Vec3,
    cyl_radius: f64,
) -> Option<AnalyticArc> {
    // Sphere centre on cylinder axis is required by caller (we
    // already filtered for that). Verify cylinder/sphere radii
    // match (rolling-ball construction guarantees this).
    if (cyl_radius - sphere_radius).abs() > sphere_radius * 1e-2 {
        return None;
    }
    let axis_unit = cyl_axis.normalized()?;
    let (perp1, perp2) = orthonormal_basis(axis_unit);

    // Circle centre = projection of sphere centre on cylinder axis.
    let to_centre = Vec3::new(
        sphere_centre.x - cyl_origin.x,
        sphere_centre.y - cyl_origin.y,
        sphere_centre.z - cyl_origin.z,
    );
    let centre_along = to_centre.dot(&axis_unit);
    let circle_centre = Point3::new(
        cyl_origin.x + axis_unit.x * centre_along,
        cyl_origin.y + axis_unit.y * centre_along,
        cyl_origin.z + axis_unit.z * centre_along,
    );

    let theta = |p: Point3| -> f64 {
        let d = Vec3::new(
            p.x - circle_centre.x,
            p.y - circle_centre.y,
            p.z - circle_centre.z,
        );
        d.dot(&perp2).atan2(d.dot(&perp1))
    };
    let theta_a = theta(p_origin);
    let theta_b = theta(p_dest);
    // Pick the shorter arc direction (delta in [-π, π]).
    let mut delta = theta_b - theta_a;
    while delta > std::f64::consts::PI {
        delta -= std::f64::consts::TAU;
    }
    while delta < -std::f64::consts::PI {
        delta += std::f64::consts::TAU;
    }
    let theta_b_adj = theta_a + delta;

    Some(AnalyticArc {
        circle_centre,
        u_basis: perp1 * cyl_radius,
        v_basis: perp2 * cyl_radius,
        theta_a,
        theta_b: theta_b_adj,
    })
}

/// Sample the analytic arc at uniform parameters and refit a
/// degree-3 NurbsCurve pcurve in the cap surface's UV space via
/// closest-point projection. Returns `None` if any projection
/// fails or the fit raises an error.
fn refit_arc_pcurve_on_surface(
    arc: &AnalyticArc,
    surface: &SurfaceGeometry,
    t0: f64,
    t1: f64,
) -> Option<NurbsCurve> {
    const N: usize = 12;
    let mut points: Vec<Point3> = Vec::with_capacity(N);
    let mut params: Vec<f64> = Vec::with_capacity(N);
    for i in 0..N {
        let s = i as f64 / (N - 1) as f64;
        let theta = t0 + s * (t1 - t0);
        // 3D point on arc: circle_centre + cos·u_basis + sin·v_basis.
        let p = Point3::new(
            arc.circle_centre.x + theta.cos() * arc.u_basis.x + theta.sin() * arc.v_basis.x,
            arc.circle_centre.y + theta.cos() * arc.u_basis.y + theta.sin() * arc.v_basis.y,
            arc.circle_centre.z + theta.cos() * arc.u_basis.z + theta.sin() * arc.v_basis.z,
        );
        let (u, v) = project_point_to_uv(surface, &p)?;
        points.push(Point3::new(u, v, 0.0));
        params.push(s);
    }
    let degree = 3usize;
    let n_cps = (degree + 2).min(points.len());
    fitting::fit_curve(&points, &params, degree, n_cps).ok()
}

/// Re-derive the per-fillet sphere for `cap_face` by walking its
/// outer loop and finding the `per_fillet_caps[i]` whose centre
/// lies on at least one of the loop's cylinder-twin axes (within
/// `r * 5%`). Same logic as `var_radius_snap::re_derive_cap_sphere`,
/// duplicated here to avoid pub-promoting that helper.
fn re_derive_cap_sphere(
    brep: &BRep,
    cap_face: FaceHandle,
    per_fillet_caps: &[(Point3, f64)],
) -> Option<(Point3, f64)> {
    let face_data = brep.faces.get(cap_face)?;
    let loop_data = brep.loops.get(face_data.outer_loop)?;
    let start = loop_data.half_edge;
    let hes = collect_loop_half_edges(brep, start);

    for &he_h in &hes {
        let Some(he) = brep.half_edges.get(he_h) else {
            continue;
        };
        let twin_face = brep
            .half_edges
            .get(he.twin)
            .and_then(|t| brep.loops.get(t.loop_))
            .map(|l| l.face);
        let Some(twin_face) = twin_face else { continue };
        let Some(SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        }) = brep
            .faces
            .get(twin_face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| s.geometry.clone())
        else {
            continue;
        };
        let Some(axis_unit) = axis.normalized() else {
            continue;
        };
        for &(centre, sphere_radius) in per_fillet_caps {
            if (sphere_radius - radius).abs() > radius * 0.05 {
                continue;
            }
            let v = Vec3::new(
                centre.x - origin.x,
                centre.y - origin.y,
                centre.z - origin.z,
            );
            let along = v.dot(&axis_unit);
            let proj = Point3::new(
                origin.x + axis_unit.x * along,
                origin.y + axis_unit.y * along,
                origin.z + axis_unit.z * along,
            );
            if centre.distance_to(&proj) < radius * 0.05 {
                return Some((centre, sphere_radius));
            }
        }
    }
    None
}

fn collect_loop_half_edges(brep: &BRep, start: HalfEdgeHandle) -> Vec<HalfEdgeHandle> {
    let mut hes = Vec::new();
    let mut cursor = start;
    let mut safety = 64usize;
    loop {
        if safety == 0 {
            break;
        }
        safety -= 1;
        hes.push(cursor);
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    hes
}
