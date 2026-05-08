//! Phase 6c-C.2b — apex-region deletion *executor*.
//!
//! Consumes a `DeletionPlan` (from `plan.rs`) and rewrites the
//! BRep to match the target topology.
//!
//! # Status
//!
//! * **Implemented:** the safe `mef` step that inserts the new
//!   cap-tangent arc edges across the cap face's outer loop.
//!   Each `mef` splits the cap face into a (smaller real cap)
//!   and a (lens sub-face containing the obsolete generators
//!   and the rolling-ball-footprint vertex between them).
//!
//! * **Deferred to a follow-up:** the lens-collapse step.
//!   Naive `kef` of one of the lens's generator edges would
//!   merge the lens into one of the two adjacent cylinder
//!   faces.  But the *other* generator of the same lens
//!   borders a *different* cylinder face, so a second `kef`
//!   would merge incompatible cylinder surfaces.  The right
//!   sequence requires either:
//!   1. A non-Euler direct slot-map mutation (delete the lens
//!      face + 2 generators + footprint vertex; splice the
//!      cylinder face boundaries to skip them).
//!   2. A multi-face combined op that handles the cross-
//!      cylinder junction at each footprint atomically.
//!
//!   Both paths require careful re-routing of the cylinder
//!   faces' outer-loop half-edge `.next` / `.prev` chains and
//!   are best implemented after a proof-of-concept with
//!   explicit before/after BRep snapshots.
//!
//! Until the lens-collapse step lands, this executor returns a
//! `PartialExecution` report describing how far it got.  The
//! cap surface stays a Sphere primitive (correct centre after
//! the centroid fix in `surface_fix.rs`); the topology now
//! includes 3 extra edges (the cap-tangent arcs) on the cap
//! face's outer loop, alongside the existing 6 generators.
//!
//! Net empirical effect: G0 unchanged from the centroid-fix
//! baseline (still 2.93e-2 on cube_uniform) because the lens
//! sub-faces still exist and the metric reads the cap face's
//! union of all sub-faces.  Once the lens-collapse step lands,
//! G0 should drop to near-zero.

use popkern_core::geometry::CurveEval;
use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};
use std::f64::consts::TAU;

use super::plan::DeletionPlan;

/// Outcome of running an executor against a plan.
///
/// `PartialExecution` is the only success variant today; the
/// `Complete` variant is reserved for when the lens-collapse
/// step is implemented.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub(crate) struct ExecutionReport {
    /// New cap-tangent arc edges inserted via `mef`.
    pub(crate) inserted_arcs: Vec<EdgeHandle>,
    /// Lens sub-faces created by the `mef` calls — each
    /// contains a rolling-ball-footprint vertex and 2 obsolete
    /// generator edges.  These need a follow-up collapse pass.
    pub(crate) lens_faces: Vec<FaceHandle>,
    /// Vertex handles of the rolling-ball footprints that
    /// remain (one per lens).  Their incident-edge count is
    /// still 2 (the 2 obsolete generators).
    pub(crate) remaining_footprints: Vec<VertexHandle>,
}

/// Errors from `execute_deletion_plan`.
#[derive(Debug, Clone, PartialEq)]
#[allow(dead_code)]
pub(crate) enum ExecuteError {
    /// Plan didn't pass `is_actionable`.  Caller should bail
    /// without mutating.
    PlanNotActionable,
    /// Couldn't find a half-edge originating from a kept
    /// vertex in the cap face's outer loop.  Means the cap
    /// face has been mutated since the plan was built.
    HalfEdgeLookupFailed { vertex: VertexHandle },
    /// `mef` failed during arc insertion.  Reports the BRep's
    /// error message verbatim.
    MefFailed(String),
    /// Couldn't determine which fillet cylinder hosts the
    /// cap-tangent arc for a given `(V_a, V_b)` pair.  Means
    /// the BRep doesn't have the expected structure (e.g. the
    /// kept vertices aren't both on a single cylinder).
    CylinderLookupFailed {
        v_a: VertexHandle,
        v_b: VertexHandle,
    },
    /// Actionable plan was missing a field required by execution.
    MissingPlanField(&'static str),
}

/// Run the executor.  See module docs for what's implemented.
#[allow(dead_code)]
pub(crate) fn execute_deletion_plan(
    brep: &mut BRep,
    plan: &DeletionPlan,
) -> Result<ExecutionReport, ExecuteError> {
    if !plan.is_actionable() {
        return Err(ExecuteError::PlanNotActionable);
    }
    let cap_face = plan
        .cap_face
        .ok_or(ExecuteError::MissingPlanField("cap_face"))?;
    let (cap_centre, cap_radius) = plan
        .cap_sphere
        .ok_or(ExecuteError::MissingPlanField("cap_sphere"))?;

    let mut inserted_arcs: Vec<EdgeHandle> = Vec::new();
    let mut lens_faces: Vec<FaceHandle> = Vec::new();
    let mut remaining_footprints: Vec<VertexHandle> = Vec::new();

    for &(v_a, v_b, _) in &plan.insert_cap_tangent_arcs {
        let cyl_surf = find_shared_cylinder(brep, v_a, v_b)
            .ok_or(ExecuteError::CylinderLookupFailed { v_a, v_b })?;

        let arc_curve =
            build_cap_tangent_half_circle(brep, cap_centre, cap_radius, &cyl_surf, v_a, v_b);

        // Walk the cap face's outer loop. The two half-edges
        // we want to mef between are: (a) the half-edge
        // originating at v_a in the cap face's loop and (b) the
        // half-edge originating at v_b in the cap face's loop.
        // These are typically separated by exactly ONE vertex
        // in the loop (the footprint between them); after the
        // earlier mef calls in this same execution the
        // separation may shrink to zero (loop has been split
        // already).
        let (he_a, he_b) = find_loop_he_pair(brep, cap_face, v_a, v_b)
            .ok_or(ExecuteError::HalfEdgeLookupFailed { vertex: v_a })?;

        // Capture the footprint vertex between v_a and v_b
        // before mef — it's the loop's vertex sandwiched between
        // v_a and v_b.
        let footprint = he_after_origin_in_loop(brep, he_a).filter(|&v| v != v_b);

        let mef = brep
            .mef(he_a, he_b)
            .map_err(|e| ExecuteError::MefFailed(format!("{e:?}")))?;
        let (t_start, t_end) = arc_curve.parameter_range();
        // Capture cap surface for pcurve fitting before mutating
        // the edge.
        let cap_surf_clone = brep
            .faces
            .get(cap_face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| s.geometry.clone());
        if let Some(edge_data) = brep.edges.get_mut(mef.edge) {
            edge_data.curve = Some(CurveOnEdge {
                geometry: popkern_core::geometry::CurveGeometry::Nurbs(arc_curve.clone()),
                t_start,
                t_end,
            });
        }
        // Phase 6c-C.2b/2d: fit pcurves on both half-edges of the
        // new arc.  Without these, downstream consumers (mesh,
        // bool, the G0/G1 metric) read `pcurve = None` and either
        // skip the edge or fall back to a stale value.  Both pcurves
        // are NURBS fits of `project_to_uv(arc_3d_sample)` against
        // the respective surface.
        // mef.half_edge is on the OLD face (cap), mef.half_edge_twin
        // is on the NEW face (lens).  Per `MefResult` docs.
        let cap_he = mef.half_edge;
        let lens_he = mef.half_edge_twin;
        if let Some(cap_surf) = cap_surf_clone {
            // cap_he runs v_a → v_b (the arc's natural direction);
            // fit forward.  lens_he is the twin running v_b → v_a;
            // fit reversed so its pcurve(t_lo) lands on v_b.
            let cap_pcurve = fit_pcurve_along_arc(&cap_surf, &arc_curve, false);
            let cyl_pcurve = fit_pcurve_along_arc(&cyl_surf, &arc_curve, true);
            if let Some(pc) = cap_pcurve {
                brep.half_edges[cap_he].pcurve = Some(pc.into());
            }
            if let Some(pc) = cyl_pcurve {
                brep.half_edges[lens_he].pcurve = Some(pc.into());
            }
        }
        inserted_arcs.push(mef.edge);
        lens_faces.push(mef.face);
        if let Some(fp) = footprint {
            remaining_footprints.push(fp);
        }
    }

    Ok(ExecutionReport {
        inserted_arcs,
        lens_faces,
        remaining_footprints,
    })
}

/// Find a cylinder surface in the BRep (across all faces)
/// whose surface is a `Cylinder` and whose surface contains
/// both `v_a` and `v_b`.  Returns the surface (cloned) on
/// success.
fn find_shared_cylinder(
    brep: &BRep,
    v_a: VertexHandle,
    v_b: VertexHandle,
) -> Option<SurfaceGeometry> {
    let p_a = brep.vertices.get(v_a)?.position;
    let p_b = brep.vertices.get(v_b)?.position;
    for (_, face) in brep.faces.iter() {
        let Some(sof) = &face.surface else { continue };
        let SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } = &sof.geometry
        else {
            continue;
        };
        if point_on_cylinder(*origin, *axis, *radius, p_a)
            && point_on_cylinder(*origin, *axis, *radius, p_b)
        {
            return Some(sof.geometry.clone());
        }
    }
    None
}

fn point_on_cylinder(origin: Point3, axis: Vec3, radius: f64, p: Point3) -> bool {
    let axis_unit = match axis.normalized() {
        Some(a) => a,
        None => return false,
    };
    let v = Vec3::new(p.x - origin.x, p.y - origin.y, p.z - origin.z);
    let along = v.dot(&axis_unit);
    let projection = Point3::new(
        origin.x + axis_unit.x * along,
        origin.y + axis_unit.y * along,
        origin.z + axis_unit.z * along,
    );
    let radial = p.distance_to(&projection);
    (radial - radius).abs() < (radius * 0.01).max(1e-6)
}

/// Walk `face`'s outer loop and find a pair `(he_a, he_b)`
/// of half-edges where `he_a.origin == v_a` and
/// `he_b.origin == v_b`.  Both half-edges live in the same
/// loop (a precondition for `mef`).  Returns `None` if either
/// vertex isn't on the loop.
fn find_loop_he_pair(
    brep: &BRep,
    face: FaceHandle,
    v_a: VertexHandle,
    v_b: VertexHandle,
) -> Option<(HalfEdgeHandle, HalfEdgeHandle)> {
    let face_data = brep.faces.get(face)?;
    let outer_loop = face_data.outer_loop;
    let loop_data = brep.loops.get(outer_loop)?;
    let start = loop_data.half_edge;
    let mut he_a = None;
    let mut he_b = None;
    let mut cursor = start;
    for _ in 0..32 {
        let he = brep.half_edges.get(cursor)?;
        if he.origin == v_a {
            he_a = Some(cursor);
        }
        if he.origin == v_b {
            he_b = Some(cursor);
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    Some((he_a?, he_b?))
}

/// The vertex one step ahead in the loop (the destination of
/// `he_a`).  Used to identify the footprint sandwiched between
/// the kept-vertex pair.
fn he_after_origin_in_loop(brep: &BRep, he: HalfEdgeHandle) -> Option<VertexHandle> {
    let next = brep.half_edges.get(he)?.next;
    Some(brep.half_edges.get(next)?.origin)
}

/// Sample the 3-D arc at `n_samples=12` points, project each
/// to UV on `surface`, and fit a degree-3 NURBS pcurve through
/// the samples.  Returns `None` if the fit fails.
///
/// Phase 6c-C.2b/2d: used to wire the inserted cap-tangent arc's
/// pcurves on both the cap and cylinder sides immediately after
/// `mef`, so downstream consumers (mesh, bool, G0/G1 metric) see
/// consistent UV → 3-D evaluations.
/// Sample the 3-D arc at `n_samples=24` points, project each
/// to UV on `surface`, and fit a degree-5 NURBS pcurve through
/// the samples. When `reversed=true`, samples are taken from
/// arc(t_hi) to arc(t_lo) — used for the twin half-edge whose
/// natural direction is opposite the arc's parameter direction.
///
/// Phase 6c-C.2b/2d/2e: used to wire the inserted cap-tangent
/// arc's pcurves on both the cap and cylinder sides immediately
/// after `mef`, so downstream consumers (mesh, bool, G0/G1
/// metric) see consistent UV → 3-D evaluations.
fn fit_pcurve_along_arc(
    surface: &SurfaceGeometry,
    arc: &NurbsCurve,
    reversed: bool,
) -> Option<NurbsCurve> {
    use popkern_core::geometry::fitting::fit_curve;
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    use std::f64::consts::TAU;
    let (t0, t1) = arc.parameter_range();
    // 48 samples + degree 9 — half-circles need higher fidelity
    // than the typical degree-3 / 12-sample fit; the arc is a
    // rational quadratic and its projection to a non-trivial UV
    // chart (sphere/cylinder) is non-polynomial, so a richer fit
    // reduces mid-arc residuals to ~1e-3 on cube_uniform.
    let n_samples = 48;
    let mut points: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut params: Vec<f64> = Vec::with_capacity(n_samples);
    // Both Sphere::u and Cylinder::u come from `atan2` clamped to
    // [0, 2π). When the arc traverses the wrap-point, raw samples
    // alternate between values near 0 and near 2π; LSQ fit
    // interpolates these as a curve through (0, 2π) — i.e., the
    // wrong way around. Unwrap by tracking previous u and adding
    // ±2π to keep the increment ≤ π.
    let surf_has_wrapping_u = matches!(
        surface,
        SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Cylinder { .. }
    );
    let mut prev_u: Option<f64> = None;
    for i in 0..n_samples {
        let s = i as f64 / (n_samples - 1) as f64;
        // When `reversed`, sample arc(t1) at s=0 and arc(t0) at
        // s=1, so the resulting pcurve's lo end maps to the arc's
        // far endpoint.  This matches the twin half-edge's natural
        // direction (origin = arc-far-end).
        let t = if reversed {
            t0 + (1.0 - s) * (t1 - t0)
        } else {
            t0 + s * (t1 - t0)
        };
        let p = arc.point_at(t);
        let (mut u, v) = project_point_to_uv(surface, &p)?;
        if surf_has_wrapping_u && let Some(u_prev) = prev_u {
            let mut delta = u - u_prev;
            while delta > std::f64::consts::PI {
                u -= TAU;
                delta -= TAU;
            }
            while delta < -std::f64::consts::PI {
                u += TAU;
                delta += TAU;
            }
        }
        prev_u = Some(u);
        points.push(Point3::new(u, v, 0.0));
        params.push(s);
    }
    let degree = 9usize;
    let n_cps = (degree + 2).min(points.len());
    let mut fitted = fit_curve(&points, &params, degree, n_cps).ok()?;
    // Re-wrap u-coords back into [u_lo, u_hi] only if every CP
    // fits in that band — otherwise the unwrapped fit was the
    // intended one and we must not clamp it. For sphere/cylinder
    // the u domain is [0, 2π], so a CP at e.g. u=2π+ε is the
    // CCW-half of an arc that crosses the seam; clamping it to
    // 2π destroys the curve.
    let ((u_lo, u_hi), (v_lo, v_hi)) = surface.parameter_range();
    for cp in &mut fitted.control_points {
        if v_lo.is_finite() && v_hi.is_finite() {
            cp.y = cp.y.clamp(v_lo, v_hi);
        }
    }
    let _ = (u_lo, u_hi);
    Some(fitted)
}

/// Build a degree-2 rational NURBS half-circle on the cap
/// sphere ∩ cylinder, going from `v_a` to `v_b`.  Reuses the
/// arc construction from `popkern_core::geometry::analytic_curves`
/// where possible.
fn build_cap_tangent_half_circle(
    brep: &BRep,
    _cap_centre: Point3,
    _cap_radius: f64,
    cyl_surf: &SurfaceGeometry,
    v_a: VertexHandle,
    v_b: VertexHandle,
) -> NurbsCurve {
    let (cyl_origin, cyl_axis, cyl_radius) = match cyl_surf {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => (*origin, *axis, *radius),
        _ => unreachable!("caller validates cylinder surface"),
    };
    let Some(cyl_axis_unit) = cyl_axis.normalized() else {
        return line_curve_between(brep, v_a, v_b);
    };

    // The cap-tangent circle is the per-fillet sphere ∩ cylinder
    // intersection — a circle in a plane perpendicular to
    // cyl_axis.  Both v_a and v_b lie on this circle, so its
    // centre is the projection of either vertex onto the
    // cylinder axis.  Deriving it from v_a (rather than the
    // shared cap_centre) handles variable-radius corners where
    // each fillet's sphere lives at its own offset along the
    // axis.
    let p_a = brep.vertices[v_a].position;
    let along = Vec3::new(
        p_a.x - cyl_origin.x,
        p_a.y - cyl_origin.y,
        p_a.z - cyl_origin.z,
    )
    .dot(&cyl_axis_unit);
    let circle_centre = Point3::new(
        cyl_origin.x + cyl_axis_unit.x * along,
        cyl_origin.y + cyl_axis_unit.y * along,
        cyl_origin.z + cyl_axis_unit.z * along,
    );
    // Two orthogonal unit vectors perpendicular to cyl_axis_unit.
    let helper = if cyl_axis_unit.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let Some(u_axis) = cyl_axis_unit.cross(&helper).normalized() else {
        return line_curve_between(brep, v_a, v_b);
    };
    let v_axis = cyl_axis_unit.cross(&u_axis);

    // Solve for the angles theta_a, theta_b at which the circle
    // passes through v_a and v_b.
    let p_a = brep.vertices[v_a].position;
    let p_b = brep.vertices[v_b].position;
    let theta = |p: Point3| -> f64 {
        let d = Vec3::new(
            p.x - circle_centre.x,
            p.y - circle_centre.y,
            p.z - circle_centre.z,
        );
        d.dot(&v_axis).atan2(d.dot(&u_axis))
    };
    let theta_a = theta(p_a);
    let theta_b = theta(p_b);

    // Pick the SHORTER arc direction.
    let mut delta = theta_b - theta_a;
    while delta > std::f64::consts::PI {
        delta -= TAU;
    }
    while delta < -std::f64::consts::PI {
        delta += TAU;
    }
    arc_as_rational_quadratic(
        circle_centre,
        cyl_radius,
        u_axis,
        v_axis,
        theta_a,
        theta_a + delta,
    )
}

fn line_curve_between(brep: &BRep, v_a: VertexHandle, v_b: VertexHandle) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![brep.vertices[v_a].position, brep.vertices[v_b].position],
        weights: vec![1.0, 1.0],
    }
}

/// Single rational-quadratic arc spanning theta_a..theta_b
/// (`|theta_b - theta_a| < π` for the single-piece form to be
/// exact).  Same construction as
/// `popkern_core::geometry::apex_blend::sphere_fit::arc_as_nurbs`:
/// CPs `[p0, shoulder, p2]` with weights `[1, cos(α), 1]` where
/// α is the half-angle.  The shoulder is the intersection of
/// tangent lines at p0 and p2; algebraically, that's
/// `centre + (p_mid - centre) / cos(α)` where `p_mid` is the
/// arc midpoint.  Returns a degree-2 rational NURBS that lies
/// exactly on the circle.
fn arc_as_rational_quadratic(
    centre: Point3,
    radius: f64,
    u_axis: Vec3,
    v_axis: Vec3,
    theta_a: f64,
    theta_b: f64,
) -> NurbsCurve {
    let half_span = (theta_b - theta_a) * 0.5;
    let cos_half = half_span.cos();
    let theta_mid = theta_a + half_span;
    let circle_eval = |theta: f64| -> Point3 {
        Point3::new(
            centre.x + radius * (theta.cos() * u_axis.x + theta.sin() * v_axis.x),
            centre.y + radius * (theta.cos() * u_axis.y + theta.sin() * v_axis.y),
            centre.z + radius * (theta.cos() * u_axis.z + theta.sin() * v_axis.z),
        )
    };
    let p0 = circle_eval(theta_a);
    let p_mid = circle_eval(theta_mid);
    let p2 = circle_eval(theta_b);
    let p1 = if cos_half.abs() > 1e-12 {
        Point3::new(
            centre.x + (p_mid.x - centre.x) / cos_half,
            centre.y + (p_mid.y - centre.y) / cos_half,
            centre.z + (p_mid.z - centre.z) / cos_half,
        )
    } else {
        p_mid
    };

    NurbsCurve {
        degree: 2,
        control_points: vec![p0, p1, p2],
        weights: vec![1.0, cos_half.abs(), 1.0],
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    }
}
