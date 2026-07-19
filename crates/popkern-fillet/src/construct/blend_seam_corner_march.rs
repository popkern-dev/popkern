//! PK-152 stage 2: mismatched-radius pair-corner seam via the
//! marching intersector.
//!
//! With unequal radii the strip×strip intersection leaves the
//! bisector-plane family (a space quartic), and the corner topology
//! differs structurally from the equal-radius case: the near anchor
//! sits on ONE strip only, joined to the larger strip's terminal
//! circle by a Δr mismatch stub on the shared support edge. The true
//! closure:
//!
//! - the seam runs from the far anchor (rail crossing, on both
//!   strips) to the point W where the intersection curve crosses the
//!   stub's support plane — W lies on the smaller strip's rail line
//!   by construction (rail = strip ∩ support plane);
//! - the near anchor MOVES to W (its rail edge shortens along its
//!   own line, no topology change);
//! - the mismatch stub becomes the larger strip's exact
//!   cylinder×plane ellipse arc from W to its retained endpoint.
//!
//! The march is exact per sample: parameterized by the near strip's
//! chart angle, the far strip's constraint is a quadratic in the
//! axial coordinate, so there is no integrator drift; the exit is
//! bisected onto the support plane. The seam 3-D curve and chart
//! pcurves interpolate the marched samples (degree 3, endpoints
//! pinned); stub and rail re-authoring is analytic. Same
//! transactional/decline contract as stage 1.

use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
use popkern_core::geometry::fitting::interpolate_curve;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveEval, CurveGeometry, CurveOnEdge, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::types::BlendSeamCornerDeclineReason as Reason;

use super::blend_seam_corner::{PairCornerCandidate, off_cylinder, strip_cylinder};
use super::blend_seam_corner_repair::{
    commit_edge_authoring, neighbour_chart_point, perform_flap_surgery,
};
use super::terminal_arc_trim::StripCylinder;

/// Support plane extracted from a face.
struct SupportPlane {
    origin: Point3,
    normal: Vec3,
    u_axis: Vec3,
    v_axis: Vec3,
}

fn face_plane(brep: &BRep, face: FaceHandle) -> Option<SupportPlane> {
    match brep.faces.get(face)?.surface.as_ref().map(|s| &s.geometry) {
        Some(SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        }) => Some(SupportPlane {
            origin: *origin,
            normal: *normal,
            u_axis: *u_axis,
            v_axis: *v_axis,
        }),
        _ => None,
    }
}

/// One edge incident to the near anchor, classified by its flanks.
struct AnchorEdge {
    edge: EdgeHandle,
    /// The flank that is NOT a strip of this corner.
    support: FaceHandle,
    /// The endpoint that is not the near anchor.
    far_vertex: VertexHandle,
}

/// March the strip×strip intersection from `start` (on both
/// cylinders) until it crosses the support plane; returns exact
/// per-sample points from `start` to the bisected exit point.
fn march_seam(
    cyl_on: &StripCylinder,
    cyl_off: &StripCylinder,
    start: Point3,
    toward: Point3,
    plane: &SupportPlane,
    tolerance: &Tolerance,
) -> Result<Vec<Point3>, Reason> {
    let (p1, p2) = orthonormal_basis(cyl_on.axis);
    let r = cyl_on.radius;
    let a_on = cyl_on.axis;
    let a_off = cyl_off.axis;
    let c = a_on.dot(&a_off);
    let cross_len = a_on.cross(&a_off).length();
    if cross_len < 1e-6 {
        return Err(Reason::ParallelAxes);
    }
    let point_on = |u: f64, t: f64| -> Point3 {
        cyl_on.origin + p1 * (r * u.cos()) + p2 * (r * u.sin()) + a_on * t
    };
    // The strip_off membership constraint is quadratic in the axial
    // coordinate t at fixed chart angle u.
    let solve_roots = |u: f64| -> Option<(f64, f64)> {
        let base = cyl_on.origin + p1 * (r * u.cos()) + p2 * (r * u.sin());
        let d0 = base - cyl_off.origin;
        let a_ = 1.0 - c * c;
        let b_ = d0.dot(&a_on) - c * d0.dot(&a_off);
        let c0 = d0.dot(&d0) - d0.dot(&a_off).powi(2) - cyl_off.radius * cyl_off.radius;
        let disc = b_ * b_ - a_ * c0;
        if disc < 0.0 {
            return None;
        }
        let sq = disc.sqrt();
        Some(((-b_ + sq) / a_, (-b_ - sq) / a_))
    };
    let solve_t = |u: f64, t_prev: f64| -> Option<f64> {
        let (t1, t2) = solve_roots(u)?;
        Some(if (t1 - t_prev).abs() <= (t2 - t_prev).abs() {
            t1
        } else {
            t2
        })
    };
    // The start anchor is the SINGULAR point of the intersection
    // quartic (both strips tangent to the shared support there): the
    // first step's two roots are symmetric around the start station
    // and continuity cannot disambiguate them. Resolve the first
    // step toward the corner region instead.
    let solve_t_first = |u: f64| -> Option<f64> {
        let (t1, t2) = solve_roots(u)?;
        let base = |t: f64| cyl_on.origin + p1 * (r * u.cos()) + p2 * (r * u.sin()) + a_on * t;
        Some(
            if (base(t1) - toward).length() <= (base(t2) - toward).length() {
                t1
            } else {
                t2
            },
        )
    };
    let d_start = start - cyl_on.origin;
    let radial_start = d_start - a_on * d_start.dot(&a_on);
    let u0 = radial_start.dot(&p2).atan2(radial_start.dot(&p1));
    let t0 = d_start.dot(&a_on);
    let s_of = |p: Point3| -> f64 { (p - plane.origin).dot(&plane.normal) };
    let s_start = s_of(start);
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    if s_start.abs() <= pos_tol {
        return Err(Reason::MarchUnresolved);
    }
    let side = s_start.signum();
    // Pick the march direction that closes on the support plane.
    let step = 0.02_f64;
    let mut dir = 0.0;
    for cand in [1.0_f64, -1.0] {
        if let Some(t) = solve_t(u0 + cand * step, t0) {
            let p = point_on(u0 + cand * step, t);
            if s_of(p) * side < s_start * side {
                dir = cand;
                break;
            }
        }
    }
    if dir == 0.0 {
        return Err(Reason::MarchUnresolved);
    }
    // Exit angle. A fillet SUPPORT plane is tangent to the strip
    // along the rail generator, so the seam meets it tangentially —
    // there is no sign change to bisect. The tangency angle is
    // analytic: the radial direction from the axis toward the plane.
    // A transversal plane (axis closer than r) keeps the
    // sign-crossing exit instead.
    let n_ax = plane.normal - a_on * plane.normal.dot(&a_on);
    let s_axis = (cyl_on.origin - plane.origin).dot(&plane.normal);
    let n_perp = n_ax.normalized();
    let ue = match n_perp {
        Some(n_perp) if (s_axis.abs() - r).abs() <= pos_tol.max(r * 1e-6) => {
            // Tangent support: exit at the rail generator's angle,
            // unwrapped in the march direction.
            let d_t = -n_perp * s_axis.signum();
            let u_rail = d_t.dot(&p2).atan2(d_t.dot(&p1));
            let mut du = (u_rail - u0) % std::f64::consts::TAU;
            if du < 0.0 {
                du += std::f64::consts::TAU;
            }
            if dir < 0.0 {
                du -= std::f64::consts::TAU;
            }
            if du.abs() < 1e-9 {
                return Err(Reason::MarchUnresolved);
            }
            u0 + du
        }
        _ => {
            // Transversal plane: step until the seam crosses it,
            // then bisect.
            let mut u_prev = u0;
            let mut t_prev = t0;
            let mut u_exit = None;
            for k in 1..=400 {
                let u = u0 + dir * step * k as f64;
                let Some(t) = solve_t(u, t_prev) else {
                    return Err(Reason::MarchUnresolved);
                };
                if s_of(point_on(u, t)) * side <= 0.0 {
                    let (mut lo, mut hi) = (u_prev, u);
                    let mut t_seed = t_prev;
                    for _ in 0..80 {
                        let mid = 0.5 * (lo + hi);
                        let Some(tm) = solve_t(mid, t_seed) else {
                            return Err(Reason::MarchUnresolved);
                        };
                        t_seed = tm;
                        if s_of(point_on(mid, tm)) * side <= 0.0 {
                            hi = mid;
                        } else {
                            lo = mid;
                        }
                    }
                    u_exit = Some(0.5 * (lo + hi));
                    break;
                }
                u_prev = u;
                t_prev = t;
            }
            u_exit.ok_or(Reason::MarchUnresolved)?
        }
    };
    // Uniform exact resample from start to exit; the seam must stay
    // on the body side of the plane the whole way.
    let n = 65;
    let mut pts = Vec::with_capacity(n);
    pts.push(start);
    let mut t_seed = t0;
    for k in 1..n {
        let u = u0 + (ue - u0) * k as f64 / (n - 1) as f64;
        let solved = if k == 1 {
            solve_t_first(u)
        } else {
            solve_t(u, t_seed)
        };
        let Some(t) = solved else {
            return Err(Reason::MarchUnresolved);
        };
        t_seed = t;
        let p = point_on(u, t);
        if s_of(p) * side < -pos_tol {
            return Err(Reason::MarchUnresolved);
        }
        pts.push(p);
    }
    // The exit point must land on the plane.
    let w = *pts.last().ok_or(Reason::MarchUnresolved)?;
    if s_of(w).abs() > pos_tol {
        return Err(Reason::MarchUnresolved);
    }
    Ok(pts)
}

/// Chart polyline of the marched points in one strip's chart,
/// sequentially period-unwrapped from an optional anchor.
fn chart_polyline(cyl: &StripCylinder, pts: &[Point3], anchor_u: Option<f64>) -> Vec<(f64, f64)> {
    let (p1, p2) = orthonormal_basis(cyl.axis);
    let mut out = Vec::with_capacity(pts.len());
    let mut prev_u: Option<f64> = anchor_u;
    for p in pts {
        let d = *p - cyl.origin;
        let radial = d - cyl.axis * d.dot(&cyl.axis);
        let mut u = radial.dot(&p2).atan2(radial.dot(&p1));
        if let Some(pu) = prev_u {
            while u - pu > std::f64::consts::PI {
                u -= std::f64::consts::TAU;
            }
            while pu - u > std::f64::consts::PI {
                u += std::f64::consts::TAU;
            }
        }
        prev_u = Some(u);
        out.push((u, d.dot(&cyl.axis)));
    }
    out
}

/// Degree-3 interpolating pcurve through chart samples with the
/// given normalized parameters.
fn chart_pcurve(uv: &[(f64, f64)], params: &[f64]) -> Result<Pcurve, Reason> {
    let pts: Vec<Point3> = uv.iter().map(|(u, v)| Point3::new(*u, *v, 0.0)).collect();
    let curve = interpolate_curve(&pts, params, 3).map_err(|_| Reason::PcurveClosureFailed)?;
    Ok(Pcurve::from_nurbs(curve))
}

/// Chord-length parameters normalized to [0, 1].
fn chord_params(pts: &[Point3]) -> Vec<f64> {
    let mut acc = vec![0.0_f64];
    for w in pts.windows(2) {
        let d = (w[1] - w[0]).length();
        acc.push(acc.last().copied().unwrap_or(0.0) + d.max(1e-15));
    }
    let total = *acc.last().unwrap_or(&1.0);
    acc.iter().map(|a| a / total).collect()
}

/// Exact cylinder×plane ellipse arc between two points on both
/// surfaces, using the conjugate-diameter form (chart-angle
/// parameterized, short branch).
fn cylinder_plane_arc(
    cyl: &StripCylinder,
    plane: &SupportPlane,
    from: Point3,
    to: Point3,
    tolerance: &Tolerance,
) -> Result<(CurveOnEdge, f64, f64), Reason> {
    let (p1, p2) = orthonormal_basis(cyl.axis);
    let denom = plane.normal.dot(&cyl.axis);
    if denom.abs() < 1e-9 {
        return Err(Reason::StubUnresolved);
    }
    let t_c = plane.normal.dot(&(plane.origin - cyl.origin)) / denom;
    let center = cyl.origin + cyl.axis * t_c;
    let u_basis = (p1 - cyl.axis * (plane.normal.dot(&p1) / denom)) * cyl.radius;
    let v_basis = (p2 - cyl.axis * (plane.normal.dot(&p2) / denom)) * cyl.radius;
    let angle_of = |p: Point3| -> f64 {
        let d = p - cyl.origin;
        let radial = d - cyl.axis * d.dot(&cyl.axis);
        radial.dot(&p2).atan2(radial.dot(&p1))
    };
    let th_from = angle_of(from);
    let mut th_to = angle_of(to);
    while th_to - th_from > std::f64::consts::PI {
        th_to -= std::f64::consts::TAU;
    }
    while th_from - th_to > std::f64::consts::PI {
        th_to += std::f64::consts::TAU;
    }
    if (th_to - th_from).abs() < 1e-9 {
        return Err(Reason::StubUnresolved);
    }
    let arc = CurveOnEdge {
        geometry: CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            start_angle: th_from,
            end_angle: th_to,
        },
        t_start: th_from,
        t_end: th_to,
    };
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    if (arc.geometry.point_at(th_from) - from).length() > pos_tol
        || (arc.geometry.point_at(th_to) - to).length() > pos_tol
    {
        return Err(Reason::StubUnresolved);
    }
    Ok((arc, th_from, th_to))
}

/// Sinusoidal pcurve of a plane in a cylinder chart over
/// `[u_from, u_to]` (see stage 1: canonical plane×cylinder form).
fn plane_sinusoid(
    cyl: &StripCylinder,
    plane: &SupportPlane,
    u_from: f64,
    u_to: f64,
) -> Result<Pcurve, Reason> {
    let (p1, p2) = orthonormal_basis(cyl.axis);
    let denom = plane.normal.dot(&cyl.axis);
    if denom.abs() < 1e-9 {
        return Err(Reason::StubUnresolved);
    }
    let n = plane.normal;
    let v_const = -n.dot(&(cyl.origin - plane.origin)) / denom;
    let v_cos = -cyl.radius * n.dot(&p1) / denom;
    let v_sin = -cyl.radius * n.dot(&p2) / denom;
    Ok(Pcurve::sinusoidal(0.0, v_const, v_cos, v_sin, u_from, u_to))
}

/// Mismatched-radius pair-corner repair (stage 2).
pub(super) fn repair_pair_corner_mixed(
    brep: &mut BRep,
    c: &PairCornerCandidate,
    tolerance: &Tolerance,
) -> Result<(), Reason> {
    let cyl_n = strip_cylinder(brep, c.strip_n).ok_or(Reason::EulerRejected)?;
    let cyl_p = strip_cylinder(brep, c.strip_p).ok_or(Reason::EulerRejected)?;
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let pos_of =
        |b: &BRep, v: VertexHandle| -> Option<Point3> { Some(b.vertices.get(v)?.position) };
    let a_pos = pos_of(brep, c.anchor_near).ok_or(Reason::EulerRejected)?;
    let b_pos = pos_of(brep, c.anchor_far).ok_or(Reason::EulerRejected)?;

    // The near anchor must sit on exactly one strip.
    let on_n = off_cylinder(a_pos, &cyl_n) <= pos_tol;
    let on_p = off_cylinder(a_pos, &cyl_p) <= pos_tol;
    let (strip_on, cyl_on, strip_off, cyl_off) = match (on_n, on_p) {
        (true, false) => (c.strip_n, cyl_n.clone(), c.strip_p, cyl_p.clone()),
        (false, true) => (c.strip_p, cyl_p.clone(), c.strip_n, cyl_n.clone()),
        _ => return Err(Reason::StubUnresolved),
    };

    // Classify the near anchor's other incident edges: the rail
    // (strip_on × support) and the Δr mismatch stub (strip_off ×
    // support), sharing one support face.
    let mut rail: Option<AnchorEdge> = None;
    let mut stub: Option<AnchorEdge> = None;
    let mut extra = 0usize;
    for (eh, e) in brep.edges.iter() {
        if eh == c.seam_edge {
            continue;
        }
        let Some(he) = brep.half_edges.get(e.half_edge) else {
            continue;
        };
        let Some(tw) = brep.half_edges.get(he.twin) else {
            continue;
        };
        if he.origin != c.anchor_near && tw.origin != c.anchor_near {
            continue;
        }
        let far_vertex = if he.origin == c.anchor_near {
            tw.origin
        } else {
            he.origin
        };
        let (Some(lf), Some(lt)) = (brep.loops.get(he.loop_), brep.loops.get(tw.loop_)) else {
            continue;
        };
        let flanks = [lf.face, lt.face];
        let entry = |support: FaceHandle| AnchorEdge {
            edge: eh,
            support,
            far_vertex,
        };
        if flanks.contains(&strip_on) && !flanks.contains(&strip_off) {
            let support = if flanks[0] == strip_on {
                flanks[1]
            } else {
                flanks[0]
            };
            rail = Some(entry(support));
        } else if flanks.contains(&strip_off) && !flanks.contains(&strip_on) {
            let support = if flanks[0] == strip_off {
                flanks[1]
            } else {
                flanks[0]
            };
            stub = Some(entry(support));
        } else {
            extra += 1;
        }
    }
    let rail = rail.ok_or(Reason::StubUnresolved)?;
    let stub = stub.ok_or(Reason::StubUnresolved)?;
    if extra != 0 || rail.support != stub.support {
        return Err(Reason::StubUnresolved);
    }
    let plane = face_plane(brep, stub.support).ok_or(Reason::StubUnresolved)?;
    let stub_far_pos = pos_of(brep, stub.far_vertex).ok_or(Reason::StubUnresolved)?;
    if off_cylinder(stub_far_pos, &cyl_off) > pos_tol
        || (stub_far_pos - plane.origin).dot(&plane.normal).abs() > pos_tol
    {
        return Err(Reason::StubUnresolved);
    }
    let rail_far_pos = pos_of(brep, rail.far_vertex).ok_or(Reason::StubUnresolved)?;

    // March the seam from the far anchor to the support-plane exit.
    let corner_pos = pos_of(brep, c.corner).unwrap_or(a_pos);
    let pts = march_seam(&cyl_on, &cyl_off, b_pos, corner_pos, &plane, tolerance)?;
    let w = *pts.last().ok_or(Reason::MarchUnresolved)?;
    // W replaces the near anchor; sanity: it stays within the corner
    // neighbourhood.
    if (w - a_pos).length() > cyl_off.radius + cyl_on.radius {
        return Err(Reason::MarchUnresolved);
    }

    // Author everything BEFORE surgery where possible.
    let params = chord_params(&pts);
    let seam_nurbs = interpolate_curve(&pts, &params, 3).map_err(|_| Reason::MarchUnresolved)?;
    let primary_seam = brep
        .edges
        .get(c.seam_edge)
        .ok_or(Reason::EulerRejected)?
        .half_edge;
    // After surgery he_n runs anchor_near(=W) → anchor_far; the
    // marched points run B → W.
    let seam_curve = {
        let (t0, t1) = {
            use popkern_core::geometry::CurveEval as _;
            seam_nurbs.parameter_range()
        };
        if primary_seam == c.he_n {
            // Curve must run W → B: reverse the interpolation input.
            let mut rev = pts.clone();
            rev.reverse();
            let curve = interpolate_curve(&rev, &params_reversed(&params), 3)
                .map_err(|_| Reason::MarchUnresolved)?;
            let (r0, r1) = curve.parameter_range();
            CurveOnEdge {
                geometry: CurveGeometry::Nurbs(curve),
                t_start: r0,
                t_end: r1,
            }
        } else {
            CurveOnEdge {
                geometry: CurveGeometry::Nurbs(seam_nurbs.clone()),
                t_start: t0,
                t_end: t1,
            }
        }
    };

    // Stub arc: strip_off × support plane, W → retained endpoint.
    let (stub_arc_from, stub_arc_to) = (w, stub_far_pos);
    let primary_stub = brep
        .edges
        .get(stub.edge)
        .ok_or(Reason::StubUnresolved)?
        .half_edge;
    let stub_primary_he = brep
        .half_edges
        .get(primary_stub)
        .ok_or(Reason::StubUnresolved)?;
    let stub_primary_from_anchor = stub_primary_he.origin == c.anchor_near;
    let (sp_from, sp_to) = if stub_primary_from_anchor {
        (stub_arc_from, stub_arc_to)
    } else {
        (stub_arc_to, stub_arc_from)
    };
    let (stub_curve, _, _) = cylinder_plane_arc(&cyl_off, &plane, sp_from, sp_to, tolerance)?;

    perform_flap_surgery(brep, c)?;
    if let Some(v) = brep.vertices.get_mut(c.anchor_near) {
        v.position = w;
    }

    // Seam pcurves: chart polylines anchored at the far-anchor end.
    let seam_staged = {
        let mut staged: Vec<(FaceHandle, HalfEdgeHandle, Pcurve)> = Vec::with_capacity(2);
        for (face, cyl, he) in [(c.strip_n, &cyl_n, c.he_n), (c.strip_p, &cyl_p, c.he_p)] {
            let h = brep.half_edges.get(he).ok_or(Reason::EulerRejected)?;
            // The far anchor is he_n's destination and he_p's origin.
            let b_end_anchor = if he == c.he_n {
                neighbour_chart_point(brep, face, h.next, false)
            } else {
                neighbour_chart_point(brep, face, h.prev, true)
            };
            let uv = chart_polyline(cyl, &pts, b_end_anchor.map(|(u, _)| u));
            if let Some((au, av)) = b_end_anchor {
                let (u0c, v0c) = uv[0];
                let chart_tol = (tolerance.linear() * 100.0).max(1e-6);
                if ((u0c - au).abs() * cyl.radius).max((v0c - av).abs()) > chart_tol {
                    return Err(Reason::PcurveClosureFailed);
                }
            }
            // Orient along the half-edge: he_n runs W→B (reverse of
            // march order), he_p runs B→W (march order).
            let pcurve = if he == c.he_n {
                let mut rev = uv.clone();
                rev.reverse();
                chart_pcurve(&rev, &params_reversed(&params))?
            } else {
                chart_pcurve(&uv, &params)?
            };
            staged.push((face, he, pcurve));
        }
        staged
    };
    commit_edge_authoring(brep, c.seam_edge, seam_curve, seam_staged);

    // Stub: strip_off-side sinusoid + support-side pcurve.
    let stub_staged = author_stub_pcurves(
        brep,
        &stub,
        &stub_curve,
        strip_off,
        &cyl_off,
        &plane,
        tolerance,
    )?;
    commit_edge_authoring(brep, stub.edge, stub_curve, stub_staged);

    // Rail: same line, shortened to W.
    author_rail(
        brep,
        &rail,
        c.anchor_near,
        rail_far_pos,
        w,
        strip_on,
        &cyl_on,
        &plane,
        tolerance,
    )?;

    // Postconditions.
    for eh in [c.seam_edge, stub.edge, rail.edge] {
        if brep.measure_edge_deviation(eh) > (tolerance.linear() * 100.0).max(1e-7) {
            return Err(Reason::ValidationFailed);
        }
    }
    brep.validate().map_err(|_| Reason::ValidationFailed)?;
    Ok(())
}

fn params_reversed(params: &[f64]) -> Vec<f64> {
    let mut rev: Vec<f64> = params.iter().map(|t| 1.0 - t).collect();
    rev.reverse();
    rev
}

/// Author the stub edge's pcurves on the larger strip (sinusoid) and
/// the support plane (projected endpoints; elliptical when the chart
/// image is axis-orthogonal, sampled otherwise).
#[allow(clippy::too_many_arguments)]
fn author_stub_pcurves(
    brep: &BRep,
    stub: &AnchorEdge,
    arc: &CurveOnEdge,
    strip_off: FaceHandle,
    cyl_off: &StripCylinder,
    plane: &SupportPlane,
    tolerance: &Tolerance,
) -> Result<Vec<(FaceHandle, HalfEdgeHandle, Pcurve)>, Reason> {
    let e = brep.edges.get(stub.edge).ok_or(Reason::StubUnresolved)?;
    let he = brep
        .half_edges
        .get(e.half_edge)
        .ok_or(Reason::StubUnresolved)?;
    let tw = brep.half_edges.get(he.twin).ok_or(Reason::StubUnresolved)?;
    let curve = arc;
    let mut staged = Vec::with_capacity(2);
    for (heh, hh) in [(e.half_edge, he), (he.twin, tw)] {
        let face = brep
            .loops
            .get(hh.loop_)
            .map(|l| l.face)
            .ok_or(Reason::StubUnresolved)?;
        // Sample the committed arc curve along this half-edge's
        // direction and project into the face chart.
        let (t0, t1) = (curve.t_start, curve.t_end);
        let forward = e.half_edge == heh;
        let n = 17;
        let mut samples3d = Vec::with_capacity(n);
        for k in 0..n {
            let f = k as f64 / (n - 1) as f64;
            let f = if forward { f } else { 1.0 - f };
            samples3d.push(curve.geometry.point_at(t0 + (t1 - t0) * f));
        }
        let pcurve = if face == strip_off {
            let uv = chart_polyline(cyl_off, &samples3d, None);
            plane_sinusoid(
                cyl_off,
                plane,
                uv.first().ok_or(Reason::StubUnresolved)?.0,
                uv.last().ok_or(Reason::StubUnresolved)?.0,
            )?
        } else {
            // Support-plane chart: project the samples exactly.
            let uv: Vec<Point3> = samples3d
                .iter()
                .map(|p| {
                    let d = *p - plane.origin;
                    Point3::new(d.dot(&plane.u_axis), d.dot(&plane.v_axis), 0.0)
                })
                .collect();
            let params: Vec<f64> = (0..n).map(|k| k as f64 / (n - 1) as f64).collect();
            let fit = interpolate_curve(&uv, &params, 3).map_err(|_| Reason::StubUnresolved)?;
            Pcurve::from_nurbs(fit)
        };
        // Sample gate: the pcurve trace must reproduce the arc.
        let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
        let surface = brep
            .faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| s.geometry.clone())
            .ok_or(Reason::StubUnresolved)?;
        use popkern_core::geometry::SurfaceEval;
        let (pt0, pt1) = pcurve.parameter_range();
        for k in 0..=8 {
            let t = pt0 + (pt1 - pt0) * k as f64 / 8.0;
            let uvp = pcurve.point_at(t);
            let p3 = surface.point_at(uvp.x, uvp.y);
            let f = k as f64 / 8.0;
            let expect = samples3d[((f * (n - 1) as f64).round() as usize).min(n - 1)];
            if (p3 - expect).length() > pos_tol.max(1e-7) {
                return Err(Reason::StubUnresolved);
            }
        }
        staged.push((face, heh, pcurve));
    }
    Ok(staged)
}

/// Re-author the rail edge after the near anchor moved along its own
/// line: same line curve, Linear pcurves refreshed on both flanks.
#[allow(clippy::too_many_arguments)]
fn author_rail(
    brep: &mut BRep,
    rail: &AnchorEdge,
    anchor: VertexHandle,
    far_pos: Point3,
    w: Point3,
    strip_on: FaceHandle,
    cyl_on: &StripCylinder,
    plane: &SupportPlane,
    tolerance: &Tolerance,
) -> Result<(), Reason> {
    let e = brep.edges.get(rail.edge).ok_or(Reason::StubUnresolved)?;
    let primary = e.half_edge;
    let he = brep.half_edges.get(primary).ok_or(Reason::StubUnresolved)?;
    let tw = brep.half_edges.get(he.twin).ok_or(Reason::StubUnresolved)?;
    let primary_from_anchor = he.origin == anchor;
    let (p_start, p_end) = if primary_from_anchor {
        (w, far_pos)
    } else {
        (far_pos, w)
    };
    let line = CurveOnEdge::line_segment(p_start, p_end, tolerance.linear())
        .ok_or(Reason::StubUnresolved)?;
    let mut staged: Vec<(FaceHandle, HalfEdgeHandle, Pcurve)> = Vec::with_capacity(2);
    for (heh, hh) in [(primary, he), (he.twin, tw)] {
        let face = brep
            .loops
            .get(hh.loop_)
            .map(|l| l.face)
            .ok_or(Reason::StubUnresolved)?;
        let (from3, to3) = if hh.origin == anchor {
            (w, far_pos)
        } else {
            (far_pos, w)
        };
        let uv_of = |p: Point3| -> (f64, f64) {
            if face == strip_on {
                let d = p - cyl_on.origin;
                let radial = d - cyl_on.axis * d.dot(&cyl_on.axis);
                (
                    radial
                        .dot(&orthonormal_basis(cyl_on.axis).1)
                        .atan2(radial.dot(&orthonormal_basis(cyl_on.axis).0)),
                    d.dot(&cyl_on.axis),
                )
            } else {
                let d = p - plane.origin;
                (d.dot(&plane.u_axis), d.dot(&plane.v_axis))
            }
        };
        let mut from_uv = uv_of(from3);
        let to_uv = uv_of(to3);
        if face == strip_on {
            // Rail is an axial generator: both endpoints share one
            // chart angle; unwrap the moved end onto the far end's
            // branch and pin it exactly (iso-u by construction).
            let mut u = from_uv.0;
            while u - to_uv.0 > std::f64::consts::PI {
                u -= std::f64::consts::TAU;
            }
            while to_uv.0 - u > std::f64::consts::PI {
                u += std::f64::consts::TAU;
            }
            if (u - to_uv.0).abs() * cyl_on.radius > (tolerance.linear() * 100.0).max(1e-6) {
                return Err(Reason::StubUnresolved);
            }
            from_uv = (to_uv.0, from_uv.1);
        }
        staged.push((face, heh, Pcurve::linear(from_uv, to_uv, 0.0, 1.0)));
    }
    commit_edge_authoring(brep, rail.edge, line, staged);
    Ok(())
}
