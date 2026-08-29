//! PK-152 blend×blend corner-seam analytic geometry.
//!
//! Two equal-radius cylinder strips with intersecting, non-parallel
//! axes meet in a pair of planar ellipses (the bisector planes of
//! the axes). The corner seam is the short arc of one ellipse
//! between the two trusted rail-crossing anchors. Everything here is
//! exact by construction: the 3-D curve is a `CurveGeometry::Ellipse`
//! in the resolved bisector plane, and each strip-side pcurve is a
//! `PcurveGeometry::Sinusoidal` — the canonical exact form of a
//! plane×cylinder intersection in cylinder UV. Pure computation;
//! nothing here mutates the BRep.

use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveEval, CurveGeometry, CurveOnEdge, Point3, SurfaceEval, Vec3};
use popkern_core::tolerance::Tolerance;

use crate::types::BlendSeamCornerDeclineReason as Reason;

use super::terminal_arc_trim::StripCylinder;

/// Resolved bisector-plane frame of an equal-radius cylinder pair.
pub(super) struct SeamFrame {
    /// Ellipse centre: the axes' intersection point.
    pub(super) center: Point3,
    /// Unit normal of the resolved bisector plane.
    pub(super) normal: Vec3,
    /// Unit in-plane direction of the major semi-axis.
    pub(super) major_dir: Vec3,
    /// Unit in-plane direction of the minor semi-axis
    /// (perpendicular to both cylinder axes).
    pub(super) minor_dir: Vec3,
    /// Major semi-radius `r / √(1 − c²)` with `c = major_dir·axis`.
    pub(super) semi_major: f64,
    /// Minor semi-radius = the shared strip radius.
    pub(super) radius: f64,
}

/// Solve the seam frame for an equal-radius pair, resolving the
/// bisector branch that contains both trusted anchors.
pub(super) fn solve_equal_radius_seam_frame(
    cyl_a: &StripCylinder,
    cyl_b: &StripCylinder,
    anchor_1: Point3,
    anchor_2: Point3,
    tolerance: &Tolerance,
) -> Result<SeamFrame, Reason> {
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let radius = cyl_a.radius;
    if (cyl_a.radius - cyl_b.radius).abs() > pos_tol.max(radius * 1e-9) {
        return Err(Reason::MismatchedRadii);
    }
    let d1 = cyl_a.axis.normalized().ok_or(Reason::ParallelAxes)?;
    let d2 = cyl_b.axis.normalized().ok_or(Reason::ParallelAxes)?;
    let cross = d1.cross(&d2);
    if cross.length() < 1e-6 {
        return Err(Reason::ParallelAxes);
    }
    // Closest points between the two axis lines; they must meet.
    let w0 = cyl_a.origin - cyl_b.origin;
    let b = d1.dot(&d2);
    let denom = 1.0 - b * b;
    let d = d1.dot(&w0);
    let e = d2.dot(&w0);
    let s = (b * e - d) / denom;
    let t = (e - b * d) / denom;
    let p1 = cyl_a.origin + d1 * s;
    let p2 = cyl_b.origin + d2 * t;
    if (p1 - p2).length() > pos_tol {
        return Err(Reason::OffsetAxes);
    }
    let center = p1 + (p2 - p1) * 0.5;

    // The intersection satisfies (q·d1)² = (q·d2)², splitting into
    // the two bisector planes with normals d1∓d2. Resolve the branch
    // by requiring both trusted anchors on it — and not on the other.
    let mut chosen: Option<Vec3> = None;
    for candidate in [d1 - d2, d1 + d2] {
        let Some(n) = candidate.normalized() else {
            continue;
        };
        let resid = (anchor_1 - center)
            .dot(&n)
            .abs()
            .max((anchor_2 - center).dot(&n).abs());
        if resid <= pos_tol {
            if chosen.is_some() {
                return Err(Reason::BranchUnresolved);
            }
            chosen = Some(n);
        }
    }
    let normal = chosen.ok_or(Reason::BranchUnresolved)?;
    let minor_dir = cross.normalized().ok_or(Reason::ParallelAxes)?;
    let major_dir = normal
        .cross(&minor_dir)
        .normalized()
        .ok_or(Reason::BranchUnresolved)?;
    let c = major_dir.dot(&d1);
    let denom_sq = 1.0 - c * c;
    if denom_sq < 1e-12 {
        return Err(Reason::BranchUnresolved);
    }
    Ok(SeamFrame {
        center,
        normal,
        major_dir,
        minor_dir,
        semi_major: radius / denom_sq.sqrt(),
        radius,
    })
}

/// Ellipse angle of a point that must lie on the seam. Gates both
/// the in-plane residual and the on-ellipse residual.
pub(super) fn seam_angle(
    frame: &SeamFrame,
    p: Point3,
    tolerance: &Tolerance,
) -> Result<f64, Reason> {
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let q = p - frame.center;
    if q.dot(&frame.normal).abs() > pos_tol {
        return Err(Reason::AnchorsOffSeam);
    }
    let alpha = q.dot(&frame.major_dir) / frame.semi_major;
    let beta = q.dot(&frame.minor_dir) / frame.radius;
    if (alpha * alpha + beta * beta - 1.0).abs() > 1e-6 {
        return Err(Reason::AnchorsOffSeam);
    }
    Ok(beta.atan2(alpha))
}

/// Author the seam 3-D curve from `p_start` to `p_end` along the
/// short ellipse branch nearest the stale corner.
pub(super) fn author_seam_curve(
    frame: &SeamFrame,
    p_start: Point3,
    p_end: Point3,
    stale_corner: Point3,
    tolerance: &Tolerance,
) -> Result<CurveOnEdge, Reason> {
    let mut psi_s = seam_angle(frame, p_start, tolerance)?;
    let psi_e = seam_angle(frame, p_end, tolerance)?;
    let mut sweep = psi_e - psi_s;
    while sweep > std::f64::consts::PI {
        sweep -= std::f64::consts::TAU;
    }
    while sweep < -std::f64::consts::PI {
        sweep += std::f64::consts::TAU;
    }
    if !(1e-9..=std::f64::consts::PI - 1e-9).contains(&sweep.abs()) {
        return Err(Reason::SweepUnresolved);
    }
    // Canonical ascending sweep: mirror the minor axis if needed.
    let mut minor = frame.minor_dir;
    if sweep < 0.0 {
        minor = -minor;
        psi_s = -psi_s;
        sweep = -sweep;
    }
    let mid = frame.center
        + frame.major_dir * (frame.semi_major * (psi_s + 0.5 * sweep).cos())
        + minor * (frame.radius * (psi_s + 0.5 * sweep).sin());
    // Branch witness: the short seam arc passes the stale corner
    // region; the complementary arc swings to the far side of both
    // strips (≈ 2·semi_major away).
    if (mid - stale_corner).length() > frame.semi_major {
        return Err(Reason::SweepUnresolved);
    }
    let arc = CurveOnEdge {
        geometry: CurveGeometry::Ellipse {
            center: frame.center,
            u_basis: frame.major_dir * frame.semi_major,
            v_basis: minor * frame.radius,
            start_angle: psi_s,
            end_angle: psi_s + sweep,
        },
        t_start: psi_s,
        t_end: psi_s + sweep,
    };
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    if (arc.geometry.point_at(arc.t_start) - p_start).length() > pos_tol
        || (arc.geometry.point_at(arc.t_end) - p_end).length() > pos_tol
    {
        return Err(Reason::SweepUnresolved);
    }
    Ok(arc)
}

/// Chart anchor for one seam endpoint: the adjacent half-edge's
/// pcurve endpoint, when available, is the authoritative chart
/// position of the shared vertex (period branch included).
pub(super) struct ChartAnchor {
    pub(super) uv: Option<(f64, f64)>,
}

/// Author one strip-side seam pcurve as the exact sinusoid of the
/// bisector plane in this cylinder's chart, anchored to the loop
/// neighbours' chart branch and sample-verified against the seam.
pub(super) fn author_side_pcurve(
    cylinder: &SurfaceGeometry,
    frame: &SeamFrame,
    origin_pos: Point3,
    dest_pos: Point3,
    origin_anchor: &ChartAnchor,
    dest_anchor: &ChartAnchor,
    tolerance: &Tolerance,
) -> Result<Pcurve, Reason> {
    let SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    } = cylinder
    else {
        return Err(Reason::PcurveClosureFailed);
    };
    let (perp1, perp2) = orthonormal_basis(*axis);
    let axial = axis.dot(&frame.normal);
    if axial.abs() < 1e-9 {
        return Err(Reason::PcurveClosureFailed);
    }
    // Plane n·(p − C0) = 0 through the cylinder chart
    // p(u, v) = o + r·cos u·perp1 + r·sin u·perp2 + v·axis:
    //   v(u) = −(n·(o − C0) + r cos u (n·perp1) + r sin u (n·perp2)) / (n·axis)
    let n = frame.normal;
    let v_const = -n.dot(&(*origin - frame.center)) / axial;
    let v_cos = -radius * n.dot(&perp1) / axial;
    let v_sin = -radius * n.dot(&perp2) / axial;

    let chart_tol = (tolerance.linear() * 100.0).max(1e-6);
    let u_of = |p: Point3| -> f64 {
        let d = p - *origin;
        let radial = d - *axis * d.dot(axis);
        radial.dot(&perp2).atan2(radial.dot(&perp1))
    };
    let unwrap_near = |mut u: f64, target: f64| -> f64 {
        while u - target > std::f64::consts::PI {
            u -= std::f64::consts::TAU;
        }
        while target - u > std::f64::consts::PI {
            u += std::f64::consts::TAU;
        }
        u
    };
    let v_at = |u: f64| -> f64 { v_const + v_cos * u.cos() + v_sin * u.sin() };
    // Endpoint chart stations, unwrapped to the neighbour's branch
    // and closure-checked against it when a neighbour chart anchor
    // exists (u compared as arc length, v as axial length).
    let resolve = |raw_pos: Point3, anchor: &ChartAnchor, fallback: Option<f64>| -> Option<f64> {
        let raw = u_of(raw_pos);
        match (anchor.uv, fallback) {
            (Some((au, av)), _) => {
                let u = unwrap_near(raw, au);
                if ((u - au).abs() * radius).max((v_at(u) - av).abs()) > chart_tol {
                    return None;
                }
                Some(u)
            }
            (None, Some(near)) => Some(unwrap_near(raw, near)),
            (None, None) => Some(raw),
        }
    };
    let u_origin = resolve(origin_pos, origin_anchor, None).ok_or(Reason::PcurveClosureFailed)?;
    let u_dest =
        resolve(dest_pos, dest_anchor, Some(u_origin)).ok_or(Reason::PcurveClosureFailed)?;
    let du = (u_dest - u_origin).abs();
    if !(1e-9..std::f64::consts::PI).contains(&du) {
        return Err(Reason::PcurveClosureFailed);
    }
    let pcurve = Pcurve::sinusoidal(0.0, v_const, v_cos, v_sin, u_origin, u_dest);
    // Sample gate: the pcurve's surface trace must stay on the seam
    // ellipse and hit both endpoint positions.
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    for k in 0..=8 {
        let t = u_origin + (u_dest - u_origin) * (k as f64 / 8.0);
        let uv = pcurve.point_at(t);
        let p3 = cylinder.point_at(uv.x, uv.y);
        seam_angle(frame, p3, tolerance)?;
        if k == 0 && (p3 - origin_pos).length() > pos_tol {
            return Err(Reason::PcurveClosureFailed);
        }
        if k == 8 && (p3 - dest_pos).length() > pos_tol {
            return Err(Reason::PcurveClosureFailed);
        }
    }
    Ok(pcurve)
}
