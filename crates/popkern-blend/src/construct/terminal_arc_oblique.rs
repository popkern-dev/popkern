//! PK-530: terminal arcs on OBLIQUE planar end faces.
//!
//! A cylinder strip ending on a planar face that is not perpendicular
//! to its axis has an elliptical section: the two anchor feet sit at
//! different axial stations. The circular settlement
//! (`terminal_arc_trim`) declined those ends `AnchorsOffSurface`, and
//! the pipeline served the stale L of straight stubs through the
//! original corner (square-pyramid base edge, r = 0.07: 7.1e-2 off the
//! strip, `fillet_edges` returning `Ok`).
//!
//! The section is exact. In the strip's own chart (angle `u` around
//! `orthonormal_basis(axis)`, `v` along the unit axis), the plane
//! `n·(p − P0) = 0` meets the cylinder where
//! `v(u) = v0 + v_cos·cos u + v_sin·sin u`, so:
//! - the 3-D curve is a `CurveGeometry::Ellipse` in `u` (center on the
//!   axis, bases `r·e_i − axis·r(n·e_i)/(n·axis)`);
//! - the strip pcurve is the `Sinusoidal` pcurve over the same `u`
//!   span;
//! - the cap pcurve is the ellipse's image in the plane chart,
//!   converted from conjugate to principal axes.
//!
//! The corner vertex moves onto the section at the stale corner's
//! angle (the circular rule), and must fall strictly between the feet
//! on the minor branch. Everything is sample-gated before the caller
//! commits. The topology does not change.

use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::{CurveEval, CurveGeometry, CurveOnEdge, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, HalfEdgeHandle};

use crate::types::TerminalArcTrimDeclineReason as Reason;

use super::terminal_arc_author::StagedStubArc;
use super::terminal_arc_trim::StripCylinder;

/// Smallest `|n·axis|` for a closed elliptical section. Below it the
/// end plane runs nearly along the strip and the section degenerates
/// toward a pair of lines.
const MIN_AXIAL_COSINE: f64 = 0.05;

/// The cylinder ∩ end-plane section in the strip chart's angle.
pub(super) struct ObliqueSection {
    origin: Point3,
    axis: Vec3,
    e1: Vec3,
    e2: Vec3,
    radius: f64,
    v0: f64,
    v_cos: f64,
    v_sin: f64,
}

impl ObliqueSection {
    /// The section of `cyl` by the plane through `plane_origin` with
    /// unit normal `normal`.
    pub(super) fn new(cyl: &StripCylinder, plane_origin: Point3, normal: Vec3) -> Option<Self> {
        let axis = cyl.axis.normalized()?;
        let k = normal.dot(&axis);
        if k.abs() < MIN_AXIAL_COSINE {
            return None;
        }
        let (e1, e2) = orthonormal_basis(axis);
        Some(Self {
            origin: cyl.origin,
            axis,
            e1,
            e2,
            radius: cyl.radius,
            v0: normal.dot(&(plane_origin - cyl.origin)) / k,
            v_cos: -cyl.radius * normal.dot(&e1) / k,
            v_sin: -cyl.radius * normal.dot(&e2) / k,
        })
    }

    fn v_at(&self, u: f64) -> f64 {
        self.v0 + self.v_cos * u.cos() + self.v_sin * u.sin()
    }

    fn point(&self, u: f64) -> Point3 {
        self.origin
            + self.e1 * (self.radius * u.cos())
            + self.e2 * (self.radius * u.sin())
            + self.axis * self.v_at(u)
    }

    /// Chart angle of `p`'s radial direction.
    pub(super) fn angle_of(&self, p: Point3) -> f64 {
        let d = p - self.origin;
        d.dot(&self.e2).atan2(d.dot(&self.e1))
    }

    /// Distance from `p` to the section point at its own angle: zero
    /// exactly when `p` lies on both the cylinder and the plane.
    pub(super) fn miss(&self, p: Point3) -> f64 {
        (self.point(self.angle_of(p)) - p).length()
    }

    /// The section point at `p`'s angle.
    pub(super) fn settle(&self, p: Point3) -> Point3 {
        self.point(self.angle_of(p))
    }
}

/// `to` unwrapped onto the 2π branch nearest `from`.
pub(super) fn unwrap_near(to: f64, from: f64) -> f64 {
    to - ((to - from) / std::f64::consts::TAU).round() * std::f64::consts::TAU
}

/// One stub's exact section arc from chart angle `u_start` to `u_end`
/// (the strip half-edge's own direction): its 3-D curve and strip
/// pcurve, sample-gated against the section. Pure computation.
pub(super) struct SectionArc {
    pub(super) edge: popkern_core::topology::EdgeHandle,
    pub(super) twin: HalfEdgeHandle,
    pub(super) curve: CurveOnEdge,
    pub(super) strip_pcurve: Pcurve,
}

pub(super) fn section_arc(
    brep: &BRep,
    he_strip: HalfEdgeHandle,
    section: &ObliqueSection,
    u_start: f64,
    u_end: f64,
    tolerance: &Tolerance,
) -> Result<SectionArc, Reason> {
    let he = brep
        .half_edges
        .get(he_strip)
        .ok_or(Reason::ArcBranchUnresolved)?;
    let (twin, edge) = (he.twin, he.edge);
    let primary = brep
        .edges
        .get(edge)
        .ok_or(Reason::ArcBranchUnresolved)?
        .half_edge;
    let sweep = u_end - u_start;
    if !(1e-9..std::f64::consts::PI).contains(&sweep.abs()) {
        return Err(Reason::ArcBranchUnresolved);
    }
    let center = section.origin + section.axis * section.v0;
    let u_basis = section.e1 * section.radius + section.axis * section.v_cos;
    let v_basis = section.e2 * section.radius + section.axis * section.v_sin;
    // The stored curve follows the primary half-edge; an Ellipse
    // needs an ascending angle, so a descending span mirrors v_basis.
    let (a0, a1) = if primary == he_strip {
        (u_start, u_end)
    } else {
        (u_end, u_start)
    };
    let (v_basis, start_angle, end_angle) = if a1 >= a0 {
        (v_basis, a0, a1)
    } else {
        (-v_basis, -a0, -a1)
    };
    let curve = CurveOnEdge {
        geometry: CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            start_angle,
            end_angle,
        },
        t_start: start_angle,
        t_end: end_angle,
    };
    let strip_pcurve = Pcurve::sinusoidal(
        0.0,
        section.v0,
        section.v_cos,
        section.v_sin,
        u_start,
        u_end,
    );
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let strip_surface = popkern_core::geometry::surface::SurfaceGeometry::Cylinder {
        origin: section.origin,
        axis: section.axis,
        radius: section.radius,
    };
    use popkern_core::geometry::SurfaceEval;
    for k in 0..=8 {
        let f = k as f64 / 8.0;
        let u = u_start + (u_end - u_start) * f;
        let exact = section.point(u);
        let t = if primary == he_strip {
            curve.t_start + (curve.t_end - curve.t_start) * f
        } else {
            curve.t_end - (curve.t_end - curve.t_start) * f
        };
        let sp = strip_pcurve.point_at(u);
        let worst = (curve.geometry.point_at(t) - exact)
            .length()
            .max((strip_surface.point_at(sp.x, sp.y) - exact).length());
        if worst > pos_tol {
            return Err(Reason::ArcBranchUnresolved);
        }
    }
    Ok(SectionArc {
        edge,
        twin,
        curve,
        strip_pcurve,
    })
}

/// Author one stub as an exact section arc (see [`section_arc`]) with
/// its cap pcurve: the ellipse's image in the cap chart `plane`
/// (origin, unit u axis, unit v axis). Pure computation.
pub(super) fn compute_oblique_stub_arc(
    brep: &BRep,
    he_strip: HalfEdgeHandle,
    section: &ObliqueSection,
    u_start: f64,
    u_end: f64,
    plane: (Point3, Vec3, Vec3),
    tolerance: &Tolerance,
) -> Result<StagedStubArc, Reason> {
    let arc = section_arc(brep, he_strip, section, u_start, u_end, tolerance)?;
    let center = section.origin + section.axis * section.v0;
    // Cap pcurve: the ellipse in the plane chart, conjugate vectors
    // U, V at centre c, converted to principal axes. The twin runs
    // opposite, from u_end back to u_start.
    let (p0, pu, pv) = plane;
    let uv = |d: Vec3| (d.dot(&pu), d.dot(&pv));
    let c = uv(center - p0);
    let big_u = uv(section.e1 * section.radius + section.axis * section.v_cos);
    let big_v = uv(section.e2 * section.radius + section.axis * section.v_sin);
    let dot = big_u.0 * big_v.0 + big_u.1 * big_v.1;
    let nu = big_u.0 * big_u.0 + big_u.1 * big_u.1;
    let nv = big_v.0 * big_v.0 + big_v.1 * big_v.1;
    let t0 = 0.5 * (2.0 * dot).atan2(nu - nv);
    let (c0, s0) = (t0.cos(), t0.sin());
    let major = (big_u.0 * c0 + big_v.0 * s0, big_u.1 * c0 + big_v.1 * s0);
    let minor = (-big_u.0 * s0 + big_v.0 * c0, -big_u.1 * s0 + big_v.1 * c0);
    let a = major.0.hypot(major.1);
    let b = minor.0.hypot(minor.1);
    if a < 1e-12 || b < 1e-12 {
        return Err(Reason::ArcBranchUnresolved);
    }
    let rotation = major.1.atan2(major.0);
    // Point(u) = c + major·cos(u − t0) + minor·sin(u − t0); `minor`
    // is ±b along the rotated v axis, and a negative orientation
    // mirrors the angle.
    let orientation = if major.0 * minor.1 - major.1 * minor.0 >= 0.0 {
        1.0
    } else {
        -1.0
    };
    let alpha = |u: f64| orientation * (u - t0);
    let cap_pcurve = Pcurve::elliptical(c, a, b, rotation, alpha(u_end), alpha(u_start), 0.0, 1.0);
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    for k in 0..=8 {
        let f = k as f64 / 8.0;
        let exact = section.point(u_start + (u_end - u_start) * f);
        let cp = cap_pcurve.point_at(1.0 - f);
        if (p0 + pu * cp.x + pv * cp.y - exact).length() > pos_tol {
            return Err(Reason::ArcBranchUnresolved);
        }
    }
    Ok(StagedStubArc {
        edge: arc.edge,
        curve: arc.curve,
        strip_half_edge: he_strip,
        strip_pcurve: arc.strip_pcurve,
        cap_half_edge: arc.twin,
        cap_pcurve,
    })
}

/// Stage an oblique end's settlement: the settled corner and both stub
/// arcs, or why not. `points` are `[anchor_in, stale corner,
/// anchor_out]`; `cap` is the cap chart plus its unit normal.
pub(super) fn stage_oblique_end(
    brep: &BRep,
    c: &super::terminal_arc_trim::CornerCandidate,
    cyl: &StripCylinder,
    points: [Point3; 3],
    cap: ((Point3, Vec3, Vec3), Vec3),
    tolerance: &Tolerance,
) -> Result<(Point3, Vec<StagedStubArc>), Reason> {
    let [anchor_in, old_corner, anchor_out] = points;
    let ((plane_origin, pu, pv), normal) = cap;
    let section =
        ObliqueSection::new(cyl, plane_origin, normal).ok_or(Reason::AnchorsOffSurface)?;
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    if section.miss(anchor_in) > pos_tol || section.miss(anchor_out) > pos_tol {
        return Err(Reason::AnchorsOffSurface);
    }
    // Chart branch of the strip's existing pcurves at the in-anchor.
    let strip = popkern_core::geometry::surface::SurfaceGeometry::Cylinder {
        origin: cyl.origin,
        axis: cyl.axis,
        radius: cyl.radius,
    };
    let (chart_u, _) =
        popkern_core::geometry::uv_projection::project_point_to_uv(&strip, &anchor_in)
            .ok_or(Reason::ArcBranchUnresolved)?;
    let u_in = unwrap_near(section.angle_of(anchor_in), chart_u);
    let u_corner = unwrap_near(section.angle_of(old_corner), u_in);
    let u_out = unwrap_near(section.angle_of(anchor_out), u_corner);
    // The settled corner must sit strictly between the feet on the
    // minor branch.
    let (first, second) = (u_corner - u_in, u_out - u_corner);
    if first * second <= 0.0 || (u_out - u_in).abs() >= std::f64::consts::PI {
        return Err(Reason::ArcBranchUnresolved);
    }
    let new_corner = section.settle(old_corner);
    let semi_major = cyl.radius
        / normal
            .normalized()
            .ok_or(Reason::AnchorsOffSurface)?
            .dot(&section.axis)
            .abs();
    if (new_corner - old_corner).length() > semi_major {
        return Err(Reason::ArcBranchUnresolved);
    }
    let plane = (plane_origin, pu, pv);
    let staged = vec![
        compute_oblique_stub_arc(brep, c.he_in, &section, u_in, u_corner, plane, tolerance)?,
        compute_oblique_stub_arc(brep, c.he_out, &section, u_corner, u_out, plane, tolerance)?,
    ];
    Ok((new_corner, staged))
}

/// Settle a FREE end: the strip runs out at a sheet boundary, and its
/// two stubs through the stale corner are boundary edges with no cap
/// face (open-shell fillet with `CapKind::Open`: one stub 4.1e-2 off
/// the strip, the other left without a curve). The stubs lie along
/// the two supports' rim edges, so the plane through foot, corner and
/// foot is the rim plane. The strip's section by it is the exact end,
/// and it is settled like an oblique cap with strip pcurves only.
/// Topology is unchanged; the caller holds no snapshot, so everything
/// is staged before committing.
pub(super) fn repair_free_end(
    brep: &mut BRep,
    c: &super::terminal_arc_trim::CornerCandidate,
    cyl: &StripCylinder,
    tolerance: &Tolerance,
) -> Result<(), Reason> {
    // Incident edges counted through loop successors: a sheet-boundary
    // half-edge has no twin to name its destination.
    let mut incident: Vec<popkern_core::topology::EdgeHandle> = brep
        .half_edges
        .values()
        .filter(|h| {
            h.origin == c.corner
                || brep
                    .half_edges
                    .get(h.next)
                    .is_some_and(|n| n.origin == c.corner)
        })
        .map(|h| h.edge)
        .collect();
    incident.sort();
    incident.dedup();
    if incident.len() != 2 {
        return Err(Reason::VertexNotIsolated);
    }
    let pos = |v| {
        brep.vertices
            .get(v)
            .map(|x: &popkern_core::topology::Vertex| x.position)
    };
    let he_in = brep
        .half_edges
        .get(c.he_in)
        .ok_or(Reason::ArcBranchUnresolved)?;
    let out_twin = brep
        .half_edges
        .get(c.he_out)
        .and_then(|h| brep.half_edges.get(h.next))
        .ok_or(Reason::ArcBranchUnresolved)?;
    let foot_a = pos(he_in.origin).ok_or(Reason::ArcBranchUnresolved)?;
    let foot_b = pos(out_twin.origin).ok_or(Reason::ArcBranchUnresolved)?;
    let stale = pos(c.corner).ok_or(Reason::ArcBranchUnresolved)?;
    let normal = (stale - foot_a)
        .cross(&(foot_b - stale))
        .normalized()
        .ok_or(Reason::AnchorsOffSurface)?;
    let section = ObliqueSection::new(cyl, stale, normal).ok_or(Reason::AnchorsOffSurface)?;
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    if section.miss(foot_a) > pos_tol || section.miss(foot_b) > pos_tol {
        return Err(Reason::AnchorsOffSurface);
    }
    let strip = popkern_core::geometry::surface::SurfaceGeometry::Cylinder {
        origin: cyl.origin,
        axis: cyl.axis,
        radius: cyl.radius,
    };
    let (chart_u, _) = popkern_core::geometry::uv_projection::project_point_to_uv(&strip, &foot_a)
        .ok_or(Reason::ArcBranchUnresolved)?;
    let u_in = unwrap_near(section.angle_of(foot_a), chart_u);
    let u_corner = unwrap_near(section.angle_of(stale), u_in);
    let u_out = unwrap_near(section.angle_of(foot_b), u_corner);
    if (u_corner - u_in) * (u_out - u_corner) <= 0.0 || (u_out - u_in).abs() >= std::f64::consts::PI
    {
        return Err(Reason::ArcBranchUnresolved);
    }
    let settled = section.settle(stale);
    let arcs = [
        section_arc(brep, c.he_in, &section, u_in, u_corner, tolerance)?,
        section_arc(brep, c.he_out, &section, u_corner, u_out, tolerance)?,
    ];
    if let Some(v) = brep.vertices.get_mut(c.corner) {
        v.position = settled;
    }
    for (arc, he) in arcs.into_iter().zip([c.he_in, c.he_out]) {
        if let Some(e) = brep.edges.get_mut(arc.edge) {
            e.curve = Some(arc.curve);
        }
        if let Some(h) = brep.half_edges.get_mut(he) {
            h.pcurve = Some(arc.strip_pcurve.clone());
        }
        brep.set_face_use_pcurve(c.strip_face, he, arc.strip_pcurve);
    }
    Ok(())
}
