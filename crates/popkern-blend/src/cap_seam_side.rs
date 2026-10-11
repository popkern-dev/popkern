//! One side of a measured seam (split from `cap_seam_measure` for the
//! 500-line cap): the pcurve-true evaluation, its normal, and its normal
//! curvature. PK-579: curvature comes from exact jets. Position finite
//! differences declined exactly where curvature varies fast, and the
//! surfaces' own second-partial evaluators (differences of analytic
//! first partials) err by ~1e-3 near dual-jet corners, the target's
//! resolution. The exact-form conversion is costly, so one jet source per
//! FACE is shared through a [`SurfaceJetCache`] by every seam of a
//! measurement call.

use std::cell::RefCell;
use std::collections::HashMap;

use popkern_core::geometry::NurbsSurface;
use popkern_core::geometry::certify::surface_jet::gregory_exact_form;
use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

/// Exact spline forms per face, built lazily and shared by every seam a
/// measurement call samples (a cap face bounds four or more seams). A
/// NURBS face is its own form; a Gregory face's exact rational form is
/// the one the certifier encloses. Point-exact second partials from it
/// cost a few microseconds; an order-6 interval jet cost ~380 µs, and
/// the 12-edge G2 cube spent 9.6 s in 25k of them.
#[derive(Default)]
pub(crate) struct SurfaceJetCache {
    forms: RefCell<HashMap<FaceHandle, Option<NurbsSurface>>>,
}

impl SurfaceJetCache {
    /// Exact `(S_u, S_v, S_uu, S_uv, S_vv)` of `face`'s `surface` at
    /// `uv`, or `None` when the kind has no spline form (the analytic
    /// kinds, whose own evaluators are exact enough) or the read is
    /// non-finite.
    fn jet(
        &self,
        face: FaceHandle,
        surface: &SurfaceGeometry,
        uv: (f64, f64),
    ) -> Option<[Vec3; 5]> {
        // A Gregory patch with closed-form seconds is read directly:
        // its exact rational form is ill-conditioned next to activated
        // dual-jet corners (see `GregoryPatch::analytic_second_partials`).
        if let SurfaceGeometry::Gregory(patch) = surface
            && let Some((suu, suv, svv)) = patch.analytic_second_partials(uv.0, uv.1)
        {
            let jet = [
                patch.partial_u(uv.0, uv.1),
                patch.partial_v(uv.0, uv.1),
                suu,
                suv,
                svv,
            ];
            return jet
                .iter()
                .all(|d| d.x.is_finite() && d.y.is_finite() && d.z.is_finite())
                .then_some(jet);
        }
        let mut map = self.forms.borrow_mut();
        let form = map
            .entry(face)
            .or_insert_with(|| match surface {
                SurfaceGeometry::Nurbs(n) => Some(n.clone()),
                SurfaceGeometry::Gregory(patch) => gregory_exact_form(patch).ok(),
                _ => None,
            })
            .as_ref()?;
        let (_, su, sv, suu, suv, svv) = form.point_and_two_partials(uv.0, uv.1);
        let jet = [su, sv, suu, suv, svv];
        jet.iter()
            .all(|d| d.x.is_finite() && d.y.is_finite() && d.z.is_finite())
            .then_some(jet)
    }
}

pub(crate) struct SeamSide<'a> {
    surface: &'a SurfaceGeometry,
    face: FaceHandle,
    jets: &'a SurfaceJetCache,
    pcurve: &'a Pcurve,
    /// Pcurve runs opposite to the edge curve (twin half-edges trace
    /// the loop direction, not the curve direction) — resolved once
    /// per side by matching pcurve-endpoint surface images to the
    /// curve endpoints.
    reversed: bool,
}

pub(crate) fn seam_side<'a>(
    brep: &'a BRep,
    he: HalfEdgeHandle,
    curve: &CurveOnEdge,
    jets: &'a SurfaceJetCache,
) -> Option<SeamSide<'a>> {
    let he_data = brep.half_edges.get(he)?;
    let pcurve = he_data.pcurve.as_ref()?;
    let face = brep.loops.get(he_data.loop_)?.face;
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    let curve_start = (&curve.geometry as &dyn CurveEval).point_at(curve.t_start);
    let uv0 = pcurve.point_at(pcurve.t_start);
    let uv1 = pcurve.point_at(pcurve.t_end);
    let s0 = surface.point_at(uv0.x, uv0.y);
    let s1 = surface.point_at(uv1.x, uv1.y);
    let reversed = s1.distance_to(&curve_start) < s0.distance_to(&curve_start);
    Some(SeamSide {
        surface,
        face,
        jets,
        pcurve,
        reversed,
    })
}

impl SeamSide<'_> {
    /// Evaluate at fraction `f` of the seam in CURVE direction.
    pub(crate) fn eval_fraction(&self, f: f64) -> ((f64, f64), Point3) {
        let g = if self.reversed { 1.0 - f } else { f };
        let uv = self
            .pcurve
            .point_at(self.pcurve.t_start + g * (self.pcurve.t_end - self.pcurve.t_start));
        ((uv.x, uv.y), self.surface.point_at(uv.x, uv.y))
    }

    /// Finite-difference unit normal, or `None` where the chart is
    /// too degenerate to trust (see module docs).
    pub(crate) fn normal(&self, uv: (f64, f64)) -> Option<Vec3> {
        // Exact partials where the kind has an exact form: position
        // differences cannot separate the partials on a nearly
        // degenerate strip chart (pentagon5_variable's small-radius rims
        // left every sample without a normal, and the seam unmeasured).
        if let Some(jet) = self.jets.jet(self.face, self.surface, uv) {
            return chart_normal(jet[0], jet[1]);
        }
        let h = 1e-6;
        let su = self.surface.point_at(uv.0 + h, uv.1) - self.surface.point_at(uv.0 - h, uv.1);
        let sv = self.surface.point_at(uv.0, uv.1 + h) - self.surface.point_at(uv.0, uv.1 - h);
        chart_normal(su, sv)
    }

    /// Normal curvature along `dir`, signed against the caller's
    /// `normal` (hemisphere-aligned across a seam so signed
    /// curvatures compare — a chart-orientation flip must not read a
    /// perfect seam as a `2κ` jump), with a two-step convergence
    /// guard: near a degenerate polar chart the FD truncation error
    /// scales quadratically with the step (measured 1.7 → 0.16 →
    /// 1.9e-3 across steps 1e-3/1e-4/1e-5 on the same pentagon fan
    /// seam), so a single-step estimate reads chart noise as a
    /// curvature crease. Samples whose two-step estimates disagree
    /// are skipped; agreeing ones are Richardson-extrapolated.
    pub(crate) fn normal_curvature_signed(
        &self,
        uv: (f64, f64),
        dir: Vec3,
        normal: Vec3,
    ) -> Option<f64> {
        if let Some(jet) = self.local_jet(uv) {
            return normal_curvature_from_jet(&jet, dir, normal);
        }
        let k_h = self.normal_curvature_at_step(uv, dir, normal, 1e-4)?;
        let k_h2 = self.normal_curvature_at_step(uv, dir, normal, 5e-5)?;
        if (k_h - k_h2).abs() > 0.1 * (1.0 + k_h2.abs()) {
            return None;
        }
        Some((4.0 * k_h2 - k_h) / 3.0)
    }

    /// First and second partials `(S_u, S_v, S_uu, S_uv, S_vv)` at `uv`:
    /// exact jets from the shared cache, else the surface's own
    /// derivative evaluators (for kinds with no exact form). `None` only
    /// for a non-finite read.
    fn local_jet(&self, uv: (f64, f64)) -> Option<[Vec3; 5]> {
        if let Some(jet) = self.jets.jet(self.face, self.surface, uv) {
            return Some(jet);
        }
        let (u, v) = uv;
        let jet = [
            self.surface.partial_u(u, v),
            self.surface.partial_v(u, v),
            self.surface.second_partial_uu(u, v),
            self.surface.second_partial_uv(u, v),
            self.surface.second_partial_vv(u, v),
        ];
        jet.iter()
            .all(|d| d.x.is_finite() && d.y.is_finite() && d.z.is_finite())
            .then_some(jet)
    }

    fn normal_curvature_at_step(
        &self,
        uv: (f64, f64),
        dir: Vec3,
        normal: Vec3,
        step: f64,
    ) -> Option<f64> {
        let h = 1e-6;
        let du = self.surface.point_at(uv.0 + h, uv.1) - self.surface.point_at(uv.0 - h, uv.1);
        let dv = self.surface.point_at(uv.0, uv.1 + h) - self.surface.point_at(uv.0, uv.1 - h);
        let su = du * (0.5 / h);
        let sv = dv * (0.5 / h);
        let (e, f_, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
        let (r1, r2) = (dir.dot(&su), dir.dot(&sv));
        let det = e * g - f_ * f_;
        if det.abs() < 1e-30 {
            return None;
        }
        let a = (r1 * g - r2 * f_) / det;
        let b = (r2 * e - r1 * f_) / det;
        let p0 = self.surface.point_at(uv.0, uv.1);
        let p1 = self.surface.point_at(uv.0 + a * step, uv.1 + b * step);
        let p2 = self.surface.point_at(uv.0 - a * step, uv.1 - b * step);
        let second = (p1 - p0) + (p2 - p0);
        let first_len = (p1 - p2).length() * 0.5;
        if first_len < 1e-30 {
            return None;
        }
        Some(second.dot(&normal) / (first_len * first_len))
    }
}

/// Relative size of `|S_u × S_v|` against `|S_u|·|S_v|` below which
/// the partials are (nearly) parallel: a genuinely singular chart
/// point (a pole row, a collapsed edge).
const SINGULAR_SINE: f64 = 1.0e-6;

/// Unit normal of a chart with partials `su`, `sv`, or `None` at a
/// singular point. The test is the ANGLE between the partials, not
/// their speed ratio: a speed-ratio test (`min < 1e-3·max`) declined
/// every sample of pentagon5_variable's small-radius rims, whose strip
/// chart is merely anisotropic (PK-579), leaving two seams unmeasured.
fn chart_normal(su: Vec3, sv: Vec3) -> Option<Vec3> {
    let n = su.cross(&sv);
    let scale = su.length() * sv.length();
    if scale.is_nan() || scale <= 1e-300 || n.length() < SINGULAR_SINE * scale {
        return None;
    }
    n.normalized()
}

/// Normal curvature along `dir` from first and second partials: `dir`'s tangent
/// components `(a, b)` in the `(S_u, S_v)` frame (least squares through
/// the first fundamental form), then `II(a, b) / I(a, b)` signed
/// against `normal`. `None` only where the chart is singular (the
/// first fundamental form degenerates), never for convergence.
fn normal_curvature_from_jet(jet: &[Vec3; 5], dir: Vec3, normal: Vec3) -> Option<f64> {
    let [su, sv, suu, suv, svv] = *jet;
    let (e, f_, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
    let det = e * g - f_ * f_;
    if det <= 1e-24 * e * g || det.is_nan() {
        return None;
    }
    let (r1, r2) = (dir.dot(&su), dir.dot(&sv));
    let a = (r1 * g - r2 * f_) / det;
    let b = (r2 * e - r1 * f_) / det;
    let first = e * a * a + 2.0 * f_ * a * b + g * b * b;
    if first <= 0.0 || first.is_nan() {
        return None;
    }
    // The jet's OWN normal, sign-aligned to the caller's: mixing a
    // normal from one representation with second partials of another
    // moved readings by ~1e-3 near rational dual-jet corners.
    let own = su.cross(&sv).normalized()?;
    let n = if own.dot(&normal) < 0.0 {
        own * -1.0
    } else {
        own
    };
    let second = a * a * suu.dot(&n) + 2.0 * a * b * suv.dot(&n) + b * b * svv.dot(&n);
    Some(second / first)
}
