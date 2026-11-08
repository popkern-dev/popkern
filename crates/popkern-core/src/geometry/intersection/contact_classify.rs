//! What two surfaces do at a point where they are tangent (PK-1125).
//!
//! Where the normals of two surfaces are parallel, the intersection is
//! not a transversal curve, and what it is depends on how the surfaces
//! bend away from their common tangent plane. Write both as heights
//! over that plane, along the first surface's unit normal `n`, in one
//! orthonormal frame `(e1, e2)` of it: near the point,
//!
//! `h₁(t) − h₂(t) = ½ tᵀ D t + O(|t|³)`, `D = II₁ − II₂`,
//!
//! where `IIᵢ` is surface `i`'s second fundamental form *in that common
//! frame* (each surface's own `(L, M, N)` are in its own parameter
//! basis, and differ from the frame's by its Jacobian: subtracting them
//! directly, as the first classifier did, compares different bases).
//! The relative curvature `D` decides (the relative Dupin indicatrix):
//!
//! - both eigenvalues of one sign: the surfaces touch at an isolated
//!   point, the second wholly on one side of the first near it;
//! - opposite signs: two branches cross, leaving along the null
//!   directions of `D`;
//! - an eigenvalue indistinguishable from zero: a tangent curve, a
//!   cusp, contact of higher order, or coincidence; second order cannot
//!   tell, so the gap `h₁ − h₂` itself is sampled on a ring about the
//!   point. Its sign changes count the branches (each crosses the ring
//!   twice), arcs where it stays within the tolerance are a tangent
//!   curve, and a ring wholly within it is coincidence.
//!
//! "Indistinguishable from zero" is a tolerance statement, not a
//! floating-point one: an eigenvalue `λ` moves the gap by `½ λ L²` over
//! the contact's extent `L`, so `|λ| ≤ 2ε / L²` (ε the linear
//! tolerance) cannot be told from zero anywhere the contact matters.
//! Both are explicit inputs ([`ContactScale`]); nothing here is global.
//! The ring is sampled at two radii and must agree; when it does not,
//! the contact is [`SurfaceContactKind::Undecided`] and callers decline
//! rather than guess.
//!
//! A pole or apex of a surface's chart (`|S_u × S_v| ≈ 0`) is not a
//! contact of the pair (design R2.9); at one, the frame comes from the
//! normal alone and the ring test decides, its heights found by Newton
//! from a regular parameter beside the pole.

use crate::error::CoreError;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{PortableMath, Vec3};

/// The tolerance and length scale a contact is judged at.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContactScale {
    /// Gaps within this are no gap: the faces' linear tolerance.
    pub linear: f64,
    /// The extent over which the contact matters: the size of the faces
    /// meeting there (a contact's window never exceeds them).
    pub extent: f64,
}

impl ContactScale {
    /// The relative curvature below which a direction bends no more
    /// than `linear` over `extent`: `2ε / L²`.
    pub fn curvature_floor(&self) -> f64 {
        2.0 * self.linear / (self.extent * self.extent)
    }
}

/// Which side of the first surface the second lies, near a touch: along
/// the first surface's normal as evaluated (`normal_at`), not its face's
/// material sense, which the caller applies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContactSide {
    /// The second surface is on the side the first's normal points to.
    AlongNormal,
    /// The second surface is on the side opposite the first's normal.
    AgainstNormal,
}

/// What two tangent surfaces do at a point.
#[derive(Debug, Clone, PartialEq)]
pub enum SurfaceContactKind {
    /// They touch at an isolated point.
    TangentPoint {
        /// Where the second surface lies, near the point.
        side: ContactSide,
    },
    /// Branches of their intersection cross at the point.
    CrossingPoint {
        /// The unit tangent of each branch at the point (a branch passes
        /// through, so `d` and `−d` are one branch).
        branches: Vec<Vec3>,
    },
    /// They are tangent along a curve through the point.
    TangentCurve {
        /// The curve's unit tangent at the point.
        direction: Vec3,
        /// Where the second surface lies on both sides of the curve;
        /// `None` when the sides differ (they cross along it).
        side: Option<ContactSide>,
    },
    /// They coincide about the point, within the tolerance over the ring.
    Coincident,
    /// The classification is not decided at the given scale.
    Undecided {
        /// Why.
        reason: String,
    },
}

/// The classification and the evidence behind it.
#[derive(Debug, Clone, PartialEq)]
pub struct ContactClassification {
    /// What the surfaces do.
    pub kind: SurfaceContactKind,
    /// The unit normal the heights were measured along (surface 1's).
    pub normal: Vec3,
    /// Eigenvalues of the relative curvature `D` (ascending), when the
    /// point is regular on both charts.
    pub eigenvalues: Option<[f64; 2]>,
    /// The curvature floor `2ε / L²` they were judged against.
    pub floor: f64,
    /// The ring radius and sign changes, when the ring decided.
    pub ring: Option<(f64, usize)>,
}

/// Classify the contact of `s1` at `uv1` with `s2` at `uv2`, two
/// parameters of one point where the surfaces are tangent.
///
/// # Errors
///
/// [`CoreError::PreconditionFailed`] when the normals are not parallel
/// (no tangency: the surfaces cross transversally there) or the two
/// parameters are not one point, within `scale.linear`.
pub fn classify_contact(
    s1: &SurfaceGeometry,
    uv1: (f64, f64),
    s2: &SurfaceGeometry,
    uv2: (f64, f64),
    scale: ContactScale,
) -> Result<ContactClassification, CoreError> {
    let p1 = s1.point_at(uv1.0, uv1.1);
    let p2 = s2.point_at(uv2.0, uv2.1);
    if p1.distance_to(&p2) > scale.linear {
        return Err(CoreError::PreconditionFailed(format!(
            "contact parameters are {} apart, past the tolerance {}",
            p1.distance_to(&p2),
            scale.linear
        )));
    }
    let n = unit_normal(s1, uv1)
        .or_else(|| unit_normal(s2, uv2))
        .ok_or_else(|| CoreError::PreconditionFailed("no normal at the contact".into()))?;
    if let Some(n2) = unit_normal(s2, uv2) {
        // Parallel within the angle a gap of ε subtends over the extent.
        let sin = n.cross(&n2).length();
        if sin > 4.0 * scale.linear / scale.extent + 1e-9 {
            return Err(CoreError::PreconditionFailed(format!(
                "normals are not parallel at the contact (sin {sin:.3e})"
            )));
        }
    }
    let (e1, e2) = frame(s1, uv1, n);
    let floor = scale.curvature_floor();
    let eigen = match (
        frame_form(s1, uv1, n, e1, e2),
        frame_form(s2, uv2, n, e1, e2),
    ) {
        (Some(f1), Some(f2)) => Some([
            [f1[0][0] - f2[0][0], f1[0][1] - f2[0][1]],
            [f1[1][0] - f2[1][0], f1[1][1] - f2[1][1]],
        ]),
        _ => None,
    }
    .map(symmetric_eigen);
    let base = ContactClassification {
        kind: SurfaceContactKind::Coincident,
        normal: n,
        eigenvalues: eigen.map(|(l, _)| l),
        floor,
        ring: None,
    };
    if let Some(([lo, hi], [d_lo, d_hi])) = eigen
        && lo.abs() > floor
        && hi.abs() > floor
    {
        let to_world = |d: [f64; 2]| e1 * d[0] + e2 * d[1];
        let kind = if lo.signum() == hi.signum() {
            // h₁ − h₂ > 0 everywhere: the second surface is below.
            SurfaceContactKind::TangentPoint {
                side: if lo > 0.0 {
                    ContactSide::AgainstNormal
                } else {
                    ContactSide::AlongNormal
                },
            }
        } else {
            // λ_lo c² + λ_hi s² = 0 in the eigenbasis.
            let t = (-lo / hi).sqrt();
            let (c, s) = (1.0 / (1.0 + t * t).sqrt(), t / (1.0 + t * t).sqrt());
            let a = [d_lo[0] * c + d_hi[0] * s, d_lo[1] * c + d_hi[1] * s];
            let b = [d_lo[0] * c - d_hi[0] * s, d_lo[1] * c - d_hi[1] * s];
            SurfaceContactKind::CrossingPoint {
                branches: vec![to_world(a), to_world(b)],
            }
        };
        return Ok(ContactClassification { kind, ..base });
    }
    // Second order cannot tell: the ring decides, at two radii.
    let largest = eigen.map_or(0.0, |([lo, hi], _)| lo.abs().max(hi.abs()));
    let rho = if largest > floor {
        // Where the bending direction's gap is 8ε: well above the band.
        (16.0 * scale.linear / largest)
            .sqrt()
            .min(0.25 * scale.extent)
    } else {
        0.25 * scale.extent
    };
    let outer = ring_pattern(s1, uv1, s2, uv2, (n, e1, e2), p1, rho, scale.linear);
    let inner = ring_pattern(s1, uv1, s2, uv2, (n, e1, e2), p1, 0.5 * rho, scale.linear);
    let (Some(outer), Some(inner)) = (outer, inner) else {
        return Ok(ContactClassification {
            kind: SurfaceContactKind::Undecided {
                reason: "a ring height did not converge".into(),
            },
            ..base
        });
    };
    if outer.summary() != inner.summary() {
        return Ok(ContactClassification {
            kind: SurfaceContactKind::Undecided {
                reason: format!(
                    "the ring at {rho:.3e} reads {:?}, at half that {:?}",
                    outer.summary(),
                    inner.summary()
                ),
            },
            ring: Some((rho, outer.sign_changes)),
            ..base
        });
    }
    let kind = outer.kind(e1, e2);
    Ok(ContactClassification {
        kind,
        ring: Some((rho, outer.sign_changes)),
        ..base
    })
}

/// Whether the chart is regular at `uv`: `|S_u × S_v|` against the
/// larger partial's square, so a pole (one partial vanishing, as
/// `sin v` does at a sphere's) reads singular at any radius.
pub(super) fn regular(s: &SurfaceGeometry, uv: (f64, f64)) -> bool {
    let su = s.partial_u(uv.0, uv.1);
    let sv = s.partial_v(uv.0, uv.1);
    let big = su.length().max(sv.length());
    big > 0.0 && su.cross(&sv).length() > 1e-8 * big * big
}

fn unit_normal(s: &SurfaceGeometry, uv: (f64, f64)) -> Option<Vec3> {
    let su = s.partial_u(uv.0, uv.1);
    let sv = s.partial_v(uv.0, uv.1);
    let n = su.cross(&sv);
    let evaluated = s.normal_at(uv.0, uv.1);
    if !regular(s, uv) {
        // A chart pole: the normal from the surface's own evaluator.
        return evaluated.normalized();
    }
    // The cross product is the precise direction; its sense is the
    // evaluator's (a chart's S_u × S_v may point inward).
    let n = n.normalized()?;
    Some(if n.dot(&evaluated) < 0.0 { n * -1.0 } else { n })
}

/// An orthonormal frame of the tangent plane normal to `n`, its first
/// axis along `S_u` where the chart is regular.
fn frame(s: &SurfaceGeometry, uv: (f64, f64), n: Vec3) -> (Vec3, Vec3) {
    let su = s.partial_u(uv.0, uv.1);
    let along = su - n * su.dot(&n);
    let e1 = if regular(s, uv) {
        along.normalized().unwrap_or_else(|| any_perpendicular(n))
    } else {
        any_perpendicular(n)
    };
    (e1, n.cross(&e1))
}

fn any_perpendicular(n: Vec3) -> Vec3 {
    let pick = if n.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    (pick - n * pick.dot(&n))
        .normalized()
        .unwrap_or(Vec3::new(0.0, 0.0, 1.0))
}

/// Surface `s`'s height form `½ tᵀ II t` over the tangent plane, in the
/// frame `(e1, e2)`, heights along `n`; `None` at a chart pole.
fn frame_form(
    s: &SurfaceGeometry,
    uv: (f64, f64),
    n: Vec3,
    e1: Vec3,
    e2: Vec3,
) -> Option<[[f64; 2]; 2]> {
    if !regular(s, uv) {
        return None;
    }
    let su = s.partial_u(uv.0, uv.1);
    let sv = s.partial_v(uv.0, uv.1);
    // Frame coordinates of the parameter directions: [a b; c d] maps
    // (du, dv) to (t₁, t₂).
    let (a, b, c, d) = (su.dot(&e1), sv.dot(&e1), su.dot(&e2), sv.dot(&e2));
    let det = a * d - b * c;
    // Its inverse maps frame coordinates to parameter ones.
    let inv = [[d / det, -b / det], [-c / det, a / det]];
    let l = s.second_partial_uu(uv.0, uv.1).dot(&n);
    let m = s.second_partial_uv(uv.0, uv.1).dot(&n);
    let nn = s.second_partial_vv(uv.0, uv.1).dot(&n);
    let ii = [[l, m], [m, nn]];
    // invᵀ · II · inv.
    let mut out = [[0.0; 2]; 2];
    for (i, row) in out.iter_mut().enumerate() {
        for (j, cell) in row.iter_mut().enumerate() {
            let mut acc = 0.0;
            for (k, ii_k) in ii.iter().enumerate() {
                for (l, ii_kl) in ii_k.iter().enumerate() {
                    acc += inv[k][i] * ii_kl * inv[l][j];
                }
            }
            *cell = acc;
        }
    }
    Some(out)
}

/// Eigenvalues (ascending) and unit eigenvectors of a symmetric 2 × 2.
fn symmetric_eigen(m: [[f64; 2]; 2]) -> ([f64; 2], [[f64; 2]; 2]) {
    let (a, b, d) = (m[0][0], 0.5 * (m[0][1] + m[1][0]), m[1][1]);
    let mean = 0.5 * (a + d);
    let radius = (0.25 * (a - d) * (a - d) + b * b).sqrt();
    let (lo, hi) = (mean - radius, mean + radius);
    // The angle of the high eigenvector: tan 2θ = 2b / (a − d).
    let theta = 0.5 * (2.0 * b).patan2(a - d);
    let v_hi = [theta.pcos(), theta.psin()];
    let v_lo = [-v_hi[1], v_hi[0]];
    ([lo, hi], [v_lo, v_hi])
}

#[path = "contact_classify_ring.rs"]
mod ring;
use ring::ring_pattern;

#[cfg(test)]
#[path = "contact_classify_tests.rs"]
mod tests;
