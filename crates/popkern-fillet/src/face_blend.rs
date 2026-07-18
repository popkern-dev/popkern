//! PK-142 stage 1 — face–face blends: rolling-ball surface between
//! two faces that share no edge. Design:
//! `docs/architecture/face-face-blend-design.md`.
//!
//! Stage 1 scope: constant radius, `SheetOnly` attach (the strip is
//! returned as an open sheet body; wall trimming is stage 2). The
//! gap-medial seed places the ball at the closest-pair midpoint —
//! near the two-contact solution by construction — and the march
//! pins the solver's free DOF with a plane through the stepped
//! centre (`spine_point`/`spine_plane_normal`), warm-starting each
//! step. Sections are arcs between the converged contacts, skinned
//! exactly as the edge-fillet sweep; the centre polyline IS the
//! marched span, so no carrier-overhang class exists here.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::contact::{
    BallSide, ContactOptions, ContactRequest, ContactSolution, solve_rolling_ball_contact,
};
use crate::error::BlendError;

/// Attach policy for a face–face blend.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceBlendAttach {
    /// Emit the blend strip as a separate open sheet body; the input
    /// body is not modified.
    SheetOnly,
    /// Trim both walls back to the contact curves, remove the
    /// enclosed gap region, and seam the strip in as a face of the
    /// (copied) body. Transactional: declines with
    /// [`BlendError::AttachDeclined`] and the input unchanged when
    /// the surgery preconditions do not hold.
    TrimBoth,
}

/// Specification of one face–face blend.
#[derive(Debug, Clone)]
pub struct FaceBlendSpec {
    /// First support face.
    pub face_a: FaceHandle,
    /// Second support face.
    pub face_b: FaceHandle,
    /// Rolling-ball radius (constant; variable laws are stage 3).
    pub radius: f64,
    /// Attach policy.
    pub attach: FaceBlendAttach,
}

/// Result of a face–face blend: the strip sheet body plus the
/// marched contact polylines (for callers that trim manually).
#[derive(Debug, Clone)]
pub struct FaceBlendResult {
    /// Open sheet body carrying the blend strip.
    pub sheet: BRep,
    /// For [`FaceBlendAttach::TrimBoth`]: the modified body with the
    /// walls trimmed and the strip seamed in. `None` for SheetOnly.
    pub attached: Option<BRep>,
    /// For TrimBoth: the strip face inside `attached`.
    pub strip_face: Option<FaceHandle>,
    /// Contact points on face A, in march order.
    pub contacts_a: Vec<Point3>,
    /// Contact points on face B, in march order.
    pub contacts_b: Vec<Point3>,
}

/// A support face's surface plus its solver and seeding extents.
pub(crate) struct FaceDomain<'a> {
    pub(crate) surface: &'a SurfaceGeometry,
    /// Padded range for solver headroom (tight-bounded, 5%).
    pub(crate) uv_range: ((f64, f64), (f64, f64)),
    /// Exact face extent (0% padding) — seeds must start ON the
    /// face, not in the padding band, or the first contact solve can
    /// pin the march outside the trim.
    pub(crate) uv_core: ((f64, f64), (f64, f64)),
}

fn face_domain<'a>(brep: &'a BRep, face: FaceHandle) -> Result<FaceDomain<'a>, BlendError> {
    let surface = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "face".into(),
            value: format!("{face:?}"),
            reason: "face is missing or has no surface".into(),
        })?;
    // Period-aware face box: vertex-projected ranges collapse to the
    // seam angle on full-wrap cylinder faces (measured: a rod blend
    // seeded on the rod's FAR side). The sampled exact box is the
    // core; the solver range pads it 5 % (min 0.05), capped at one
    // period for periodic u.
    let (cu, cv) =
        popkern_core::geometry::uv_projection::face_exact_parameter_range(brep, face, surface);
    let (u_period, _) = surface.parameter_periodicity();
    let pad = |(lo, hi): (f64, f64)| -> (f64, f64) {
        let p = (0.05 * (hi - lo).abs()).max(0.05);
        (lo - p, hi + p)
    };
    let u_range = match u_period {
        Some(period) if (cu.1 - cu.0) >= period - 1e-9 => cu,
        Some(period) => {
            let p = (0.05 * (cu.1 - cu.0).abs())
                .max(0.05)
                .min(0.5 * (period - (cu.1 - cu.0)));
            (cu.0 - p, cu.1 + p)
        }
        None => pad(cu),
    };
    let v_range = pad(cv);
    Ok(FaceDomain {
        surface,
        uv_range: (u_range, v_range),
        uv_core: (cu, cv),
    })
}

/// Closest pair between two bounded surface domains: coarse grid
/// cross-product, then alternating closest-point polish.
fn closest_pair(a: &FaceDomain, b: &FaceDomain) -> ((f64, f64), (f64, f64)) {
    const N: usize = 16;
    let sample = |d: &FaceDomain, i: usize, j: usize| -> (f64, f64) {
        let ((u0, u1), (v0, v1)) = d.uv_core;
        (
            u0 + (u1 - u0) * i as f64 / (N - 1) as f64,
            v0 + (v1 - v0) * j as f64 / (N - 1) as f64,
        )
    };
    let mut best = ((0.0, 0.0), (0.0, 0.0), f64::INFINITY);
    for i in 0..N {
        for j in 0..N {
            let uva = sample(a, i, j);
            let pa = a.surface.point_at(uva.0, uva.1);
            for k in 0..N {
                for l in 0..N {
                    let uvb = sample(b, k, l);
                    let pb = b.surface.point_at(uvb.0, uvb.1);
                    let d = pa.distance_to(&pb);
                    if d < best.2 {
                        best = (uva, uvb, d);
                    }
                }
            }
        }
    }
    let (mut uva, mut uvb, _) = best;
    // Alternating polish: project each point onto the other surface,
    // clamped to the CORE extent so the seed stays on the faces.
    for _ in 0..12 {
        let pa = a.surface.point_at(uva.0, uva.1);
        if let Some(nb) = b.surface.closest_point(&pa) {
            uvb = clamp_uv(nb, b.uv_core);
        }
        let pb = b.surface.point_at(uvb.0, uvb.1);
        if let Some(na) = a.surface.closest_point(&pb) {
            uva = clamp_uv(na, a.uv_core);
        }
    }
    (uva, uvb)
}

pub(crate) fn clamp_uv(uv: (f64, f64), range: ((f64, f64), (f64, f64))) -> (f64, f64) {
    (
        uv.0.clamp(range.0.0, range.0.1),
        uv.1.clamp(range.1.0, range.1.1),
    )
}

/// Blend two faces that share no edge with a constant-radius rolling
/// ball. See the module docs for scope; declines are structured.
pub fn face_blend(
    brep: &BRep,
    spec: &FaceBlendSpec,
    tolerance: Tolerance,
) -> Result<FaceBlendResult, BlendError> {
    if !(spec.radius.is_finite() && spec.radius > 0.0) {
        return Err(BlendError::InvalidParameter {
            parameter: "radius".into(),
            value: format!("{}", spec.radius),
            reason: "radius must be finite and positive".into(),
        });
    }
    let a = face_domain(brep, spec.face_a)?;
    let b = face_domain(brep, spec.face_b)?;

    // Gap-medial seed.
    let (uva, uvb) = closest_pair(&a, &b);
    let pa = a.surface.point_at(uva.0, uva.1);
    let pb = b.surface.point_at(uvb.0, uvb.1);
    let gap = pa.distance_to(&pb);
    let slack = tolerance.linear().max(1.0e-9) * 10.0;
    if gap > 2.0 * spec.radius + slack {
        return Err(BlendError::InvalidParameter {
            parameter: "radius".into(),
            value: format!("{}", spec.radius),
            reason: format!(
                "gap between faces ({gap:.6}) exceeds the ball diameter ({:.6}) — \
                 no two-contact solution exists",
                2.0 * spec.radius
            ),
        });
    }
    if gap < slack {
        return Err(BlendError::InvalidParameter {
            parameter: "face_b".into(),
            value: format!("{:?}", spec.face_b),
            reason: "faces touch or intersect — use fillet_edges on the shared edge".into(),
        });
    }
    let centre = Point3::new(
        (pa.x + pb.x) * 0.5,
        (pa.y + pb.y) * 0.5,
        (pa.z + pb.z) * 0.5,
    );
    let side = |surface: &SurfaceGeometry, uv: (f64, f64), p: Point3| -> BallSide {
        let n = surface.normal_at(uv.0, uv.1);
        if n.dot(&(centre - p)) >= 0.0 {
            BallSide::Outside
        } else {
            BallSide::Inside
        }
    };
    let side_a = side(a.surface, uva, pa);
    let side_b = side(b.surface, uvb, pb);

    // Parallel-support fast path: antiparallel normals make the
    // contact Jacobian singular (NormalsParallel), and tangency to
    // BOTH walls exists only at radius = gap/2 exactly (line-contact
    // family). Construct the march in closed form: centres on the
    // medial line, contacts as perpendicular feet.
    let n_seed_a = a.surface.normal_at(uva.0, uva.1);
    let n_seed_b = b.surface.normal_at(uvb.0, uvb.1);
    if n_seed_a.cross(&n_seed_b).length() < 1.0e-8 {
        if (gap - 2.0 * spec.radius).abs() > slack {
            return Err(BlendError::InvalidParameter {
                parameter: "radius".into(),
                value: format!("{}", spec.radius),
                reason: format!(
                    "parallel supports admit tangency only at radius = gap/2                      (gap {gap:.6}; requested diameter {:.6})",
                    2.0 * spec.radius
                ),
            });
        }
        return crate::face_blend_parallel::parallel_gap_blend(
            brep, spec, &a, &b, uva, uvb, tolerance,
        );
    }

    // March tangent: the common tangent of the two contact planes;
    // for (anti)parallel supports any support-tangent direction
    // perpendicular to the gap works — take A's u-derivative
    // projected off the gap normal.
    let march_tangent = |sol_uv0: (f64, f64), sol_uv1: (f64, f64)| -> Option<Vec3> {
        let na = a.surface.normal_at(sol_uv0.0, sol_uv0.1);
        let nb = b.surface.normal_at(sol_uv1.0, sol_uv1.1);
        let cross = na.cross(&nb);
        if cross.length() > 1.0e-8 {
            return cross.normalized();
        }
        let du = a.surface.partial_u(sol_uv0.0, sol_uv0.1);
        (du - na * du.dot(&na)).normalized()
    };

    let options = ContactOptions {
        clamp_initial_uv: true,
        ..ContactOptions::default()
    };
    let seed_tangent = march_tangent(uva, uvb).ok_or_else(|| BlendError::InvalidParameter {
        parameter: "face_a".into(),
        value: format!("{:?}", spec.face_a),
        reason: "degenerate support tangent at the gap-medial seed".into(),
    })?;
    let solve = |seed_uv0: (f64, f64),
                 seed_uv1: (f64, f64),
                 spine_point: Point3,
                 tangent: Vec3|
     -> Result<ContactSolution, String> {
        let request = ContactRequest {
            surface_0: a.surface as &dyn SurfaceEval,
            surface_1: b.surface as &dyn SurfaceEval,
            side_0: side_a,
            side_1: side_b,
            radius: spec.radius,
            initial_uv0: seed_uv0,
            uv0_range: a.uv_range,
            initial_uv1: seed_uv1,
            uv1_range: b.uv_range,
            spine_point,
            spine_plane_normal: tangent,
        };
        solve_rolling_ball_contact(&request, &options, &tolerance, None)
            .map_err(|e| format!("{e:?}"))
    };
    let first =
        solve(uva, uvb, centre, seed_tangent).map_err(|reason| BlendError::InvalidParameter {
            parameter: "radius".into(),
            value: format!("{}", spec.radius),
            reason: format!("gap-medial contact failed: {reason}"),
        })?;

    // March both directions until a contact leaves its trim domain.
    let step = (spec.radius * 0.5).max(tolerance.linear() * 100.0);
    let max_steps = 400usize;
    let mut forward: Vec<ContactSolution> = Vec::new();
    let mut backward: Vec<ContactSolution> = Vec::new();
    for (sign, out) in [(1.0, &mut forward), (-1.0, &mut backward)] {
        let mut prev = first.clone();
        for _ in 0..max_steps {
            let Some(tangent) = march_tangent(prev.uv0, prev.uv1) else {
                break;
            };
            // A full step near the trim border can push the solve out
            // of the padded chart; halve toward the border before
            // giving up so the chain reaches the trim.
            let mut advanced = None;
            let mut try_step = step * sign;
            for _ in 0..4 {
                let target = prev.center + tangent * try_step;
                if let Ok(sol) = solve(prev.uv0, prev.uv1, target, tangent) {
                    advanced = Some(sol);
                    break;
                }
                try_step *= 0.5;
            }
            let Some(sol) = advanced else {
                break;
            };
            // Trim-domain exit or stagnation ends the march.
            if sol.center.distance_to(&prev.center) < step * 0.1
                || on_range_border(sol.uv0, a.uv_range)
                || on_range_border(sol.uv1, b.uv_range)
            {
                out.push(sol);
                break;
            }
            out.push(sol.clone());
            prev = sol;
        }
    }
    let mut chain: Vec<ContactSolution> = backward.into_iter().rev().collect();
    chain.push(first);
    chain.extend(forward);
    if chain.len() < 2 {
        return Err(BlendError::InvalidParameter {
            parameter: "radius".into(),
            value: format!("{}", spec.radius),
            reason: "gap-medial march produced fewer than two sections".into(),
        });
    }

    let resolver = |target: Point3, warm: &ContactSolution| -> Option<ContactSolution> {
        let tangent = march_tangent(warm.uv0, warm.uv1)?;
        solve(warm.uv0, warm.uv1, target, tangent).ok()
    };
    crate::face_blend_attach::finish_face_blend(brep, spec, chain, &resolver, None, tolerance)
}

fn on_range_border(uv: (f64, f64), range: ((f64, f64), (f64, f64))) -> bool {
    let eps_u = (range.0.1 - range.0.0) * 1.0e-3;
    let eps_v = (range.1.1 - range.1.0) * 1.0e-3;
    uv.0 - range.0.0 < eps_u
        || range.0.1 - uv.0 < eps_u
        || uv.1 - range.1.0 < eps_v
        || range.1.1 - uv.1 < eps_v
}

/// Fit helper shared with the sheet builder.
pub(crate) fn polyline_nurbs(points: &[Point3]) -> Result<NurbsCurve, String> {
    crate::center_curve::fit_nurbs_through_points(points)
}
