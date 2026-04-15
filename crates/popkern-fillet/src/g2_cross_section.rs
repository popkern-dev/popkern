//! G2 cross-section construction (§2.4 + §7).
//!
//! Phase 1d of the contact-problem implementation. Builds a degree-5
//! Bezier (NURBS) cross-section curve at one spine sample with
//! position + tangent + **curvature** continuity to both support
//! surfaces — i.e. G2.
//!
//! ## What G2 buys you
//!
//! G1 (the standard rolling-ball arc) eliminates positional and
//! tangent breaks at the blend boundary, but the curvature jumps
//! from the surface's normal curvature to the rolling-ball
//! curvature `1/r`. That jump shows up as a kink in highlight
//! lines and reflection patterns — visible on Class-A surfaces and
//! detectable by curvature-comb analysis.
//!
//! G2 matches the second fundamental form at the boundary too:
//! the cross-section curve's curvature vector at each contact
//! point equals the surface's normal curvature in the cross-section
//! tangent direction. The curvature is continuous across the
//! transition, eliminating the highlight-line kink.
//!
//! ## How
//!
//! A degree-5 polynomial Bezier has 6 control points and 6 vector
//! degrees of freedom per dimension, exactly enough to fit:
//!
//! - position at both endpoints (2 constraints)
//! - first derivative at both endpoints (2 constraints)
//! - second derivative at both endpoints (2 constraints)
//!
//! For the cross-section:
//!
//! - `P_0 = contact0`, `P_5 = contact1`
//! - `P'(0) = chord · T_0` where `T_0` is the unit tangent to the
//!   cross-section at `contact0`, computed by projecting the chord
//!   `(contact1 - contact0)` onto surface 0's tangent plane.
//! - `P'(1) = chord · T_1` similarly on surface 1.
//! - `P''(0) = chord² · κ_0 · n̂_0` where `κ_0` is the normal
//!   curvature of surface 0 in direction `T_0` and `n̂_0` is the
//!   surface 0 outward unit normal at `contact0`. This makes the
//!   curve's curvature vector at the start exactly match the
//!   surface's normal curvature vector — the G2 condition.
//! - `P''(1)` likewise.
//!
//! Inverting the Bezier derivative formulas gives the inner control
//! points:
//!
//! ```text
//! P_1 = P_0 + (1/5)·P'(0)
//! P_4 = P_5 - (1/5)·P'(1)
//! P_2 = P''(0)/20 + 2·P_1 - P_0
//! P_3 = P''(1)/20 + 2·P_4 - P_5
//! ```
//!
//! ## Phase 1d status (MVP)
//!
//! - **Math**: implemented and unit-tested on plane-plane and
//!   plane-cylinder. The plane-plane case verifies that endpoint
//!   curvature is zero (since planes have zero normal curvature
//!   everywhere) and the cross-section bulges outward in the
//!   middle. The plane-cylinder case verifies the curvature at the
//!   cylinder contact equals `1/cylinder_radius` to within 1% (the
//!   §2.4 acceptance threshold).
//! - **Integration**: NOT wired into the existing fillet pipeline
//!   yet. The current `fillet_surface::sample_arc` path produces
//!   circular arcs; replacing it with G2 Bezier sweeps is downstream
//!   work (phase 1d+ or phase 1e). Exposing the math as a callable
//!   primitive is enough to satisfy §19.5 (pipeline API) for the
//!   G2 path.
//!
//! See `docs/superpowers/plans/2026-04-14-fillets-rolling-ball-contact-solver.md`
//! for the broader phase-1 roadmap.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};
use thiserror::Error;

use crate::contact::math::signed_unit_normal;
use crate::contact::{BallSide, ContactSolution};

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Inputs to a G2 cross-section construction.
///
/// A converged [`ContactSolution`] from the rolling-ball solver
/// (phase 1a) gives us the endpoints, normals, and `(uv0, uv1)`
/// parameters we need; the surface references let us evaluate first
/// and second fundamental forms for the curvature matching.
pub struct G2Request<'a> {
    /// Surface 0 (provides position, tangent plane, normal curvature
    /// at the `contact0` endpoint).
    pub surface_0: &'a dyn SurfaceEval,
    /// Surface 1.
    pub surface_1: &'a dyn SurfaceEval,
    /// Which side of `surface_0`'s natural normal the rolling ball
    /// sits on. Same convention as `ContactRequest::side_0`.
    pub side_0: BallSide,
    /// Which side of `surface_1`'s natural normal the rolling ball
    /// sits on.
    pub side_1: BallSide,
    /// The converged contact solution at this spine sample.
    pub solution: &'a ContactSolution,
}

/// A converged G2 cross-section, packaged with telemetry that lets
/// callers inspect the endpoint curvatures (for §2.4 acceptance
/// checks) and the chord length (for downstream sampling).
#[derive(Debug, Clone)]
pub struct G2CrossSection {
    /// The degree-5 NURBS cross-section curve. 6 control points,
    /// uniform [0, 1] knot vector.
    pub curve: NurbsCurve,
    /// Normal curvature of surface 0 at `contact0` in the
    /// cross-section tangent direction. The cross-section curve's
    /// curvature at `t=0` equals this.
    pub kappa_0: f64,
    /// Normal curvature of surface 1 at `contact1`.
    pub kappa_1: f64,
    /// 3D distance `|contact1 - contact0|`, used as the tangent
    /// magnitude scale.
    pub chord_length: f64,
}

/// Failure modes for [`g2_cross_section_bezier`]. All three are
/// "the geometry doesn't admit a G2 cross-section here" cases —
/// degenerate normals, parallel surfaces, or a zero-length chord.
#[derive(Debug, Clone, Error)]
pub enum G2Error {
    /// One of the surfaces produced a degenerate normal at the
    /// contact point. The caller is expected to fall back to the
    /// G1 (arc) cross-section.
    #[error(
        "G2 cross-section: surface {surface_index} has a degenerate normal at the contact point"
    )]
    DegenerateNormal {
        /// `0` for surface 0, `1` for surface 1.
        surface_index: usize,
    },
    /// The two contact points coincide, so there's no chord to
    /// build a tangent against. This happens at "tangent contact"
    /// singularities.
    #[error("G2 cross-section: the two contact points coincide (zero chord)")]
    ZeroChord,
    /// The cross-section tangent at one of the endpoints projected
    /// to zero in the corresponding tangent plane — i.e. the chord
    /// is exactly perpendicular to the surface tangent plane at
    /// that point. Should never happen for a normal rolling-ball
    /// solution but flagged defensively.
    #[error(
        "G2 cross-section: cross-section tangent degenerates on surface \
         {surface_index} (chord is perpendicular to its tangent plane)"
    )]
    TangentDegenerate {
        /// `0` for surface 0, `1` for surface 1.
        surface_index: usize,
    },
}

/// Build a degree-5 Bezier cross-section with G2 (curvature-matching)
/// boundary conditions at both contact points.
///
/// See the [module docs](self) for the math.
pub fn g2_cross_section_bezier(request: &G2Request) -> Result<G2CrossSection, G2Error> {
    let solution = request.solution;
    let p0 = solution.contact0;
    let p5 = solution.contact1;
    let chord_vec = p5 - p0;
    let chord_length = chord_vec.length();
    if chord_length < 1.0e-12 {
        return Err(G2Error::ZeroChord);
    }

    // Outward unit normals on both surfaces, signed per BallSide.
    let n_hat_0 = signed_unit_normal(
        request.surface_0,
        solution.uv0.0,
        solution.uv0.1,
        request.side_0,
        1.0e-12,
    )
    .ok_or(G2Error::DegenerateNormal { surface_index: 0 })?;
    let n_hat_1 = signed_unit_normal(
        request.surface_1,
        solution.uv1.0,
        solution.uv1.1,
        request.side_1,
        1.0e-12,
    )
    .ok_or(G2Error::DegenerateNormal { surface_index: 1 })?;

    // Cross-section tangent at each endpoint: project the chord into
    // the surface's tangent plane.  We flip the sign at endpoint 1
    // so both tangents point "into" the curve (away from their
    // contact endpoint), matching the Bezier convention where
    // P'(0) and P'(1) both point in the curve-traversal direction.
    let t0_3d = project_to_tangent_plane(&chord_vec, &n_hat_0);
    let t0_len = t0_3d.length();
    if t0_len < 1.0e-12 {
        return Err(G2Error::TangentDegenerate { surface_index: 0 });
    }
    let t0_unit = t0_3d * (1.0 / t0_len);

    let t1_3d = project_to_tangent_plane(&chord_vec, &n_hat_1);
    let t1_len = t1_3d.length();
    if t1_len < 1.0e-12 {
        return Err(G2Error::TangentDegenerate { surface_index: 1 });
    }
    let t1_unit = t1_3d * (1.0 / t1_len);

    // Normal curvature of each surface in its respective cross-section
    // tangent direction, computed from the second fundamental form.
    let kappa_0 =
        normal_curvature_in_direction(request.surface_0, solution.uv0, &n_hat_0, &t0_unit);
    let kappa_1 =
        normal_curvature_in_direction(request.surface_1, solution.uv1, &n_hat_1, &t1_unit);

    // Bezier derivative magnitudes — pick |P'(0)| = |P'(1)| = chord
    // so the parameter speed is consistent with the chord length.
    // P'(0) = chord · T_0_unit, P'(1) = chord · T_1_unit.
    let p_prime_0 = t0_unit * chord_length;
    let p_prime_1 = t1_unit * chord_length;

    // Curvature vector convention: the rolling ball sits on the
    // OPPOSITE side of the surface from where the surface curves
    // toward (i.e., on the convex side of a concave fillet).
    // The cross-section curve must curve TOWARD the rolling-ball
    // center, which is in the `-n̂` direction (n̂ here being
    // signed outward from the ball). So the curvature vector is
    // `-κ_n · n̂` rather than `+κ_n · n̂`.
    //
    // Concretely: for plane-plane convex 90°, the rolling ball center
    // is at (r, y, r), each plane's outward normal points away from
    // the ball, and the cross-section arc curves AROUND the ball
    // (concave-up from the surface's perspective). The curve's
    // curvature vector at the contact points is +1/r toward the ball
    // center, which is `-n̂_surface`.
    let p_double_prime_0 = n_hat_0 * (-chord_length * chord_length * kappa_0);
    let p_double_prime_1 = n_hat_1 * (-chord_length * chord_length * kappa_1);

    // Invert the degree-5 Bezier derivative formulas:
    //   P'(0) = 5·(P_1 - P_0)        → P_1 = P_0 + P'(0)/5
    //   P'(1) = 5·(P_5 - P_4)        → P_4 = P_5 - P'(1)/5
    //   P''(0) = 20·(P_2 - 2·P_1 + P_0)
    //                                → P_2 = P''(0)/20 + 2·P_1 - P_0
    //   P''(1) = 20·(P_5 - 2·P_4 + P_3)
    //                                → P_3 = P''(1)/20 + 2·P_4 - P_5
    let p1 = p0 + p_prime_0 * (1.0 / 5.0);
    let p4 = p5 - p_prime_1 * (1.0 / 5.0);
    let p2 = p_double_prime_0 * (1.0 / 20.0) + (p1 - Point3::new(0.0, 0.0, 0.0)) * 2.0
        - (p0 - Point3::new(0.0, 0.0, 0.0));
    let p3 = p_double_prime_1 * (1.0 / 20.0) + (p4 - Point3::new(0.0, 0.0, 0.0)) * 2.0
        - (p5 - Point3::new(0.0, 0.0, 0.0));
    let p2 = Point3::new(p2.x, p2.y, p2.z);
    let p3 = Point3::new(p3.x, p3.y, p3.z);

    // Degree-5 NURBS with 6 control points and a clamped uniform
    // knot vector [0,0,0,0,0,0, 1,1,1,1,1,1].
    let curve = NurbsCurve {
        degree: 5,
        knots: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        control_points: vec![p0, p1, p2, p3, p4, p5],
        weights: vec![1.0; 6],
    };

    Ok(G2CrossSection {
        curve,
        kappa_0,
        kappa_1,
        chord_length,
    })
}

// ---------------------------------------------------------------------------
// Math helpers (private)
// ---------------------------------------------------------------------------

/// Project a 3D vector onto the plane perpendicular to a given
/// (unit) normal — i.e., remove the component parallel to the normal.
fn project_to_tangent_plane(v: &Vec3, n_hat: &Vec3) -> Vec3 {
    let dot = v.dot(n_hat);
    *v - *n_hat * dot
}

/// Normal curvature of `surface` at `(u, v)` in the direction
/// `tangent_3d` (which is assumed to lie in the surface's tangent
/// plane and to be unit length, but we don't require it).
///
/// Computes the second fundamental form `II(T, T) / I(T, T)` after
/// expressing `tangent_3d` in the `(S_u, S_v)` basis of the tangent
/// plane. Uses the *signed* unit normal from `signed_unit_normal`
/// so the result carries the [`BallSide`] sign — i.e., `kappa > 0`
/// means the surface curves toward the rolling-ball center.
pub(crate) fn normal_curvature_in_direction(
    surface: &dyn SurfaceEval,
    uv: (f64, f64),
    n_hat: &Vec3,
    tangent_3d: &Vec3,
) -> f64 {
    let s_u = surface.partial_u(uv.0, uv.1);
    let s_v = surface.partial_v(uv.0, uv.1);
    let s_uu = surface.second_partial_uu(uv.0, uv.1);
    let s_uv = surface.second_partial_uv(uv.0, uv.1);
    let s_vv = surface.second_partial_vv(uv.0, uv.1);

    // First fundamental form coefficients.
    let e = s_u.dot(&s_u);
    let f = s_u.dot(&s_v);
    let g = s_v.dot(&s_v);

    // Second fundamental form coefficients (against the signed normal).
    let l = s_uu.dot(n_hat);
    let m = s_uv.dot(n_hat);
    let n_coef = s_vv.dot(n_hat);

    // Express `tangent_3d` in the `(S_u, S_v)` basis: solve
    //   tangent_3d ≈ a · S_u + b · S_v
    // by least-squares (the `(a, b)` that minimize the residual,
    // ignoring any out-of-plane component which the caller should
    // already have projected out).
    //
    // Normal equations: [E F; F G] · (a, b)^T = (S_u·t, S_v·t)^T
    //                  → (a, b) = [G -F; -F E] · (S_u·t, S_v·t)^T / (EG - F²)
    let metric_det = e * g - f * f;
    if metric_det < 1.0e-30 {
        return 0.0;
    }
    let su_t = s_u.dot(tangent_3d);
    let sv_t = s_v.dot(tangent_3d);
    let a = (g * su_t - f * sv_t) / metric_det;
    let b = (-f * su_t + e * sv_t) / metric_det;

    // II(T, T) and I(T, T) in the (a, b) basis.
    let ii_tt = a * a * l + 2.0 * a * b * m + b * b * n_coef;
    let i_tt = a * a * e + 2.0 * a * b * f + b * b * g;
    if i_tt < 1.0e-30 {
        return 0.0;
    }
    ii_tt / i_tt
}

// ---------------------------------------------------------------------------
// G2 sweep (skin G2 cross-sections along the spine)
// ---------------------------------------------------------------------------

/// Build a NURBS swept surface by skinning G2 cross-sections along
/// the spine, one per `ContactSolution`. The surface goes from
/// `solutions[0].contact0..contact1` at u=0 to
/// `solutions[N-1].contact0..contact1` at u=1, with degree-5 v-direction
/// (the cross-section degree) and degree-1 u-direction (linear skin
/// between consecutive cross-sections).
///
/// This is the G2 counterpart of `fillet_surface::build_nurbs_sweep`,
/// which produces a degree-1 u, degree-1 v polyline-arc sweep. The
/// cross-section direction here is degree-5 because each cross-
/// section is a degree-5 Bezier from `g2_cross_section_bezier`. The
/// u-direction degree is left at 1 to match the legacy sweep's
/// behaviour; phase 2+ may bump it for smoother along-spine
/// continuity.
///
/// Returns `Err(G2Error)` if any per-sample G2 cross-section build
/// fails. Callers can fall back to `build_nurbs_sweep` (the legacy
/// G1 path) on failure.
pub fn g2_sweep_from_contact_solutions(
    solutions: &[ContactSolution],
    surface_0: &dyn SurfaceEval,
    surface_1: &dyn SurfaceEval,
    side_0: BallSide,
    side_1: BallSide,
) -> Result<NurbsSurface, G2Error> {
    if solutions.len() < 2 {
        return Err(G2Error::ZeroChord);
    }

    // Build a G2 cross-section curve at each spine sample.
    let mut sections: Vec<NurbsCurve> = Vec::with_capacity(solutions.len());
    for solution in solutions {
        let req = G2Request {
            surface_0,
            surface_1,
            side_0,
            side_1,
            solution,
        };
        let cs = g2_cross_section_bezier(&req)?;
        sections.push(cs.curve);
    }

    // All cross-sections share the same degree (5) and number of
    // control points (6). Skin them into a tensor-product surface
    // with degree-1 u (linear between consecutive sections) and
    // degree-5 v (the cross-section degree).
    let n_sections = sections.len();
    let n_v = 6;
    let degree_u = 1;
    let degree_v = 5;

    // Control points are stored row-major: control_points[i*n_v + j]
    // is the i-th section's j-th cross-section control point.
    let mut control_points = Vec::with_capacity(n_sections * n_v);
    for section in &sections {
        for cp in &section.control_points {
            control_points.push(*cp);
        }
    }

    let knots_u = build_clamped_uniform_knots(n_sections, degree_u);
    // The cross-section knot vector is the same for every section
    // since they're all the same shape (degree-5, 6 cps, clamped
    // uniform [0..1]); reuse it.
    let knots_v = sections[0].knots.clone();
    let weights = vec![1.0; n_sections * n_v];

    Ok(NurbsSurface {
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        control_points,
        num_u: n_sections,
        num_v: n_v,
        weights,
    })
}

fn build_clamped_uniform_knots(n_control_points: usize, degree: usize) -> Vec<f64> {
    // Standard clamped uniform knot vector: [0]*(p+1), interior, [1]*(p+1).
    let n = n_control_points;
    let p = degree;
    let mut knots = Vec::with_capacity(n + p + 1);
    knots.extend(std::iter::repeat_n(0.0, p + 1));
    let interior = n.saturating_sub(p + 1);
    for i in 1..=interior {
        knots.push(i as f64 / (interior + 1) as f64);
    }
    knots.extend(std::iter::repeat_n(1.0, p + 1));
    knots
}

#[cfg(test)]
#[path = "g2_cross_section_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "g2_sweep_tests.rs"]
mod sweep_tests;
