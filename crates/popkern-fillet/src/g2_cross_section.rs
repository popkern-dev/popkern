//! G2 cross-section construction (§2.4 + §7).
//!
//! The G2 cross-section stage of the contact-problem implementation. Builds a degree-5
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
//! ## Status (MVP)
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

use crate::contact::ContactSolution;

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
    /// The converged contact solution at this spine sample. The
    /// contact normals are derived from the converged ball center
    /// (`n̂ᵢ = (contactᵢ − center)/r`), so no `BallSide` translation
    /// is needed — the cross-section is always consistent with the
    /// solution that produced it (PK-091).
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
    /// The built cross-section crosses to the non-material side of a
    /// support tangent plane — the mirror-image ("weld bead")
    /// configuration. This is the postcondition gate that catches
    /// wrong-side contact solutions (PK-091); the blend interior must
    /// stay on the ball-center side of both contact tangent planes.
    #[error(
        "G2 cross-section: interior crosses {height:.3e} beyond the support          tangent plane (non-material side) — wrong-side contact solution"
    )]
    MaterialSide {
        /// Worst signed height beyond the tangent plane (positive =
        /// non-material side).
        height: f64,
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

    // Outward-from-ball unit normals, derived from the CONVERGED
    // solution geometry (n̂ᵢ = (contactᵢ − center)/r). This is exact
    // for any converged contact and immune to orientation/side
    // convention slips — the class of defect behind PK-091.
    let radial_0 = p0 - solution.center;
    let radial_0_len = radial_0.length();
    if radial_0_len < 1.0e-12 {
        return Err(G2Error::DegenerateNormal { surface_index: 0 });
    }
    let n_hat_0 = radial_0 * (1.0 / radial_0_len);
    let radial_1 = p5 - solution.center;
    let radial_1_len = radial_1.length();
    if radial_1_len < 1.0e-12 {
        return Err(G2Error::DegenerateNormal { surface_index: 1 });
    }
    let n_hat_1 = radial_1 * (1.0 / radial_1_len);

    // Cross-section tangent at each endpoint: project the chord into
    // the surface's tangent plane (sign-invariant in n̂).
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

    // End-derivative magnitude from the tangent turn angle: a curve
    // turning by θ over the span has arc length ≈ chord·(θ/2)/sin(θ/2);
    // using the bare chord (the previous choice) under-speeds the ends
    // and concentrates curvature mid-span (measured 22% over the
    // circular value on a pentagonal-pyramid edge).
    let cos_turn = t0_unit.dot(&t1_unit).clamp(-1.0, 1.0);
    let half_turn = 0.5 * cos_turn.acos();
    let magnitude = if half_turn < 1.0e-9 {
        chord_length
    } else {
        chord_length * half_turn / half_turn.sin()
    };
    let p_prime_0 = t0_unit * magnitude;
    let p_prime_1 = t1_unit * magnitude;

    // G2 boundary condition: the cross-section's curvature VECTOR at
    // each contact must EQUAL the support's normal curvature vector
    // κ·n̂ (which is invariant under flipping n̂, so the sign is
    // unambiguous). With |P'| = magnitude, the normal component of
    // P'' must be magnitude²·κ·n̂.
    //
    // NOTE the + sign: the previous construction used −κ·n̂, i.e. the
    // NEGATION of the support's curvature vector, rationalized as
    // "the curve must bend toward the ball". That is the mid-span
    // behaviour, not the boundary condition; at the contact the curve
    // must match the SUPPORT. The error was masked on planar supports
    // (κ = 0) and by a magnitude-only acceptance test.
    let p_double_prime_0 = n_hat_0 * (magnitude * magnitude * kappa_0);
    let p_double_prime_1 = n_hat_1 * (magnitude * magnitude * kappa_1);

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
    let cps = [p0, p1, p2, p3, p4, p5];

    // Material-side postcondition gate (PK-091). The blend interior
    // must stay on the ball-center side of each contact tangent
    // plane. The midpoint check holds for every valid rolling-ball
    // blend (mid-span the curve wraps the ball); the full-curve check
    // applies only when both supports are locally flat, where the
    // support cannot legitimately carry the curve above its own
    // tangent plane.
    let tol = 1.0e-6 * chord_length;
    let mid = de_casteljau_deg5(&cps, 0.5);
    let mid_h = ((mid - p0).dot(&n_hat_0)).max((mid - p5).dot(&n_hat_1));
    if mid_h > tol {
        return Err(G2Error::MaterialSide { height: mid_h });
    }
    let flat_scale = 1.0e-9 / chord_length.max(1.0e-12);
    if kappa_0.abs() < flat_scale && kappa_1.abs() < flat_scale {
        let mut worst = f64::NEG_INFINITY;
        for step in 1..8 {
            let t = step as f64 / 8.0;
            let p = de_casteljau_deg5(&cps, t);
            worst = worst
                .max((p - p0).dot(&n_hat_0))
                .max((p - p5).dot(&n_hat_1));
        }
        if worst > tol {
            return Err(G2Error::MaterialSide { height: worst });
        }
    }

    // Degree-5 NURBS with 6 control points and a clamped uniform
    // knot vector [0,0,0,0,0,0, 1,1,1,1,1,1].
    let curve = NurbsCurve {
        degree: 5,
        knots: vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
        control_points: cps.to_vec(),
        weights: vec![1.0; 6],
    };

    Ok(G2CrossSection {
        curve,
        kappa_0,
        kappa_1,
        chord_length,
    })
}

/// Evaluate a degree-5 Bezier at `t` by de Casteljau.
fn de_casteljau_deg5(cps: &[Point3; 6], t: f64) -> Point3 {
    let mut pts: Vec<Vec3> = cps
        .iter()
        .map(|p| *p - Point3::new(0.0, 0.0, 0.0))
        .collect();
    for level in (1..6).rev() {
        for i in 0..level {
            pts[i] = pts[i] * (1.0 - t) + pts[i + 1] * t;
        }
    }
    Point3::new(pts[0].x, pts[0].y, pts[0].z)
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
/// (the cross-section degree) and a degree-3 interpolating loft in u
/// (C² along the spine). The previous degree-1 skin left C⁰ facets at
/// every spine sample, which downstream consumers (tessellation
/// normals, the KP20 vertex-blend border marches) picked up as facet
/// noise (PK-091).
///
/// Falls back to a degree-1 skin below 4 sections.
///
/// Returns `Err(G2Error)` if any per-sample G2 cross-section build
/// fails (including the material-side postcondition). Callers can
/// fall back to `build_nurbs_sweep` (the legacy G1 path) on failure.
pub fn g2_sweep_from_contact_solutions(
    solutions: &[ContactSolution],
    surface_0: &dyn SurfaceEval,
    surface_1: &dyn SurfaceEval,
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
            solution,
        };
        let cs = g2_cross_section_bezier(&req)?;
        sections.push(cs.curve);
    }

    let n_sections = sections.len();
    let n_v = 6;
    let degree_v = 5;
    let knots_v = sections[0].knots.clone();

    if n_sections < 4 {
        // Not enough sections for a cubic loft: linear skin.
        let mut control_points = Vec::with_capacity(n_sections * n_v);
        for section in &sections {
            control_points.extend(section.control_points.iter().copied());
        }
        return Ok(NurbsSurface {
            degree_u: 1,
            degree_v,
            knots_u: build_clamped_uniform_knots(n_sections, 1),
            knots_v,
            control_points,
            num_u: n_sections,
            num_v: n_v,
            weights: vec![1.0; n_sections * n_v],
        });
    }

    // Cubic interpolating loft: for each v-column of section control
    // points, solve for u-direction control points so the surface
    // passes through every section at its parameter. Uniform
    // parameters with de Boor knot averaging give a shared,
    // well-conditioned basis for all columns.
    let n = n_sections;
    let degree_u = 3;
    let params: Vec<f64> = (0..n).map(|i| i as f64 / (n - 1) as f64).collect();
    let mut knots_u = vec![0.0; 4];
    for j in 1..=(n - 4) {
        knots_u.push((params[j] + params[j + 1] + params[j + 2]) / 3.0);
    }
    knots_u.extend(std::iter::repeat_n(1.0, 4));

    // Basis matrix rows: first/last rows are identities (clamped).
    let mut matrix = vec![vec![0.0_f64; n]; n];
    matrix[0][0] = 1.0;
    matrix[n - 1][n - 1] = 1.0;
    for (k, row) in matrix.iter_mut().enumerate().take(n - 1).skip(1) {
        for (i, entry) in row.iter_mut().enumerate() {
            *entry = bspline_basis(&knots_u, i, degree_u, params[k]);
        }
    }

    // LU-free dense Gaussian elimination with partial pivoting on the
    // shared matrix, solving all 3·n_v right-hand sides at once.
    let n_rhs = 3 * n_v;
    let mut rhs = vec![vec![0.0_f64; n_rhs]; n];
    for (k, section) in sections.iter().enumerate() {
        for (j, cp) in section.control_points.iter().enumerate() {
            rhs[k][3 * j] = cp.x;
            rhs[k][3 * j + 1] = cp.y;
            rhs[k][3 * j + 2] = cp.z;
        }
    }
    for col in 0..n {
        let mut pivot = col;
        for row in (col + 1)..n {
            if matrix[row][col].abs() > matrix[pivot][col].abs() {
                pivot = row;
            }
        }
        if matrix[pivot][col].abs() < 1.0e-14 {
            return Err(G2Error::ZeroChord);
        }
        matrix.swap(col, pivot);
        rhs.swap(col, pivot);
        for row in (col + 1)..n {
            let factor = matrix[row][col] / matrix[col][col];
            if factor == 0.0 {
                continue;
            }
            let (upper, lower) = matrix.split_at_mut(row);
            for (c, entry) in lower[0].iter_mut().enumerate().skip(col) {
                *entry -= factor * upper[col][c];
            }
            let (rhs_upper, rhs_lower) = rhs.split_at_mut(row);
            for (c, entry) in rhs_lower[0].iter_mut().enumerate() {
                *entry -= factor * rhs_upper[col][c];
            }
        }
    }
    for col in (0..n).rev() {
        for c in 0..n_rhs {
            let mut acc = rhs[col][c];
            for (k, rhs_row) in rhs.iter().enumerate().take(n).skip(col + 1) {
                acc -= matrix[col][k] * rhs_row[c];
            }
            rhs[col][c] = acc / matrix[col][col];
        }
    }

    let mut control_points = Vec::with_capacity(n * n_v);
    for row in rhs.iter().take(n) {
        for j in 0..n_v {
            control_points.push(Point3::new(row[3 * j], row[3 * j + 1], row[3 * j + 2]));
        }
    }

    Ok(NurbsSurface {
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        control_points,
        num_u: n,
        num_v: n_v,
        weights: vec![1.0; n * n_v],
    })
}

/// Cox–de Boor evaluation of the `i`-th B-spline basis of degree `p`
/// over `knots` at parameter `t` (used only for interior loft rows,
/// so the half-open-interval endpoint convention never bites).
fn bspline_basis(knots: &[f64], i: usize, p: usize, t: f64) -> f64 {
    if p == 0 {
        return if knots[i] <= t && t < knots[i + 1] {
            1.0
        } else {
            0.0
        };
    }
    let mut acc = 0.0;
    let denom_left = knots[i + p] - knots[i];
    if denom_left > 1.0e-14 {
        acc += (t - knots[i]) / denom_left * bspline_basis(knots, i, p - 1, t);
    }
    let denom_right = knots[i + p + 1] - knots[i + 1];
    if denom_right > 1.0e-14 {
        acc += (knots[i + p + 1] - t) / denom_right * bspline_basis(knots, i + 1, p - 1, t);
    }
    acc
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
