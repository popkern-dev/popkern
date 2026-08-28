//! B-spline derivative evaluation with analytical Jacobians.
//!
//! Tangent, second derivative, curvature, and arc-length evaluation for
//! the [`Spline`] element.  Point evaluation and shared infrastructure
//! live in the sibling [`spline_eval`](super::spline_eval) module.

use super::spline_eval::{
    SplineEvalResult, accumulate_jac, accumulate_scalar_jac, basis_functions, find_span,
    periodic_knot_vector,
};
use crate::Sketch;
use crate::elements::Spline;

// ── Gauss-Legendre 8-point quadrature constants ──────────────────────────

/// Gauss-Legendre nodes on [-1, 1].
const GL_NODES: [f64; 8] = [
    -0.960_289_856_497_536_3,
    -0.796_666_477_413_626_7,
    -0.525_532_409_916_329,
    -0.183_434_642_495_649_8,
    0.183_434_642_495_649_8,
    0.525_532_409_916_329,
    0.796_666_477_413_626_7,
    0.960_289_856_497_536_3,
];

/// Gauss-Legendre weights on [-1, 1].
const GL_WEIGHTS: [f64; 8] = [
    0.101_228_536_290_376_3,
    0.222_381_034_453_374_5,
    0.313_706_645_877_887_3,
    0.362_683_783_378_362,
    0.362_683_783_378_362,
    0.313_706_645_877_887_3,
    0.222_381_034_453_374_5,
    0.101_228_536_290_376_3,
];

// ── First derivative (tangent) ───────────────────────────────────────────

/// Evaluate the first derivative (tangent) of the B-spline at `t`.
///
/// The derivative of a degree-`d` B-spline is a degree-`(d-1)` B-spline
/// with derived control points
/// `Q_j = d * (P_{j+1} - P_j) / (knots[j+d+1] - knots[j+1])`.
pub(crate) fn eval_tangent(spline: &Spline, t: f64, sketch: &Sketch) -> SplineEvalResult {
    let n = spline.control_points.len();
    let degree = spline.degree;

    if degree == 0 {
        return SplineEvalResult {
            x: 0.0,
            y: 0.0,
            jacobian: vec![],
        };
    }

    let knots = if spline.periodic {
        periodic_knot_vector(&spline.knots, n, degree)
    } else {
        spline.knots.clone()
    };
    let n_eff = if spline.periodic { n + degree } else { n };

    // Derivative knot vector: knots[1..n_eff+degree] (interior of original).
    let deriv_knots: Vec<f64> = knots[1..n_eff + degree].to_vec();
    let n_deriv = n_eff - 1;

    let (span, t) = find_span(&deriv_knots, degree - 1, n_deriv, t);
    let basis = basis_functions(&deriv_knots, span, degree - 1, t);

    let d = degree as f64;
    let mut tx = 0.0;
    let mut ty = 0.0;
    let mut jac_map: Vec<(usize, f64, f64)> = Vec::new();

    for (k, &nk) in basis.iter().enumerate().take(degree) {
        let q_idx = span - (degree - 1) + k;
        let orig_j = q_idx;
        let denom = knots[orig_j + degree + 1] - knots[orig_j + 1];

        let coeff = if denom.abs() > f64::EPSILON {
            d / denom
        } else {
            0.0
        };

        let cp_idx_a = if spline.periodic { orig_j % n } else { orig_j };
        let cp_idx_b = if spline.periodic {
            (orig_j + 1) % n
        } else {
            orig_j + 1
        };

        let pt_a = &sketch.points[spline.control_points[cp_idx_a]];
        let pt_b = &sketch.points[spline.control_points[cp_idx_b]];
        let [iax, iay] = pt_a.param_indices();
        let [ibx, iby] = pt_b.param_indices();
        let ax = sketch.params[iax];
        let ay = sketch.params[iay];
        let bx = sketch.params[ibx];
        let by = sketch.params[iby];

        let qx = coeff * (bx - ax);
        let qy = coeff * (by - ay);
        tx += nk * qx;
        ty += nk * qy;

        let dc_a = -nk * coeff;
        let dc_b = nk * coeff;

        accumulate_jac(&mut jac_map, iax, dc_a, 0.0);
        accumulate_jac(&mut jac_map, iay, 0.0, dc_a);
        accumulate_jac(&mut jac_map, ibx, dc_b, 0.0);
        accumulate_jac(&mut jac_map, iby, 0.0, dc_b);
    }

    jac_map.retain(|(_, dx, dy)| dx.abs() > f64::EPSILON || dy.abs() > f64::EPSILON);

    SplineEvalResult {
        x: tx,
        y: ty,
        jacobian: jac_map,
    }
}

// ── Second derivative ────────────────────────────────────────────────────

/// Evaluate the second derivative of the B-spline at `t`.
///
/// Applies the derivative formula twice: the second derivative is a
/// degree-`(d-2)` B-spline with control points derived from the first
/// derivative control points.
pub(crate) fn eval_second_deriv(spline: &Spline, t: f64, sketch: &Sketch) -> SplineEvalResult {
    let n = spline.control_points.len();
    let degree = spline.degree;

    if degree <= 1 {
        return SplineEvalResult {
            x: 0.0,
            y: 0.0,
            jacobian: vec![],
        };
    }

    let knots = if spline.periodic {
        periodic_knot_vector(&spline.knots, n, degree)
    } else {
        spline.knots.clone()
    };
    let n_eff = if spline.periodic { n + degree } else { n };
    let d = degree as f64;

    // Second-derivative knot vector: knots[2..n_eff+degree-1].
    let deriv2_knots: Vec<f64> = knots[2..n_eff + degree - 1].to_vec();
    let n_deriv2 = n_eff - 2;

    if n_deriv2 == 0 {
        return SplineEvalResult {
            x: 0.0,
            y: 0.0,
            jacobian: vec![],
        };
    }

    let (span, t) = find_span(&deriv2_knots, degree - 2, n_deriv2, t);
    let basis = basis_functions(&deriv2_knots, span, degree - 2, t);

    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut jac_map: Vec<(usize, f64, f64)> = Vec::new();

    for (k, &nk) in basis.iter().enumerate().take(degree - 1) {
        let r_idx = span - (degree - 2) + k;
        let orig_r = r_idx;

        let denom2 = knots[orig_r + degree + 1] - knots[orig_r + 2];
        let coeff2 = if denom2.abs() > f64::EPSILON {
            (d - 1.0) / denom2
        } else {
            0.0
        };

        let cp_indices: [usize; 3] = if spline.periodic {
            [orig_r % n, (orig_r + 1) % n, (orig_r + 2) % n]
        } else {
            [orig_r, orig_r + 1, orig_r + 2]
        };

        let mut px = [0.0; 3];
        let mut py = [0.0; 3];
        let mut ix = [0usize; 3];
        let mut iy = [0usize; 3];
        for (m, &ci) in cp_indices.iter().enumerate() {
            let pt = &sketch.points[spline.control_points[ci]];
            let [i_x, i_y] = pt.param_indices();
            ix[m] = i_x;
            iy[m] = i_y;
            px[m] = sketch.params[i_x];
            py[m] = sketch.params[i_y];
        }

        let denom_q0 = knots[orig_r + degree + 1] - knots[orig_r + 1];
        let denom_q1 = knots[orig_r + degree + 2] - knots[orig_r + 2];

        let coeff_q0 = if denom_q0.abs() > f64::EPSILON {
            d / denom_q0
        } else {
            0.0
        };
        let coeff_q1 = if denom_q1.abs() > f64::EPSILON {
            d / denom_q1
        } else {
            0.0
        };

        // R = coeff2 * (Q1 - Q0)
        //   = coeff2 * (coeff_q1*(P2-P1) - coeff_q0*(P1-P0))
        //   = coeff2 * (coeff_q0*P0 - (coeff_q0+coeff_q1)*P1 + coeff_q1*P2)
        let c0 = coeff2 * coeff_q0;
        let c1 = -coeff2 * (coeff_q0 + coeff_q1);
        let c2 = coeff2 * coeff_q1;

        let rx = c0 * px[0] + c1 * px[1] + c2 * px[2];
        let ry = c0 * py[0] + c1 * py[1] + c2 * py[2];
        sx += nk * rx;
        sy += nk * ry;

        let coeffs = [nk * c0, nk * c1, nk * c2];
        for (m, &c) in coeffs.iter().enumerate() {
            accumulate_jac(&mut jac_map, ix[m], c, 0.0);
            accumulate_jac(&mut jac_map, iy[m], 0.0, c);
        }
    }

    jac_map.retain(|(_, dx, dy)| dx.abs() > f64::EPSILON || dy.abs() > f64::EPSILON);

    SplineEvalResult {
        x: sx,
        y: sy,
        jacobian: jac_map,
    }
}

// ── Curvature ────────────────────────────────────────────────────────────

/// Compute signed curvature at parameter `t`.
///
/// `kappa = (x'y'' - y'x'') / (x'^2 + y'^2)^(3/2)`
///
/// Returns `(curvature, sparse_jacobian)` where the Jacobian maps
/// `param_index -> d(kappa)/d(param)`.
pub(crate) fn eval_curvature(spline: &Spline, t: f64, sketch: &Sketch) -> (f64, Vec<(usize, f64)>) {
    let tang = eval_tangent(spline, t, sketch);
    let accel = eval_second_deriv(spline, t, sketch);

    let xp = tang.x;
    let yp = tang.y;
    let xpp = accel.x;
    let ypp = accel.y;

    let cross = xp * ypp - yp * xpp;
    let speed_sq = xp * xp + yp * yp;
    let speed = speed_sq.sqrt();
    let speed_cubed = speed * speed_sq;

    if speed_cubed < f64::EPSILON {
        return (0.0, vec![]);
    }

    let kappa = cross / speed_cubed;
    let inv_s3 = 1.0 / speed_cubed;
    let coeff_speed = 3.0 * kappa / speed_sq;

    let mut jac_map: Vec<(usize, f64)> = Vec::new();

    // Contributions from tangent derivatives.
    for &(idx, dxp_dp, dyp_dp) in &tang.jacobian {
        let dcross = ypp * dxp_dp - xpp * dyp_dp;
        let dspeed_sq = xp * dxp_dp + yp * dyp_dp;
        let dk = inv_s3 * dcross - coeff_speed * dspeed_sq;
        accumulate_scalar_jac(&mut jac_map, idx, dk);
    }

    // Contributions from second-derivative derivatives.
    for &(idx, dxpp_dp, dypp_dp) in &accel.jacobian {
        let dcross = xp * dypp_dp - yp * dxpp_dp;
        let dk = inv_s3 * dcross;
        accumulate_scalar_jac(&mut jac_map, idx, dk);
    }

    jac_map.retain(|(_, v)| v.abs() > f64::EPSILON);
    (kappa, jac_map)
}

// ── Arc length ───────────────────────────────────────────────────────────

/// Compute arc length of the spline from `t0` to `t1` using 8-point
/// Gauss-Legendre quadrature.
///
/// Returns `(length, sparse_jacobian)` where the Jacobian maps
/// `param_index -> d(length)/d(param)`.
pub(crate) fn eval_arc_length(
    spline: &Spline,
    t0: f64,
    t1: f64,
    sketch: &Sketch,
) -> (f64, Vec<(usize, f64)>) {
    let half_dt = (t1 - t0) / 2.0;
    let mid_t = (t1 + t0) / 2.0;
    let mut length = 0.0;
    let mut jac_map: Vec<(usize, f64)> = Vec::new();

    for i in 0..8 {
        let ti = half_dt * GL_NODES[i] + mid_t;
        let tang = eval_tangent(spline, ti, sketch);
        let speed = (tang.x * tang.x + tang.y * tang.y).sqrt();
        let wi = GL_WEIGHTS[i] * half_dt;
        length += wi * speed;

        if speed > f64::EPSILON {
            let inv_speed = 1.0 / speed;
            for &(idx, dxp, dyp) in &tang.jacobian {
                let d_speed = (tang.x * dxp + tang.y * dyp) * inv_speed;
                accumulate_scalar_jac(&mut jac_map, idx, wi * d_speed);
            }
        }
    }

    jac_map.retain(|(_, v)| v.abs() > f64::EPSILON);
    (length, jac_map)
}
