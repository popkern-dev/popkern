//! Exact fillet station geometry (PK-177).
//!
//! Shared solvers for EXACT station rows: the alternating-offset
//! centre solve (equidistant `r` from both supports, ~1e-10) and
//! the two-segment piecewise-rational arc row for sweeps beyond the
//! single-segment rational range. The chart extension and the loft
//! builder both consume these so every station — interior,
//! extension, wide-arc — is exact station geometry rather than
//! walker-accuracy samples.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_math::PortableMath;

/// Re-solve one station's support contacts exactly: alternating
/// offset projection to the centre equidistant `radius` from both
/// supports (seeded between the walker footprints), then the exact
/// foot points. `None` when the solve does not converge or wanders
/// more than half a radius from the seeds.
///
/// Returns the contacts with the solved ball centre: the row
/// built on these three is exact (both contacts `radius` from the
/// centre). A centre sampled from the fitted centre curve instead is
/// another station's — on a cylinder tee's saddle the rows' centres
/// sat 6.03–6.12 from their r = 6 contacts, the exact rational row was
/// refused, and the sampled skin left the rails 6.5e-3 off (PK-734).
pub(crate) fn refine_station(
    s1: &SurfaceGeometry,
    s2: &SurfaceGeometry,
    fp1: &Point3,
    fp2: &Point3,
    radius: f64,
) -> Option<(Point3, Point3, Point3)> {
    use popkern_core::geometry::SurfaceEval;
    let project = |sg: &SurfaceGeometry, p: Point3| -> Option<Point3> {
        let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(sg, &p)?;
        Some(sg.point_at(u, v))
    };
    let seed = Point3::new(
        0.5 * (fp1.x + fp2.x),
        0.5 * (fp1.y + fp2.y),
        0.5 * (fp1.z + fp2.z),
    );
    let mut centre = seed;
    let mut converged = false;
    for _ in 0..200 {
        let q1 = project(s1, centre)?;
        let d1 = centre.distance_to(&q1);
        if d1 > 1e-12 {
            centre = q1 + (centre - q1) * (radius / d1);
        }
        let q2 = project(s2, centre)?;
        let d2 = centre.distance_to(&q2);
        if d2 > 1e-12 {
            centre = q2 + (centre - q2) * (radius / d2);
        }
        let r1 = (centre.distance_to(&project(s1, centre)?) - radius).abs();
        let r2 = (centre.distance_to(&project(s2, centre)?) - radius).abs();
        if r1.max(r2) < 1e-10 {
            converged = true;
            break;
        }
    }
    if !converged {
        return None;
    }
    let c1 = project(s1, centre)?;
    let c2 = project(s2, centre)?;
    // Containment: an exact solve should land within footprint
    // accuracy of the walker's contacts; a slide along the support
    // means a different station solved.
    if c1.distance_to(fp1) > radius * 0.5 || c2.distance_to(fp2) > radius * 0.5 {
        return None;
    }
    Some((c1, c2, centre))
}

/// One EXACT two-segment piecewise-rational arc row from `start` to
/// `end` about `center`: control points [start, M1, mid, M2, end]
/// with segment-mid weights cos(θ/4) and a doubled interior knot at
/// 0.5. Exact for sweeps the single-segment row declines
/// (~115°..170°); `None` outside its own well-conditioned range.
pub(crate) fn bi_rational_arc_row(
    center: &Point3,
    start: &Point3,
    end: &Point3,
) -> Option<([Point3; 5], f64)> {
    let v1 = *start - *center;
    let v2 = *end - *center;
    let r1 = v1.length();
    let r2 = v2.length();
    if r1 < 1e-12 || r2 < 1e-12 {
        return None;
    }
    // Row is exact for the MEAN circle; a footprint radius mismatch
    // of ε leaves an ε/2 residual. 1e-3 relative admits walker
    // footprints (~1e-5 class) at a sub-precision residual — still
    // an order under the sampled skin's chordal error.
    if (r1 - r2).abs() > 1e-3 * r1.max(r2) {
        return None;
    }
    let u1 = v1 * (1.0 / r1);
    let u2 = v2 * (1.0 / r2);
    let cos_theta = u1.dot(&u2).clamp(-1.0, 1.0);
    let theta = cos_theta.pacos();
    if !(1e-9..=3.0).contains(&theta) {
        return None;
    }
    let r = 0.5 * (r1 + r2);
    let bis = (u1 + u2).normalized()?;
    let mid = *center + bis * r;
    let w = (theta * 0.25).pcos();
    let m1 = *center + (u1 + bis).normalized()? * (r / w);
    let m2 = *center + (u2 + bis).normalized()? * (r / w);
    Some(([*start, m1, mid, m2, *end], w))
}

/// Degree-1 polyline NURBS through the station points (clamped
/// uniform knots) — the exact rail of a degree-1 loft through the
/// same stations.
pub(crate) fn station_polyline(
    pts: &[Point3],
) -> Result<popkern_core::geometry::NurbsCurve, String> {
    if pts.len() < 2 {
        return Err("polyline trim needs at least 2 stations".into());
    }
    Ok(popkern_core::geometry::NurbsCurve {
        degree: 1,
        knots: crate::fillet_surface::build_clamped_uniform_knots(pts.len(), 1),
        control_points: pts.to_vec(),
        weights: vec![1.0; pts.len()],
    })
}

/// Degree-2 B-spline basis functions at `t` for clamped knots `u`
/// (n control points): returns the (at most 3) nonzero values as
/// (first_index, [n0, n1, n2]).
fn basis2(u: &[f64], n: usize, t: f64) -> (usize, [f64; 3]) {
    // Find knot span k with u[k] <= t < u[k+1], within [2, n-1].
    let mut k = 2;
    while k + 1 < n && t >= u[k + 1] {
        k += 1;
    }
    let mut left = [0.0_f64; 3];
    let mut right = [0.0_f64; 3];
    let mut vals = [1.0_f64, 0.0, 0.0];
    for j in 1..=2 {
        left[j] = t - u[k + 1 - j];
        right[j] = u[k + j] - t;
        let mut saved = 0.0;
        for r in 0..j {
            let denom = right[r + 1] + left[j - r];
            let temp = if denom.abs() > 1e-300 {
                vals[r] / denom
            } else {
                0.0
            };
            vals[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        vals[j] = saved;
    }
    (k - 2, vals)
}

/// Interpolating degree-2 loft through exact station rows (by chord
/// length along the rails): solves,
/// per v-column in HOMOGENEOUS space, for control points whose
/// degree-2 B-spline passes through every row (a plain degree-2
/// B-spline does not interpolate its control points — using rows as
/// CPs would lose G0 at the stations). Rows must share one arc
/// layout (m points, m weights per row).
pub(crate) fn interpolating_loft(
    rows: &[Vec<Point3>],
    row_weights: &[Vec<f64>],
) -> Option<popkern_core::geometry::NurbsSurface> {
    let n = rows.len();
    if n < 3 || rows.len() != row_weights.len() {
        return None;
    }
    let m = rows[0].len();
    if m < 3 || rows.iter().any(|r| r.len() != m) || row_weights.iter().any(|w| w.len() != m) {
        return None;
    }
    // Chord-length interpolation parameters (the mean of the first and
    // last columns' chords: the two rails) + averaged clamped knots.
    // Uniform parameters on unevenly spaced stations bend the strip
    // between them: an equal tee's r 8 saddle strips, whose stations
    // bunch toward their vanishing tips, parted 9.6e-5 from their
    // supports there (PK-767), where chord length keeps them under the
    // 1e-5 seam gate.
    let mut params = Vec::with_capacity(n);
    let mut arc = 0.0;
    for i in 0..n {
        if i > 0 {
            let chord = |j: usize| rows[i][j].distance_to(&rows[i - 1][j]);
            arc += 0.5 * (chord(0) + chord(m - 1));
        }
        params.push(arc);
    }
    // A repeated station has no chord to measure: uniform, as before.
    if arc > 0.0 && params.windows(2).all(|w| w[1] > w[0]) {
        for p in params.iter_mut() {
            *p /= arc;
        }
    } else {
        params = (0..n).map(|i| (i as f64) / ((n - 1) as f64)).collect();
    }
    let mut knots = vec![0.0; n + 3];
    for j in 1..(n - 2) {
        knots[j + 2] = 0.5 * (params[j] + params[j + 1]);
    }
    for k in knots.iter_mut().skip(n) {
        *k = 1.0;
    }
    // Dense n x n basis matrix (n is tens at most).
    let mut a = vec![0.0_f64; n * n];
    for (i, &t) in params.iter().enumerate() {
        let (first, vals) = basis2(&knots, n, t.min(1.0 - 1e-12));
        for (o, v) in vals.iter().enumerate() {
            if first + o < n {
                a[i * n + first + o] = *v;
            }
        }
    }
    // Homogeneous RHS per column: (wx, wy, wz, w).
    let mut rhs = vec![[0.0_f64; 4]; n * m];
    for i in 0..n {
        for j in 0..m {
            let w = row_weights[i][j];
            let p = rows[i][j];
            rhs[i * m + j] = [p.x * w, p.y * w, p.z * w, w];
        }
    }
    // Gaussian elimination with partial pivoting on the shared
    // matrix, applied to all m columns x 4 components.
    let mut perm: Vec<usize> = (0..n).collect();
    for col in 0..n {
        let (piv, _) = (col..n)
            .map(|r| (r, a[perm[r] * n + col].abs()))
            .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap_or(std::cmp::Ordering::Equal))?;
        perm.swap(col, piv);
        let p_row = perm[col];
        let diag = a[p_row * n + col];
        if diag.abs() < 1e-14 {
            return None;
        }
        for &row in perm.iter().skip(col + 1) {
            let f = a[row * n + col] / diag;
            if f == 0.0 {
                continue;
            }
            for c in col..n {
                a[row * n + c] -= f * a[p_row * n + c];
            }
            for j in 0..m {
                let pivot_vals = rhs[p_row * m + j];
                let target = &mut rhs[row * m + j];
                for (t, pv) in target.iter_mut().zip(pivot_vals.iter()) {
                    *t -= f * pv;
                }
            }
        }
    }
    let mut solution = vec![[0.0_f64; 4]; n * m];
    for col in (0..n).rev() {
        let p_row = perm[col];
        for j in 0..m {
            let mut acc = rhs[p_row * m + j];
            for c in (col + 1)..n {
                let sol = solution[c * m + j];
                for (av, sv) in acc.iter_mut().zip(sol.iter()) {
                    *av -= a[p_row * n + c] * sv;
                }
            }
            let diag = a[p_row * n + col];
            for comp in acc.iter_mut() {
                *comp /= diag;
            }
            solution[col * m + j] = acc;
        }
    }
    // Dehomogenize.
    let mut control_points = Vec::with_capacity(n * m);
    let mut weights = Vec::with_capacity(n * m);
    for h in &solution {
        let w = h[3];
        if !(w.is_finite() && w > 1e-9) {
            return None;
        }
        control_points.push(Point3::new(h[0] / w, h[1] / w, h[2] / w));
        weights.push(w);
    }
    Some(popkern_core::geometry::NurbsSurface {
        degree_u: 2,
        degree_v: if m == 5 { 2 } else { 2.min(m - 1) },
        knots_u: knots,
        knots_v: if m == 5 {
            vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0]
        } else if m == 3 {
            vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]
        } else {
            crate::fillet_surface::build_clamped_uniform_knots(m, 2.min(m - 1))
        },
        control_points,
        num_u: n,
        num_v: m,
        weights,
    })
}

/// Extract the degree-1 loft's exact station rows and re-loft them
/// through the interpolating degree-2 solve.
pub(crate) fn smooth_reloft_with_rails(
    loft: &popkern_core::geometry::NurbsSurface,
) -> Option<(
    popkern_core::geometry::NurbsSurface,
    popkern_core::geometry::NurbsCurve,
    popkern_core::geometry::NurbsCurve,
)> {
    let smooth = smooth_reloft(loft)?;
    let rail = |col: usize| -> popkern_core::geometry::NurbsCurve {
        let n = smooth.num_u;
        let m = smooth.num_v;
        popkern_core::geometry::NurbsCurve {
            degree: 2,
            knots: smooth.knots_u.clone(),
            control_points: (0..n).map(|i| smooth.control_points[i * m + col]).collect(),
            weights: (0..n).map(|i| smooth.weights[i * m + col]).collect(),
        }
    };
    let t1 = rail(0);
    let t2 = rail(smooth.num_v - 1);
    Some((smooth, t1, t2))
}

fn smooth_reloft(
    loft: &popkern_core::geometry::NurbsSurface,
) -> Option<popkern_core::geometry::NurbsSurface> {
    let n = loft.num_u;
    let m = loft.num_v;
    if n < 3 || loft.control_points.len() != n * m {
        return None;
    }
    let rows: Vec<Vec<Point3>> = (0..n)
        .map(|i| loft.control_points[i * m..(i + 1) * m].to_vec())
        .collect();
    let row_weights: Vec<Vec<f64>> = (0..n)
        .map(|i| loft.weights[i * m..(i + 1) * m].to_vec())
        .collect();
    crate::station_geometry::interpolating_loft(&rows, &row_weights)
}
