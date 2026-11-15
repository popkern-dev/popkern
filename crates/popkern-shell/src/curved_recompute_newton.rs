//! Newton refinement of a vertex onto the intersection of surfaces.
//!
//! Sequential projection onto each surface in turn converges linearly,
//! by about cos² of the angle between the surfaces each round: where an
//! opening's plane cuts the offset sphere near its top the two meet at
//! a shallow angle and 200 rounds stop well short (a ball cut at 3.9 of
//! radius 5 and shelled 1 kept its inner rim at the outer rim's radius
//! 3.13, where the offset sphere meets the plane at 0.89; at 3.5 it
//! converged). This takes the projection's last iterate and steps on the
//! surfaces' distances with the minimal-norm Gauss–Newton correction
//! `p ← p − Jᵀ (J Jᵀ)⁻¹ f`, the rows of `J` the unit normals at each
//! surface's nearest point, quadratically convergent wherever those
//! normals are independent; dependent normals (a tangency, where the
//! intersection is not a transversal point or curve) decline.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};

use crate::curved_recompute::project_to_surface;

/// Newton rounds before giving up.
const ROUNDS: usize = 50;

/// `pos` moved onto every surface within `accept`, by Gauss–Newton on
/// the signed distances; `None` when the normals are dependent or the
/// iteration does not settle.
pub(crate) fn refine_onto_surfaces(
    surfaces: &[&SurfaceGeometry],
    mut pos: Point3,
    accept: f64,
) -> Option<Point3> {
    if surfaces.is_empty() || surfaces.len() > 3 {
        return None;
    }
    for _ in 0..ROUNDS {
        let mut rows: Vec<(Vec3, f64)> = Vec::with_capacity(surfaces.len());
        for s in surfaces {
            let q = project_to_surface(s, &pos)?;
            let n = unit_normal(s, &q)?;
            // Signed distance along the surface normal at the foot.
            rows.push((n, (pos - q).dot(&n)));
        }
        if rows.iter().all(|(_, f)| f.abs() <= 0.25 * accept) {
            return Some(pos);
        }
        let step = minimal_norm_step(&rows)?;
        pos = pos - step;
    }
    None
}

/// The unit normal of `s` at its point `q`.
fn unit_normal(s: &SurfaceGeometry, q: &Point3) -> Option<Vec3> {
    let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(s, q)?;
    s.normal_at(u, v).normalized()
}

/// `Jᵀ (J Jᵀ)⁻¹ f` for unit rows `J` and residuals `f`; `None` when the
/// Gram matrix is singular (dependent normals).
fn minimal_norm_step(rows: &[(Vec3, f64)]) -> Option<Vec3> {
    let k = rows.len();
    let mut gram = [[0.0_f64; 3]; 3];
    for i in 0..k {
        for j in 0..k {
            gram[i][j] = rows[i].0.dot(&rows[j].0);
        }
    }
    let rhs: Vec<f64> = rows.iter().map(|r| r.1).collect();
    let lambda = solve_small(&gram, &rhs, k)?;
    let mut step = Vec3::new(0.0, 0.0, 0.0);
    for (i, l) in lambda.iter().enumerate() {
        step = step + rows[i].0 * *l;
    }
    Some(step)
}

/// Solve the leading `k × k` block of `a` against `b` by Gaussian
/// elimination with partial pivoting; `None` when a pivot vanishes.
fn solve_small(a: &[[f64; 3]; 3], b: &[f64], k: usize) -> Option<Vec<f64>> {
    let mut m = [[0.0_f64; 4]; 3];
    for i in 0..k {
        m[i][..k].copy_from_slice(&a[i][..k]);
        m[i][3] = b[i];
    }
    for col in 0..k {
        let piv = (col..k).max_by(|&x, &y| m[x][col].abs().total_cmp(&m[y][col].abs()))?;
        // Normals within ~2e-3° of dependent (sin² under 1e-9): a
        // tangency, no transversal point. A pinhole opening meets its
        // offset sphere at 0.4° (sin² 5e-5) and solves cleanly.
        if m[piv][col].abs() < 1e-9 {
            return None;
        }
        m.swap(col, piv);
        let pivot_row = m[col];
        for (r, row) in m.iter_mut().enumerate().take(k) {
            if r != col {
                let f = row[col] / pivot_row[col];
                for (x, p) in row.iter_mut().zip(pivot_row.iter()).skip(col) {
                    *x -= f * p;
                }
            }
        }
    }
    Some((0..k).map(|i| m[i][3] / m[i][i]).collect())
}
