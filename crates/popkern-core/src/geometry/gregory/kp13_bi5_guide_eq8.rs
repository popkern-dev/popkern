use super::super::bernstein_tri5::{TRI5_CP_COUNT, flat_index};
use super::super::direct_linalg::solve_dense;
use super::super::kp13_bi5_chi::ChiJet;
use super::super::kp13_bi5_jet::{TaylorJet2, compose_chain_rule};
use super::Bi5GuideSector;
use crate::geometry::Point3;

// ---- tri-5 → bi-5 degree-raise ----

/// Closed-form tri-5 → bi-5 degree raise.
///
/// Converts 21 tri-5 BB coefs (one sector of `g△^k`) into the 6 × 6
/// bi-5 BB coefs of the same polynomial expressed on `[0, 1]²`.
///
/// The tri-5 polynomial
/// `P(u, v) = Σ_{i+j ≤ 5} C_{ij} · B^5_{5−i−j, i, j}(u, v)`
/// (treating `(u, v)` as the natural-simplex coords `(u_s, v_s)`
/// via the linear shear `L`) equals the bi-5 polynomial
/// `Q(u, v) = Σ_{p, q ∈ 0..6} D_{pq} · B^5_p(u) · B^5_q(v)`
/// where each `D_{pq}` is a specific linear combination of the 21
/// `C_{ij}`. The coefficients are derived from the multinomial
/// expansion of `(1−u−v)^{5−i−j}` followed by the closed-form
/// monomial-to-Bernstein conversion.
pub fn tri5_to_bi5_degree_raise(cps_tri: &[Point3; TRI5_CP_COUNT]) -> [[Point3; 6]; 6] {
    let zero = Point3::new(0.0, 0.0, 0.0);

    // Step 1: tri-5 BB → bivariate monomials.
    // monomials[r][s] (for r + s ≤ 5) is the coefficient of u^r v^s.
    let mut monomials = [[zero; 6]; 6];
    for i in 0..=5 {
        for j in 0..=(5 - i) {
            let cp = cps_tri[flat_index(5, i, j)];
            // B^5_{a,b,c} with (a, b, c) = (5 − i − j, i, j):
            //   K · (1 − u − v)^a · u^b · v^c
            // = K · Σ_{r ≥ b, s ≥ c, r+s ≤ 5}
            //       [a! / ((5 − r − s)! (r − b)! (s − c)!) · (−1)^{r+s-b-c}]
            //       · u^r · v^s
            // with K = 5! / (a! b! c!), so K · a! = 5! / (b! c!).
            for r in i..=5 {
                for s in j..=(5 - r) {
                    let r_b = r - i;
                    let s_c = s - j;
                    let i_exp = 5 - r - s;
                    let sign = if (r_b + s_c) % 2 == 0 { 1.0 } else { -1.0 };
                    let numer = factorial(5) as f64;
                    let denom = (factorial(i)
                        * factorial(j)
                        * factorial(i_exp)
                        * factorial(r_b)
                        * factorial(s_c)) as f64;
                    let coef = sign * numer / denom;
                    monomials[r][s].x += coef * cp.x;
                    monomials[r][s].y += coef * cp.y;
                    monomials[r][s].z += coef * cp.z;
                }
            }
        }
    }

    // Step 2: bivariate monomials → bi-5 BB.
    // u^r = Σ_{p=r}^{5} [C(5 − r, p − r) / C(5, p)] · B^5_p(u).
    // Bivariate: u^r v^s = Σ_{p ≥ r, q ≥ s} w(r, p) · w(s, q) · B^5_p(u) B^5_q(v).
    let mut bi5 = [[zero; 6]; 6];
    for p in 0..6 {
        for q in 0..6 {
            let denom = (binomial(5, p) * binomial(5, q)) as f64;
            for r in 0..=p {
                for s in 0..=q {
                    if r + s > 5 {
                        continue;
                    }
                    let numer = (binomial(5 - r, p - r) * binomial(5 - s, q - s)) as f64;
                    let w = numer / denom;
                    bi5[p][q].x += w * monomials[r][s].x;
                    bi5[p][q].y += w * monomials[r][s].y;
                    bi5[p][q].z += w * monomials[r][s].z;
                }
            }
        }
    }
    bi5
}

/// Error from the CC-χ eq (8) solve.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Eq8CcError {
    /// The 6 × 6 linear system at sector `k` is singular to machine
    /// precision. Not expected for any `n ≥ 3` via
    /// [`CatmullClarkChi`] but surfaced defensively.
    SingularSystem(usize),
}

impl std::fmt::Display for Eq8CcError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Eq8CcError::SingularSystem(k) => {
                write!(f, "eq (8) CC 6×6 system singular at sector {k}")
            }
        }
    }
}

impl std::error::Error for Eq8CcError {}

/// Flat list of the six bi-5 CPs at corner `(u, v) = (1, 1)` that
/// determine the 2nd-order Taylor jet at that corner. All other 30
/// CPs stay pinned at their degree-raised values during eq (8).
const EQ8_CORNER_CPS: [(usize, usize); 6] = [(5, 5), (4, 5), (5, 4), (3, 5), (5, 3), (4, 4)];

/// Apply the CC-χ eq (8) override to one sector's bi-5 CPs.
///
/// Solves a 6 × 6 linear system per axis: each row is a jet-
/// component constraint and each column is a unit perturbation of
/// one `EQ8_CORNER_CPS` entry. The same matrix applies to all three
/// axes (linearity); only the RHS per axis differs.
///
/// On success writes the six CPs; on singular system leaves them as
/// their degree-raise values (caller should discard the partially-
/// mutated sector).
pub(super) fn apply_eq8_cc_override(
    sector: &mut Bi5GuideSector,
    target_jet: &TaylorJet2,
    chi_jet_simplex: &ChiJet,
) -> Result<(), ()> {
    let (u, v) = chi_jet_simplex.pos;

    // Save current values and zero the 6 corner CPs to isolate the
    // "base" contribution (composite jet without the unknowns).
    let mut saved: [Point3; 6] = [Point3::new(0.0, 0.0, 0.0); 6];
    for (i, &(p, q)) in EQ8_CORNER_CPS.iter().enumerate() {
        saved[i] = sector.cps[p][q];
        sector.cps[p][q] = Point3::new(0.0, 0.0, 0.0);
    }

    let base_outer = sector.outer_hessian(u, v);
    let base_jet = compose_chain_rule(&base_outer, chi_jet_simplex);

    // Build 6 × 6 matrix A where A[i][j] is the i-th jet component
    // contribution (in axis 0) when CP j is set to (1, 0, 0).
    let mut a_mat = vec![vec![0.0_f64; 6]; 6];
    for (j, &(p, q)) in EQ8_CORNER_CPS.iter().enumerate() {
        sector.cps[p][q] = Point3::new(1.0, 0.0, 0.0);
        let outer = sector.outer_hessian(u, v);
        let jet = compose_chain_rule(&outer, chi_jet_simplex);
        a_mat[0][j] = jet.value.x - base_jet.value.x;
        a_mat[1][j] = jet.d_u.x - base_jet.d_u.x;
        a_mat[2][j] = jet.d_v.x - base_jet.d_v.x;
        a_mat[3][j] = jet.d_uu.x - base_jet.d_uu.x;
        a_mat[4][j] = jet.d_uv.x - base_jet.d_uv.x;
        a_mat[5][j] = jet.d_vv.x - base_jet.d_vv.x;
        sector.cps[p][q] = Point3::new(0.0, 0.0, 0.0);
    }

    // Solve per axis.
    let target_components = [
        target_jet.value,
        target_jet.d_u,
        target_jet.d_v,
        target_jet.d_uu,
        target_jet.d_uv,
        target_jet.d_vv,
    ];
    let base_components = [
        base_jet.value,
        base_jet.d_u,
        base_jet.d_v,
        base_jet.d_uu,
        base_jet.d_uv,
        base_jet.d_vv,
    ];

    let mut solutions = [[0.0_f64; 6]; 3];
    for axis in 0..3 {
        let mut rhs = vec![0.0_f64; 6];
        for i in 0..6 {
            let t = axis_component(target_components[i], axis);
            let b = axis_component(base_components[i], axis);
            rhs[i] = t - b;
        }
        let x = solve_dense(a_mat.clone(), rhs).ok_or(())?;
        solutions[axis].copy_from_slice(&x[..6]);
    }

    // Write solutions into the 6 corner CPs.
    for (i, &(p, q)) in EQ8_CORNER_CPS.iter().enumerate() {
        sector.cps[p][q] = Point3::new(solutions[0][i], solutions[1][i], solutions[2][i]);
    }

    // `saved` discarded on success (we overwrote the 6 CPs with
    // new solutions).
    let _ = saved;
    Ok(())
}

#[inline]
fn axis_component(p: Point3, axis: usize) -> f64 {
    match axis {
        0 => p.x,
        1 => p.y,
        2 => p.z,
        _ => unreachable!("axis {axis} out of range"),
    }
}

/// Apply the bilinear-χ eq (8) override to one sector's bi-5 CPs.
///
/// Writes `P_{5,5}, P_{4,5}, P_{5,4}, P_{3,5}, P_{5,3}, P_{4,4}` so
/// that the 2nd-order bi-5 jet at corner `(1, 1)` equals `jet`.
/// The other 30 CPs (including `P_{3,3}, P_{3,4}, P_{4,3}`) are
/// untouched.
pub(super) fn apply_eq8_bilinear_override(sector: &mut Bi5GuideSector, jet: &TaylorJet2) {
    let value = jet.value;
    let d_u = jet.d_u;
    let d_v = jet.d_v;
    let d_uu = jet.d_uu;
    let d_uv = jet.d_uv;
    let d_vv = jet.d_vv;

    let p55 = value;
    let p45 = Point3::new(
        p55.x - d_u.x / 5.0,
        p55.y - d_u.y / 5.0,
        p55.z - d_u.z / 5.0,
    );
    let p54 = Point3::new(
        p55.x - d_v.x / 5.0,
        p55.y - d_v.y / 5.0,
        p55.z - d_v.z / 5.0,
    );
    let p35 = Point3::new(
        2.0 * p45.x - p55.x + d_uu.x / 20.0,
        2.0 * p45.y - p55.y + d_uu.y / 20.0,
        2.0 * p45.z - p55.z + d_uu.z / 20.0,
    );
    let p53 = Point3::new(
        2.0 * p54.x - p55.x + d_vv.x / 20.0,
        2.0 * p54.y - p55.y + d_vv.y / 20.0,
        2.0 * p54.z - p55.z + d_vv.z / 20.0,
    );
    let p44 = Point3::new(
        p45.x + p54.x - p55.x + d_uv.x / 25.0,
        p45.y + p54.y - p55.y + d_uv.y / 25.0,
        p45.z + p54.z - p55.z + d_uv.z / 25.0,
    );

    sector.cps[5][5] = p55;
    sector.cps[4][5] = p45;
    sector.cps[5][4] = p54;
    sector.cps[3][5] = p35;
    sector.cps[5][3] = p53;
    sector.cps[4][4] = p44;
}

// ---- small-integer helpers ----

#[inline]
const fn factorial(n: usize) -> u64 {
    let mut acc: u64 = 1;
    let mut k = 1;
    while k <= n {
        acc *= k as u64;
        k += 1;
    }
    acc
}

#[inline]
const fn binomial(n: usize, k: usize) -> u64 {
    if k > n {
        return 0;
    }
    factorial(n) / (factorial(k) * factorial(n - k))
}
