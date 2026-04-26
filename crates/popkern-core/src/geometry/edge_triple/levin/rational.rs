//! Exact-rational foundation for the symbolic Levin pipeline (Q[√d]).
//!
//! The f64 path in `dispatch.rs` is fast and good enough for typical
//! engineering input, but accumulates floating-point noise on the
//! pencil's characteristic polynomial — making it impossible to
//! *prove* a fibre is exactly rank 2 (vs rank 3 with a tiny non-zero
//! eigenvalue), or to detect a true tangent intersection symbolically.
//!
//! This module gives the rational-arithmetic primitives that let us
//! work out the pencil exactly when the input quadric matrices have
//! rational entries:
//!
//! - `Rat = BigRational` (arbitrary-precision rationals; no overflow).
//! - `Mat4SymRat` — 4×4 symmetric matrix of `Rat`s.
//! - `surface_to_quadric_rat` for `Plane | Sphere | Cylinder | Cone`
//!   inputs whose parameters round-trip exactly to rationals (most
//!   typical CAD inputs do).
//! - `det4_rat` (cofactor expansion), `linear_comb_rat` (`λA + μB`).
//! - `char_poly_rat` — coefficients of `det(λA + B)` as a degree-4
//!   polynomial in `λ`, computed by Lagrange interpolation at five
//!   rational `λ` values (0, ±1, ±2). All-integer arithmetic.
//! - `rational_roots` — rational-root theorem to enumerate exact
//!   rational λ at which the pencil is degenerate.
//! - `find_rational_degenerate_fiber` — exact rank-deficient fibre
//!   when one exists in `Q`.
//!
//! Q[√d] adjunction (the actual symbolic eigendecomposition for
//! rank-2 plane extraction) is layered on top in a follow-up module.
//! With this foundation, the disambiguation between rank-2 and rank-3
//! fibres is *provable* — when `det(λA + B) = 0` has a rational root
//! AND the resulting matrix has a rational kernel of dimension 2, we
//! know we're in the rank-2 case symbolically.

#![allow(dead_code)] // foundation module — symbolic path layered on top

use num_rational::BigRational;
use num_traits::{One, Signed, Zero};

use crate::geometry::Vec3;
use crate::geometry::surface::SurfaceGeometry;

/// Arbitrary-precision rational. `BigRational` over `BigInt` — no
/// overflow even for pathologically scaled quadric coefficients.
pub(crate) type Rat = BigRational;

/// 4×4 symmetric matrix of exact rationals. Stored densely; symmetry
/// is maintained by every constructor in this module.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Mat4SymRat(pub(crate) [[Rat; 4]; 4]);

impl Mat4SymRat {
    pub(crate) fn zero() -> Self {
        let z = Rat::zero();
        Mat4SymRat([
            [z.clone(), z.clone(), z.clone(), z.clone()],
            [z.clone(), z.clone(), z.clone(), z.clone()],
            [z.clone(), z.clone(), z.clone(), z.clone()],
            [z.clone(), z.clone(), z.clone(), z],
        ])
    }
}

/// Convert an `f64` to an exact rational using `BigRational::from_float`
/// (which represents the float bit-exactly as `m / 2^k`).
fn rat_from_f64(x: f64) -> Option<Rat> {
    BigRational::from_float(x)
}

fn rat_zero() -> Rat {
    Rat::zero()
}

fn rat_one() -> Rat {
    Rat::one()
}

/// Convert an analytic surface to its exact rational quadric form.
/// Returns `None` for surfaces with non-rational parameters (NaN,
/// infinite, or the `Cone` half-angle's `cos²θ` doesn't round-trip
/// exactly to rational — handled below).
pub(crate) fn surface_to_quadric_rat(surf: &SurfaceGeometry) -> Option<Mat4SymRat> {
    match surf {
        SurfaceGeometry::Plane { origin, normal, .. } => Some(plane_quadric_rat(*origin, *normal)?),
        SurfaceGeometry::Sphere { center, radius } => Some(sphere_quadric_rat(*center, *radius)?),
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => Some(cylinder_quadric_rat(*origin, *axis, *radius)?),
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => Some(cone_quadric_rat(*apex, *axis, *half_angle)?),
        _ => None,
    }
}

fn plane_quadric_rat(origin: crate::geometry::Point3, normal: Vec3) -> Option<Mat4SymRat> {
    let nx = rat_from_f64(normal.x)?;
    let ny = rat_from_f64(normal.y)?;
    let nz = rat_from_f64(normal.z)?;
    let ox = rat_from_f64(origin.x)?;
    let oy = rat_from_f64(origin.y)?;
    let oz = rat_from_f64(origin.z)?;
    let np = &nx * &ox + &ny * &oy + &nz * &oz;
    let v = [nx, ny, nz, -np];
    let mut a = mat4_zero();
    for (i, row) in a.iter_mut().enumerate() {
        for (j, slot) in row.iter_mut().enumerate() {
            *slot = &v[i] * &v[j];
        }
    }
    Some(Mat4SymRat(a))
}

fn sphere_quadric_rat(center: crate::geometry::Point3, radius: f64) -> Option<Mat4SymRat> {
    let cx = rat_from_f64(center.x)?;
    let cy = rat_from_f64(center.y)?;
    let cz = rat_from_f64(center.z)?;
    let r = rat_from_f64(radius)?;
    let cdot = &cx * &cx + &cy * &cy + &cz * &cz;
    let r_sq = &r * &r;
    let mut a = mat4_zero();
    a[0][0] = rat_one();
    a[1][1] = rat_one();
    a[2][2] = rat_one();
    a[0][3] = -&cx;
    a[3][0] = -cx;
    a[1][3] = -&cy;
    a[3][1] = -cy;
    a[2][3] = -&cz;
    a[3][2] = -cz;
    a[3][3] = cdot - r_sq;
    Some(Mat4SymRat(a))
}

fn cylinder_quadric_rat(
    origin: crate::geometry::Point3,
    axis: Vec3,
    radius: f64,
) -> Option<Mat4SymRat> {
    let ax = rat_from_f64(axis.x)?;
    let ay = rat_from_f64(axis.y)?;
    let az = rat_from_f64(axis.z)?;
    let ox = rat_from_f64(origin.x)?;
    let oy = rat_from_f64(origin.y)?;
    let oz = rat_from_f64(origin.z)?;
    let r = rat_from_f64(radius)?;
    // M = I − a·aᵀ
    let m00 = rat_one() - &ax * &ax;
    let m01 = -(&ax * &ay);
    let m02 = -(&ax * &az);
    let m11 = rat_one() - &ay * &ay;
    let m12 = -(&ay * &az);
    let m22 = rat_one() - &az * &az;
    let m = [
        [m00.clone(), m01.clone(), m02.clone()],
        [m01.clone(), m11.clone(), m12.clone()],
        [m02.clone(), m12.clone(), m22.clone()],
    ];
    let o = [ox.clone(), oy.clone(), oz.clone()];
    let mo = mat3_vec3_rat(&m, &o);
    let omo = &o[0] * &mo[0] + &o[1] * &mo[1] + &o[2] * &mo[2];
    let r_sq = &r * &r;
    let mut a = mat4_zero();
    a[0][0] = m00;
    a[0][1] = m01.clone();
    a[1][0] = m01;
    a[0][2] = m02.clone();
    a[2][0] = m02;
    a[1][1] = m11;
    a[1][2] = m12.clone();
    a[2][1] = m12;
    a[2][2] = m22;
    a[0][3] = -mo[0].clone();
    a[3][0] = -mo[0].clone();
    a[1][3] = -mo[1].clone();
    a[3][1] = -mo[1].clone();
    a[2][3] = -mo[2].clone();
    a[3][2] = -mo[2].clone();
    a[3][3] = omo - r_sq;
    Some(Mat4SymRat(a))
}

fn cone_quadric_rat(
    apex: crate::geometry::Point3,
    axis: Vec3,
    half_angle: f64,
) -> Option<Mat4SymRat> {
    // cos²θ comes from f64 trig. The result generally won't be a
    // *short* rational (it's an irrational like cos²(π/4) = 0.5,
    // which IS rational, but cos²(0.3) is not). We accept whatever
    // float `cos²θ` gives and `from_float` represents it bit-exactly
    // as m / 2^k.
    let cos2 = rat_from_f64(half_angle.cos().powi(2))?;
    let ax = rat_from_f64(axis.x)?;
    let ay = rat_from_f64(axis.y)?;
    let az = rat_from_f64(axis.z)?;
    let vx = rat_from_f64(apex.x)?;
    let vy = rat_from_f64(apex.y)?;
    let vz = rat_from_f64(apex.z)?;
    // Q3 = cos²θ · I − a·aᵀ
    let q00 = &cos2 - &ax * &ax;
    let q01 = -(&ax * &ay);
    let q02 = -(&ax * &az);
    let q11 = &cos2 - &ay * &ay;
    let q12 = -(&ay * &az);
    let q22 = &cos2 - &az * &az;
    let q = [
        [q00.clone(), q01.clone(), q02.clone()],
        [q01.clone(), q11.clone(), q12.clone()],
        [q02.clone(), q12.clone(), q22.clone()],
    ];
    let v = [vx, vy, vz];
    let qv = mat3_vec3_rat(&q, &v);
    let vqv = &v[0] * &qv[0] + &v[1] * &qv[1] + &v[2] * &qv[2];
    let mut a = mat4_zero();
    a[0][0] = q00;
    a[0][1] = q01.clone();
    a[1][0] = q01;
    a[0][2] = q02.clone();
    a[2][0] = q02;
    a[1][1] = q11;
    a[1][2] = q12.clone();
    a[2][1] = q12;
    a[2][2] = q22;
    a[0][3] = -qv[0].clone();
    a[3][0] = -qv[0].clone();
    a[1][3] = -qv[1].clone();
    a[3][1] = -qv[1].clone();
    a[2][3] = -qv[2].clone();
    a[3][2] = -qv[2].clone();
    a[3][3] = vqv;
    Some(Mat4SymRat(a))
}

fn mat3_vec3_rat(m: &[[Rat; 3]; 3], v: &[Rat; 3]) -> [Rat; 3] {
    [
        &m[0][0] * &v[0] + &m[0][1] * &v[1] + &m[0][2] * &v[2],
        &m[1][0] * &v[0] + &m[1][1] * &v[1] + &m[1][2] * &v[2],
        &m[2][0] * &v[0] + &m[2][1] * &v[1] + &m[2][2] * &v[2],
    ]
}

fn mat4_zero() -> [[Rat; 4]; 4] {
    let z = Rat::zero();
    [
        [z.clone(), z.clone(), z.clone(), z.clone()],
        [z.clone(), z.clone(), z.clone(), z.clone()],
        [z.clone(), z.clone(), z.clone(), z.clone()],
        [z.clone(), z.clone(), z.clone(), z],
    ]
}

/// `λ·A + μ·B` (exact).
pub(crate) fn linear_comb_rat(
    a: &Mat4SymRat,
    b: &Mat4SymRat,
    lambda: &Rat,
    mu: &Rat,
) -> Mat4SymRat {
    let mut out = mat4_zero();
    for (i, row) in out.iter_mut().enumerate() {
        for (j, slot) in row.iter_mut().enumerate() {
            *slot = lambda * &a.0[i][j] + mu * &b.0[i][j];
        }
    }
    Mat4SymRat(out)
}

/// Determinant of a 4×4 rational matrix via row-0 cofactor expansion.
pub(crate) fn det4_rat(m: &[[Rat; 4]; 4]) -> Rat {
    let mut det = Rat::zero();
    for j in 0..4 {
        let cofactor = minor3x3_det_rat(m, 0, j);
        if j % 2 == 0 {
            det += &m[0][j] * &cofactor;
        } else {
            det -= &m[0][j] * &cofactor;
        }
    }
    det
}

fn minor3x3_det_rat(m: &[[Rat; 4]; 4], skip_row: usize, skip_col: usize) -> Rat {
    let z = Rat::zero();
    let mut sub = [
        [z.clone(), z.clone(), z.clone()],
        [z.clone(), z.clone(), z.clone()],
        [z.clone(), z.clone(), z],
    ];
    let mut ri = 0;
    for (i, row) in m.iter().enumerate() {
        if i == skip_row {
            continue;
        }
        let mut cj = 0;
        for (j, val) in row.iter().enumerate() {
            if j == skip_col {
                continue;
            }
            sub[ri][cj] = val.clone();
            cj += 1;
        }
        ri += 1;
    }
    &sub[0][0] * (&sub[1][1] * &sub[2][2] - &sub[1][2] * &sub[2][1])
        - &sub[0][1] * (&sub[1][0] * &sub[2][2] - &sub[1][2] * &sub[2][0])
        + &sub[0][2] * (&sub[1][0] * &sub[2][1] - &sub[1][1] * &sub[2][0])
}

/// Coefficients of `det(λA + B)` as a degree-4 univariate polynomial
/// in `λ`, returned in monomial-descending order:
/// `c[0]·λ⁴ + c[1]·λ³ + c[2]·λ² + c[3]·λ + c[4]`.
///
/// Computed exactly via Lagrange interpolation at λ ∈ {0, 1, -1, 2, -2}.
pub(crate) fn char_poly_rat(a: &Mat4SymRat, b: &Mat4SymRat) -> [Rat; 5] {
    let zero = Rat::zero();
    let one = Rat::one();
    let two = &one + &one;
    let neg_one = -one.clone();
    let neg_two = -two.clone();
    let lambdas = [
        zero.clone(),
        one.clone(),
        neg_one.clone(),
        two.clone(),
        neg_two.clone(),
    ];
    let s: [Rat; 5] = [
        det4_rat(&linear_comb_rat(a, b, &lambdas[0], &one).0),
        det4_rat(&linear_comb_rat(a, b, &lambdas[1], &one).0),
        det4_rat(&linear_comb_rat(a, b, &lambdas[2], &one).0),
        det4_rat(&linear_comb_rat(a, b, &lambdas[3], &one).0),
        det4_rat(&linear_comb_rat(a, b, &lambdas[4], &one).0),
    ];
    // Closed-form inverse Vandermonde for {0, 1, -1, 2, -2} (same as
    // the f64 path in `pencil::char_poly_in_lambda`):
    //   c4 = s[0]
    //   2c0 + 2c2 = s[1] + s[2] - 2·c4
    //   32c0 + 8c2 = s[3] + s[4] - 2·c4
    //   2c1 + 2c3 = s[1] - s[2]
    //   16c1 + 4c3 = s[3] - s[4]
    let c4 = s[0].clone();
    let aa = &s[1] + &s[2] - &two * &c4;
    let bb = &s[3] + &s[4] - &two * &c4;
    // c0 = (bb − 4·aa) / 24
    let four = &two + &two;
    let twenty_four = &four * Rat::from(BigRational::from_integer(num_bigint::BigInt::from(6)));
    let c0 = (&bb - &four * &aa) / twenty_four;
    // c2 = aa/2 − c0
    let c2 = &aa / &two - &c0;
    let cc = &s[1] - &s[2];
    let dd = &s[3] - &s[4];
    // c1 = (dd − 2·cc) / 12
    let twelve = Rat::from(BigRational::from_integer(num_bigint::BigInt::from(12)));
    let c1 = (&dd - &two * &cc) / twelve;
    // c3 = cc/2 − c1
    let c3 = cc / &two - &c1;
    [c0, c1, c2, c3, c4]
}

/// Find all rational roots of a polynomial with rational coefficients
/// `c[0]·x^n + ... + c[n]` via the rational-root theorem.
///
/// Strategy: clear denominators to get an integer-coefficient
/// polynomial, then enumerate all `±p/q` with `p` dividing the
/// constant term and `q` dividing the leading coefficient.
/// Returns roots in arbitrary order.
pub(crate) fn rational_roots(coefs: &[Rat]) -> Vec<Rat> {
    if coefs.is_empty() {
        return Vec::new();
    }
    // Find common denominator.
    let mut common_denom = num_bigint::BigInt::from(1);
    for c in coefs {
        common_denom = lcm_bigint(&common_denom, c.denom());
    }
    // Integer coefficients (still in monomial-descending order).
    let int_coefs: Vec<num_bigint::BigInt> = coefs
        .iter()
        .map(|c| (c * Rat::from(BigRational::from_integer(common_denom.clone()))).to_integer())
        .collect();
    // Skip leading zeros so the leading coef is non-zero.
    let mut start = 0;
    while start < int_coefs.len() - 1 && int_coefs[start].is_zero() {
        start += 1;
    }
    let coefs = &int_coefs[start..];
    if coefs.is_empty() {
        return Vec::new();
    }
    let leading = &coefs[0];
    let constant = &coefs[coefs.len() - 1];
    // Special case: constant is zero ⇒ x=0 is a root. Factor it out
    // and recurse on the deflated polynomial: c[0]·x^n + … + c[n−1]·x
    //   = x · (c[0]·x^(n−1) + … + c[n−1])
    // so the deflated coefs are `coefs[..coefs.len() − 1]`.
    if constant.is_zero() {
        let deflated: Vec<Rat> = coefs[..coefs.len() - 1]
            .iter()
            .map(|i| Rat::from(BigRational::from_integer(i.clone())))
            .collect();
        let mut roots = vec![Rat::zero()];
        for r in rational_roots(&deflated) {
            if !roots.iter().any(|x| x == &r) {
                roots.push(r);
            }
        }
        return roots;
    }
    let p_divisors = divisors_bigint(constant);
    let q_divisors = divisors_bigint(leading);
    let mut roots = Vec::new();
    for p in &p_divisors {
        for q in &q_divisors {
            for sign in [1, -1] {
                let p_signed = if sign == -1 { -p } else { p.clone() };
                let candidate = Rat::new(p_signed.clone(), q.clone());
                if eval_int_poly(coefs, &candidate).is_zero()
                    && !roots.iter().any(|x| x == &candidate)
                {
                    roots.push(candidate);
                }
            }
        }
    }
    roots
}

fn eval_int_poly(coefs: &[num_bigint::BigInt], x: &Rat) -> Rat {
    let mut acc = Rat::zero();
    for c in coefs {
        acc = acc * x + Rat::from(BigRational::from_integer(c.clone()));
    }
    acc
}

fn lcm_bigint(a: &num_bigint::BigInt, b: &num_bigint::BigInt) -> num_bigint::BigInt {
    use num_integer::Integer;
    if a.is_zero() || b.is_zero() {
        return num_bigint::BigInt::from(1);
    }
    a.lcm(b)
}

fn divisors_bigint(n: &num_bigint::BigInt) -> Vec<num_bigint::BigInt> {
    let abs_n = n.abs();
    if abs_n.is_zero() {
        return vec![num_bigint::BigInt::from(1)];
    }
    let one = num_bigint::BigInt::from(1);
    let mut out = vec![one.clone()];
    let mut d = num_bigint::BigInt::from(2);
    while &d * &d <= abs_n {
        if (&abs_n % &d).is_zero() {
            out.push(d.clone());
            let q = &abs_n / &d;
            if q != d {
                out.push(q);
            }
        }
        d = &d + &one;
    }
    if abs_n != one {
        out.push(abs_n);
    }
    out
}

/// True iff the pencil `λA + μB` is *identically* rank-deficient —
/// i.e., `det(λA + μB) ≡ 0` as a polynomial in `λ`. This happens
/// exactly when `A` and `B` are scalar multiples of each other (so
/// the surfaces they encode are *coincident*, modulo scaling), or
/// when both `A` and `B` are themselves rank-deficient and their
/// pencil shares a common kernel.
///
/// CAD-relevant case: two duplicated surfaces in a BRep model. The
/// f64 dispatch chokes on this (every λ is "almost" a degenerate
/// root, with FP noise making it impossible to pick one). The
/// symbolic test is exact: every coefficient of the characteristic
/// polynomial is rational zero.
pub(crate) fn pencil_is_identically_degenerate(a: &Mat4SymRat, b: &Mat4SymRat) -> bool {
    let coefs = char_poly_rat(a, b);
    coefs.iter().all(|c| c.is_zero())
}

/// True iff the two quadric matrices are *scalar multiples* of each
/// other — i.e., `A = c · B` for some non-zero rational `c`. This is
/// the strict "coincident surface" relation: same algebraic surface,
/// possibly with different homogeneous-coordinate scaling. Note
/// that `pencil_is_identically_degenerate` is *weaker* — coaxial
/// cylinders with different radii have a degenerate pencil (shared
/// kernel along the axis) but are NOT coincident.
pub(crate) fn matrices_are_proportional(a: &Mat4SymRat, b: &Mat4SymRat) -> bool {
    // Find the first non-zero entry of B; use it to derive the ratio
    // c = a_ij / b_ij. Then verify A = c·B everywhere.
    let mut ratio: Option<Rat> = None;
    for i in 0..4 {
        for j in 0..4 {
            if !b.0[i][j].is_zero() {
                ratio = Some(&a.0[i][j] / &b.0[i][j]);
                break;
            }
        }
        if ratio.is_some() {
            break;
        }
    }
    let c = match ratio {
        Some(c) => c,
        // B is the zero matrix; A coincides with B iff A is also zero.
        None => return a.0.iter().flat_map(|r| r.iter()).all(|x| x.is_zero()),
    };
    for i in 0..4 {
        for j in 0..4 {
            if a.0[i][j] != &c * &b.0[i][j] {
                return false;
            }
        }
    }
    true
}

/// True iff `surf_a` and `surf_b` describe the same algebraic
/// surface (modulo scaling) — exactly. This is the user-facing
/// CAD "duplicated surface" query. Symbolic answer for any pair
/// of analytic surfaces with rational parameters.
pub(crate) fn surfaces_are_coincident(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
) -> Option<bool> {
    let a = surface_to_quadric_rat(surf_a)?;
    let b = surface_to_quadric_rat(surf_b)?;
    Some(matrices_are_proportional(&a, &b))
}

/// Search the pencil for an exact rational degenerate fibre. Returns
/// the (λ, μ=1) parameter and the resulting rank-deficient matrix.
/// Returns `None` if the characteristic polynomial has no rational
/// roots (the pencil's degenerate fibres all have irrational λ — Q[√d]
/// adjunction would catch them).
pub(crate) fn find_rational_degenerate_fiber(
    a: &Mat4SymRat,
    b: &Mat4SymRat,
) -> Option<(Rat, Mat4SymRat)> {
    let coefs = char_poly_rat(a, b);
    let roots = rational_roots(&coefs);
    let one = Rat::one();
    for lam in roots {
        let m = linear_comb_rat(a, b, &lam, &one);
        if det4_rat(&m.0).is_zero() {
            return Some((lam, m));
        }
    }
    // Fallback: A or B itself rank-deficient.
    if det4_rat(&a.0).is_zero() {
        return Some((Rat::one(), a.clone()));
    }
    if det4_rat(&b.0).is_zero() {
        return Some((Rat::zero(), b.clone()));
    }
    None
}
