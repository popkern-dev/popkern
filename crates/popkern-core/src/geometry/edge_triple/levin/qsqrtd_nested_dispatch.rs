use num_bigint::BigInt;
use num_traits::{One, ToPrimitive, Zero};

use super::qsqrtd::QSqrtD;
use super::qsqrtd_nested::QSqrtDSqrtL;
use super::rational::{Mat4SymRat, Rat, det4_rat, linear_comb_rat};

/// True iff every matrix entry has numerator and denominator below
/// `max_bits`; used to keep the symbolic path bounded on inputs
/// round-tripped from f64 irrationals.
pub(crate) fn matrix_entries_within_bit_budget(m: &Mat4SymRat, max_bits: u64) -> bool {
    for row in &m.0 {
        for entry in row {
            if entry.numer().bits() > max_bits || entry.denom().bits() > max_bits {
                return false;
            }
        }
    }
    true
}

/// Walk a bounded set of small rational `lambda` candidates and
/// return the first rank-deficient fiber `lambda A + B`.
pub(crate) fn find_small_rational_fiber(a: &Mat4SymRat, b: &Mat4SymRat) -> Option<Mat4SymRat> {
    let one = Rat::one();
    let candidates: [(BigInt, BigInt); 19] = [
        (BigInt::from(-3), BigInt::from(1)),
        (BigInt::from(-2), BigInt::from(1)),
        (BigInt::from(-1), BigInt::from(1)),
        (BigInt::from(0), BigInt::from(1)),
        (BigInt::from(1), BigInt::from(1)),
        (BigInt::from(2), BigInt::from(1)),
        (BigInt::from(3), BigInt::from(1)),
        (BigInt::from(-1), BigInt::from(2)),
        (BigInt::from(1), BigInt::from(2)),
        (BigInt::from(-3), BigInt::from(2)),
        (BigInt::from(3), BigInt::from(2)),
        (BigInt::from(-1), BigInt::from(3)),
        (BigInt::from(1), BigInt::from(3)),
        (BigInt::from(-2), BigInt::from(3)),
        (BigInt::from(2), BigInt::from(3)),
        (BigInt::from(-1), BigInt::from(4)),
        (BigInt::from(1), BigInt::from(4)),
        (BigInt::from(-3), BigInt::from(4)),
        (BigInt::from(3), BigInt::from(4)),
    ];
    for (p, q) in candidates {
        let lam = Rat::new(p, q);
        let m = linear_comb_rat(a, b, &lam, &one);
        if det4_rat(&m.0).is_zero() {
            return Some(m);
        }
    }
    None
}

/// Convert a `QSqrtD` value `a + b*sqrt(d)` to its real numeric
/// value.
pub(crate) fn qsqrtd_to_f64(v: &QSqrtD) -> Option<f64> {
    let a_f = v.a.to_f64()?;
    if v.b.is_zero() {
        return Some(a_f);
    }
    let b_f = v.b.to_f64()?;
    let d_f = v.d.to_f64()?;
    if d_f < 0.0 {
        return None;
    }
    Some(a_f + b_f * d_f.sqrt())
}

/// Convert a nested-radical value `p0 + p1*sqrt(alpha)` to f64.
pub(crate) fn qsqrtd_sqrt_l_to_f64(v: &QSqrtDSqrtL) -> Option<f64> {
    let p0_f = qsqrtd_to_f64(&v.p0)?;
    if v.p1.is_zero() {
        return Some(p0_f);
    }
    let p1_f = qsqrtd_to_f64(&v.p1)?;
    let alpha_f = qsqrtd_to_f64(&v.alpha)?;
    if alpha_f < 0.0 {
        return None;
    }
    Some(p0_f + p1_f * alpha_f.sqrt())
}

/// Convert a 4-vector of nested-radical plane coefficients to f64.
pub(crate) fn qsqrtd_sqrt_l_array_to_f64(v: &[QSqrtDSqrtL; 4]) -> Option<[f64; 4]> {
    Some([
        qsqrtd_sqrt_l_to_f64(&v[0])?,
        qsqrtd_sqrt_l_to_f64(&v[1])?,
        qsqrtd_sqrt_l_to_f64(&v[2])?,
        qsqrtd_sqrt_l_to_f64(&v[3])?,
    ])
}
