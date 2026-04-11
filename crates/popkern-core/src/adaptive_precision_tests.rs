use super::*;

#[test]
fn adaptive_sign_clear() {
    assert_eq!(adaptive_sign(1.0, 1e-10), Some(std::cmp::Ordering::Greater));
    assert_eq!(adaptive_sign(-1.0, 1e-10), Some(std::cmp::Ordering::Less));
}

#[test]
fn adaptive_sign_ambiguous() {
    assert_eq!(adaptive_sign(1e-15, 1e-10), None);
}

#[test]
fn double_precision_add() {
    let a = DoublePrecision::from_f64(1.0);
    let b = DoublePrecision::from_f64(1e-16);
    let sum = a.dd_add(b);
    // Should preserve the small value that f64 would lose
    assert!((sum.to_f64() - 1.0).abs() < 1e-14);
}

#[test]
fn double_precision_mul() {
    let a = DoublePrecision::from_f64(3.0);
    let b = DoublePrecision::from_f64(7.0);
    let product = a.dd_mul(b);
    assert!((product.to_f64() - 21.0).abs() < 1e-14);
}

#[test]
fn double_precision_sign() {
    let val = DoublePrecision::from_f64(0.5);
    assert_eq!(val.sign(1e-10), Some(std::cmp::Ordering::Greater));

    let neg = DoublePrecision::from_f64(-0.5);
    assert_eq!(neg.sign(1e-10), Some(std::cmp::Ordering::Less));

    let tiny = DoublePrecision::from_f64(1e-15);
    assert_eq!(tiny.sign(1e-10), None);
}
