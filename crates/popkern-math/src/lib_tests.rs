use super::PortableMath;

/// The portable functions agree with the platform ones to a few ulps
/// (they compute the same functions; only the last bits may differ).
#[test]
fn values_track_the_inherent_functions() {
    for i in 0..2000 {
        let x = -10.0 + 0.01 * f64::from(i);
        let close = |a: f64, b: f64| (a - b).abs() <= 4.0 * f64::EPSILON * (1.0 + b.abs());
        assert!(close(x.psin(), x.sin()));
        assert!(close(x.pcos(), x.cos()));
        assert!(close(x.patan2(0.7), x.atan2(0.7)));
        assert!(close(x.phypot(0.3), x.hypot(0.3)));
        let (s, c) = x.psin_cos();
        assert_eq!((s, c), (x.psin(), x.pcos()));
        let xf = x as f32;
        assert!((xf.psin() - xf.sin()).abs() <= 4.0 * f32::EPSILON);
    }
    assert!((2.0f64.ppowf(0.5) - std::f64::consts::SQRT_2).abs() < 1e-15);
    assert!((1.0f64.pexp().pln() - 1.0).abs() < 1e-15);
}

#[test]
fn edge_values_are_exact() {
    assert_eq!(0.0f64.psin().to_bits(), 0.0f64.to_bits());
    assert_eq!(0.0f64.pcos(), 1.0);
    assert_eq!(1.0f64.pacos(), 0.0);
    assert_eq!(0.0f64.patan2(1.0), 0.0);
    assert!(f64::NAN.psin().is_nan());
    assert_eq!(1.0f64.pln(), 0.0);
}

/// Pinned bits. The portable functions are the kernel's numeric
/// contract: the same input gives these bits on every target (the
/// wasm leg is `gallery_parity`). A changed pin means the `libm`
/// version changed — a kernel-wide numeric change to review, not a
/// test to refresh.
#[test]
fn bits_are_pinned() {
    let x = 0.7310585786300049f64;
    assert_eq!(x.psin().to_bits(), 0x3fe5_5d74_7e06_2c92);
    assert_eq!(x.pcos().to_bits(), 0x3fe7_d2ae_a980_dbbe);
    assert_eq!(x.patan2(-0.3).to_bits(), 0x3fff_5cff_e8cd_141b);
    assert_eq!(x.ppowf(2.7).to_bits(), 0x3fdb_7831_905a_ab76);
    assert_eq!(x.pexp().to_bits(), 0x4000_9e44_242c_a509);
    assert_eq!(x.pacos().to_bits(), 0x3fe8_0792_35d8_aaa8);
}
