use super::*;

// ---- Determinism ----

#[test]
fn from_gtri5_deterministic() {
    let params = Bi5FreeParams::zero(5);
    let g_tri = apply_free_params(&params, Point3::new(1.0, 2.0, 3.0)).unwrap();
    let a = Bi5Guide::from_gtri5(&g_tri);
    let b = Bi5Guide::from_gtri5(&g_tri);
    assert_eq!(a, b);
}
