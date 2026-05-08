use super::*;

pub(super) fn approx_p(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

// Evaluate a tri-5 polynomial at (u, v) using the codebase's
// extended basis (no simplex restriction — needed because bi-5 eval
// probes u + v > 1).
pub(super) fn eval_tri5_extended(cps: &[Point3; TRI5_CP_COUNT], u: f64, v: f64) -> Point3 {
    let b = bernstein_tri5_all_extended(u, v);
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for i in 0..TRI5_CP_COUNT {
        out.x += b[i] * cps[i].x;
        out.y += b[i] * cps[i].y;
        out.z += b[i] * cps[i].z;
    }
    out
}
