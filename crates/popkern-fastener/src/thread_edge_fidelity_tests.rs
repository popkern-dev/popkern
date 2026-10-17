//! A cut thread's edge curves lie on both of their faces (PK-568).
//! On the M10x1.5 stud, the flank x shank curve once left the flank by
//! 2.81e-3 near the runout start; a boundary sample placed there
//! tripped the mesher's off-surface gate and sent the whole flank to
//! the uniform-grid fallback. Every thread edge is now checked at
//! 1000 uniform fractions against both adjacent faces, at the body
//! tolerance (measured worst 9.7e-7 after PK-675; 1.9e-6 on the
//! quadratic-arc helix).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::{
    project_point_to_uv_multisheet, project_point_to_uv_near,
};
use popkern_core::tolerance::Tolerance;

use crate::builders::{RodParams, threaded_rod};
use crate::spec::ThreadSpec;

#[test]
fn stud_thread_edges_lie_on_both_faces() {
    let tolerance = Tolerance::default();
    let spec = match ThreadSpec::parse("M10x1.5") {
        Ok(s) => s,
        Err(e) => panic!("spec: {e}"),
    };
    let params = RodParams {
        spec,
        length: 8.0,
        tip_chamfer: false,
        cosmetic_thread: false,
    };
    let rod = match threaded_rod(&params, tolerance) {
        Ok(r) => r,
        Err(e) => panic!("stud: {e}"),
    };
    let mut worst = (0.0f64, String::new());
    for (edge_handle, edge) in rod.edges.iter() {
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let Some(first) = rod.half_edges.get(edge.half_edge) else {
            continue;
        };
        for side in [edge.half_edge, first.twin] {
            let face = rod
                .half_edges
                .get(side)
                .and_then(|d| rod.loops.get(d.loop_))
                .map(|l| l.face);
            let Some(surface) = face
                .and_then(|f| rod.faces.get(f))
                .and_then(|f| f.surface.as_ref())
            else {
                continue;
            };
            let mut previous: Option<(f64, f64)> = None;
            for i in 0..=1000 {
                let fraction = f64::from(i) / 1000.0;
                let t = curve.t_start + (curve.t_end - curve.t_start) * fraction;
                let point = curve.geometry.point_at(t);
                let distance =
                    |uv: (f64, f64)| point.distance_to(&surface.geometry.point_at(uv.0, uv.1));
                // Continue from the previous witness; solve cold only
                // when the continuation leaves the validator's band
                // (check 9: ten times the tolerance).
                let near = previous
                    .and_then(|(u, v)| project_point_to_uv_near(&surface.geometry, &point, u, v))
                    .filter(|&uv| distance(uv) <= 10.0 * tolerance.linear());
                let Some(uv) =
                    near.or_else(|| project_point_to_uv_multisheet(&surface.geometry, &point))
                else {
                    break;
                };
                previous = Some(uv);
                let d = distance(uv);
                if d > worst.0 {
                    worst = (d, format!("{edge_handle:?} on {face:?} at f = {fraction}"));
                }
            }
        }
    }
    assert!(
        worst.0 <= tolerance.linear(),
        "thread edge leaves its face by {:.3e}: {}",
        worst.0,
        worst.1
    );
}
