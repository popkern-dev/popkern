//! Pcurve quality tests for fillet results.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn default_tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn fillet_result_100_percent_pcurve_coverage() {
    let brep = BRep::unit_cube(default_tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let result = popkern_fillet::fillet_edge(&brep, edge, 0.1, default_tol(), None).unwrap();
    for (he_handle, he) in result.brep.half_edges.iter() {
        assert!(
            he.pcurve.is_some(),
            "half-edge {:?} missing pcurve",
            he_handle
        );
    }
}

#[test]
fn fillet_pcurves_have_uv_within_surface_range() {
    use popkern_core::geometry::surface::SurfaceEval;

    let brep = BRep::unit_cube(default_tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let result = popkern_fillet::fillet_edge(&brep, edge, 0.1, default_tol(), None).unwrap();

    for fh in result.brep.faces.keys() {
        let face = &result.brep.faces[fh];
        let Some(ref surf) = face.surface else {
            continue;
        };
        let ((u_min, u_max), (v_min, v_max)) = surf.geometry.parameter_range();

        // Walk outer loop half-edges.
        let loop_ = &result.brep.loops[face.outer_loop];
        let start = loop_.half_edge;
        let mut cursor = start;
        loop {
            let he = &result.brep.half_edges[cursor];
            if let Some(ref pcurve) = he.pcurve {
                let (t0, t1) = pcurve.parameter_range();
                for i in 0..=10 {
                    let t = t0 + (t1 - t0) * (i as f64 / 10.0);
                    let pt = pcurve.point_at(t);
                    // Allow pcurves to extend up to one period beyond
                    // the native range for seam-unwrapped coordinates.
                    let u_period = u_max - u_min;
                    let v_period = v_max - v_min;
                    if u_min.is_finite() {
                        assert!(
                            pt.x >= u_min - u_period - 0.01 && pt.x <= u_max + u_period + 0.01,
                            "pcurve u={} far outside [{}, {}]",
                            pt.x,
                            u_min,
                            u_max
                        );
                    }
                    if v_min.is_finite() {
                        assert!(
                            pt.y >= v_min - v_period - 0.01 && pt.y <= v_max + v_period + 0.01,
                            "pcurve v={} far outside [{}, {}]",
                            pt.y,
                            v_min,
                            v_max
                        );
                    }
                }
            }
            cursor = he.next;
            if cursor == start {
                break;
            }
        }
    }
}
