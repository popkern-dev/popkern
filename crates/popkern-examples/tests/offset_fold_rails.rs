//! Outward offset of the wavy plate up to its reported max-safe
//! distance (PK-883): every edge's pcurve images trace their edge.
//!
//! The offset measured 1.1e-3..1.9e-3 off at EVERY distance, not just
//! near the fold the max-safe report is bounded by: nine-station rail
//! polylines, a corner-sphere chart pole on the top-corner caps, and —
//! past half the vanish distance — rails, corner vertices and the
//! rolling-ball radius read by closest-point projection onto the offset
//! surface, which lands on the wrong point once the spine lies beyond
//! the offset surface's focal point. Rails, corners and the join
//! classifier now take the original face's normal-field image.

use popkern_core::tolerance::Tolerance;
use popkern_shell::{ShellDirection, ShellOptions, max_safe_offset, offset::offset_solid};

#[test]
fn the_wavy_plate_offsets_honestly_up_to_its_max_safe_distance() {
    let tol = Tolerance::default();
    let fixture = popkern_examples::offset_corpus::wavy_plate(tol).expect("wavy plate");
    let max = max_safe_offset(&fixture.brep, ShellDirection::Outward, tol)
        .expect("max-safe report")
        .max_distance;
    for fraction in [0.2, 0.5, 0.8, 0.9, 1.0] {
        let d = max * fraction;
        let body = offset_solid(&fixture.brep, d, &ShellOptions::default(), tol)
            .unwrap_or_else(|e| panic!("offset {d}: {e}"))
            .brep;
        let worst = body
            .half_edges
            .keys()
            .filter_map(|he| body.measure_half_edge_pcurve_deviation(he))
            .map(|m| m.max_deviation)
            .fold(0.0_f64, f64::max);
        assert!(
            worst <= tol.linear(),
            "offset {d}: worst pcurve image {worst:.3e}"
        );
        let errors = body.validate();
        assert!(errors.is_ok(), "offset {d}: {errors:?}");
    }
}
