//! Tests for the KP20 paper-border extractors.

use super::*;
use crate::{FilletRadius, FilletSpec, preview_fillet};
use popkern_core::geometry::gregory::kp20_paper_border_corner_mismatch;
use popkern_core::tolerance::Tolerance;

/// End-to-end extraction gate on real rolling-ball data: the loop
/// closes exactly, every ring march is exact (cylinder/plane),
/// and the corner mismatch is finite and bounded — the measured
/// G¹-input incompatibility at c, recorded rather than hidden.
#[test]
fn pyramid_apex_extraction_is_closed_and_exact() {
    let tol = Tolerance::default();
    let pyramid = popkern_core::topology::BRep::n_gonal_pyramid(5, 1.0, 1.2, tol).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .find(|(_, v)| v.position.x.abs() < 1e-9 && v.position.y.abs() < 1e-9 && v.position.z > 0.0)
        .map(|(h, _)| h)
        .expect("apex");
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = pyramid.edges[eh].half_edge;
            let origin = pyramid.half_edges[he].origin;
            let twin_origin = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
            (origin == apex || twin_origin == apex).then_some(eh)
        })
        .collect();
    let specs: Vec<FilletSpec> = apex_edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.12),
        })
        .collect();
    let preview = preview_fillet(&pyramid, &specs, tol).expect("preview");
    let fillet_surfaces: Vec<FilletSurfaceResult> = preview
        .per_edge
        .iter()
        .map(|eg| FilletSurfaceResult {
            surface: eg.blend_surface.clone(),
            trim_curve_1: eg.trim_curve_1.clone(),
            trim_curve_2: eg.trim_curve_2.clone(),
            continuity_downgrade: None,
        })
        .collect();
    let indices: Vec<usize> = (0..fillet_surfaces.len()).collect();
    let boundary =
        crate::blend::compute_blend_boundary(apex, &indices, &fillet_surfaces, &pyramid, &tol)
            .expect("boundary");
    let (border, report) =
        kp20_paper_border_from_blend_boundary(&boundary, &fillet_surfaces, &pyramid)
            .expect("extract");
    assert_eq!(report.n, 5);
    assert!(
        report.max_loop_gap < 1.0e-9,
        "loop gap {}",
        report.max_loop_gap
    );
    assert_eq!(
        report.first_order_marches, 0,
        "cylinder/plane marches must be exact"
    );
    let c = (core::f64::consts::TAU / 5.0).cos();
    let mut mismatch = 0.0_f64;
    for sector in &border.sectors {
        mismatch = mismatch.max(kp20_paper_border_corner_mismatch(
            5,
            &sector[0].rows,
            &sector[1].rows,
            2.0 * c / 3.0,
        ));
    }
    // Rolling-ball input is G1-class: the corner mismatch is real
    // and must be REPORTED, not zero. Bound it loosely so gross
    // extraction bugs (orders of magnitude off) still fail.
    assert!(
        mismatch > 1.0e-4 && mismatch < 1.0,
        "corner mismatch out of expected G1-class band: {mismatch:e}",
    );
}
