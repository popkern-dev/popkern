//! Reachability + acceptance tests for the KP15 TBdata rung (PK-090)
//! on REAL fillet fixtures — the backstop must be measured, not
//! merely wired.

use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn fixture(n: usize, height: f64, radii: &[f64]) -> (Vec<FilletSurfaceResult>, Point3) {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(n, 1.0, height, tol).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| {
            a.1.position
                .z
                .partial_cmp(&b.1.position.z)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|(h, v)| (h, v.position))
        .expect("apex");
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = pyramid.edges[eh].half_edge;
            let origin = pyramid.half_edges[he].origin;
            let twin_origin = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
            (origin == apex.0 || twin_origin == apex.0).then_some(eh)
        })
        .collect();
    let specs: Vec<crate::FilletSpec> = apex_edges
        .iter()
        .zip(radii.iter())
        .map(|(&edge, &radius)| crate::FilletSpec {
            edge,
            radius: crate::FilletRadius::Constant(radius),
        })
        .collect();
    let preview =
        crate::preview_fillet_with_continuity(&pyramid, &specs, tol, crate::ContinuityLevel::G2)
            .expect("G2 preview");
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
    (fillet_surfaces, apex.1)
}

/// Ring corners for the rung, taken from the extractor itself
/// (identity rotation) — exercises the alignment path end-to-end.
fn extractor_corners(fillet_surfaces: &[FilletSurfaceResult], apex: Point3) -> Vec<Point3> {
    let (border, _) =
        crate::construct::kp20_paper_border_extract::kp20_paper_border_from_fillet_setback(
            fillet_surfaces,
            apex,
        )
        .expect("border");
    border
        .sectors
        .iter()
        .map(|sector| sector[0].rows[0][6])
        .collect()
}

#[test]
fn tbdata_rung_commits_on_pentagon5_uniform() {
    let (fillet_surfaces, apex) = fixture(5, 0.8, &[0.05; 5]);
    let corners = extractor_corners(&fillet_surfaces, apex);
    let patches = try_kp15_setback_tbdata_rung(apex, &corners, &fillet_surfaces)
        .expect("TBdata rung must commit on the exact-corner fixture");
    assert_eq!(patches.len(), 5);
    for (i, p) in patches.iter().enumerate() {
        assert!(
            !p.self_intersects_with_tolerance(SELF_INT_TOL),
            "patch {i} folded"
        );
    }
}

/// PIN → GATE (PK-131, flipped 2026-07-18): the TBdata rung COMMITS
/// six unfolded patches on the skew-corner variable-radius hexagon.
/// The pin originally accepted a fold-descent/eq-gates decline as the
/// measured all-or-nothing outcome; the PK-087/PK-090 fold-descent
/// hardening (damped rim LSQ, outlier repair-not-drop, seam-curve
/// refresh) closed it silently — the either-way pin masked the win.
/// A decline is a regression now.
#[test]
fn tbdata_rung_outcome_is_pinned_on_hexagon6_variable() {
    let (fillet_surfaces, apex) = fixture(6, 0.8, &[0.02, 0.04, 0.06, 0.06, 0.04, 0.02]);
    let corners = extractor_corners(&fillet_surfaces, apex);
    let patches = try_kp15_setback_tbdata_rung(apex, &corners, &fillet_surfaces)
        .expect("TBdata rung must commit on hexagon6_variable (PK-131 flipped pin)");
    assert_eq!(patches.len(), 6);
    for (i, p) in patches.iter().enumerate() {
        assert!(
            !p.self_intersects_with_tolerance(SELF_INT_TOL),
            "patch {i} folded"
        );
    }
}
