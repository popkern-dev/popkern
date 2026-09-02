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

/// PK-399: on the blunt pyramids (h = 3, contact lines meeting at
/// ~21° on the shared faces) the setback corners are EXACT trim
/// crossings — on both strips to roundoff, with coherent strip
/// normals — and no corner is served by the gapped-corner
/// containment. The former window-scan minimiser stalled at 3.4e-3
/// here and every corner read as gapped. The residual corner-block
/// mismatch is the recorded in-plane compromise of the shared corner
/// block (1.2e-4 on the pentagon G2; the mean block without the
/// projection measured rim tangency 6.0e-2).
#[test]
fn blunt_pyramid_setback_corners_are_exact_trim_crossings() {
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    let tol = Tolerance::default();
    for (n, continuity) in [
        (5_usize, crate::types::ContinuityLevel::G2),
        (6, crate::types::ContinuityLevel::G2),
        (5, crate::types::ContinuityLevel::G1),
    ] {
        let pyramid =
            popkern_core::topology::BRep::n_gonal_pyramid(n, 1.0, 3.0, tol).expect("pyramid");
        let apex = pyramid
            .vertices
            .iter()
            .find(|(_, v)| {
                v.position.x.abs() < 1e-9 && v.position.y.abs() < 1e-9 && v.position.z > 0.0
            })
            .map(|(h, _)| h)
            .expect("apex");
        let apex_position = pyramid.vertices[apex].position;
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
                radius: FilletRadius::Constant(0.05),
            })
            .collect();
        let preview = crate::preview_fillet_with_continuity(&pyramid, &specs, tol, continuity)
            .expect("preview");
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
        let (_, report) = kp20_paper_border_from_fillet_setback(&fillet_surfaces, apex_position)
            .expect("setback extract");
        let label = format!("n={n} {continuity:?}");
        assert_eq!(report.gapped_corners, 0, "[{label}] gapped corners");
        assert_eq!(report.reflex_corners, 0, "[{label}] reflex corners");
        assert!(
            report.max_loop_gap < 1.0e-12,
            "[{label}] corner gap {:.3e}",
            report.max_loop_gap
        );
        assert_eq!(report.trim_corners.len(), n);
        for (k, arc) in report.setback_arcs.iter().enumerate() {
            let next = &report.setback_arcs[(k + 1) % n];
            let corner = arc.corner_next;
            assert!(
                corner.distance_to(&next.corner_prev) < 1.0e-12,
                "[{label}] corner {k} not shared"
            );
            let mut normals = Vec::with_capacity(2);
            for fi in [arc.fillet_index, next.fillet_index] {
                let surface = &fillet_surfaces[fi].surface;
                let (u, v) = project_point_to_uv(surface, &corner).expect("corner projects");
                let foot = surface.point_at(u, v).distance_to(&corner);
                assert!(
                    foot < 1.0e-9,
                    "[{label}] corner {k} is {foot:.3e} off strip {fi}"
                );
                normals.push(surface.normal_at(u, v).normalized().expect("normal"));
            }
            let sin_angle = normals[0].cross(&normals[1]).length();
            assert!(
                sin_angle < 1.0e-9,
                "[{label}] strips disagree at corner {k}: sin {sin_angle:.3e}"
            );
        }
        let worst_block = report
            .corner_block_mismatch
            .iter()
            .copied()
            .fold(0.0_f64, f64::max);
        // PK-400: both pieces at every corner carry ONE block by
        // construction (ring-wide fit); the assembly's overlap average
        // is a no-op to roundoff (was 2e-5..1.2e-4 with pinned fits).
        assert!(
            worst_block < 1.0e-8,
            "[{label}] corner-block mismatch {worst_block:.3e}"
        );
        assert!(!report.ring_fit_fallback, "[{label}] ring fit fell back");
        assert_eq!(report.section_axis_fallbacks, 0, "[{label}] axis fallbacks");
        assert_turns_are_the_face_angle(&label, n, 3.0, &report.corner_turns);
    }
}

/// Every cap corner turns by the shared face's apex angle
/// (`2·asin(sin(π/n) / √(1 + h²))` for a regular pyramid of unit base
/// radius and height `h`): the rim arcs' ends are cross-section cuts
/// whatever the radii (PK-401).
fn assert_turns_are_the_face_angle(label: &str, n: usize, h: f64, turns: &[f64]) {
    let face_angle = 2.0 * ((core::f64::consts::PI / n as f64).sin() / (1.0 + h * h).sqrt()).asin();
    assert_eq!(turns.len(), n, "[{label}] one turn per corner");
    for (k, turn) in turns.iter().enumerate() {
        assert!(
            (turn - face_angle).abs() < 1.0e-6,
            "[{label}] corner {k} turns {:.4}° against the face angle {:.4}°",
            turn.to_degrees(),
            face_angle.to_degrees()
        );
    }
}

/// Mixed-radius layouts (hexagon6_variable 0.02/0.04/0.06/0.06/0.04/
/// 0.02, pentagon5_variable 0.02/0.04/0.06/0.04/0.02, h = 0.8): the
/// exact trim-intersection corners of one strip sit at different
/// stations along its edge, and the straight-uv rim between them
/// turned the cap boundary FLAT or REFLEX at the corner between the
/// two small strips (−25.9° on the hexagon, −13.3° on the pentagon;
/// PK-399 served it through the χ̃-calibrated containment at live g1
/// 8.8e-4). With cross-sectional rim ends (PK-401) every corner turns
/// by the face angle, no corner is reflex or gapped, and both pieces
/// at every corner carry one block.
#[test]
fn mixed_radius_setback_corners_are_convex_face_angle_turns() {
    let tol = Tolerance::default();
    for (n, radii) in [
        (6_usize, vec![0.02, 0.04, 0.06, 0.06, 0.04, 0.02]),
        (5, vec![0.02, 0.04, 0.06, 0.04, 0.02]),
    ] {
        let pyramid =
            popkern_core::topology::BRep::n_gonal_pyramid(n, 1.0, 0.8, tol).expect("pyramid");
        let apex = pyramid
            .vertices
            .iter()
            .find(|(_, v)| {
                v.position.x.abs() < 1e-9 && v.position.y.abs() < 1e-9 && v.position.z > 0.0
            })
            .map(|(h, _)| h)
            .expect("apex");
        let apex_position = pyramid.vertices[apex].position;
        let mut apex_edges: Vec<_> = pyramid
            .edges
            .iter()
            .filter_map(|(eh, _)| {
                let he = pyramid.edges[eh].half_edge;
                let origin = pyramid.half_edges[he].origin;
                let twin_origin = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
                (origin == apex || twin_origin == apex).then_some(eh)
            })
            .collect();
        apex_edges.sort();
        let specs: Vec<FilletSpec> = apex_edges
            .iter()
            .zip(radii.iter().copied())
            .map(|(&edge, radius)| FilletSpec {
                edge,
                radius: FilletRadius::Constant(radius),
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
        let (_, report) = kp20_paper_border_from_fillet_setback(&fillet_surfaces, apex_position)
            .expect("setback extract");
        let label = format!("n={n} variable");
        assert_eq!(report.gapped_corners, 0, "[{label}] gapped corners");
        assert_eq!(report.reflex_corners, 0, "[{label}] reflex corners");
        assert_eq!(report.section_axis_fallbacks, 0, "[{label}] axis fallbacks");
        assert!(!report.ring_fit_fallback, "[{label}] ring fit fell back");
        assert!(
            report.max_loop_gap < 1.0e-12,
            "[{label}] corner gap {:.3e}",
            report.max_loop_gap
        );
        assert_turns_are_the_face_angle(&label, n, 0.8, &report.corner_turns);
        let worst_block = report
            .corner_block_mismatch
            .iter()
            .copied()
            .fold(0.0_f64, f64::max);
        assert!(
            worst_block < 1.0e-8,
            "[{label}] corner-block mismatch {worst_block:.3e}"
        );
    }
}
