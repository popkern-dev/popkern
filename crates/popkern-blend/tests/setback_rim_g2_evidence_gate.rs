//! PK-519 diagnostic probe: per-seam rim transverse-G2 residuals on
//! the pentagon flagship (base 1.0, height 1.2, r 0.16, G2 request).
//! Run with `cargo test -p popkern-blend --test probe_pk519_rim_quarters -- --ignored --nocapture`.
//! Removed or promoted into a gate before PK-519 closes.

use popkern_blend::blend_continuity_evidence::{
    ContinuityGates, measure_blend_continuity_evidence,
};
use popkern_blend::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
#[ignore = "PK-519 diagnostic probe"]
fn rim_quarter_residuals() {
    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol).expect("pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .max_by(|a, b| {
            a.1.position
                .z
                .partial_cmp(&b.1.position.z)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
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
            radius: FilletRadius::Constant(0.16),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..Default::default()
    };
    let pre_faces: std::collections::BTreeSet<_> = pyramid.faces.keys().collect();
    let result = fillet_edges(&pyramid, &specs, &options, tol, None).expect("pipeline");
    let brep = &result.brep;
    brep.validate().expect("valid");
    let blend_faces: Vec<_> = brep
        .faces
        .keys()
        .filter(|f| !pre_faces.contains(f))
        .collect();
    let evidence = measure_blend_continuity_evidence(
        brep,
        &blend_faces,
        ContinuityLevel::G2,
        &ContinuityGates::commercial(),
        100,
    );
    eprintln!(
        "achieved={:?} requested={:?} worst_g0={:.3e} worst_g1={:.3e} worst_g2={:.3e} seams={} sharp={} unmeasured={}",
        evidence.achieved,
        evidence.requested,
        evidence.worst_g0_gap,
        evidence.worst_g1_normal,
        evidence.worst_g2_residual,
        evidence.seams.len(),
        evidence.sharp_seams,
        evidence.unmeasured_seams
    );
    let mut seams = evidence.seams.clone();
    seams.sort_by(|a, b| {
        b.worst_g2_residual
            .partial_cmp(&a.worst_g2_residual)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    for s in seams.iter().take(30) {
        // Surface kinds on both sides for classing (rim = Gregory|Nurbs).
        let kinds = seam_kinds(brep, s.edge);
        eprintln!(
            "  seam {:?} kind={:?} gated={} kinds={} g0={:.3e} g1={:.3e} g2={:.3e} angle={:.2}",
            s.edge,
            s.kind,
            s.gated,
            kinds,
            s.worst_g0_gap,
            s.worst_g1_normal,
            s.worst_g2_residual,
            s.worst_normal_angle_deg
        );
    }
}

fn seam_kinds(brep: &BRep, edge: popkern_core::topology::EdgeHandle) -> String {
    let e = &brep.edges[edge];
    let he = e.half_edge;
    let twin = brep.half_edges[he].twin;
    let mut out = Vec::new();
    for h in [he, twin] {
        let face = brep.loops[brep.half_edges[h].loop_].face;
        let kind = match brep.faces[face].surface.as_ref().map(|s| &s.geometry) {
            Some(popkern_core::geometry::SurfaceGeometry::Gregory(_)) => "Gregory",
            Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) => "Nurbs",
            Some(popkern_core::geometry::SurfaceGeometry::Plane { .. }) => "Plane",
            Some(_) => "Other",
            None => "None",
        };
        out.push(kind);
    }
    out.join("|")
}
