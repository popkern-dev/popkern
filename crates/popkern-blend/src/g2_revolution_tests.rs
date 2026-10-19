//! PK-389: the G2 fillet of a circular plane×cylinder rim is an exact
//! surface of revolution whose rows are the rails.

use crate::blend_continuity_evidence::{ContinuityGates, measure_blend_continuity_evidence};
use crate::types::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec};
use crate::{BlendResult, fillet_edges};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn rim_edges(brep: &BRep, z: f64) -> Vec<popkern_core::topology::EdgeHandle> {
    let mut rim = Vec::new();
    for (e, d) in brep.edges.iter() {
        let he = &brep.half_edges[d.half_edge];
        let a = brep.vertices[he.origin].position;
        let b = brep.vertices[brep.half_edges[he.twin].origin].position;
        let f1 = brep.loops[he.loop_].face;
        let f2 = brep.loops[brep.half_edges[he.twin].loop_].face;
        if f1 != f2 && (a.z - z).abs() < 1e-9 && (b.z - z).abs() < 1e-9 {
            rim.push(e);
        }
    }
    rim
}

fn g2_rim(radius: f64) -> BlendResult {
    let brep = BRep::cylinder(1.0, 2.0, tol()).expect("cylinder");
    let specs: Vec<FilletSpec> = rim_edges(&brep, 2.0)
        .into_iter()
        .map(|edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        })
        .collect();
    assert!(!specs.is_empty(), "cylinder has a top rim");
    let opts = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..FilletOptions::default()
    };
    match fillet_edges(&brep, &specs, &opts, tol(), None) {
        Ok(r) => r,
        Err(e) => panic!("G2 rim fillet declined: {e}"),
    }
}

fn blend_faces(result: &BlendResult) -> Vec<(FaceHandle, &popkern_core::geometry::NurbsSurface)> {
    result
        .brep
        .faces
        .iter()
        .filter_map(|(h, f)| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Nurbs(n)) => Some((h, n)),
            _ => None,
        })
        .collect()
}

#[test]
fn circular_rim_g2_is_an_exact_revolution_that_validates() {
    let result = g2_rim(0.2);
    let issues = result.brep.validate_all();
    assert!(issues.is_empty(), "served body must validate: {issues:?}");
    let blends = blend_faces(&result);
    assert_eq!(blends.len(), 1, "one blend strip");
    let (face, nurbs) = blends[0];
    assert_eq!(nurbs.degree_u, 2, "rational-quadratic revolution in u");
    // Four 90° segments, written as a ring in periodic form (PK-893):
    // eight distinct rows and the first two repeated.
    assert_eq!(nurbs.num_u, 10, "four 90° segments");
    assert_eq!(
        nurbs.u_period(),
        Some(1.0),
        "a closed rim's strip is a ring"
    );
    assert_eq!(nurbs.degree_v, 5, "degree-5 G2 section in v");
    assert!(
        nurbs
            .weights
            .iter()
            .any(|w| (w - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-12),
        "rational arc weights present"
    );
    let evidence = measure_blend_continuity_evidence(
        &result.brep,
        &[face],
        ContinuityLevel::G2,
        &ContinuityGates::commercial(),
        33,
    );
    assert!(evidence.measured(), "{evidence:?}");
    assert!(
        !evidence.downgraded(),
        "measured G2 must hold: {evidence:?}"
    );
    assert!(
        !result
            .edge_warnings
            .iter()
            .any(|(_, w)| matches!(w, crate::error::BlendWarning::ContinuityDowngraded { .. })),
        "no downgrade: {:?}",
        result.edge_warnings
    );
}

#[test]
fn revolution_rows_sit_on_the_rails_for_every_radius() {
    for radius in [0.05, 0.2, 0.45] {
        let result = g2_rim(radius);
        let deviations = result.brep.validate_curve_surface_all();
        assert!(deviations.is_empty(), "r={radius}: {deviations:?}");
    }
}
