//! Feature edges and inferred blend roles (PH-011, gate FI1) on the
//! frozen STEP fixtures.

use std::f64::consts::FRAC_PI_2;
use std::path::Path;

use popkern_analysis::{AnalysisModel, FaceRole, Provenance, RoleSource, Surface};
use popkern_analysis_brep::{FeatureTolerances, detect_features, model_from_breps};
use popkern_analysis_cli::step_model::{
    breps_from_step, model_from_step_file, tessellation_options,
};

fn load(name: &str) -> AnalysisModel {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(format!("../../tests/fixtures/analysis/{name}.step"));
    model_from_step_file(&path).unwrap_or_else(|e| panic!("{name}: {e}"))
}

fn fillet_radii(m: &AnalysisModel) -> Vec<f64> {
    m.faces
        .iter()
        .filter_map(|f| match f.role {
            Some(FaceRole::Fillet { radius }) => {
                assert_eq!(f.role_source, RoleSource::Inferred);
                Some(radius)
            }
            _ => None,
        })
        .collect()
}

/// Re-runs detection on solid 0 for its diagnostics.
fn rerun(
    m: &AnalysisModel,
    tolerances: &FeatureTolerances,
) -> popkern_analysis_brep::SolidFeatures {
    let s = &m.solids[0];
    let mut faces: Vec<_> = m
        .faces
        .iter()
        .filter(|f| f.solid == s.id)
        .cloned()
        .collect();
    for f in &mut faces {
        f.role = None;
    }
    detect_features(
        &s.vertices,
        &s.triangles,
        &s.triangle_face,
        &mut faces,
        tolerances,
        &|_, _| None,
    )
}

#[test]
fn fi1_stepped_bar_fillets_are_inferred_with_their_radius() {
    for r in [1.0, 2.0, 3.0, 5.0] {
        let name = format!("stepped_bar_r{r}");
        let m = load(&name);
        let radii = fillet_radii(&m);
        assert_eq!(radii.len(), 2, "{name}: {radii:?}");
        for x in &radii {
            assert!((x - r).abs() <= 1e-9, "{name}: radius {x}");
        }
        // The fillets replace the shoulder's concave corners.
        assert!(m.concave_edges.is_empty(), "{name}: {:?}", m.concave_edges);
        // Measured margin under the G1 threshold.
        let found = rerun(&m, &FeatureTolerances::default());
        assert!(
            found.max_tangent_angle < 1e-9,
            "{name}: {}",
            found.max_tangent_angle
        );
    }
}

#[test]
fn fi1_negative_controls_bores_and_dimple_are_not_fillets() {
    for name in ["block_two_bores", "block_dimple"] {
        let m = load(name);
        assert!(fillet_radii(&m).is_empty(), "{name}");
        assert!(m.concave_edges.is_empty(), "{name}: {:?}", m.concave_edges);
    }
    // With G1 switched off, the stepped bar's fillets are not inferred:
    // the inference rests on measured tangency, not on surface kind.
    let m = load("stepped_bar_r2");
    let mut faces: Vec<_> = m.faces.clone();
    for f in &mut faces {
        f.role = None;
    }
    let s = &m.solids[0];
    let strict = FeatureTolerances {
        tangent_angle: -1.0,
        ..FeatureTolerances::default()
    };
    detect_features(
        &s.vertices,
        &s.triangles,
        &s.triangle_face,
        &mut faces,
        &strict,
        &|_, _| None,
    );
    assert!(faces.iter().all(|f| f.role.is_none()));
}

#[test]
fn l_bracket_inner_corner_is_the_one_concave_edge() {
    let m = load("l_bracket");
    assert!(fillet_radii(&m).is_empty());
    assert_eq!(m.concave_edges.len(), 1, "{:?}", m.concave_edges);
    let e = &m.concave_edges[0];
    assert!((e.dihedral - FRAC_PI_2).abs() <= 1e-9, "{}", e.dihedral);
    let mut z: Vec<f64> = Vec::new();
    for p in &e.polyline {
        assert!(
            (p[0] - 5.0).abs() < 1e-9 && (p[1] - 5.0).abs() < 1e-9,
            "{p:?}"
        );
        z.push(p[2]);
    }
    z.sort_by(f64::total_cmp);
    assert!(
        z[0].abs() < 1e-9 && (z[z.len() - 1] - 20.0).abs() < 1e-9,
        "{z:?}"
    );
}

/// A frozen gallery part's model with `tolerances`.
fn gallery(name: &str, tolerances: &FeatureTolerances) -> AnalysisModel {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(format!("../../tests/fixtures/analysis/gallery/{name}.step"));
    let bytes = std::fs::read(&path).unwrap_or_else(|e| panic!("{name}: {e}"));
    let bodies = breps_from_step(&bytes).unwrap_or_else(|e| panic!("{name}: {e}"));
    let breps: Vec<_> = bodies.iter().map(|b| b.brep.clone()).collect();
    model_from_breps(
        &breps,
        &tessellation_options(),
        tolerances,
        Provenance::Step {
            file_sha256: String::new(),
        },
    )
    .unwrap_or_else(|e| panic!("{name}: {e}"))
}

#[test]
fn variable_radius_blends_are_inferred_with_exact_freeform_normals() {
    // tapered_fillet: a round growing from 0.04 to 0.18 along the front
    // edge and a constant 0.04 round on the back, both B-spline faces.
    // Exact normals from popkern's surfaces measure their G1 joins
    // below 1e-6 rad, so a 1e-6 tangent tolerance still infers both;
    // the radius is the smallest along the face (measured 0.0398).
    let strict = FeatureTolerances {
        tangent_angle: 1e-6,
        ..FeatureTolerances::default()
    };
    for tolerances in [FeatureTolerances::default(), strict] {
        let m = gallery("tapered_fillet", &tolerances);
        let blends: Vec<(u32, f64)> = m
            .faces
            .iter()
            .filter_map(|f| match (&f.surface, f.role) {
                (Surface::Freeform, Some(FaceRole::Fillet { radius })) => Some((f.id.0, radius)),
                _ => None,
            })
            .collect();
        assert_eq!(blends.len(), 2, "{blends:?}");
        for (_, r) in blends {
            assert!((r - 0.04).abs() <= 0.02 * 0.04, "radius {r}");
        }
    }
}

#[test]
fn chamfers_are_inferred_on_exactly_the_chamfer_faces() {
    // chamfered_bracket: an L section with equal chamfers on the two
    // bottom edges, distance-angle chamfers on the two top edges and
    // two-distance chamfers on the step, one of them on the concave
    // inner corner: six chamfer faces. The narrow top and end faces
    // between two chamfers pass the strip tests too and are excluded
    // because chamfers never share an edge.
    let m = gallery("chamfered_bracket", &FeatureTolerances::default());
    let chamfers: Vec<u32> = m
        .faces
        .iter()
        .filter(|f| f.role == Some(FaceRole::Chamfer))
        .map(|f| f.id.0)
        .collect();
    assert_eq!(chamfers, [1, 2, 7, 9, 10, 13]);
    for f in m.faces.iter().filter(|f| f.role.is_some()) {
        assert_eq!(f.role_source, RoleSource::Inferred);
        let Surface::Plane { normal, .. } = f.surface else {
            panic!("chamfer {} is not planar", f.id.0);
        };
        // Every chamfer is oblique to the axes the bracket is built on.
        assert!(
            normal.iter().filter(|c| c.abs() > 1e-9).count() >= 2,
            "{normal:?}"
        );
    }
    // Negative controls: no chamfer or fillet roles on parts without.
    for name in ["l_bracket", "block_two_bores"] {
        let m = load(name);
        assert!(m.faces.iter().all(|f| f.role.is_none()), "{name}");
    }
}
