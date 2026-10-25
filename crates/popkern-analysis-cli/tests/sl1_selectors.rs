//! Gate SL1 (docs/notes/analysis/12-verification-benchmarks.md):
//! selectors on the frozen STEP fixtures resolve to the exact faces;
//! a 0-match and a 2-match reference decline with typed errors. Also
//! the adapter postconditions on the same fixtures.

use std::f64::consts::PI;
use std::path::{Path, PathBuf};

use popkern_analysis::{
    AnalysisError, AnalysisModel, Expect, FaceDef, FaceId, FaceRef, Selector, Tolerance,
};
use popkern_analysis_brep::{AdapterError, FeatureTolerances, model_from_breps, resolve_faces};
use popkern_analysis_cli::setup_file::parse_setup;
use popkern_analysis_cli::step_model::{
    breps_from_step, model_from_step_file, tessellation_options,
};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("../../tests/fixtures/analysis/{name}.step"))
}

fn model(name: &str) -> AnalysisModel {
    model_from_step_file(&fixture(name)).unwrap_or_else(|e| panic!("{name}: {e}"))
}

fn def(name: &str, selector: Selector, expect: Expect) -> FaceDef {
    FaceDef {
        name: name.to_string(),
        reference: FaceRef::Select(selector),
        expect,
    }
}

fn plane(normal: [f64; 3], through: [f64; 3]) -> Selector {
    Selector::Plane {
        normal,
        through,
        angle_tol: None,
        dist_tol: None,
    }
}

fn cylinder(axis: [f64; 3], radius: f64, through: Option<[f64; 3]>) -> Selector {
    Selector::Cylinder {
        axis,
        radius,
        through,
        angle_tol: None,
        dist_tol: None,
    }
}

fn resolve_one(m: &AnalysisModel, selector: Selector) -> Result<Vec<FaceId>, AnalysisError> {
    let defs = [def("f", selector, Expect::default())];
    resolve_faces(m, &defs, &Tolerance::default()).map(|r| r.faces[0].faces.clone())
}

fn ids(ids: &[u32]) -> Vec<FaceId> {
    ids.iter().map(|&i| FaceId(i)).collect()
}

#[test]
fn plane_cylinder_and_sphere_selectors_resolve_exactly() {
    let a = model("block_two_bores");
    let cases = [
        (plane([0.0, 0.0, -1.0], [0.0, 0.0, 0.0]), 0),
        (plane([-1.0, 0.0, 0.0], [0.0, 3.0, 3.0]), 1),
        (plane([0.0, 1.0, 0.0], [1.0, 10.0, 1.0]), 2),
        (plane([1.0, 0.0, 0.0], [10.0, 0.0, 0.0]), 3),
        (plane([0.0, -1.0, 0.0], [4.0, 0.0, 9.0]), 4),
        (plane([0.0, 0.0, 1.0], [5.0, 5.0, 10.0]), 5),
        (cylinder([0.0, 0.0, 1.0], 1.5, Some([2.5, 5.0, 0.0])), 6),
        // either axis sense, any point on the axis
        (cylinder([0.0, 0.0, -1.0], 1.5, Some([7.5, 5.0, 3.0])), 7),
        (Selector::NearestTo([2.5, 5.0, 5.0]), 6),
    ];
    for (selector, face) in cases {
        assert_eq!(
            resolve_one(&a, selector.clone()),
            Ok(ids(&[face])),
            "{selector:?}"
        );
    }
    let b = model("block_dimple");
    let sphere = Selector::Sphere {
        centre: [5.0, 5.0, 11.0],
        radius: 3.0,
        dist_tol: None,
    };
    assert_eq!(resolve_one(&b, sphere), Ok(ids(&[6])));
    assert_eq!(
        resolve_one(&b, plane([0.0, 0.0, 1.0], [0.0, 0.0, 10.0])),
        Ok(ids(&[5]))
    );
}

#[test]
fn a_zero_match_reference_is_unresolved() {
    let a = model("block_two_bores");
    let b = model("block_dimple");
    let misses = [
        (&a, plane([0.0, 0.0, 1.0], [5.0, 5.0, 7.0])),
        // the bottom face's inward normal: the sense must be outward
        (&a, plane([0.0, 0.0, 1.0], [0.0, 0.0, 0.0])),
        (&a, cylinder([0.0, 0.0, 1.0], 1.6, None)),
        (&a, cylinder([1.0, 0.0, 0.0], 1.5, None)),
        (
            &b,
            Selector::Sphere {
                centre: [5.0, 5.0, 11.0],
                radius: 2.9,
                dist_tol: None,
            },
        ),
    ];
    for (m, selector) in misses {
        assert_eq!(
            resolve_one(m, selector.clone()),
            Err(AnalysisError::SelectorUnresolved {
                name: "f".to_string()
            }),
            "{selector:?}"
        );
    }
}

#[test]
fn a_two_match_reference_is_ambiguous_with_candidates_in_id_order() {
    let a = model("block_two_bores");
    let both_bores = cylinder([0.0, 0.0, 1.0], 1.5, None);
    assert_eq!(
        resolve_one(&a, both_bores.clone()),
        Err(AnalysisError::SelectorAmbiguous {
            name: "f".to_string(),
            matched: 2,
            candidates: ids(&[6, 7]),
        })
    );
    // The same selector is fine when the setup expects two.
    for expect in [Expect::Exactly(2), Expect::AtLeast(2)] {
        let defs = [def("bores", both_bores.clone(), expect)];
        let r = resolve_faces(&a, &defs, &Tolerance::default()).unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(r.faces[0].faces, ids(&[6, 7]));
    }
    let caps = Selector::AnyOf(vec![
        plane([0.0, 0.0, 1.0], [0.0, 0.0, 10.0]),
        plane([0.0, 0.0, -1.0], [0.0, 0.0, 0.0]),
    ]);
    assert!(matches!(
        resolve_one(&a, caps),
        Err(AnalysisError::SelectorAmbiguous { matched: 2, .. })
    ));
}

#[test]
fn a_setup_file_resolves_end_to_end() {
    let setup = parse_setup(
        r#"
units = "mm"
material = "PETG (FDM)"
[faces]
base = { plane = { normal = [0, 0, -1], through = [0, 0, 0] } }
bore = { cylinder = { axis = [0, 0, 1], radius = 1.5, through = [7.5, 5, 0] } }
side = { id = 3 }
[[support]]
faces = ["base"]
kind = "fixed"
[[load]]
faces = ["bore"]
force_n = 20
direction = "-x"
"#,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let r = resolve_faces(&model("block_two_bores"), &setup.faces, &setup.tolerance)
        .unwrap_or_else(|e| panic!("{e}"));
    let resolved: Vec<(&str, Vec<FaceId>)> = r
        .faces
        .iter()
        .map(|f| (f.name.as_str(), f.faces.clone()))
        .collect();
    assert_eq!(
        resolved,
        [
            ("base", ids(&[0])),
            ("bore", ids(&[7])),
            ("side", ids(&[3]))
        ]
    );
    assert_eq!(r.warnings.len(), 1, "a raw id reference is warned");
}

#[test]
fn the_adapter_output_is_closed_outward_and_volume_consistent() {
    // (fixture, analytic volume, curved area): mesh volume is within
    // chord · curved area of the exact volume.
    let h = 2.0;
    let r = 3.0;
    let cases = [
        (
            "block_two_bores",
            1000.0 - 2.0 * PI * 1.5 * 1.5 * 10.0,
            2.0 * 2.0 * PI * 1.5 * 10.0,
        ),
        (
            "block_dimple",
            1000.0 - PI * h * h * (3.0 * r - h) / 3.0,
            2.0 * PI * r * h,
        ),
    ];
    for (name, exact, curved_area) in cases {
        let m = model(name);
        assert_eq!(m.validate(), Ok(()), "{name}");
        assert_eq!(m.solids.len(), 1);
        let s = &m.solids[0];
        let volume = popkern_analysis_brep::orient::signed_volume(&s.vertices, &s.triangles);
        let bound = m.chord_tolerance * curved_area;
        assert!(
            (volume - exact).abs() <= bound,
            "{name}: {volume} vs {exact} (bound {bound})"
        );
        let body = s.exact_volume.unwrap_or(f64::NAN);
        assert!((body - exact).abs() <= bound, "{name}: body_volume {body}");
        assert_eq!(popkern_analysis_brep::weld::open_edges(&s.triangles), 0);
        // flat side faces measure exactly; the bored top within chord
        assert!((m.faces[1].area - 100.0).abs() < 1e-9, "{name}");
    }
    let a = model("block_two_bores");
    let bored = 100.0 - 2.0 * PI * 1.5 * 1.5;
    for face in [0, 5] {
        let area = a.faces[face].area;
        assert!(
            (area - bored).abs() <= 1e-3 * bored,
            "face {face}: {area} vs {bored}"
        );
    }
}

#[test]
fn an_open_shell_declines_not_watertight() {
    let bytes = std::fs::read(fixture("block_two_bores")).unwrap_or_else(|e| panic!("{e}"));
    let mut breps: Vec<_> = breps_from_step(&bytes)
        .unwrap_or_else(|e| panic!("{e}"))
        .into_iter()
        .map(|b| b.brep)
        .collect();
    let top = breps[0]
        .faces
        .keys()
        .nth(5)
        .unwrap_or_else(|| panic!("6 faces"));
    breps[0].faces.remove(top);
    let provenance = popkern_analysis::Provenance::Step {
        file_sha256: String::new(),
    };
    let result = model_from_breps(
        &breps,
        &tessellation_options(),
        &FeatureTolerances::default(),
        provenance,
    );
    assert!(
        matches!(result, Err(AdapterError::NotWatertight { open_edges, .. }) if open_edges > 0),
        "{result:?}"
    );
}

#[test]
fn two_imports_give_identical_models_and_resolutions() {
    for name in ["block_two_bores", "block_dimple"] {
        let first = model(name);
        let second = model(name);
        assert_eq!(first, second, "{name}");
        let defs = [def(
            "n",
            Selector::NearestTo([5.0, 5.0, 9.0]),
            Expect::default(),
        )];
        let r1 = resolve_faces(&first, &defs, &Tolerance::default());
        let r2 = resolve_faces(&second, &defs, &Tolerance::default());
        assert_eq!(r1, r2);
    }
}
