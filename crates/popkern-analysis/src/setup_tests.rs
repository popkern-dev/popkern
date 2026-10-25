use crate::face_ref::{Expect, FaceRef, Selector, Tolerance};
use crate::options::AnalysisOptions;
use crate::report::SolveReport;
use crate::setup::*;

fn plane(normal: [f64; 3]) -> FaceRef {
    FaceRef::Select(Selector::Plane {
        normal,
        through: [0.0; 3],
        angle_tol: None,
        dist_tol: None,
    })
}

fn face(name: &str, reference: FaceRef) -> FaceDef {
    FaceDef {
        name: name.to_string(),
        reference,
        expect: Expect::default(),
    }
}

fn valid() -> StructuralSetup {
    StructuralSetup {
        material: "PLA (FDM)".to_string(),
        build_direction: Some([0.0, 0.0, 1.0]),
        tolerance: Tolerance::default(),
        faces: vec![
            face("base", plane([0.0, 0.0, -1.0])),
            face(
                "hook",
                FaceRef::Select(Selector::NearestTo([40.0, 0.0, 12.0])),
            ),
        ],
        supports: vec![SupportSpec {
            faces: vec!["base".to_string()],
            kind: SupportKind::Fixed,
        }],
        loads: vec![LoadSpec {
            faces: vec!["hook".to_string()],
            total_force_n: 50.0,
            direction: LoadDirection::Vector([0.0, 0.0, -1.0]),
        }],
        thermal: None,
        modal: None,
        infill: None,
        generative: None,
        options: AnalysisOptions::default(),
    }
}

#[test]
fn a_complete_setup_validates() {
    assert_eq!(valid().validate(), Ok(()));
}

#[test]
fn an_empty_material_is_rejected() {
    let mut s = valid();
    s.material = "  ".to_string();
    assert_eq!(s.validate(), Err(SetupError::MissingMaterial));
}

#[test]
fn an_undefined_face_name_is_rejected() {
    let mut s = valid();
    s.loads[0].faces.push("lid".to_string());
    assert_eq!(
        s.validate(),
        Err(SetupError::UndefinedFace {
            used_by: "load 0".to_string(),
            name: "lid".to_string()
        })
    );
}

#[test]
fn a_duplicate_face_name_is_rejected() {
    let mut s = valid();
    s.faces.push(face("base", plane([1.0, 0.0, 0.0])));
    assert_eq!(
        s.validate(),
        Err(SetupError::DuplicateFaceName("base".to_string()))
    );
}

#[test]
fn zero_directions_are_rejected() {
    let mut s = valid();
    s.loads[0].direction = LoadDirection::Vector([0.0; 3]);
    assert_eq!(
        s.validate(),
        Err(SetupError::ZeroDirection("load 0".to_string()))
    );
    let mut s = valid();
    s.supports[0].kind = SupportKind::Roller {
        normal: [0.0, f64::NAN, 0.0],
    };
    assert_eq!(
        s.validate(),
        Err(SetupError::ZeroDirection("support 0".to_string()))
    );
    let mut s = valid();
    s.build_direction = Some([0.0; 3]);
    assert!(matches!(s.validate(), Err(SetupError::ZeroDirection(_))));
}

#[test]
fn non_positive_tolerances_are_rejected() {
    let mut s = valid();
    s.tolerance.distance = 0.0;
    assert!(matches!(
        s.validate(),
        Err(SetupError::NonPositiveTolerance { .. })
    ));
    let mut s = valid();
    s.faces[0].reference = FaceRef::Select(Selector::Plane {
        normal: [0.0, 0.0, 1.0],
        through: [0.0; 3],
        angle_tol: Some(-1.0),
        dist_tol: None,
    });
    assert_eq!(
        s.validate(),
        Err(SetupError::BadSelector {
            face: "base".to_string(),
            what: "angle tolerance".to_string(),
            value: -1.0
        })
    );
}

#[test]
fn malformed_selectors_are_rejected() {
    let cases = [
        Selector::Plane {
            normal: [0.0; 3],
            through: [0.0; 3],
            angle_tol: None,
            dist_tol: None,
        },
        Selector::Cylinder {
            axis: [0.0, 0.0, 1.0],
            radius: 0.0,
            through: None,
            angle_tol: None,
            dist_tol: None,
        },
        Selector::Sphere {
            centre: [f64::INFINITY, 0.0, 0.0],
            radius: 1.0,
            dist_tol: None,
        },
        Selector::AnyOf(Vec::new()),
    ];
    for selector in cases {
        let mut s = valid();
        s.faces[0].reference = FaceRef::Select(selector.clone());
        assert!(
            matches!(s.validate(), Err(SetupError::BadSelector { .. })),
            "{selector:?} was accepted"
        );
    }
}

#[test]
fn bad_forces_and_empty_face_lists_are_rejected() {
    let mut s = valid();
    s.loads[0].total_force_n = -5.0;
    assert!(matches!(s.validate(), Err(SetupError::BadForce { .. })));
    let mut s = valid();
    s.supports[0].faces.clear();
    assert_eq!(
        s.validate(),
        Err(SetupError::NoFaces("support 0".to_string()))
    );
    let mut s = valid();
    s.faces[0].expect = Expect::Exactly(0);
    assert_eq!(
        s.validate(),
        Err(SetupError::ExpectsNothing("base".to_string()))
    );
}

#[test]
fn out_of_range_options_are_rejected() {
    let mut s = valid();
    s.options.aggregation_threshold = 1.0;
    let Err(SetupError::Option(e)) = s.validate() else {
        panic!("aggregation_threshold = 1 was accepted");
    };
    assert_eq!(e.field, "aggregation_threshold");
    let o = AnalysisOptions {
        cg_rel_tol: 0.0,
        ..AnalysisOptions::default()
    };
    assert!(o.validate().is_err());
    let o = AnalysisOptions {
        subcell_depth: 7,
        ..AnalysisOptions::default()
    };
    assert!(o.validate().is_err());
    assert_eq!(AnalysisOptions::default().validate(), Ok(()));
}

#[test]
fn the_report_echoes_the_options_it_ran_with() {
    let options = AnalysisOptions {
        target_cell_size: Some(1.5),
        cg_rel_tol: 1.0e-6,
        ..AnalysisOptions::default()
    };
    let report = SolveReport::new(options);
    assert_eq!(report.options, options);
}

#[test]
fn an_unknown_or_excluded_material_is_rejected() {
    let mut s = valid();
    s.material = "Unobtainium".to_string();
    assert!(matches!(
        s.validate(),
        Err(SetupError::Material(
            crate::material::MaterialError::Unknown { .. }
        ))
    ));
    s.material = "TPU 95A".to_string();
    assert!(matches!(
        s.validate(),
        Err(SetupError::Material(
            crate::material::MaterialError::Excluded { .. }
        ))
    ));
}
