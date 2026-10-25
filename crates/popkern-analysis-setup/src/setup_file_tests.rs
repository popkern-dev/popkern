use popkern_analysis::{
    AnalysisOptions, Axis, Expect, FaceRef, LoadDirection, Selector, SetupError, SupportKind,
    Tolerance,
};

use crate::setup_file::{SetupFileError, parse_setup, write_setup};

/// The example from docs/architecture/analysis-library-boundary.md,
/// verbatim.
const BOUNDARY_DOC_EXAMPLE: &str = r#"
units = "mm"
material = "PLA (FDM)"
build_direction = [0, 0, 1]         # optional

[tolerance]
angle_deg = 0.5
distance = 0.01

[faces]
base   = { plane = { normal = [0, 0, -1], through = [0, 0, 0] } }
bolt_a = { cylinder = { axis = [0, 0, 1], radius = 2.5, through = [10, 5, 0] } }
hook   = { nearest_to = [40, 0, 12] }

[[support]]
faces = ["base"]
kind = "fixed"

[[load]]
faces = ["hook"]
force_n = 50
direction = [0, 0, -1]

[options]
max_active_cells = 30000
"#;

#[test]
fn the_boundary_doc_example_parses_field_by_field() {
    let s = match parse_setup(BOUNDARY_DOC_EXAMPLE) {
        Ok(s) => s,
        Err(e) => panic!("example rejected: {e}"),
    };
    assert_eq!(s.material, "PLA (FDM)");
    assert_eq!(s.build_direction, Some([0.0, 0.0, 1.0]));
    assert_eq!(
        s.tolerance,
        Tolerance {
            angle_deg: 0.5,
            distance: 0.01
        }
    );
    let names: Vec<&str> = s.faces.iter().map(|f| f.name.as_str()).collect();
    assert_eq!(names, ["base", "bolt_a", "hook"]);
    assert_eq!(
        s.faces[0].reference,
        FaceRef::Select(Selector::Plane {
            normal: [0.0, 0.0, -1.0],
            through: [0.0, 0.0, 0.0],
            angle_tol: None,
            dist_tol: None
        })
    );
    assert_eq!(
        s.faces[1].reference,
        FaceRef::Select(Selector::Cylinder {
            axis: [0.0, 0.0, 1.0],
            radius: 2.5,
            through: Some([10.0, 5.0, 0.0]),
            angle_tol: None,
            dist_tol: None
        })
    );
    assert_eq!(
        s.faces[2].reference,
        FaceRef::Select(Selector::NearestTo([40.0, 0.0, 12.0]))
    );
    assert!(s.faces.iter().all(|f| f.expect == Expect::Exactly(1)));
    assert_eq!(s.supports.len(), 1);
    assert_eq!(s.supports[0].faces, ["base"]);
    assert_eq!(s.supports[0].kind, SupportKind::Fixed);
    assert_eq!(s.loads.len(), 1);
    assert_eq!(s.loads[0].faces, ["hook"]);
    assert_eq!(s.loads[0].total_force_n, 50.0);
    assert_eq!(
        s.loads[0].direction,
        LoadDirection::Vector([0.0, 0.0, -1.0])
    );
    assert_eq!(s.options, AnalysisOptions::default());
}

#[test]
fn a_written_setup_parses_back_equal() {
    let mut text = BOUNDARY_DOC_EXAMPLE.to_string();
    text.push_str("cg_rel_tol = 1e-6\n");
    text = text.replace(
        "hook   = { nearest_to = [40, 0, 12] }",
        "hook   = { nearest_to = [40, 0, 12] }\n\
         lugs = { any_of = [{ sphere = { centre = [1, 2, 3], radius = 4 } }, { name = \"x\" }], expect = { at_least = 2 } }",
    );
    // `name` inside a selector list is not a selector key.
    assert!(matches!(
        parse_setup(&text),
        Err(SetupFileError::UnknownField { .. })
    ));
    let text = text.replace("{ name = \"x\" }", "{ nearest_to = [0, 0, 0] }");
    let mut first = parse_setup(&text).unwrap_or_else(|e| panic!("{e}"));
    first.loads[0].direction = LoadDirection::Axis(Axis::MinusZ);
    first.supports[0].kind = SupportKind::Roller {
        normal: [1.0, 0.0, 0.0],
    };
    let written = write_setup(&first).unwrap_or_else(|e| panic!("{e}"));
    let second = parse_setup(&written).unwrap_or_else(|e| panic!("{e}\n{written}"));
    assert_eq!(first, second, "round trip changed the setup:\n{written}");
    assert_eq!(second.options.cg_rel_tol, 1e-6);
    assert_eq!(second.faces[3].name, "lugs");
    assert_eq!(second.faces[3].expect, Expect::AtLeast(2));
    for kind in [
        SupportKind::SurfaceRoller,
        SupportKind::DirectionalRoller {
            direction: [0.5, 0.25, 0.0],
        },
    ] {
        first.supports[0].kind = kind;
        let written = write_setup(&first).unwrap_or_else(|e| panic!("{e}"));
        let again = parse_setup(&written).unwrap_or_else(|e| panic!("{e}\n{written}"));
        assert_eq!(first, again, "round trip changed {kind:?}:\n{written}");
    }
}

fn rejects(text: &str) -> SetupFileError {
    match parse_setup(text) {
        Ok(_) => panic!("accepted:\n{text}"),
        Err(e) => e,
    }
}

#[test]
fn an_unknown_key_is_its_own_error() {
    let e = rejects(&BOUNDARY_DOC_EXAMPLE.replace("force_n = 50", "force_n = 50\nforce_lbf = 11"));
    assert!(matches!(e, SetupFileError::UnknownField { .. }), "{e:?}");
    let e = rejects(&BOUNDARY_DOC_EXAMPLE.replace("max_active_cells", "max_cells"));
    assert!(matches!(e, SetupFileError::UnknownField { .. }), "{e:?}");
}

#[test]
fn a_missing_material_is_its_own_error() {
    let e = rejects(&BOUNDARY_DOC_EXAMPLE.replace("material = \"PLA (FDM)\"", ""));
    assert_eq!(
        e,
        SetupFileError::MissingField {
            field: "material".to_string()
        }
    );
    let e = rejects(&BOUNDARY_DOC_EXAMPLE.replace("force_n = 50", ""));
    assert_eq!(
        e,
        SetupFileError::MissingField {
            field: "force_n".to_string()
        }
    );
}

#[test]
fn an_undefined_face_name_is_its_own_error() {
    let e = rejects(&BOUNDARY_DOC_EXAMPLE.replace("faces = [\"hook\"]", "faces = [\"hoook\"]"));
    assert_eq!(
        e,
        SetupFileError::Invalid(SetupError::UndefinedFace {
            used_by: "load 0".to_string(),
            name: "hoook".to_string()
        })
    );
}

#[test]
fn a_zero_direction_is_its_own_error() {
    let e =
        rejects(&BOUNDARY_DOC_EXAMPLE.replace("direction = [0, 0, -1]", "direction = [0, 0, 0]"));
    assert_eq!(
        e,
        SetupFileError::Invalid(SetupError::ZeroDirection("load 0".to_string()))
    );
}

#[test]
fn a_non_positive_tolerance_is_its_own_error() {
    let e = rejects(&BOUNDARY_DOC_EXAMPLE.replace("distance = 0.01", "distance = -0.01"));
    assert!(
        matches!(
            e,
            SetupFileError::Invalid(SetupError::NonPositiveTolerance { .. })
        ),
        "{e:?}"
    );
}

#[test]
fn other_file_mistakes_have_typed_errors() {
    let e = rejects(&BOUNDARY_DOC_EXAMPLE.replace("units = \"mm\"", "units = \"in\""));
    assert_eq!(e, SetupFileError::UnsupportedUnits("in".to_string()));
    let e = rejects(&BOUNDARY_DOC_EXAMPLE.replace(
        "hook   = { nearest_to = [40, 0, 12] }",
        "hook   = { nearest_to = [40, 0, 12], id = 3 }",
    ));
    assert_eq!(
        e,
        SetupFileError::FaceReferenceCount {
            face: "hook".to_string(),
            found: 2
        }
    );
    let e = rejects(&BOUNDARY_DOC_EXAMPLE.replace(
        "kind = \"fixed\"",
        "kind = \"roller\"\nnormal = [0, 0, -1]\ndirection = [1, 0, 0]",
    ));
    assert_eq!(e, SetupFileError::RollerNormalAndDirection(0));
    let e = rejects(&BOUNDARY_DOC_EXAMPLE.replace(
        "kind = \"fixed\"",
        "kind = \"fixed\"\ndirection = [1, 0, 0]",
    ));
    assert_eq!(e, SetupFileError::NormalOnFixed(0));
    let e =
        rejects(&BOUNDARY_DOC_EXAMPLE.replace("direction = [0, 0, -1]", "direction = \"down\""));
    assert!(
        matches!(e, SetupFileError::UnknownDirection { load: 0, .. }),
        "{e:?}"
    );
    let e = rejects("units = \"mm\"\nmaterial = ");
    assert!(matches!(e, SetupFileError::Syntax { .. }), "{e:?}");
}

#[test]
fn a_thermal_table_round_trips_and_needs_one_source() {
    use popkern_analysis::setup::TemperatureSource;
    let base = "units = \"mm\"\nmaterial = \"Aluminium 6061-T6\"\n[faces]\nbase = { plane = { normal = [0, 0, -1], through = [0, 0, 0] } }\n[[support]]\nfaces = [\"base\"]\nkind = \"fixed\"\n";
    let with = |t: &str| format!("{base}[thermal]\nreference_c = 20\n{t}\n");
    let s = parse_setup(&with("uniform_c = 80")).unwrap_or_else(|e| panic!("{e}"));
    let t = s.thermal.clone().unwrap_or_else(|| panic!("no thermal"));
    assert_eq!(t.temperature, TemperatureSource::Uniform(80.0));
    let text = write_setup(&s).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(parse_setup(&text), Ok(s));
    let f = parse_setup(&with("setup = \"heat.toml\"")).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(
        f.thermal.map(|t| t.temperature),
        Some(TemperatureSource::ThermalSetup("heat.toml".to_string()))
    );
    for bad in ["", "uniform_c = 80\nsetup = \"heat.toml\""] {
        assert_eq!(
            parse_setup(&with(bad)),
            Err(SetupFileError::ThermalSourceCount)
        );
    }
}
