use popkern_analysis::thermal_setup::{
    ConvectionFaces, ConvectionPreset, HeatTransfer, ThermalSetupError,
};

use crate::setup_file::{SetupFileError, parse_setup};
use crate::thermal_file::{ThermalFileError, parse_thermal_setup, write_thermal_setup};

const EXAMPLE: &str = r#"
analysis = "thermal"
units = "mm"
material = "Aluminium 6061-T6"

[faces]
base = { plane = { normal = [0, 0, -1], through = [0, 0, 0] } }
chip = { nearest_to = [5, 5, 10] }

[[fixed_temperature]]
faces = ["base"]
temperature_c = 30

[[heat_flux]]
faces = ["chip"]
power_w = 3.5

[[convection]]
faces = "all_other"
preset = "fan"
ambient_c = 22

[[heat_source]]
power_w = 0.5
"#;

fn rejects(text: &str) -> ThermalFileError {
    match parse_thermal_setup(text) {
        Ok(_) => panic!("accepted:\n{text}"),
        Err(e) => e,
    }
}

#[test]
fn the_example_parses_field_by_field_and_round_trips() {
    let s = parse_thermal_setup(EXAMPLE).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(s.material, "Aluminium 6061-T6");
    assert_eq!(s.fixed[0].faces, ["base"]);
    assert_eq!(s.fixed[0].temperature_c, 30.0);
    assert_eq!(s.flux[0].power_w, 3.5);
    assert_eq!(s.convection[0].faces, ConvectionFaces::AllOther);
    assert_eq!(
        s.convection[0].h,
        HeatTransfer::Preset(ConvectionPreset::Fan)
    );
    assert_eq!(s.convection[0].ambient_c, 22.0);
    assert_eq!(s.sources[0].power_w, 0.5);
    assert_eq!(s.sources[0].body, None);
    let mut changed = s.clone();
    changed.sources[0].body = Some("heater".to_string());
    changed.convection[0].h = HeatTransfer::Value(12.5);
    changed.convection[0].faces = ConvectionFaces::Named(vec!["chip".to_string()]);
    changed.flux.clear();
    for setup in [s, changed] {
        let text = write_thermal_setup(&setup).unwrap_or_else(|e| panic!("{e}"));
        let back = parse_thermal_setup(&text).unwrap_or_else(|e| panic!("{e}\n{text}"));
        assert_eq!(back, setup, "{text}");
    }
}

#[test]
fn malformed_thermal_files_have_typed_errors() {
    assert_eq!(
        rejects(&EXAMPLE.replace("analysis = \"thermal\"", "")),
        ThermalFileError::NotThermal
    );
    let e = rejects(&format!(
        "{EXAMPLE}\n[[support]]\nfaces = [\"base\"]\nkind = \"fixed\"\n"
    ));
    assert!(
        matches!(
            e,
            ThermalFileError::File(SetupFileError::UnknownField { .. })
        ),
        "{e:?}"
    );
    let e = rejects(&EXAMPLE.replace("preset = \"fan\"", "preset = \"fan\"\nh = 20"));
    assert_eq!(e, ThermalFileError::CoefficientCount(0));
    let e = rejects(&EXAMPLE.replace("preset = \"fan\"", ""));
    assert_eq!(e, ThermalFileError::CoefficientCount(0));
    let e = rejects(&EXAMPLE.replace("\"fan\"", "\"hurricane\""));
    assert!(matches!(
        e,
        ThermalFileError::UnknownKeyword { what: "preset", .. }
    ));
    let e = rejects(&EXAMPLE.replace("\"all_other\"", "\"everything\""));
    assert!(matches!(
        e,
        ThermalFileError::UnknownKeyword { what: "faces", .. }
    ));
    let e = rejects(&EXAMPLE.replace("temperature_c = 30", "temperature_c = -400"));
    assert!(matches!(
        e,
        ThermalFileError::Invalid(ThermalSetupError::BadTemperature { .. })
    ));
    let flux_only = EXAMPLE
        .replace(
            "[[fixed_temperature]]\nfaces = [\"base\"]\ntemperature_c = 30\n",
            "",
        )
        .replace(
            "[[convection]]\nfaces = \"all_other\"\npreset = \"fan\"\nambient_c = 22\n",
            "",
        );
    assert_eq!(
        rejects(&flux_only),
        ThermalFileError::Invalid(ThermalSetupError::NoTemperatureReference)
    );
    // An excluded material (every preset has a cited conductivity since
    // PH-009).
    let tpu = EXAMPLE.replace("Aluminium 6061-T6", "TPU 95A");
    assert!(matches!(
        rejects(&tpu),
        ThermalFileError::Invalid(ThermalSetupError::Material(_))
    ));
}

#[test]
fn the_structural_reader_refuses_a_thermal_file() {
    let e = match parse_setup(EXAMPLE) {
        Ok(_) => panic!("structural reader accepted a thermal file"),
        Err(e) => e,
    };
    assert!(
        matches!(
            e,
            SetupFileError::UnknownField { .. } | SetupFileError::WrongAnalysis { .. }
        ),
        "{e:?}"
    );
    let minimal = "analysis = \"thermal\"\nunits = \"mm\"\nmaterial = \"ABS (FDM)\"\n";
    assert_eq!(
        parse_setup(minimal),
        Err(SetupFileError::WrongAnalysis {
            expected: "structural",
            found: "thermal".to_string()
        })
    );
}

#[test]
fn an_emittance_adds_radiation_and_round_trips() {
    use popkern_analysis::thermal_setup::Emittance;
    let with = |e: &str| {
        EXAMPLE.replace(
            "ambient_c = 22",
            &format!("ambient_c = 22\nemittance = {e}"),
        )
    };
    let s = parse_thermal_setup(&with("0.8")).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(s.convection[0].emittance, Some(Emittance::Value(0.8)));
    let m = parse_thermal_setup(&with("\"material\"")).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(m.convection[0].emittance, Some(Emittance::Material));
    let text = write_thermal_setup(&m).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(parse_thermal_setup(&text), Ok(m));
    assert!(matches!(
        rejects(&with("\"shiny\"")),
        ThermalFileError::UnknownKeyword {
            what: "emittance",
            ..
        }
    ));
    assert!(matches!(
        rejects(&with("1.5")),
        ThermalFileError::Invalid(ThermalSetupError::BadEmittance { .. })
    ));
}

#[test]
fn a_transient_table_round_trips_and_bad_values_decline() {
    let with = |t: &str| format!("{EXAMPLE}\n[transient]\n{t}\n");
    let s = parse_thermal_setup(&with("initial_c = 20\nstep_s = 0.5\nduration_s = 30"))
        .unwrap_or_else(|e| panic!("{e}"));
    let t = s.transient.unwrap_or_else(|| panic!("no transient"));
    assert_eq!(
        (t.initial_c, t.step_s, t.duration_s, t.steps()),
        (20.0, 0.5, 30.0, 60)
    );
    let text = write_thermal_setup(&s).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(parse_thermal_setup(&text), Ok(s));
    assert!(matches!(
        rejects(&with("initial_c = 20\nstep_s = 0\nduration_s = 30")),
        ThermalFileError::Invalid(ThermalSetupError::BadTransient { .. })
    ));
    assert!(matches!(
        rejects(&with("initial_c = 20\nstep_s = 1e-6\nduration_s = 30")),
        ThermalFileError::Invalid(ThermalSetupError::BadTransient { .. })
    ));
}
