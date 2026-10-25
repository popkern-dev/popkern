use crate::face_ref::{Expect, FaceRef, Selector, Tolerance};
use crate::material::MaterialError;
use crate::options::AnalysisOptions;
use crate::setup::{FaceDef, SetupError};
use crate::thermal_setup::*;

fn face(name: &str) -> FaceDef {
    FaceDef {
        name: name.to_string(),
        reference: FaceRef::Select(Selector::NearestTo([0.0, 0.0, 0.0])),
        expect: Expect::default(),
    }
}

fn valid() -> ThermalSetup {
    ThermalSetup {
        material: "ABS (FDM)".to_string(),
        tolerance: Tolerance::default(),
        faces: vec![face("base"), face("chip")],
        fixed: vec![FixedTemperature {
            faces: vec!["base".to_string()],
            temperature_c: 25.0,
        }],
        flux: vec![HeatFlux {
            faces: vec!["chip".to_string()],
            power_w: 2.0,
        }],
        convection: vec![Convection {
            faces: ConvectionFaces::AllOther,
            h: HeatTransfer::Preset(ConvectionPreset::StillAirMatte),
            ambient_c: 22.0,
            emittance: None,
        }],
        sources: Vec::new(),
        transient: None,
        options: AnalysisOptions::default(),
    }
}

#[test]
fn a_valid_setup_returns_the_conductivity() {
    assert_eq!(valid().validate(), Ok(0.14));
}

#[test]
fn a_flux_only_setup_has_no_temperature_reference() {
    let mut s = valid();
    s.fixed.clear();
    s.convection.clear();
    assert_eq!(s.validate(), Err(ThermalSetupError::NoTemperatureReference));
    // Convection alone is a reference.
    let mut s = valid();
    s.fixed.clear();
    assert!(s.validate().is_ok());
}

/// PH-009: PLA's and PETG's conductivity are cited now (printed, 100%
/// infill; Tychanicz-Kwiecień 2025), so their thermal setups validate.
#[test]
fn printed_pla_and_petg_conduct_with_their_cited_values() {
    for (name, k) in [("PLA (FDM)", 0.182), ("PETG (FDM)", 0.198)] {
        let mut s = valid();
        s.material = name.to_string();
        assert_eq!(s.validate(), Ok(k), "{name}");
    }
}

#[test]
fn malformed_entries_have_typed_errors() {
    let mut s = valid();
    s.flux[0].faces = vec!["base".to_string()];
    assert_eq!(
        s.validate(),
        Err(ThermalSetupError::ConflictingBoundary("base".to_string()))
    );
    let mut s = valid();
    s.fixed[0].temperature_c = -300.0;
    assert!(matches!(
        s.validate(),
        Err(ThermalSetupError::BadTemperature { .. })
    ));
    let mut s = valid();
    s.convection[0].h = HeatTransfer::Value(0.0);
    assert!(matches!(
        s.validate(),
        Err(ThermalSetupError::BadCoefficient { .. })
    ));
    let mut s = valid();
    s.flux[0].power_w = f64::NAN;
    assert!(matches!(
        s.validate(),
        Err(ThermalSetupError::BadPower { .. })
    ));
    let mut s = valid();
    s.convection.push(s.convection[0].clone());
    assert!(matches!(
        s.validate(),
        Err(ThermalSetupError::ConflictingBoundary(_))
    ));
    let mut s = valid();
    s.fixed[0].faces = vec!["lid".to_string()];
    assert!(matches!(
        s.validate(),
        Err(ThermalSetupError::Setup(SetupError::UndefinedFace { .. }))
    ));
    let mut s = valid();
    s.material = "Unobtainium".to_string();
    assert!(matches!(
        s.validate(),
        Err(ThermalSetupError::Material(MaterialError::Unknown { .. }))
    ));
}

#[test]
fn presets_and_the_radiation_coefficient() {
    let h = |p| HeatTransfer::Preset(p).value();
    assert_eq!(h(ConvectionPreset::StillAirMatte), 10.0);
    assert_eq!(h(ConvectionPreset::StillAirBrightMetal), 5.0);
    assert_eq!(h(ConvectionPreset::Fan), 25.0);
    assert_eq!(h(ConvectionPreset::StrongAirflow), 50.0);
    // T6: 4σT³ε = 6.1 W/m²K at 300 K, ε = 1 (Lienhard eqn 2.31).
    let h_rad = radiation_coefficient(300.0, 1.0);
    assert!((h_rad - 6.124).abs() < 1e-3, "{h_rad}");
}
