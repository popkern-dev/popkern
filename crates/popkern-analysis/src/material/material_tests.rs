use super::*;

fn get(name: &str) -> &'static Material {
    preset(name).unwrap_or_else(|e| panic!("{e}"))
}

/// name, E (flexural), stiff bound (tensile), ν, σ_inplane, behaviour,
/// σ_layer, HDT, density: notes 11 and 13, and the Ultimaker TDS.
#[allow(clippy::type_complexity)]
const EXPECTED: [(
    &str,
    f64,
    Option<f64>,
    f64,
    f64,
    Behaviour,
    Option<f64>,
    Option<f64>,
    Option<f64>,
); 6] = [
    (
        "PLA (FDM)",
        3019.0,
        Some(3250.0),
        0.36,
        52.5,
        Behaviour::Brittle,
        Some(33.1),
        Some(58.8),
        Some(1240.0),
    ),
    (
        "ABS (FDM)",
        1430.0,
        Some(1962.0),
        0.35,
        38.1,
        Behaviour::Yielding,
        Some(19.0),
        Some(86.6),
        Some(1100.0),
    ),
    (
        "PETG (FDM)",
        1882.0,
        Some(1939.0),
        0.38,
        46.2,
        Behaviour::Brittle,
        Some(19.0),
        Some(76.2),
        Some(1270.0),
    ),
    (
        "PC (FDM)",
        1619.0,
        Some(2394.0),
        0.37,
        53.3,
        Behaviour::Yielding,
        Some(25.2),
        Some(104.5),
        Some(1200.0),
    ),
    (
        "Aluminium 6061-T6",
        69_000.0,
        None,
        0.33,
        276.0,
        Behaviour::Yielding,
        None,
        None,
        Some(2700.0),
    ),
    (
        "Steel S235",
        200_000.0,
        None,
        0.30,
        235.0,
        Behaviour::Yielding,
        None,
        None,
        Some(7830.0),
    ),
];

#[test]
fn presets_match_the_notes_exactly() {
    assert_eq!(PRESETS.len(), EXPECTED.len());
    for (m, e) in PRESETS.iter().zip(EXPECTED) {
        let (name, modulus, stiff, nu, inplane, behaviour, layer, hdt, density) = e;
        assert_eq!(m.name, name);
        assert_eq!(m.youngs_modulus.value(), Some(modulus), "{name} E");
        assert_eq!(
            m.youngs_modulus_stiff_bound.value(),
            stiff,
            "{name} stiff bound"
        );
        assert_eq!(m.poisson_ratio.value(), Some(nu), "{name} nu");
        assert_eq!(m.strength.inplane.value(), Some(inplane), "{name} inplane");
        assert_eq!(m.strength.inplane_behaviour, behaviour, "{name} behaviour");
        assert_eq!(m.strength.interlayer.value(), layer, "{name} interlayer");
        assert_eq!(m.heat_deflection_temp.value(), hdt, "{name} HDT");
        assert_eq!(m.density.value(), density, "{name} density");
        assert_eq!(get(name), m);
    }
}

#[test]
fn every_value_carries_a_source() {
    for m in PRESETS.iter() {
        let all = [
            m.youngs_modulus,
            m.youngs_modulus_stiff_bound,
            m.poisson_ratio,
            m.strength.inplane,
            m.strength.interlayer,
            m.strength.interlayer_shear,
            m.strength.ultimate_compressive,
            m.strength.interlayer_shear_scatter,
            m.strength.ultimate_compressive_scatter,
            m.density,
            m.conductivity,
            m.specific_heat,
            m.emittance_default,
            m.heat_deflection_temp,
            m.thermal_expansion,
        ];
        for p in all {
            match p {
                Property::NotApplicable => assert_eq!(p.source(), None),
                _ => assert!(
                    p.source().is_some_and(|s| !s.trim().is_empty()),
                    "{}: {p:?} has no source",
                    m.name
                ),
            }
        }
    }
}

#[test]
fn polymer_e_is_flexural_with_tensile_as_the_stiff_bound() {
    // XY flexural / tensile from cross-cutting Thread 6.
    let ratios = [
        ("PLA (FDM)", 0.93),
        ("ABS (FDM)", 0.73),
        ("PETG (FDM)", 0.97),
        ("PC (FDM)", 0.68),
    ];
    for (name, ratio) in ratios {
        let m = get(name);
        let flexural = m.youngs_modulus.value().unwrap_or(f64::NAN);
        let tensile = m.youngs_modulus_stiff_bound.value().unwrap_or(f64::NAN);
        assert!(
            (flexural / tensile - ratio).abs() < 0.005,
            "{name}: {flexural}/{tensile}"
        );
        assert!(
            m.youngs_modulus
                .source()
                .is_some_and(|s| s.contains("flexural")),
            "{name}"
        );
        assert!(
            m.youngs_modulus_stiff_bound
                .source()
                .is_some_and(|s| s.contains("tensile")),
            "{name}"
        );
    }
}

#[test]
fn every_preset_validates() {
    for m in PRESETS.iter() {
        assert_eq!(m.validate(), Ok(()), "{}", m.name);
    }
}

#[test]
fn out_of_range_properties_are_rejected() {
    let base = *get("PLA (FDM)");
    let typical = |value| Property::Typical {
        value,
        note: "test",
    };
    let mut m = base;
    m.poisson_ratio = typical(0.5);
    assert!(matches!(
        m.validate(),
        Err(MaterialError::OutOfRange {
            property: "poisson_ratio",
            ..
        })
    ));
    let mut m = base;
    m.youngs_modulus = typical(0.0);
    assert!(matches!(
        m.validate(),
        Err(MaterialError::OutOfRange {
            property: "youngs_modulus",
            ..
        })
    ));
    let mut m = base;
    m.strength.inplane = typical(f64::NAN);
    assert!(matches!(
        m.validate(),
        Err(MaterialError::OutOfRange {
            property: "strength.inplane",
            ..
        })
    ));
    let mut m = base;
    m.emittance_default = typical(1.2);
    assert!(matches!(
        m.validate(),
        Err(MaterialError::OutOfRange { .. })
    ));
    let mut m = base;
    m.youngs_modulus = Property::Unsourced { note: "test" };
    assert_eq!(
        m.validate(),
        Err(MaterialError::MissingProperty {
            material: "PLA (FDM)",
            property: "youngs_modulus"
        })
    );
}

/// (name, conductivity W/(m·K), c_p J/(kg·K), CTE 1/K): every thermal
/// value as notes 11 and 13 transcribe it; `None` is unsourced (PH-009).
type ThermalRow = (&'static str, Option<f64>, Option<f64>, Option<f64>);

const THERMAL: [ThermalRow; 6] = [
    ("PLA (FDM)", Some(0.182), None, Some(96.0e-6)),
    ("ABS (FDM)", Some(0.14), None, Some(86.0e-6)),
    ("PETG (FDM)", Some(0.198), None, None),
    ("PC (FDM)", Some(0.29), Some(1250.0), None),
    ("Aluminium 6061-T6", Some(167.0), Some(896.0), Some(23.0e-6)),
    ("Steel S235", Some(64.0), Some(434.0), Some(12.0e-6)),
];

/// The thermal values match the notes exactly, and the unsourced ones
/// stay empty: a value appearing without a source fails here (and in
/// `every_value_carries_a_source`).
#[test]
fn unsourced_values_stay_empty() {
    for (m, (name, k, cp, cte)) in PRESETS.iter().zip(THERMAL) {
        assert_eq!(m.name, name);
        assert_eq!(m.conductivity.value(), k, "{name} k");
        assert_eq!(m.specific_heat.value(), cp, "{name} c_p");
        assert_eq!(m.thermal_expansion.value(), cte, "{name} CTE");
        for (p, v) in [
            (m.conductivity, k),
            (m.specific_heat, cp),
            (m.thermal_expansion, cte),
        ] {
            assert_eq!(
                v.is_none(),
                matches!(p, Property::Unsourced { .. }),
                "{name}"
            );
        }
    }
}

/// (name, S_layer, its sd, σ_uc, its sd), MPa: notes 11 (PH-080,
/// PH-084); `None` is unsourced or not applicable.
type StrengthRow = (
    &'static str,
    Option<f64>,
    Option<f64>,
    Option<f64>,
    Option<f64>,
);

// Kholil's PLA σ_c scatter is 3.14 MPa, a measurement, not π.
#[allow(clippy::approx_constant)]
const SHEAR_AND_COMPRESSION: [StrengthRow; 6] = [
    (
        "PLA (FDM)",
        Some(25.74),
        Some(1.03),
        Some(66.78),
        Some(3.14),
    ),
    (
        "ABS (FDM)",
        Some(23.01),
        Some(1.63),
        Some(36.38),
        Some(2.25),
    ),
    ("PETG (FDM)", None, None, Some(35.02), Some(2.35)),
    ("PC (FDM)", Some(30.16), Some(1.64), None, None),
    ("Aluminium 6061-T6", None, None, None, None),
    ("Steel S235", None, None, None, None),
];

/// S_layer and σ_uc match the notes exactly, a strength and its scatter
/// come together, and only printed polymers have an S_layer.
#[test]
fn shear_and_compressive_strengths_match_the_notes() {
    for (m, (name, shear, shear_sd, compressive, compressive_sd)) in
        PRESETS.iter().zip(SHEAR_AND_COMPRESSION)
    {
        let st = &m.strength;
        assert_eq!(m.name, name);
        assert_eq!(st.interlayer_shear.value(), shear, "{name} S_layer");
        assert_eq!(st.interlayer_shear_scatter.value(), shear_sd, "{name}");
        assert_eq!(st.ultimate_compressive.value(), compressive, "{name}");
        assert_eq!(
            st.ultimate_compressive_scatter.value(),
            compressive_sd,
            "{name}"
        );
        assert_eq!(
            st.interlayer.value().is_none(),
            st.interlayer_shear == Property::NotApplicable,
            "{name}"
        );
    }
}

#[test]
fn excluded_and_unknown_materials_are_declined() {
    assert!(matches!(
        preset("TPU 95A"),
        Err(MaterialError::Excluded { reason, .. }) if reason.contains("hyperelastic")
    ));
    assert!(matches!(
        preset("Nylon (FDM)"),
        Err(MaterialError::Excluded { .. })
    ));
    let Err(MaterialError::Unknown { known, .. }) = preset("Unobtainium") else {
        panic!("an unknown name resolved");
    };
    assert_eq!(known.len(), PRESETS.len());
    // exact names only: no case folding, no partial match
    assert!(preset("pla (fdm)").is_err());
    assert!(preset("PLA").is_err());
}
