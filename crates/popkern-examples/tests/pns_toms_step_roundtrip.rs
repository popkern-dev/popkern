#![cfg(toms_package)]

use popkern_core::geometry::polyhedral_spline::{
    ControlNet, PolyhedralCoverageOptions, build_control_net_surface_complex,
    weld_clean_surface_complex_boundaries,
};
use popkern_examples::step_validation::{
    StepRoundTripValidationOptions, export_and_validate_breps,
};

struct Fixture {
    stem: &'static str,
    obj: &'static str,
}

const SMOKE_FIXTURES: &[Fixture] = &[
    Fixture {
        stem: "t0",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T0.obj"),
    },
    Fixture {
        stem: "tee",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/tee.obj"),
    },
    Fixture {
        stem: "cube",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/cube.obj"),
    },
    Fixture {
        stem: "t1",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T1.obj"),
    },
    Fixture {
        stem: "t2",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T2.obj"),
    },
    Fixture {
        stem: "polar5sct",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/polar5sct.obj"),
    },
    Fixture {
        stem: "eop5sct",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/eop5sct.obj"),
    },
    Fixture {
        stem: "ngon5",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/ngon5.obj"),
    },
    Fixture {
        stem: "t012",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/T012.obj"),
    },
    Fixture {
        stem: "cylinder_partial",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/cylinder.obj"),
    },
];

const HEAVY_FIXTURES: &[Fixture] = &[
    Fixture {
        stem: "suzanne_all_configurations",
        obj: include_str!(
            "../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/suzanne_all_configurations.obj"
        ),
    },
    Fixture {
        stem: "hand_partial",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/hand.obj"),
    },
    Fixture {
        stem: "airplane_partial",
        obj: include_str!("../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/testfile/airplane.obj"),
    },
];

#[test]
fn pns_toms_step_outputs_round_trip_through_generic_import() {
    run_roundtrip_fixtures(SMOKE_FIXTURES);
}

#[test]
#[ignore = "heavy: 3 large PnS meshes (suzanne 1075 / hand 817 / airplane 1800 OBJ lines); nightly tier, PK-190 owns the CI wiring — the 10 smaller fixtures run by default above"]
fn pns_toms_heavy_step_outputs_round_trip_through_generic_import() {
    run_roundtrip_fixtures(HEAVY_FIXTURES);
}

fn run_roundtrip_fixtures(fixtures: &[Fixture]) {
    let options = PolyhedralCoverageOptions::default();
    let validation = StepRoundTripValidationOptions::default();
    for fixture in fixtures {
        let net = ControlNet::from_obj_str(fixture.obj).unwrap();
        let complex = build_control_net_surface_complex(&net, options).unwrap();
        let breps: Vec<_> = complex
            .components
            .iter()
            .map(|component| component.brep.clone())
            .collect();
        assert!(
            !breps.is_empty(),
            "{} must produce at least one BRep",
            fixture.stem
        );
        export_and_validate_breps(fixture.stem, &breps, &validation).unwrap();

        let welded = weld_clean_surface_complex_boundaries(&complex, options)
            .unwrap_or_else(|e| panic!("{} weld failed: {e}", fixture.stem));
        if welded.report.welded_pairs > 0 {
            export_and_validate_breps(
                &format!("{}_welded", fixture.stem),
                &[welded.brep],
                &validation,
            )
            .unwrap();
        }
    }
}
