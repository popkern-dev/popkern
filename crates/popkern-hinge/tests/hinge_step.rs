//! PK-695: hinge bodies exported to STEP match OpenCASCADE's volumes.
//! The OCC numbers were measured with `tools/occ-oracle/occ_solids.py`
//! on the files this test writes (2026-10-07); the test pins popkern's
//! own volumes to them so a construction change that moves the
//! geometry fails here.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fastener::ThreadSpec;
use popkern_hinge::{HingeComponent, HingeHoles, HingeKind, HingeParams, hinge_component};

/// popkern's fine-tessellation volume against OCC's exact one.
const VOLUME_REL_TOL: f64 = 1.0e-3;

fn step(name: &str, bodies: &[BRep]) -> std::path::PathBuf {
    let doc = popkern_io::StepDocument {
        products: vec![popkern_io::StepProductNode {
            name: name.to_string(),
            id: name.to_uppercase(),
            children: vec![],
            shapes: vec![std::sync::Arc::new(popkern_io::ShapeRep {
                kind: popkern_io::ShapeRepKind::AdvancedBrep,
                solids: bodies.to_vec(),
                context: popkern_io::RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: popkern_io::ProductMetadata::default(),
        }],
        roots: vec![0],
        units: popkern_io::UnitContext {
            context: popkern_io::RepresentationContext::default(),
        },
        opaque_bag: popkern_io::OpaqueBag::default(),
    };
    let bytes =
        popkern_io::export_step(&doc, &popkern_io::StepExportOptions::default()).expect("step");
    let path = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("{name}.step"));
    std::fs::write(&path, bytes).expect("write");
    path
}

fn bodies(kind: HingeKind, p: &HingeParams) -> Vec<BRep> {
    let mut components = vec![HingeComponent::LeafA, HingeComponent::LeafB];
    if kind == HingeKind::Butt {
        components.push(HingeComponent::Pin);
    }
    components
        .into_iter()
        .map(|c| hinge_component(kind, c, p, Tolerance::default()).expect("component"))
        .collect()
}

fn check(name: &str, kind: HingeKind, p: &HingeParams, occ: &[f64]) {
    let built = bodies(kind, p);
    let path = step(name, &built);
    for (i, body) in built.iter().enumerate() {
        let v = popkern_measure::body_volume(body, &Tolerance::default()).expect("volume");
        eprintln!("{name} body {i}: popkern {v:.6} ({})", path.display());
        if let Some(&o) = occ.get(i) {
            let rel = (v - o).abs() / o;
            assert!(
                rel < VOLUME_REL_TOL,
                "{name} body {i}: {v} vs OCC {o} ({rel:.2e})"
            );
        }
    }
    assert_eq!(occ.len(), built.len(), "{name}: an OCC volume per body");
}

fn holed() -> HingeParams {
    HingeParams {
        holes: Some(HingeHoles {
            per_leaf: 2,
            thread: ThreadSpec::parse("M3").expect("M3"),
            countersink: true,
        }),
        ..HingeParams::new(60.0)
    }
}

#[test]
fn default_butt_hinge_matches_occ() {
    check(
        "hinge_butt_default",
        HingeKind::Butt,
        &HingeParams::new(60.0),
        &OCC_DEFAULT,
    );
}

#[test]
fn nine_knuckle_clearance_hinge_matches_occ() {
    let p = HingeParams {
        knuckles: 9,
        clearance: 0.3,
        ..HingeParams::new(60.0)
    };
    check("hinge_butt_nine", HingeKind::Butt, &p, &OCC_NINE);
}

#[test]
fn countersunk_hinge_matches_occ() {
    check("hinge_butt_holes", HingeKind::Butt, &holed(), &OCC_HOLES);
}

/// OCC volumes (mm³) per body: leaf A, leaf B, pin.
const OCC_DEFAULT: [f64; 3] = [2_773.307_641_854, 2_288.631_761_235, 753.982_236_862];
const OCC_NINE: [f64; 3] = [2_735.730_494_520, 2_438.184_395_615, 753.982_236_862];
const OCC_HOLES: [f64; 3] = [2_698.502_393_794, 2_213.826_513_174, 753.982_236_862];
