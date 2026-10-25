//! Assemblies and multi-solid bodies (PH-010): placements applied per
//! instance, per-solid exact volumes, and non-rigid placements declined.

use std::f64::consts::PI;
use std::path::Path;
use std::sync::Arc;

use popkern_analysis_cli::step_model::{StepModelError, breps_from_document, model_from_step_file};
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_io::{
    AssemblyChild, OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind,
    StepDocument, StepProductNode, UnitContext,
};

fn fixture(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join(format!("../../tests/fixtures/analysis/{name}.step"))
}

fn bounds(v: &[[f64; 3]]) -> ([f64; 3], [f64; 3]) {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in v {
        for a in 0..3 {
            lo[a] = lo[a].min(p[a]);
            hi[a] = hi[a].max(p[a]);
        }
    }
    (lo, hi)
}

#[test]
fn both_instances_are_placed_and_measured() {
    let m =
        model_from_step_file(&fixture("assembly_two_dimples")).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(m.solids.len(), 2);
    let names: Vec<_> = m
        .solids
        .iter()
        .map(|s| s.name.clone().unwrap_or_default())
        .collect();
    assert_eq!(
        names,
        ["assembly_two_dimples/left", "assembly_two_dimples/right"]
    );
    // Dimple block: 1000 − spherical cap (h 2, R 3); curved area 2πRh.
    let exact = 1000.0 - PI * 4.0 * 7.0 / 3.0;
    let bound = m.chord_tolerance * 2.0 * PI * 3.0 * 2.0;
    for s in &m.solids {
        assert!((s.enclosed_volume() - exact).abs() <= bound, "{:?}", s.name);
        assert!(
            (s.exact_volume.unwrap_or(f64::NAN) - exact).abs() <= bound,
            "{:?}",
            s.name
        );
    }
    let (lo, hi) = bounds(&m.solids[0].vertices);
    assert!(lo.iter().all(|x| x.abs() < 1e-9) && hi.iter().all(|x| (x - 10.0).abs() < 1e-9));
    // Rotated 90° about z, then moved to (30, 5, 0): x 20..30, y 5..15.
    let (lo, hi) = bounds(&m.solids[1].vertices);
    let want_lo = [20.0, 5.0, 0.0];
    let want_hi = [30.0, 15.0, 10.0];
    for a in 0..3 {
        assert!(
            (lo[a] - want_lo[a]).abs() < 1e-9 && (hi[a] - want_hi[a]).abs() < 1e-9,
            "{lo:?} {hi:?}"
        );
    }
}

#[test]
fn each_solid_of_one_body_has_its_own_exact_volume() {
    let m = model_from_step_file(&fixture("two_boxes_one_body")).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(m.solids.len(), 2);
    let mut volumes: Vec<f64> = m
        .solids
        .iter()
        .map(|s| s.exact_volume.unwrap_or(f64::NAN))
        .collect();
    volumes.sort_by(f64::total_cmp);
    assert!((volumes[0] - 8.0).abs() <= 1e-12 * 8.0, "{volumes:?}");
    assert!((volumes[1] - 1000.0).abs() <= 1e-12 * 1000.0, "{volumes:?}");
    for s in &m.solids {
        let v = s.exact_volume.unwrap_or(f64::NAN);
        assert!((s.enclosed_volume() - v).abs() <= 1e-9 * v);
    }
}

#[test]
fn a_non_rigid_placement_is_declined_by_name() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap_or_else(|e| panic!("{e}"));
    let node = |name: &str, children, solids: Vec<BRep>| StepProductNode {
        name: name.to_string(),
        id: name.to_string(),
        children,
        shapes: if solids.is_empty() {
            Vec::new()
        } else {
            vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids,
                context: RepresentationContext::default(),
                non_brep_items: Vec::new(),
            })]
        },
        metadata: ProductMetadata::default(),
    };
    let child = |t| AssemblyChild {
        product_index: 1,
        transform: t,
        instance_name: Some("stretched".to_string()),
    };
    let doc = |t| StepDocument {
        products: vec![
            node("top", vec![child(t)], Vec::new()),
            node("cube", Vec::new(), vec![cube.clone()]),
        ],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };
    match breps_from_document(&doc(Transform3::scale(1.0, 2.0, 1.0))) {
        Err(StepModelError::NonRigidPlacement { child }) => assert_eq!(child, "top/stretched"),
        other => panic!("expected NonRigidPlacement, got {other:?}"),
    }
    // Control: a rigid placement is accepted and applied.
    let bodies = breps_from_document(&doc(Transform3::translation(5.0, 0.0, 0.0)))
        .unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(bodies.len(), 1);
    assert_eq!(bodies[0].name, "top/stretched");
}
