//! The propeller's volume against OpenCASCADE (PK-818, PK-872).
//!
//! The quadrature integrates each mesh triangle on the exact surface and
//! each boundary chord's sliver up to its trim. The propeller's hub
//! holes are bounded by blade-root edges whose fitted pcurves sit up to
//! ~1e-4 off the edges (PK-826): the sliver test located the trim on
//! the pcurve's distance and declined every hub chord, which left a
//! mesh-density bias (+4e-4 at the default chord, converging toward the
//! reference as the mesh refined). Chord ends now locate on the 3-D
//! edge and slivers end on its projection.
//!
//! Reference: OCC BRepGProp's adaptive volume (1e-9) of the exported
//! STEP, 0.030047407143187158 (2026-10-17; the 0.03006067 recorded on
//! 2026-10-13 predates later propeller changes). Measured 2.7e-6 off,
//! then 8.6e-7 after PK-877/PK-885; Richardson-extrapolated sliver fans
//! (the hub's blade-root slivers sum to 5e-6, and 16-panel fans left
//! ~2e-8 of them) bring it to 2.4e-7. The remaining residue is the five
//! blade faces' mesh dependence (5.5e-9 each between chord scales),
//! which PK-818 still tracks.

use popkern_rhai::ScriptHost;

#[test]
fn the_propeller_volume_matches_opencascade() {
    let source = include_str!("../examples/propeller.rhai");
    let mut run = ScriptHost::new()
        .run(source, "propeller.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    let parts = run
        .model
        .take()
        .expect("model")
        .to_parts_all()
        .expect("parts");
    let brep = parts[0].brep();
    let volume = popkern_measure::body_volume(brep, &brep.tolerance).expect("volume");
    let occ = 0.030047407143187158;
    let rel = (volume - occ) / occ;
    assert!(
        rel.abs() < 5e-7,
        "propeller volume {volume} vs OCC {occ} ({rel:e})"
    );
}
