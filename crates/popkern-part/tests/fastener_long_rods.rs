//! Long and slotted M2.5 threaded rods serve through the cellular
//! boolean (PK-679, PK-680), watertight, with the volume OpenCASCADE's
//! TESSELLATION measures on the exported STEP. (OCC's BRepGProp volume
//! integration is unusable on these 100+-turn NURBS bodies: one rod
//! integrated to 188.6, 196.9 or 200.4 depending on its tolerance; the
//! tessellated volume, OCC slab sums and the uniform-section estimate
//! agree within 0.15%.)

use popkern_core::tolerance::Tolerance;
use popkern_part::fastener::{ThreadSpec, ThreadedRod};
use popkern_part::{Part, Plane, Profile};

/// Signed-tetrahedron volume at a chord height fine enough for the
/// M2.5 thread (0.27 deep) and coarse enough for a 112-turn coil.
fn mesh_volume(part: &Part) -> f64 {
    let options = popkern_mesh::TessellationOptions {
        chord_height: 1.0e-3,
        max_subdivisions: 512,
        ..Default::default()
    };
    let mesh = part.mesh(&options).expect("tessellate");
    let mut v6 = 0.0;
    for t in &mesh.triangles {
        let a = mesh.vertices[t[0]].position;
        let b = mesh.vertices[t[1]].position;
        let c = mesh.vertices[t[2]].position;
        v6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    (v6 / 6.0).abs()
}

fn rod(length: f64) -> Part {
    let spec = ThreadSpec::parse("M2.5").expect("M2.5");
    let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
    ThreadedRod::new(spec, length)
        .build(Tolerance::default())
        .expect("rod serves")
}

fn assert_cellular_watertight(part: &Part, fixture: &str) {
    let fallbacks = popkern_bool::merge_fallback::take_merge_fallback_events();
    assert!(
        fallbacks.is_empty(),
        "{fixture}: {} boolean(s) fell back to graph-merge: {fallbacks:?}",
        fallbacks.len()
    );
    assert_eq!(
        part.open_edge_count().expect("mesh"),
        0,
        "{fixture}: watertight"
    );
}

fn assert_volume(part: &Part, oracle: f64, fixture: &str) {
    let volume = mesh_volume(part);
    let rel = (volume - oracle).abs() / oracle;
    assert!(
        rel < 5.0e-3,
        "{fixture}: mesh volume {volume:.4} vs OCC tessellated {oracle:.4} ({:.3}% off)",
        rel * 100.0
    );
}

/// PK-679: a 50 mm rod (112 thread turns, tip chamfer on) fell back
/// from cellular at the tip-chamfer cone — the tool walls' coarse
/// edge-hit grid missed the cone seam's pierce, the cone's record
/// crossed its chart seam without a node (ChartOverwrap 7.854); with
/// every pierce found, the rod's bottom face kept its groove fragment
/// because the ray classifier's chart-edge band spanned half a turn of
/// the 112-turn wall.
#[test]
fn m2p5_fifty_mm_rod_serves_cellular_with_occ_volume() {
    let part = rod(50.0);
    assert_cellular_watertight(&part, "M2.5 x 50 rod");
    assert_volume(&part, 195.528, "M2.5 x 50 rod");
}

/// PK-680: a narrow slot through an M2.5 x 12 rod declined at result
/// validation — consolidation averaged the slot-top junction and the
/// crest-helix crossing of the slot wall (2.0e-3 apart) onto their
/// midpoint and the crest edge left its flank by 8.7e-4.
#[test]
fn slotted_m2p5_rod_serves_cellular() {
    let part = rod(12.0);
    let slot = Profile::on(Plane::xy_at(5.6))
        .rect(4.0, 0.6)
        .extrude(0.8)
        .expect("slot");
    let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
    let slotted = part.subtract(&slot).expect("slotted rod serves");
    assert_cellular_watertight(&slotted, "M2.5 x 12 slotted rod");
    // OCC tessellated 45.815; the unslotted rod's 46.869 less the
    // material inside the slot box (~1.05 by the strip integral).
    // OCC BRepGProp read 45.29 and OCC's box-Common of the slot 0.
    assert_volume(&slotted, 45.815, "M2.5 x 12 slotted rod");
}

/// PK-679: the 50 mm rod cut at mid-length by a planar box.
#[test]
#[ignore = "heavy: two 112-turn booleans (~3 min release); PK-679 acceptance row"]
fn m2p5_fifty_mm_rod_cut_at_mid_length_serves_cellular() {
    let part = rod(50.0);
    let cut = Profile::on(Plane::xy_at(25.3))
        .rect(10.0, 10.0)
        .extrude(30.0)
        .expect("cut box");
    let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
    let half = part.subtract(&cut).expect("cut rod serves");
    assert_cellular_watertight(&half, "M2.5 x 50 rod cut at z 25.3");
    assert_volume(&half, 98.9985, "M2.5 x 50 rod cut at z 25.3");
}
