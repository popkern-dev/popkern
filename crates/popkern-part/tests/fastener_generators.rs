//! PK-460 fastener generator serves: hex bolt, socket cap,
//! countersunk screw, nut with chamfers — watertight gates plus
//! dimensional pins (the exact thread-form volume oracles live in
//! `fastener_thread.rs`; these check the assembled fasteners).

use popkern_core::tolerance::Tolerance;
use popkern_part::Part;
use popkern_part::fastener::{CountersunkScrew, HexBolt, HexNut, SocketCapScrew, ThreadSpec};

/// Bounding box over the TESSELLATED surface — the BRep vertex set
/// under-reads curved faces (a cylindrical head's only vertices sit
/// on its seam azimuth, so a vertex bbox misses the diameter).
fn bbox(part: &Part) -> ([f64; 3], [f64; 3]) {
    let mesh = part
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("tessellate");
    let mut lo = [f64::MAX; 3];
    let mut hi = [f64::MIN; 3];
    for v in &mesh.vertices {
        let p = v.position;
        for (i, c) in [p.x, p.y, p.z].into_iter().enumerate() {
            lo[i] = lo[i].min(c);
            hi[i] = hi[i].max(c);
        }
    }
    (lo, hi)
}

#[test]
fn m6_hex_bolt_serves_watertight_with_iso_dimensions() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    let bolt = HexBolt::new(spec, 8.0)
        .thread_length(5.0)
        .build(Tolerance::default())
        .expect("hex bolt");
    assert_eq!(bolt.open_edge_count().expect("mesh"), 0, "watertight");
    let (lo, hi) = bbox(&bolt);
    // Head: across flats 10 → across corners 2/√3·10 ≈ 11.55; the
    // 30° corner chamfer trims the top band but the mid-height
    // corners keep the full width. Height: head −4 → tip 8.
    assert!((hi[2] - 8.0).abs() < 1e-6, "tip at z=8: {}", hi[2]);
    assert!((lo[2] + 4.0).abs() < 1e-6, "head base at z=−4: {}", lo[2]);
    // The head hexagon has one flat facing +Y, so its corners align
    // with the ±Y axis (at 90° and 270°) and the across-corners
    // diameter shows in the Y extent; the X extent measures across
    // flats.
    let across_corners = 2.0 / 3.0_f64.sqrt() * 10.0;
    assert!(
        (hi[1] - lo[1] - across_corners).abs() < 1e-6,
        "across corners {} vs {}",
        hi[1] - lo[1],
        across_corners
    );
    assert!(
        (hi[0] - lo[0] - 10.0).abs() < 1e-6,
        "across flats {} vs 10",
        hi[0] - lo[0]
    );
}

#[test]
fn m4_socket_cap_screw_serves_watertight() {
    let spec = ThreadSpec::parse("M4").expect("M4");
    let screw = SocketCapScrew::new(spec, 8.0)
        .thread_length(3.0)
        .build(Tolerance::default())
        .expect("socket cap screw");
    assert_eq!(screw.open_edge_count().expect("mesh"), 0, "watertight");
    let (lo, hi) = bbox(&screw);
    // ISO 4762 M4: head Ø7 × height 4.
    assert!((lo[2] + 4.0).abs() < 1e-6, "head base at z=−4: {}", lo[2]);
    assert!((hi[2] - 8.0).abs() < 1e-6, "tip at z=8: {}", hi[2]);
    assert!((hi[0] - lo[0] - 7.0).abs() < 1e-6, "head diameter 7");
}

#[test]
#[ignore = "blocked on the PK-460 walk-class campaign: stranded/loop-gap/LH/two-start merge defects (fork in flight)"]
fn m5_countersunk_screw_serves_watertight() {
    let spec = ThreadSpec::parse("M5").expect("M5");
    let screw = CountersunkScrew::new(spec, 10.0)
        .thread_length(3.0)
        .build(Tolerance::default())
        .expect("countersunk screw");
    assert_eq!(screw.open_edge_count().expect("mesh"), 0, "watertight");
    let (lo, hi) = bbox(&screw);
    // Flush top at z=0, ISO head Ø11.2, overall length 10.
    assert!(lo[2].abs() < 1e-6, "flush top at z=0: {}", lo[2]);
    assert!((hi[2] - 10.0).abs() < 1e-6, "tip at z=10: {}", hi[2]);
    assert!((hi[0] - lo[0] - 11.2).abs() < 1e-6, "head diameter 11.2");
}

#[test]
#[ignore = "blocked on the PK-460 walk-class campaign: stranded/loop-gap/LH/two-start merge defects (fork in flight)"]
fn m6_nut_with_chamfers_serves_watertight() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    let nut = HexNut::new(spec)
        .clearance(0.15)
        .build(Tolerance::default())
        .expect("nut");
    assert_eq!(nut.open_edge_count().expect("mesh"), 0, "watertight");
    let (lo, hi) = bbox(&nut);
    assert!(
        lo[2].abs() < 1e-6 && (hi[2] - 5.2).abs() < 1e-6,
        "ISO thickness 5.2"
    );
}

#[test]
#[ignore = "blocked on the PK-460 walk-class campaign: stranded/loop-gap/LH/two-start merge defects (fork in flight)"]
fn left_hand_two_start_bolt_serves() {
    let spec = ThreadSpec::parse("M8")
        .expect("M8")
        .left_handed()
        .with_starts(2);
    let bolt = HexBolt::new(spec, 9.0)
        .thread_length(6.0)
        .build(Tolerance::default())
        .expect("LH two-start bolt");
    assert_eq!(bolt.open_edge_count().expect("mesh"), 0, "watertight");
}
