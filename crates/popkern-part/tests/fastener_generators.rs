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

/// Nine-turn thread pin (PK-460): thread_length 8 plans a 9-turn tool
/// whose classify-window pad rephased the surface-curve AABB grid
/// into `padded_cell`'s inflation blind spot — the groove strip's
/// seam sample voted a contested 3–2 Outside, the kept sliver
/// mis-wired a one-pitch seam slit, and the subtraction declined on a
/// `Loop not closed (gap 1.0)`. Any thread_length whose planned turn
/// count differed from the serving fixtures reproduced it (8/12/13/16
/// × tl 8 all declined; tl 6/7 all served).
#[test]
fn m6_bolt_nine_turn_thread_serves_watertight() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    let bolt = HexBolt::new(spec, 14.0)
        .thread_length(8.0)
        .build(Tolerance::default())
        .expect("nine-turn hex bolt");
    assert_eq!(bolt.open_edge_count().expect("mesh"), 0, "watertight");
}

#[test]
fn m6_nut_with_chamfers_serves_watertight() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
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
    assert_cellular_served("M6 nut, clearance 0.15");
}

/// PK-566: at clearance 0.2 both mouth-chamfer cuts fell back from
/// cellular. The chamfer cone's two threaded-band fragments wrap the
/// cone seam with S-shaped loops, and the periodic-fold interior
/// sampler placed each one's sample in the OTHER fragment's region,
/// so their Inside/Outside verdicts swapped (the tool's apex ring
/// was kept, the mouth ring dropped). Every boolean of the nut must
/// now serve through the cellular pipeline.
#[test]
fn m6_nut_clearance_0p2_mouth_chamfers_serve_cellular() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
    let nut = HexNut::new(spec)
        .clearance(0.2)
        .build(Tolerance::default())
        .expect("nut");
    assert_eq!(nut.open_edge_count().expect("mesh"), 0, "watertight");
    assert_cellular_served("M6 nut, clearance 0.2");
}

fn assert_cellular_served(fixture: &str) {
    let fallbacks = popkern_bool::merge_fallback::take_merge_fallback_events();
    assert!(
        fallbacks.is_empty(),
        "{fixture}: {} boolean(s) fell back to graph-merge: {fallbacks:?}",
        fallbacks.len()
    );
}

/// PK-491 phase-luck pins: the 2026-09-16 sweep on the corrected-AABB
/// tree measured M5 countersunk 12/5 serving with 103 OPEN MESH EDGES
/// (a validator-clean mesh lie) and 11/6 declining on 21 strands —
/// both members of the seam-conflict class the trim-boundary /
/// embedding-unfold mechanisms retired. Pinned per the close-out
/// acceptance so the class cannot silently reopen at these phases.
#[test]
fn m5_countersunk_pk491_phase_points_serve_watertight() {
    let mut failures = Vec::new();
    for (length, thread) in [(12.0, 5.0), (11.0, 6.0)] {
        let spec = ThreadSpec::parse("M5").expect("M5");
        let screw = CountersunkScrew::new(spec, length)
            .thread_length(thread)
            .build(Tolerance::default());
        let screw = match screw {
            Ok(s) => s,
            Err(e) => {
                failures.push(format!("{length}/{thread}: {e}"));
                continue;
            }
        };
        let open = screw.open_edge_count().expect("mesh");
        eprintln!("[pk491-pin] M5 countersunk {length}/{thread}: {open} open edges");
        if open != 0 {
            failures.push(format!("{length}/{thread}: {open} open edges"));
        }
    }
    assert!(
        failures.is_empty(),
        "M5 countersunk phase points: {failures:?}"
    );
}

#[test]
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

/// The fastener-kit socket cap (M4 × 10, 6 mm of thread): its tip
/// chamfer cone charts on the off-principal branch [π/2, 5π/2]. The
/// principal-cell fold cannot reach it, so the rim-riding flank ×
/// cone record dropped whole, the thread-flank face was
/// graze-filtered unsplit, and the body opened 4217 mesh edges; the
/// clip's nothing-produced retry now translates such records onto the
/// face's branch (PK-564).
#[test]
fn m4_socket_cap_with_short_thread_serves_watertight() {
    let spec = ThreadSpec::parse("M4").expect("M4");
    let screw = SocketCapScrew::new(spec, 10.0)
        .thread_length(6.0)
        .build(Tolerance::default())
        .expect("socket cap screw");
    assert_eq!(screw.open_edge_count().expect("mesh"), 0, "watertight");
}

/// Signed-tetrahedron volume of a fine tessellation. The default
/// 0.005 chord height under-reads the M2 threads' volume by 0.74–0.81%
/// (thread depth 0.25 over most of the shank); at 0.0005 every size
/// lands within the 0.5% budget of the exact OCC volume (PK-671).
fn mesh_volume(part: &Part) -> f64 {
    let options = popkern_mesh::TessellationOptions {
        chord_height: 5.0e-4,
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

/// OpenCASCADE volumes of the exported STEP (`tools/occ-oracle/
/// occ_solids.py`), mm³: M3 hex 2026-10-04; M2.5 socket re-measured
/// 2026-10-05 for the ISO 4762 head row and the full 0.55·k socket
/// depth (PK-671).
const M3_HEX_OCC_VOLUME: f64 = 128.399_742_906;
const M2P5_SOCKET_OCC_VOLUME: f64 = 87.682_763_417;
/// Tessellation chord error budget against the exact volume.
const VOLUME_REL_TOL: f64 = 5.0e-3;

/// PK-669: two of the commonest small ISO screws declined.
///
/// * M3 hex (`hex head chamfer` stage): the corner-chamfer cutter's
///   base was offset only AXIALLY off the face, tilting its cone to
///   34.1° and moving its face circle to 0.5738·s against the
///   0.5774·s circumradius — six near-graze slivers 0.013 mm high
///   instead of a 30° chamfer, which the cellular boolean read as a
///   non-manifold edge. The cone now meets the face on the 0.95·s
///   chamfer circle at a true 30°.
/// * M2.5 socket (`external thread` stage, every M2.5 × thread 6):
///   splitting the 13.9-turn root edge, the pcurve image search's
///   64-sample scan aliased onto the adjacent helix turn, the split
///   pieces lost their pcurves, the root face never split at the
///   shank's end planes and the merge stranded 16 half-claimed
///   record sides.
///
/// Both serve watertight through the cellular pipeline with the
/// volume OpenCASCADE measures on their exported STEP, within the
/// display tessellation's chord error.
#[test]
fn m3_hex_and_m2p5_socket_screws_serve_cellular_with_occ_volume() {
    assert_screws_serve(&[
        ("M3", true, M3_HEX_OCC_VOLUME),
        ("M2.5", false, M2P5_SOCKET_OCC_VOLUME),
    ]);
}

/// PK-671: M1.6, M2 and M2.5 heads come from the ISO 4017 / ISO 4762
/// tables (the derived proportions gave M2.5 hex s = 4.0 against the
/// standard 5, M2.5 socket dk = 3.75 against 4.5). Each serves
/// watertight through the cellular pipeline with the OCC volume of
/// its exported STEP (`popkern run` + `occ_solids.py`, 2026-10-05;
/// 12 mm under head, 6 mm of thread).
#[test]
fn smallest_iso_hex_and_socket_screws_serve_cellular_with_occ_volume() {
    assert_screws_serve(&[
        ("M1.6", true, 30.750_368_751),
        ("M1.6", false, 30.682_636_290),
        ("M2", true, 52.741_125_668),
        ("M2", false, 54.073_291_468),
        ("M2.5", true, 89.131_024_958),
    ]);
}

/// Build each `(size, hex?, OCC volume)` screw (12 under head, 6 of
/// thread) and require a cellular, watertight serve whose mesh volume
/// matches OCC within the tessellation budget.
fn assert_screws_serve(cases: &[(&str, bool, f64)]) {
    let mut failures = Vec::new();
    for &(size, hex, occ_volume) in cases {
        let label = format!("{size} {}", if hex { "hex" } else { "socket" });
        let spec = ThreadSpec::parse(size).expect("spec");
        let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
        let built = if hex {
            HexBolt::new(spec, 12.0)
                .thread_length(6.0)
                .build(Tolerance::default())
        } else {
            SocketCapScrew::new(spec, 12.0)
                .thread_length(6.0)
                .build(Tolerance::default())
        };
        let screw = match built {
            Ok(s) => s,
            Err(e) => {
                failures.push(format!("{label}: declined: {e}"));
                continue;
            }
        };
        let fallbacks = popkern_bool::merge_fallback::take_merge_fallback_events();
        if !fallbacks.is_empty() {
            failures.push(format!("{label}: graph-merge fallbacks {fallbacks:?}"));
        }
        let open = screw.open_edge_count().expect("mesh");
        if open != 0 {
            failures.push(format!("{label}: {open} open mesh edges"));
        }
        let volume = mesh_volume(&screw);
        let rel = (volume - occ_volume).abs() / occ_volume;
        if rel > VOLUME_REL_TOL {
            failures.push(format!(
                "{label}: mesh volume {volume:.4} vs OCC {occ_volume:.4} ({:.3}% off)",
                rel * 100.0
            ));
        }
    }
    assert!(failures.is_empty(), "small ISO screws: {failures:?}");
}

/// PK-669: the hex corner chamfer meets the exposed head face on the
/// chamfer circle (0.95 × across flats) — the face's boundary is that
/// circle, not the hexagon nicked at its corners.
#[test]
fn hex_head_chamfer_meets_the_face_on_its_circle() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    let bolt = HexBolt::new(spec, 8.0)
        .thread_length(5.0)
        .build(Tolerance::default())
        .expect("hex bolt");
    let mesh = bolt
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("tessellate");
    // ISO 4017 M6: s = 10, k = 4; the exposed face is z = −4.
    let r_face = mesh
        .vertices
        .iter()
        .map(|v| v.position)
        .filter(|p| (p.z + 4.0).abs() < 1e-9)
        .map(|p| (p.x * p.x + p.y * p.y).sqrt())
        .fold(0.0_f64, f64::max);
    assert!(
        (r_face - 0.95 * 10.0 / 2.0).abs() < 1e-6,
        "exposed face reaches r = {r_face}, chamfer circle 4.75"
    );
}

/// PK-762: full-thread socket cap screws (the arm's M4..M8 bolts)
/// serve through the cellular pipeline. With only the groove tool's
/// centre floored at the head plane, its profile half-width notched the
/// head's underside; the cellular result was non-manifold and every
/// such screw paid a graph-merge retry (M5x12 5.3 s, of it 2.7 s the
/// retry; the arm's four bolt definitions 21 s of its 75 s build).
#[test]
fn full_thread_socket_cap_screws_serve_cellular() {
    for (size, length) in [("M4", 10.0), ("M5", 12.0), ("M6", 16.0), ("M8", 20.0)] {
        let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
        let spec = ThreadSpec::parse(size).expect("spec");
        let screw = SocketCapScrew::new(spec, length)
            .build(Tolerance::default())
            .expect("socket cap screw");
        assert_cellular_served(&format!("{size}x{length} socket cap, full thread"));
        assert_eq!(
            screw.open_edge_count().expect("mesh"),
            0,
            "{size}x{length} watertight"
        );
    }
}
