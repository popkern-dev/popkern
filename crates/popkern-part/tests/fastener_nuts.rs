//! Default-option hex nuts (ISO 4032 proportions, exact basic thread
//! profile, corner and mouth chamfers): served through the cellular
//! boolean, watertight, with the volume OpenCASCADE measures on the
//! exported STEP.

use popkern_core::tolerance::Tolerance;
use popkern_part::Part;
use popkern_part::fastener::{HexNut, ThreadSpec};

/// Signed-tetrahedron volume of the display tessellation.
fn mesh_volume(part: &Part) -> f64 {
    let mesh = part
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("tessellate");
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
/// occ_solids.py` on `popkern run` exports, 2026-10-05), mm³.
const M3_NUT_OCC_VOLUME: f64 = 47.678_503_501;
const M5_NUT_OCC_VOLUME: f64 = 179.998_062_835;
/// M2: measured after PK-681 (its threaded section's OCC slab volumes
/// are uniform, 1.14110 per 0.1 mm, and sum to this whole-body value).
const M2_NUT_OCC_VOLUME: f64 = 17.695_515_551;
/// M1.6 and M2.5 (2026-10-06; OCC BRepGProp and OCC's tessellation
/// agree within 0.003%, and PK-671 measured 9.238 / 34.434).
const M1P6_NUT_OCC_VOLUME: f64 = 9.235_450_145;
const M2P5_NUT_OCC_VOLUME: f64 = 34.435_227_921;
/// Display-tessellation chord error budget against the exact volume.
const VOLUME_REL_TOL: f64 = 5.0e-3;

/// PK-673: the default M3 and M5 nuts declined at the TOP thread-mouth
/// chamfer ("stranded 3 half-claimed record sides" after a cellular
/// non-manifold fallback). Three stacked defects:
///
/// * The chamfer cone's circle at the crest radius crosses the
///   multi-turn crest band where the band's stripe straddles the
///   branch window the clip query was normalized to (u ∈ [−1.50,
///   +0.05] against [−2π, 0]); the plain crossing pass saw one of its
///   two crossings, the record read Outside, the crest band was never
///   split and dropped whole. Closed loops on multi-turn charts now
///   take crossings on every lift.
/// * The kept crest arc then runs through its closed record's
///   parameter seam; the far-tier endpoint re-trim read its ends raw
///   (decreasing window, the complement arc) and left both ends
///   6–7e-3 off their junctions. Closed-source windows are now
///   unwrapped, and a seam-crossing window re-trims exactly.
/// * The crest × flank records of the internal thread cut were fitted
///   at 204 control points for 258 exact samples, a least-squares
///   regime whose end spans ring: 1.0e-4 off the flank inside the
///   first interval. The mouth chamfer cut a short piece there and the
///   curve-surface gate declined it. The fitter's ladder now goes to
///   full interpolation past half the sample count.
#[test]
fn m3_and_m5_default_nuts_serve_cellular_with_occ_volume() {
    assert_nuts_serve(&[("M3", M3_NUT_OCC_VOLUME), ("M5", M5_NUT_OCC_VOLUME)]);
}

/// PK-681: the default M2 nut declined at BOTH mouth chamfers. Two
/// clip-substrate defects, neither M2-specific:
///
/// * Bottom mouth: the cone's circle on the crest band is one full
///   period at u ∈ [2π, 4π]; the clip query's single-period branch
///   landing folded each control point by rem_euclid, sending the
///   closing point from 2π to 0. No crossing with the band's boundary
///   was found, the record dropped, the uncut band read Inside and
///   vanished (14 single-use edges). The landing is now a rigid
///   whole-period translation.
/// * Top mouth: the cutter's seam edge pierced the crest band at
///   u = −0.0 and the edge-face trim-box prefilter rejected it against
///   [1.84, 2π] — the same azimuth. The crest arc then crossed the
///   cone's seam without a node and the cone's chart walk declined
///   (ChartOverwrap 6.53 > 2π) to a 3D-Newell walk that cannot cut the
///   winding separator. The prefilter now tests u modulo the period.
#[test]
fn smallest_default_nuts_serve_cellular_with_occ_volume() {
    assert_nuts_serve(&[
        ("M1.6", M1P6_NUT_OCC_VOLUME),
        ("M2", M2_NUT_OCC_VOLUME),
        ("M2.5", M2P5_NUT_OCC_VOLUME),
    ]);
}

/// Build each default nut and require a cellular, watertight serve
/// whose mesh volume matches OCC within the tessellation budget.
fn assert_nuts_serve(cases: &[(&str, f64)]) {
    let mut failures = Vec::new();
    for &(size, occ_volume) in cases {
        let spec = ThreadSpec::parse(size).expect("spec");
        let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
        let nut = match HexNut::new(spec).build(Tolerance::default()) {
            Ok(n) => n,
            Err(e) => {
                failures.push(format!("{size}: declined: {e}"));
                continue;
            }
        };
        let fallbacks = popkern_bool::merge_fallback::take_merge_fallback_events();
        if !fallbacks.is_empty() {
            failures.push(format!("{size}: graph-merge fallbacks {fallbacks:?}"));
        }
        let open = nut.open_edge_count().expect("mesh");
        if open != 0 {
            failures.push(format!("{size}: {open} open mesh edges"));
        }
        let volume = mesh_volume(&nut);
        let rel = (volume - occ_volume).abs() / occ_volume;
        if rel > VOLUME_REL_TOL {
            failures.push(format!(
                "{size}: mesh volume {volume:.4} vs OCC {occ_volume:.4} ({:.3}% off)",
                rel * 100.0
            ));
        }
    }
    assert!(failures.is_empty(), "default nuts: {failures:?}");
}
