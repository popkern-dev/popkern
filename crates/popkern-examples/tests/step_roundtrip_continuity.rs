//! PK-506: STEP round-trip continuity audit. Every edge seam of a
//! filleted body is certified (PK-499 whole-edge certifier, commercial
//! gates) on the native body and again on the body a reader produces
//! from its STEP export. The reader's surfaces and pcurves are
//! certified as imported, never snapped back to native geometry. A
//! seam whose certified level drops across the round trip is a
//! continuity loss the STEP validity oracles cannot see.
//!
//! Legs:
//! - `popkern`: popkern-io's reader with default options;
//! - `occ`: OpenCASCADE reads the export and writes it back
//!   (`tools/occ-oracle/occ_rewrite.py`), and popkern-io imports
//!   OCC's file. This leg SKIPS when the OCP venv is absent (the
//!   established oracle-battery pattern).
//!
//! Fixtures: the fillet mileage corpus and three KP cap vertex blends.

use popkern_blend::ContinuityGates;
use popkern_blend::blend_continuity_certification::evidence_certify_options;
use std::collections::BTreeMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::certify::SeamCertifyOptions;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::edge_seam_certificate::certify_edge_seam;
use popkern_core::topology::{BRep, EdgeHandle};
use popkern_part::prelude::*;

#[path = "step_roundtrip_continuity_support/kp_caps.rs"]
mod kp_caps;
#[path = "fillet_mileage_corpus_support/parts.rs"]
mod parts;
use parts::*;

/// The level a seam's certificate proves against the commercial gates.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Level {
    /// The bounds do not cover the seam (nothing proved).
    Unbounded,
    /// The position gap is not proved within the G0 gate.
    Open,
    G0,
    G1,
    G2,
}

/// One edge's geometric identity, stable across a round trip.
struct EdgeKey {
    ends: [Point3; 2],
    mid: Point3,
}

fn edge_key(brep: &BRep, edge: EdgeHandle) -> Option<EdgeKey> {
    let data = brep.edges.get(edge)?;
    let curve = data.curve.as_ref()?;
    let half_edge = brep.half_edges.get(data.half_edge)?;
    let twin = brep.half_edges.get(half_edge.twin)?;
    let geometry = &curve.geometry as &dyn CurveEval;
    Some(EdgeKey {
        ends: [
            brep.vertices.get(half_edge.origin)?.position,
            brep.vertices.get(twin.origin)?.position,
        ],
        mid: geometry.point_at(0.5 * (curve.t_start + curve.t_end)),
    })
}

/// The edge of `imported` carrying `key`: unordered endpoints within
/// `tol`, and the native midpoint within `tol` of its curve.
fn matched_edge(key: &EdgeKey, imported: &BRep, tol: f64) -> Option<EdgeHandle> {
    let near = |a: &Point3, b: &Point3| a.distance_to(b) <= tol;
    imported.edges.keys().find(|&candidate| {
        let Some(other) = edge_key(imported, candidate) else {
            return false;
        };
        let ends = (near(&key.ends[0], &other.ends[0]) && near(&key.ends[1], &other.ends[1]))
            || (near(&key.ends[0], &other.ends[1]) && near(&key.ends[1], &other.ends[0]));
        ends && imported.edges[candidate].curve.as_ref().is_some_and(|c| {
            let g = &c.geometry as &dyn CurveEval;
            g.closest_point(&key.mid)
                .is_some_and(|t| near(&g.point_at(t), &key.mid))
        })
    })
}

/// Refinement budget per seam. Construction evidence certifies under
/// 256 intervals; the audit affords more so that a level is decided,
/// not stopped by budget (measured: 256 -> 4096 changed two of ~900
/// verdicts; 1024 settles both).
const AUDIT_CERTIFY_BUDGET: usize = 1024;

/// A seam's certified level: the highest gate every upper bound meets.
fn certify(brep: &BRep, edge: EdgeHandle, gates: &ContinuityGates) -> Level {
    let options = SeamCertifyOptions {
        budget: AUDIT_CERTIFY_BUDGET,
        ..evidence_certify_options(brep, gates)
    };
    let c = certify_edge_seam(brep, edge, &options);
    if !c.bounds_cover_seam() {
        Level::Unbounded
    } else if c.g0_upper > gates.g0_gap {
        Level::Open
    } else if c.g1_upper > gates.g1_normal {
        Level::G0
    } else if c.g2_upper > gates.g2_curvature {
        Level::G1
    } else {
        Level::G2
    }
}

/// One fixture's audit on one leg: level drops keyed `"native->imported"`
/// (`Unmatched` when no imported edge carries the native seam).
type Drops = BTreeMap<String, usize>;

/// Every native seam with its certified level, certified once per
/// fixture and shared by both legs.
fn native_levels(native: &BRep) -> Vec<(EdgeKey, Level)> {
    let gates = ContinuityGates::commercial();
    native
        .edges
        .keys()
        .filter_map(|edge| Some((edge_key(native, edge)?, certify(native, edge, &gates))))
        .collect()
}

fn audit_leg(native: &[(EdgeKey, Level)], imported: &[BRep]) -> Drops {
    let gates = ContinuityGates::commercial();
    let mut drops = Drops::new();
    for (key, before) in native {
        let after = imported.iter().find_map(|body| {
            matched_edge(key, body, gates.g0_gap).map(|e| certify(body, e, &gates))
        });
        if after.is_none_or(|after| after < *before) {
            let to = after.map_or("Unmatched".to_string(), |l| format!("{l:?}"));
            *drops.entry(format!("{before:?}->{to}")).or_default() += 1;
        }
    }
    drops
}

fn export(brep: &BRep, path: &std::path::Path) {
    let name = path
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| "part".into());
    let doc = popkern_io::StepDocument {
        products: vec![popkern_io::StepProductNode {
            name: name.clone(),
            id: name.to_uppercase(),
            children: vec![],
            shapes: vec![std::sync::Arc::new(popkern_io::ShapeRep {
                kind: popkern_io::ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
                context: popkern_io::RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: popkern_io::ProductMetadata::default(),
        }],
        roots: vec![0],
        units: popkern_io::UnitContext::default(),
        opaque_bag: popkern_io::OpaqueBag::default(),
    };
    popkern_io::export_step_file(&doc, path, &popkern_io::StepExportOptions::default())
        .unwrap_or_else(|e| panic!("export {}: {e:?}", path.display()));
}

/// The solids popkern-io imports from `path`.
fn import(path: &std::path::Path) -> Result<Vec<BRep>, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{e}"))?;
    let result = popkern_io::import_step(&bytes, &popkern_io::StepImportOptions::default())
        .map_err(|e| format!("{e}"))?;
    Ok(result
        .document
        .products
        .into_iter()
        .flat_map(|p| p.shapes)
        .flat_map(|s| s.solids.clone())
        .collect())
}

/// OCC rewrites each `(src, dst)` pair; `None` when the venv is absent.
fn occ_rewrite(pairs: &[(std::path::PathBuf, std::path::PathBuf)]) -> Option<Vec<bool>> {
    let python = std::path::Path::new(&std::env::var("HOME").unwrap_or_default())
        .join("venvs/ocp/bin/python");
    let script = std::path::Path::new("../../tools/occ-oracle/occ_rewrite.py");
    if !python.exists() || !script.exists() {
        return None;
    }
    let mut command = std::process::Command::new(&python);
    command.arg(script);
    for (src, dst) in pairs {
        command.arg(src).arg(dst);
    }
    let output = command.output().ok()?;
    let stdout = String::from_utf8_lossy(&output.stdout);
    let verdicts: Vec<bool> = stdout
        .lines()
        .filter(|l| l.starts_with('{'))
        .map(|l| l.contains("\"ok\": true"))
        .collect();
    (verdicts.len() == pairs.len()).then_some(verdicts)
}

type Fixture = (&'static str, fn() -> Result<BRep, String>);

fn part_body(build: fn() -> Result<Part, PartError>) -> Result<BRep, String> {
    build()
        .map(|p| p.brep().clone())
        .map_err(|e| format!("{e}"))
}

macro_rules! corpus_fixture {
    ($name:literal, $build:ident) => {
        ($name, || part_body($build))
    };
}

fn fixtures() -> Vec<Fixture> {
    vec![
        corpus_fixture!("mounting_bracket", build_mounting_bracket),
        corpus_fixture!("ribbed_plate", build_ribbed_plate),
        corpus_fixture!("bossed_housing", build_bossed_housing),
        corpus_fixture!("turned_bushing", build_turned_bushing),
        corpus_fixture!("pocketed_tray", build_pocketed_tray),
        corpus_fixture!("slotted_block", build_slotted_block),
        corpus_fixture!("countersunk_plate", build_countersunk_plate),
        corpus_fixture!("oblique_countersunk_plate", build_oblique_countersunk_plate),
        corpus_fixture!("drafted_cup", build_drafted_cup),
        corpus_fixture!("cross_drilled_boss", build_cross_drilled_boss),
        corpus_fixture!("coaxial_bushing", build_coaxial_bushing),
        corpus_fixture!("tilted_hole_plate", build_tilted_hole_plate),
        corpus_fixture!("offset_bracket", build_offset_bracket),
        corpus_fixture!("rim_chamfered_disc", build_rim_chamfered_disc),
        corpus_fixture!("star_socket_plate", build_star_socket_plate),
        corpus_fixture!("mixed_corner_block", build_mixed_corner_block),
        corpus_fixture!("chamfered_peg", build_chamfered_peg),
        corpus_fixture!("skewed_dress_bar", build_skewed_dress_bar),
        corpus_fixture!("counterbored_plate", build_counterbored_plate),
        corpus_fixture!("angled_rib_plate", build_angled_rib_plate),
        corpus_fixture!("double_boss_plate", build_double_boss_plate),
        corpus_fixture!("notched_flange", build_notched_flange),
        corpus_fixture!("pentagon_tower", build_pentagon_tower),
        ("kp_pyramid_apex", kp_caps::pyramid_apex),
        ("kp_cube_corner", kp_caps::cube_corner),
        ("kp_notched_saddle", kp_caps::notched_saddle),
    ]
}

/// Recorded level drops per fixture on the popkern leg; `[]` = every
/// seam keeps its certified level. The remaining class is filed:
/// Gregory cap export (`G2->G1`) is PK-551. (PK-550 retired the
/// conic/sinusoid pcurve class: 13 `G0->Open` and 31 `G0|G1->Unbounded`
/// drops.) Update a row only with the measured delta.
const POPKERN_BASELINE: &[(&str, &[(&str, usize)])] = &[
    ("mounting_bracket", &[]),
    ("ribbed_plate", &[]),
    ("bossed_housing", &[]),
    ("turned_bushing", &[]),
    ("pocketed_tray", &[]),
    ("slotted_block", &[]),
    ("countersunk_plate", &[]),
    ("oblique_countersunk_plate", &[]),
    ("drafted_cup", &[]),
    ("cross_drilled_boss", &[]),
    ("coaxial_bushing", &[]),
    ("tilted_hole_plate", &[]),
    ("offset_bracket", &[]),
    ("rim_chamfered_disc", &[]),
    ("star_socket_plate", &[]),
    ("mixed_corner_block", &[]),
    ("chamfered_peg", &[]),
    ("skewed_dress_bar", &[]),
    ("counterbored_plate", &[]),
    ("angled_rib_plate", &[]),
    ("double_boss_plate", &[]),
    ("notched_flange", &[]),
    ("pentagon_tower", &[]),
    ("kp_pyramid_apex", &[("G2->G1", 10)]),
    ("kp_cube_corner", &[("G2->G1", 6)]),
    ("kp_notched_saddle", &[]),
];

/// Total recorded drops on the OCC leg (OpenCASCADE rewrite, then
/// popkern-io import): the popkern leg's 16 Gregory cap drops (PK-551)
/// plus tilted_hole_plate's drill-exit seam, already `Open` natively
/// (certified gap 7.7e-6..1.27e-5 against the 1e-5 gate, PK-555), which
/// OpenCASCADE's rewrite takes to `Unbounded`. PK-552 retired 55
/// foreign-pcurve drops and PK-554 the 4 exported-polyline drops.
const OCC_BASELINE_DROPS: usize = 17;

fn as_drops(row: &[(&str, usize)]) -> Drops {
    row.iter().map(|&(k, n)| (k.to_string(), n)).collect()
}

#[test]
fn certified_seam_continuity_survives_the_step_round_trip() {
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("step_roundtrip_continuity");
    std::fs::create_dir_all(&dir).unwrap_or_else(|e| panic!("{e}"));
    let fixtures = fixtures();
    assert_eq!(fixtures.len(), POPKERN_BASELINE.len());
    let mut bodies = Vec::new();
    for ((name, build), &(pinned, _)) in fixtures.into_iter().zip(POPKERN_BASELINE) {
        assert_eq!(name, pinned, "baseline rows follow the fixture order");
        let body = build().unwrap_or_else(|e| panic!("{name}: fixture declined: {e}"));
        let path = dir.join(format!("{name}.step"));
        export(&body, &path);
        bodies.push((name, native_levels(&body), path));
    }
    let mut measured = Vec::new();
    for ((name, native, path), &(_, row)) in bodies.iter().zip(POPKERN_BASELINE) {
        let imported = import(path).unwrap_or_else(|e| panic!("{name}: import: {e}"));
        let drops = audit_leg(native, &imported);
        eprintln!("[pk506] {name} popkern drops {drops:?}");
        measured.push((*name, drops, as_drops(row)));
    }
    let mismatched: Vec<_> = measured
        .iter()
        .filter(|(_, got, want)| got != want)
        .collect();
    assert!(
        mismatched.is_empty(),
        "popkern-leg drops moved off the baseline: {mismatched:?}"
    );

    let pairs: Vec<_> = bodies
        .iter()
        .map(|(name, _, path)| (path.clone(), dir.join(format!("{name}.occ.step"))))
        .collect();
    let Some(verdicts) = occ_rewrite(&pairs) else {
        eprintln!("[pk506] OCC leg skipped: OCP venv absent");
        return;
    };
    let mut occ_drops = 0;
    for (((name, native, _), (_, occ_path)), ok) in bodies.iter().zip(&pairs).zip(verdicts) {
        assert!(ok, "{name}: OpenCASCADE rewrite failed");
        let imported = import(occ_path).unwrap_or_else(|e| panic!("{name}: OCC import: {e}"));
        let drops = audit_leg(native, &imported);
        eprintln!("[pk506] {name} occ drops {drops:?}");
        occ_drops += drops.values().sum::<usize>();
    }
    assert_eq!(
        occ_drops, OCC_BASELINE_DROPS,
        "OCC-leg drops moved off the baseline"
    );
}
