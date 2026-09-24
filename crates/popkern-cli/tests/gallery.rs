//! PK-412: the examples gallery cannot rot. One table drives three
//! things — the native gate on every script (served, no failures, no
//! downgraded blends, watertight mesh, STEP through the OCC oracle
//! when present), the manifest the playground page reads
//! (`web/playground/examples.json` + the script copies beside it),
//! and the thumbnails. `POPKERN_UPDATE_GALLERY=1` regenerates the
//! manifest, copies and thumbnails; without it the test fails on any
//! drift, so a script edit that changes a feature count is a diff in
//! review.

use std::path::{Path, PathBuf};

use popkern_cli::{Outcome, execute, parse_args};
use popkern_playground::{HostOptions, Session};
use serde_json::{Value, json};

/// The gallery: file stem, title, blurb, expected feature count, and
/// the surface areas the script demonstrates (the guide cites them).
const GALLERY: &[(&str, &str, &str, usize, &str)] = &[
    (
        "boss_plate",
        "Boss plate",
        "A plate, a boss, a union and G1 rounds selected by geometry — the reference part.",
        6,
        "extrude, tool bodies, union, fillet with selectors, assert_watertight",
    ),
    (
        "revolved_knob",
        "Revolved knob",
        "A stepped section revolved about Z; exact planes, cylinders and rounded rims.",
        3,
        "path profiles, revolve, circular edge selection, fillet",
    ),
    (
        "shelled_enclosure",
        "Shelled enclosure",
        "A rounded block hollowed from the top to a uniform wall.",
        4,
        "fillet on vertical edges, shell with an open face",
    ),
    (
        "chamfered_bracket",
        "Chamfered bracket",
        "An L-bracket on the YZ plane with equal, two-distance and distance-angle chamfers.",
        5,
        "sketch planes other than XY, chamfer kinds, allow_empty",
    ),
    (
        "bored_block",
        "Bored block",
        "A through bore whose exact rims are labelled, then rounded by that label.",
        7,
        "subtract, exact boolean rims, label_edge, labeled() selection",
    ),
    (
        "tapered_fillet",
        "Tapered fillet",
        "A variable-radius round growing along one edge under a linear law.",
        3,
        "fillet radius laws",
    ),
    (
        "pyramid_apex",
        "Pyramid apex",
        "A pentagonal pyramid carved by five exact planes; its slant edges \
         rounded with curvature continuity (g2), and the five rounds closed \
         by the setback Gregory vertex cap with dual cross-jet corners.",
        18,
        "boolean carve with exact planes, chord-length edge selection, \
         vertex blend cap",
    ),
    (
        "parametric_plate",
        "Parametric plate",
        "Parameters, an expression and a named configuration drive a plate and its rib.",
        6,
        "param, param_expr, p(), expr(), configurations",
    ),
    (
        "project_box",
        "Project box",
        "A printable two-body box: shelled base with screw bosses from a loop, a bar-lip lid with loop-cut vents, exact-cone countersinks and a vents-off configuration.",
        71,
        "loops computing placements, shell, multi-body documents, configurations suppressing loop features, chord_length_between",
    ),
    (
        "ball_bearing",
        "Ball bearing",
        "A deep-groove bearing: races revolved with exact torus grooves, chamfered lead-ins, and a ring of balls placed by trigonometry — ten bodies.",
        30,
        "revolve with arc paths, exact torus and sphere surfaces, trigonometric placement loops, multi-body export",
    ),
    (
        "involute_gear",
        "Involute gear",
        "A twenty-tooth spur gear whose flanks are true involutes computed in the script, bored and keyed in one combined cut.",
        8,
        "computed polygon profiles, script functions, boolean tool unions",
    ),
    (
        "centrifugal_impeller",
        "Centrifugal impeller",
        "Seven backswept blades carved by one star-shaped cut whose spiral walls and tapered trailing edges this script computes.",
        11,
        "script functions, computed star polygons, logarithmic spirals, single-cut cavities",
    ),
    (
        "polyhedral_dice",
        "Polyhedral dice",
        "The five Platonic solids as a dice set — face normals computed from the \
         golden ratio, each die carved from a block by one exact plane cut per \
         face pair, and a filleted d6 dimpled with 21 spherical pips.",
        176,
        "plane_at custom sketch planes, intersect, computed face normals, \
         spherical dimple cuts, multi-body sets",
    ),
    (
        "coil_spring",
        "Coil spring",
        "Three springs from one verb: a ground-end compression spring, a tapered valve \
         spring, and a left-hand spring phased a quarter turn — every path an exact \
         rational helix.",
        12,
        "helical_sweep, exact helix paths, ground flat ends, taper_to, left_hand winding, \
         start_angle phasing",
    ),
    (
        "propeller",
        "Marine propeller",
        "Three lofted blades — stadium sections whose chord and pitch the script \
         computes per station — twisted, phased 120° apart and unioned onto a \
         revolved hub with a tapered nose.",
        26,
        "loft through computed sections, rational stadium profiles, rigid rotations, \
         chart-wrapping boolean unions on closed lofted walls",
    ),
    (
        "keyed_shaft",
        "Keyed shaft",
        "A shaft with a keyway cut by a moved and turned tool, and a chamfered lead-in.",
        8,
        "translate, rotate, subtract, chamfer on a circular edge",
    ),
    (
        "planetary_gearbox",
        "Planetary gearbox",
        "A 4:1 planetary set — sun, three phased planets, an internal ring gear cut by \
         a toothed disc, and a spider carrier — six meshing bodies whose involute teeth \
         and assembly phasing are all computed in the script.",
        45,
        "computed involute teeth, internal ring gears, meshing-assembly phasing, \
         the (Zs+Zr)%N constraint, single-polygon spiders, multi-body assemblies",
    ),
    (
        "threaded_stud",
        "Threaded stud",
        "The smallest real-thread part: an M10×1.5 stud from one verb — a coarse, \
         short coil that serves in seconds, with the tip left square for speed.",
        1,
        "threaded_rod, ISO thread designations, tip_chamfer option",
    ),
    (
        "hex_bolt",
        "Hex bolt",
        "A single short M6 hex bolt — head, plain shank and a real modeled ISO \
         thread from one verb, sized to build fast.",
        1,
        "bolt generator, thread_length, tip_chamfer option",
    ),
    (
        "countersunk_joint",
        "Countersunk joint",
        "A plate with a flush exact-cone clearance seat and the flat-head M5 screw \
         that fits it — one quick thread coil in the whole build.",
        5,
        "clearance_hole with countersink, countersunk head bolts, axis_at placement",
    ),
    (
        "fastener_kit",
        "Fastener kit",
        "An M6 hex bolt with its nut, a socket-head cap screw, and a countersunk \
         screw — real modeled ISO threads with heads, tip lead-ins, and runout \
         from the nominal tables.",
        9,
        "bolt and nut generators, ISO thread designations, head styles, \
         printed-fit clearance, multi-body layout",
    ),
    (
        "angled_fitting",
        "Angled fitting",
        "A wedge bracket tapped normal to its inclined face, a countersunk \
         clearance hole through the base, and the flat-head screw that fits it.",
        6,
        "tapped_hole and clearance_hole on arbitrary axes, axis_at placement, \
         countersinks, thread fit classes",
    ),
];

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("workspace root")
}

fn examples_dir() -> PathBuf {
    root().join("crates/popkern-rhai/examples")
}

fn page_dir() -> PathBuf {
    root().join("web/playground")
}

fn updating() -> bool {
    std::env::var("POPKERN_UPDATE_GALLERY").is_ok_and(|v| v == "1")
}

fn tmp(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join("gallery");
    std::fs::create_dir_all(&dir).expect("dir");
    dir.join(name)
}

fn oracle_passes(step: &Path) -> Option<bool> {
    let oracle = Path::new(&std::env::var("HOME").unwrap_or_default()).join("venvs/ocp/bin/python");
    let check = root().join("tools/occ-oracle/occ_check.py");
    if !(oracle.exists() && check.exists()) {
        return None;
    }
    let output = std::process::Command::new(&oracle)
        .arg(&check)
        .arg(step)
        .output()
        .expect("oracle");
    Some(String::from_utf8_lossy(&output.stdout).contains("\"passed\": true"))
}

#[test]
fn every_gallery_script_serves_clean_watertight_and_step_valid() {
    assert!(GALLERY.len() >= 8, "at least eight gallery scripts");
    let mut entries = Vec::new();
    for (stem, title, blurb, features, demonstrates) in GALLERY {
        let script = examples_dir().join(format!("{stem}.rhai"));
        let source = std::fs::read_to_string(&script).unwrap_or_else(|e| panic!("{stem}: {e}"));
        // Native engine: served, clean, every body watertight.
        let mut session = Session::new(HostOptions::default());
        let envelope = session.run(&source, &format!("{stem}.rhai"));
        assert!(
            envelope.outcome.is_served(),
            "[{stem}] {:?}",
            envelope.outcome
        );
        assert_eq!(envelope.features, *features, "[{stem}] feature count");
        let report = envelope.report.as_ref().expect("report");
        assert_eq!(report["failures"], 0, "[{stem}]");
        assert_eq!(
            report["downgraded_blends"], 0,
            "[{stem}] no script leans on a downgrade"
        );
        assert!(!envelope.bodies.is_empty(), "[{stem}] a body on screen");
        for body in &envelope.bodies {
            assert_eq!(body.open_edges, 0, "[{stem}] watertight mesh");
        }
        // The command: STEP + thumbnail, oracle on the STEP.
        let step = tmp(&format!("{stem}.step"));
        let png = if updating() {
            page_dir().join("examples").join(format!("{stem}.png"))
        } else {
            tmp(&format!("{stem}.png"))
        };
        let args: Vec<String> = [
            "run",
            script.to_str().expect("utf8"),
            "--step",
            step.to_str().expect("utf8"),
            "--png",
            png.to_str().expect("utf8"),
            "--width",
            "480",
            "--height",
            "360",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        let summary = execute(&parse_args(&args).expect("args"));
        assert!(
            matches!(summary.outcome, Outcome::Served),
            "[{stem}] {:?}",
            summary.outcome
        );
        match oracle_passes(&step) {
            Some(passed) => assert!(passed, "[{stem}] OCC oracle rejected the STEP"),
            None => eprintln!("[pk412] {stem}: OCC oracle absent"),
        }
        entries.push(json!({
            "file": format!("{stem}.rhai"),
            "title": title,
            "blurb": blurb,
            "features": features,
            "demonstrates": demonstrates,
            "thumbnail": format!("examples/{stem}.png"),
        }));
    }
    // The manifest and the script copies the page serves.
    let manifest = json!({ "examples": entries });
    let manifest_text = format!(
        "{}\n",
        serde_json::to_string_pretty(&manifest).expect("json")
    );
    let manifest_path = page_dir().join("examples.json");
    let copies = page_dir().join("examples");
    if updating() {
        std::fs::create_dir_all(&copies).expect("dir");
        std::fs::write(&manifest_path, &manifest_text).expect("write manifest");
        for (stem, ..) in GALLERY {
            std::fs::copy(
                examples_dir().join(format!("{stem}.rhai")),
                copies.join(format!("{stem}.rhai")),
            )
            .expect("copy");
        }
        return;
    }
    let committed = std::fs::read_to_string(&manifest_path).unwrap_or_default();
    assert_eq!(
        committed, manifest_text,
        "web/playground/examples.json drifted from the gallery table: run with \
         POPKERN_UPDATE_GALLERY=1"
    );
    let committed: Value = serde_json::from_str(&committed).expect("manifest json");
    for (stem, ..) in GALLERY {
        let source = std::fs::read(examples_dir().join(format!("{stem}.rhai"))).expect("read");
        let copy = std::fs::read(copies.join(format!("{stem}.rhai"))).unwrap_or_default();
        assert_eq!(
            source, copy,
            "[{stem}] web/playground/examples copy drifted: run with POPKERN_UPDATE_GALLERY=1"
        );
        assert!(
            copies.join(format!("{stem}.png")).exists(),
            "[{stem}] thumbnail missing: run with POPKERN_UPDATE_GALLERY=1"
        );
    }
    assert_eq!(
        committed["examples"].as_array().map(Vec::len),
        Some(GALLERY.len())
    );
}
