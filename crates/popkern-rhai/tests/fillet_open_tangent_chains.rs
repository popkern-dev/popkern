//! Open tangent-chain fillets (PK-735,
//! `docs/architecture/fillet-open-tangent-chains.md`): a run of
//! filleted edges joined end to end on one tangent line — part of a
//! rounded rectangle's rim (line, arc, line, …) or one edge alone —
//! fillets as one chain. Its interior junctions meet on the exact
//! section like a closed chain (PK-733). Its two ends are:
//! - TANGENT CAPS: the rim continues unfilleted and tangent past the
//!   end, and the strip stops on its section under a planar cap;
//! - END FACES: the end is a sharp corner, and the strip runs out on
//!   the planar face holding the corner's other edges (circle, ellipse
//!   or fitted section).
//!
//! Volumes are analytic: the strip removes (convex) or adds (concave)
//! the section `(1 − π/4) r²` along exactly the selected edges,
//! straight runs by length (an oblique cut by its centroid line) and
//! arcs by Pappus.

use popkern_rhai::ScriptHost;

fn served_volume(source: &str) -> f64 {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(source, "chain.rhai")
        .unwrap_or_else(|e| panic!("{e}\n{}", e.payload));
    let mut model = outcome.model.expect("model");
    let part = model.to_part().expect("part");
    let brep = part.brep();
    popkern_measure::body_volume(brep, &brep.tolerance).expect("volume")
}

/// Section area and its centroid offset from the sharp corner.
fn section(r: f64) -> (f64, f64) {
    let pi = std::f64::consts::PI;
    (
        (1.0 - pi / 4.0) * r * r,
        r * (10.0 - 3.0 * pi) / (3.0 * (4.0 - pi)),
    )
}

/// Section volume along `straight` mm of lines and `arcs` quarter arcs
/// of radius `corner`.
fn chain_change(straight: f64, arcs: f64, corner: f64, r: f64) -> f64 {
    let (a, x) = section(r);
    a * (straight + arcs * 0.5 * std::f64::consts::PI * (corner - x))
}

fn rounded_area(w: f64, h: f64, corner: f64) -> f64 {
    w * h - (4.0 - std::f64::consts::PI) * corner * corner
}

fn assert_rel(got: f64, want: f64, what: &str) {
    let rel = (got - want).abs() / want;
    assert!(
        rel < 1e-4,
        "{what}: volume {got} vs analytic {want} (rel {rel:.2e})"
    );
}

/// The rim edges touching any of `points` (all at one height).
fn rim_at(points: &[[f64; 3]]) -> String {
    points
        .iter()
        .map(|p| {
            format!(
                "edges().at_height({:?}).at_vertex([{:?}, {:?}, {:?}])",
                p[2], p[0], p[1], p[2]
            )
        })
        .collect::<Vec<_>>()
        .join(".or(")
        + ")".repeat(points.len() - 1).as_str()
}

/// 60 × 40 × 20 block centred on the origin, vertical corners rounded
/// 8, with `selection` filleted at `r`.
fn rounded_box(selection: &str, r: f64) -> String {
    format!(
        r#"
let m = model();
let b = m.profile(xy()).rect(60.0, 40.0).sketch("block").extrude(blind(20.0));
b = b.fillet(edges().straight().direction("z"), 8.0, #{{ label: "corners" }});
b = b.fillet({selection}, {r:?}, #{{ label: "partial_rim" }});
m.assert_watertight();
m
"#
    )
}

/// A 50 × 30 pocket 12 deep (corners rounded 6) in an 80 × 60 × 20
/// block, with `selection` on its floor rim filleted at `r`.
fn rounded_pocket(selection: &str, r: f64) -> String {
    format!(
        r#"
let m = model();
let b = m.profile(xy()).rect(80.0, 60.0).sketch("block").extrude(blind(20.0));
let pocket = m.tool(|t| {{
    let p = t.profile(xy_at(8.0)).rect(50.0, 30.0).sketch("pocket").extrude(blind(20.0));
    p.fillet(edges().straight().direction("z"), 6.0, #{{ label: "pocket_corners" }});
}});
b = b.subtract(pocket, #{{ label: "pocket" }});
b = b.fillet({selection}, {r:?}, #{{ label: "partial_floor_rim" }});
m.assert_watertight();
m
"#
    )
}

/// The rounded box's three-edge (one corner arc) and five-edge (both
/// +x corner arcs) top-rim chains, with their straight lengths and arc
/// counts.
fn box_chains() -> [(String, f64, f64); 2] {
    let z = 20.0;
    [
        (
            rim_at(&[[22.0, 20.0, z], [30.0, 12.0, z]]),
            44.0 + 24.0,
            1.0,
        ),
        (
            rim_at(&[
                [22.0, 20.0, z],
                [30.0, 12.0, z],
                [30.0, -12.0, z],
                [22.0, -20.0, z],
            ]),
            2.0 * 44.0 + 24.0,
            2.0,
        ),
    ]
}

/// The rounded pocket's three- and five-edge floor-rim chains.
fn pocket_chains() -> [(String, f64, f64); 2] {
    let z = 8.0;
    [
        (rim_at(&[[19.0, 15.0, z], [25.0, 9.0, z]]), 38.0 + 18.0, 1.0),
        (
            rim_at(&[
                [19.0, 15.0, z],
                [25.0, 9.0, z],
                [25.0, -9.0, z],
                [19.0, -15.0, z],
            ]),
            2.0 * 38.0 + 18.0,
            2.0,
        ),
    ]
}

const LONG_LINES: &str = r#"edges().at_height(20.0).straight().direction("x")"#;

/// A block with one rounded vertical corner: its top rim's line, arc,
/// line chain ends at two sharp corners.
fn one_corner_block(r: f64) -> String {
    let three = rim_at(&[[22.0, 20.0, 20.0], [30.0, 12.0, 20.0]]);
    format!(
        r#"
let m = model();
let b = m.profile(xy()).rect(60.0, 40.0).sketch("block").extrude(blind(20.0));
b = b.fillet(edges().straight().direction("z").at_vertex([30.0, 20.0, 0.0]), 8.0, #{{ label: "corner" }});
b = b.fillet({three}, {r:?}, #{{ label: "rim" }});
m.assert_watertight();
m
"#
    )
}

/// The one-corner block's chain plus the vertical edge at its
/// (30, -20) end: the chain's end meets that fillet in a mitre corner,
/// the top edge y = -20 staying sharp.
fn one_corner_block_mitre(r: f64) -> String {
    let three = rim_at(&[[22.0, 20.0, 20.0], [30.0, 12.0, 20.0]]);
    format!(
        r#"
let m = model();
let b = m.profile(xy()).rect(60.0, 40.0).sketch("block").extrude(blind(20.0));
b = b.fillet(edges().straight().direction("z").at_vertex([30.0, 20.0, 0.0]), 8.0, #{{ label: "corner" }});
b = b.fillet({three}.or(edges().straight().direction("z").at_vertex([30.0, -20.0, 0.0])), {r:?}, #{{ label: "rim" }});
m.assert_watertight();
m
"#
    )
}

/// The one-corner block's chain plus the vertical edge AND the top
/// edge y = -20 at its (30, -20) end: three fillets meet there in a
/// vertex blend (PK-770).
fn one_corner_block_vertex_blend(r: f64) -> String {
    let three = rim_at(&[[22.0, 20.0, 20.0], [30.0, 12.0, 20.0]]);
    format!(
        r#"
let m = model();
let b = m.profile(xy()).rect(60.0, 40.0).sketch("block").extrude(blind(20.0));
b = b.fillet(edges().straight().direction("z").at_vertex([30.0, 20.0, 0.0]), 8.0, #{{ label: "corner" }});
b = b.fillet({three}.or(edges().straight().direction("z").at_vertex([30.0, -20.0, 0.0])).or(edges().at_height(20.0).straight().direction("x").at_vertex([30.0, -20.0, 20.0])), {r:?}, #{{ label: "rim" }});
m.assert_watertight();
m
"#
    )
}

/// An 80 × 60 × 20 block with a 50 × 30 pocket 12 deep whose one
/// vertical corner (25, 15) is rounded 6: the floor rim's line, arc,
/// line chain plus the pocket's vertical edge at (25, -15) — a concave
/// mitre at the chain's (25, -15, 8) end.
fn one_corner_pocket_mitre(r: f64) -> String {
    let rim = rim_at(&[[19.0, 15.0, 8.0], [25.0, 9.0, 8.0]]);
    format!(
        r#"
let m = model();
let b = m.profile(xy()).rect(80.0, 60.0).sketch("block").extrude(blind(20.0));
let pocket = m.tool(|t| {{
    let p = t.profile(xy_at(8.0)).rect(50.0, 30.0).sketch("pocket").extrude(blind(20.0));
    p.fillet(edges().straight().direction("z").at_vertex([25.0, 15.0, 8.0]), 6.0, #{{ label: "pocket_corner" }});
}});
b = b.subtract(pocket, #{{ label: "pocket" }});
b = b.fillet({rim}.or(edges().straight().direction("z").at_vertex([25.0, -15.0, 8.0])), {r:?}, #{{ label: "rim" }});
m.assert_watertight();
m
"#
    )
}

/// The trapezoid block: its -x end face slants from (-30, -20) to
/// (-20, 20), x = -30 + (y + 20) / 4.
fn trapezoid_block(r: f64) -> String {
    let rim = rim_at(&[[22.0, 20.0, 20.0], [30.0, 12.0, 20.0]]);
    format!(
        r#"
let m = model();
let b = m.profile(xy()).polygon([[-30.0, -20.0], [30.0, -20.0], [30.0, 20.0], [-20.0, 20.0]]).sketch("block").extrude(blind(20.0));
b = b.fillet(edges().straight().direction("z").at_vertex([30.0, 20.0, 0.0]), 8.0, #{{ label: "corner" }});
b = b.fillet({rim}, {r:?}, #{{ label: "rim" }});
m.assert_watertight();
m
"#
    )
}

/// The hook profile: a line (-20, 0) → (0, 0), then a quarter arc
/// about (0, 20) tangent to it, to (20, 20); the top y = 20 closes it
/// square. Its top rim's line and arc are a two-edge chain whose arc
/// end runs out on the y = 20 face (a NURBS strip: fitted section);
/// the line end runs out on x = -20. 10 high.
fn hook_block(r: f64) -> String {
    format!(
        r#"
let m = model();
let b = m.profile(xy()).path(-20.0, 0.0).line_to(0.0, 0.0).arc_to([14.142135623730951, 5.857864376269049], [20.0, 20.0]).line_to(-20.0, 20.0).close().sketch("hook").extrude(blind(10.0));
b = b.fillet(edges().at_height(10.0).at_vertex([0.0, 0.0, 10.0]), {r:?}, #{{ label: "rim" }});
m.assert_watertight();
m
"#
    )
}

#[test]
fn rounded_box_partial_top_rims_fillet_as_open_chains() {
    for (rim, straight, arcs) in box_chains() {
        for r in [1.0, 3.0] {
            let want = rounded_area(60.0, 40.0, 8.0) * 20.0 - chain_change(straight, arcs, 8.0, r);
            assert_rel(
                served_volume(&rounded_box(&rim, r)),
                want,
                &format!("convex {arcs}-arc chain r = {r}"),
            );
        }
    }
}

#[test]
fn rounded_pocket_partial_floor_rims_fillet_as_open_concave_chains() {
    for (rim, straight, arcs) in pocket_chains() {
        for r in [1.0, 2.0] {
            let pocket = rounded_area(50.0, 30.0, 6.0) * 12.0;
            let want = 80.0 * 60.0 * 20.0 - pocket + chain_change(straight, arcs, 6.0, r);
            assert_rel(
                served_volume(&rounded_pocket(&rim, r)),
                want,
                &format!("concave {arcs}-arc chain r = {r}"),
            );
        }
    }
}

#[test]
fn rounded_box_single_rim_lines_fillet_with_section_end_caps() {
    // The two long top lines: one-edge chains, both ends tangent to the
    // unfilleted corner arcs.
    for r in [1.0, 3.0] {
        let want = rounded_area(60.0, 40.0, 8.0) * 20.0 - chain_change(2.0 * 44.0, 0.0, 8.0, r);
        assert_rel(
            served_volume(&rounded_box(LONG_LINES, r)),
            want,
            &format!("single lines r = {r}"),
        );
    }
}

#[test]
fn a_chain_ending_at_sharp_corners_runs_out_on_the_end_faces() {
    let c = 8.0;
    for r in [1.0, 3.0] {
        let block = 60.0 * 40.0 * 20.0 - (1.0 - std::f64::consts::PI / 4.0) * c * c * 20.0;
        // Lines 52 (y = 20, x from -30 to 22) and 32 (x = 30, y from 12 to -20).
        let want = block - chain_change(52.0 + 32.0, 1.0, c, r);
        assert_rel(
            served_volume(&one_corner_block(r)),
            want,
            &format!("sharp-ended chain r = {r}"),
        );
    }
}

#[test]
fn a_chain_ending_at_a_corner_meets_another_fillet_in_a_mitre() {
    // Each strip's removal prism is cut by the mitre plane (y + z = 0):
    // its centroid line runs to y = -(20 - x̄) (the x = 30 line) and to
    // z = 20 - x̄ (the vertical fillet, from the base).
    let c = 8.0;
    for r in [1.0, 3.0] {
        let (_, x) = section(r);
        let block = 60.0 * 40.0 * 20.0 - (1.0 - std::f64::consts::PI / 4.0) * c * c * 20.0;
        let straight = 52.0 + (32.0 - x) + (20.0 - x);
        let want = block - chain_change(straight, 1.0, c, r);
        assert_rel(
            served_volume(&one_corner_block_mitre(r)),
            want,
            &format!("mitre-ended chain r = {r}"),
        );
    }
}

#[test]
fn a_concave_chain_ending_at_a_pocket_corner_meets_the_wall_fillet_in_a_mitre() {
    // Material ADDED: the floor rim (line 19 + 6-arc + line to the
    // mitre plane) and the wall fillet from the floor up to the pocket
    // top (z = 20, a sharp run-out on the block's top face).
    let c = 6.0;
    for r in [1.0, 2.0] {
        let (_, x) = section(r);
        let pocket = (50.0 * 30.0 - (1.0 - std::f64::consts::PI / 4.0) * c * c) * 12.0;
        // Lines: y = 15 from x = -25 to 19 (44), x = 25 from y = 9 to
        // the mitre plane at y = -(15 - x̄) (24 - x̄); the wall fillet
        // from the mitre plane z = 8 + x̄ up to z = 20 (12 - x̄).
        let straight = 44.0 + (24.0 - x) + (12.0 - x);
        let want = 80.0 * 60.0 * 20.0 - pocket + chain_change(straight, 1.0, c, r);
        assert_rel(
            served_volume(&one_corner_pocket_mitre(r)),
            want,
            &format!("concave mitre r = {r}"),
        );
    }
}

#[test]
fn a_chain_ending_at_a_three_fillet_corner_serves_through_the_vertex_blend() {
    // The chain's junctions are settled on their exact sections in the
    // per-edge pipeline (PK-770), which owns the vertex blend; the
    // served volume is cross-checked against OCC in the validity test.
    // The fillets remove material, at most the mitre-class amount plus
    // the y = -20 edge's full prism.
    for r in [1.0, 3.0] {
        let got = served_volume(&one_corner_block_vertex_blend(r));
        let mitre = served_volume(&one_corner_block_mitre(r));
        let (a, _) = section(r);
        assert!(
            got < mitre && got > mitre - a * 60.0,
            "r = {r}: vertex-blend body {got} outside ({}, {mitre})",
            mitre - a * 60.0
        );
    }
}

#[test]
fn a_chain_runs_out_on_an_oblique_end_face_on_its_ellipse_section() {
    let c = 8.0;
    for r in [1.0, 3.0] {
        let (_, x) = section(r);
        // A prism cut by a plane: section area × its centroid line.
        let x_end = -30.0 + (20.0 - x + 20.0) / 4.0;
        let straight = (22.0 - x_end) + 32.0;
        let block = (2200.0 - (1.0 - std::f64::consts::PI / 4.0) * c * c) * 20.0;
        let want = block - chain_change(straight, 1.0, c, r);
        assert_rel(
            served_volume(&trapezoid_block(r)),
            want,
            &format!("oblique end r = {r}"),
        );
    }
}

#[test]
fn a_chain_with_a_curved_end_runs_out_on_its_fitted_section() {
    // Line 20 plus a quarter arc of radius 20 (convex, centre inside);
    // both ends square on their end faces.
    for r in [1.0, 3.0] {
        let (a, x) = section(r);
        let pi = std::f64::consts::PI;
        let area = 400.0 + pi * 400.0 / 4.0;
        let want = area * 10.0 - a * (20.0 + 0.5 * pi * (20.0 - x));
        assert_rel(
            served_volume(&hook_block(r)),
            want,
            &format!("hook chain r = {r}"),
        );
    }
}

/// OpenCASCADE reads every served class as one valid solid (caps
/// oriented, end faces cut, sections on both faces). Skipped when the
/// OCP venv is absent.
#[test]
fn opencascade_reads_each_open_chain_class_as_one_valid_solid() {
    let home = std::path::PathBuf::from(std::env::var("HOME").unwrap_or_default());
    let python = home.join("venvs/ocp/bin/python");
    let tools = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occ-oracle");
    if !python.exists() || !tools.exists() {
        eprintln!("[pk735] OCC oracle absent; validity not cross-checked");
        return;
    }
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("pk735");
    std::fs::create_dir_all(&dir).expect("dir");
    let [_, five] = box_chains();
    let [three, _] = pocket_chains();
    for (name, source) in [
        ("box_five", rounded_box(&five.0, 3.0)),
        ("pocket_three", rounded_pocket(&three.0, 2.0)),
        ("lines", rounded_box(LONG_LINES, 3.0)),
        ("sharp", one_corner_block(3.0)),
        ("mitre", one_corner_block_mitre(3.0)),
        ("concave_mitre", one_corner_pocket_mitre(2.0)),
        ("vertex_blend", one_corner_block_vertex_blend(3.0)),
        ("oblique", trapezoid_block(3.0)),
        ("hook", hook_block(3.0)),
    ] {
        let mut host = ScriptHost::new();
        let mut model = host
            .run(&source, "chain.rhai")
            .unwrap_or_else(|e| panic!("{name}: {e}"))
            .model
            .expect("model");
        let step = dir.join(format!("{name}.step"));
        model.to_step(&step).expect("STEP export");
        let run = |script: &str| -> serde_json::Value {
            let out = std::process::Command::new(&python)
                .arg(tools.join(script))
                .arg(&step)
                .output()
                .expect("oracle");
            serde_json::from_slice(&out.stdout).expect("oracle json")
        };
        let check = run("occ_check.py");
        assert_eq!(
            check["passed"], true,
            "{name}: OCC rejects the STEP: {check}"
        );
        let out = std::process::Command::new(&python)
            .arg(tools.join("occ_solids.py"))
            .arg(&step)
            .arg("--adaptive")
            .output()
            .expect("oracle");
        let solids: serde_json::Value = serde_json::from_slice(&out.stdout).expect("oracle json");
        let list = solids["solids"].as_array().cloned().unwrap_or_default();
        assert_eq!(list.len(), 1, "{name}: one solid: {solids}");
        // Our volume agrees with OCC's adaptive integration.
        let occ = list[0]["volume"].as_f64().expect("volume");
        let ours = served_volume(&source);
        assert!(
            (ours - occ).abs() / occ < 1e-4,
            "{name}: volume {ours} vs OCC {occ}"
        );
    }
}
