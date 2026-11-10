//! A draft leans each face so the part draws along the pull: an outer
//! wall in as it rises, a hole's or a pocket's wall out (brep-bench
//! draft/hole-tapered and draft/pocket-walls, which narrowed: the
//! rotation followed the surface's own normal, which a hole's wall
//! turns away from the material).

use popkern_rhai::ScriptHost;

fn volume(source: &str) -> f64 {
    let mut host = ScriptHost::with_options(Default::default());
    let mut model = host
        .run(source, "draft")
        .unwrap_or_else(|e| panic!("{e}"))
        .model
        .expect("model");
    let parts = model.to_parts_all().expect("parts");
    parts
        .iter()
        .map(|p| popkern_measure::body_volume(p.brep(), &p.brep().tolerance).expect("volume"))
        .sum()
}

/// A hole of radius 6 through a 30 × 30 × 20 block, drafted 3° about the
/// floor: the block less a frustum, 6 at the floor, 6 + 20 tan 3° at the
/// top.
#[test]
fn a_drafted_hole_widens_toward_the_pull() {
    let v = volume(
        "let m = model();\nlet b = m.profile(xy()).rect(30.0, 30.0).sketch(\"block\").extrude(blind(20.0));\nlet bore = m.tool(|m| m.profile(xy_at(-1.0)).circle(6.0).sketch(\"bore\").extrude(blind(22.0)));\nb = b.subtract(bore);\nb = b.draft(faces().kind(\"cylinder\"), 3.0, xy());\nm\n",
    );
    let (r0, c) = (6.0f64, 20.0f64);
    let r1 = r0 + c * (3.0f64.to_radians()).tan();
    let exact = 30.0 * 30.0 * c - std::f64::consts::PI * c * (r0 * r0 + r0 * r1 + r1 * r1) / 3.0;
    assert!((v - exact).abs() < 1e-7 * exact, "{v} vs {exact}");
}

/// A 24 × 14 pocket 10 deep in a 40 × 30 × 15 block, its walls drafted
/// 5° about its floor: abc − (A₀h + P₀th²/2 + 4t²h³/3).
#[test]
fn a_drafted_pocket_widens_toward_its_opening() {
    let v = volume(
        "let m = model();\nlet b = m.profile(xy()).rect(40.0, 30.0).at(20.0, 15.0).sketch(\"block\").extrude(blind(15.0));\nlet pocket = m.tool(|m| m.profile(xy_at(5.0)).rect(24.0, 14.0).at(20.0, 15.0).sketch(\"pocket\").extrude(blind(11.0)));\nb = b.subtract(pocket);\nb = b.draft(faces().at_height(10.0), 5.0, xy_at(5.0));\nm\n",
    );
    let t = (5.0f64.to_radians()).tan();
    let h = 10.0;
    let exact = 40.0 * 30.0 * 15.0
        - (24.0 * 14.0 * h + 76.0 * t * h * h / 2.0 + 4.0 * t * t * h * h * h / 3.0);
    assert!((v - exact).abs() < 1e-9 * exact, "{v} vs {exact}");
}

/// brep-bench draft/box-top-drafted-away: a 30 × 10 × 10 box's sides
/// drafted 30° meet at 5√3 < 10, consuming the top. That face's loop
/// turned inside out and the body served unclosed (1135 for 1154.7);
/// then it declined naming the face; now the long sides meet in a ridge
/// at height H = b / 2t and the drafted ends meet the ridge (PK-1159):
/// a hip roof, sections (a − 2zt)(b − 2zt) to H. At 20° they meet at
/// 13.7, above the top, and the draft serves the frustum.
#[test]
fn a_draft_consuming_a_face_serves_the_ridge() {
    let box_drafted = |deg: f64| {
        format!(
            "let m = model();\nlet b = m.profile(xy()).rect(30.0, 10.0).at(15.0, 5.0).sketch(\"box\").extrude(blind(10.0));\nb = b.draft(faces().except(faces().facing(\"z\").or(faces().facing_neg(\"z\"))), {deg:?}, xy());\nm\n"
        )
    };
    let (a, b) = (30.0, 10.0);
    let t = (30.0f64.to_radians()).tan();
    let ridge = b / (2.0 * t);
    let roof = a * b * ridge - (a + b) * t * ridge * ridge + 4.0 / 3.0 * t * t * ridge.powi(3);
    let v = volume(&box_drafted(30.0));
    assert!((v - roof).abs() < 1e-9 * roof, "{v} vs {roof}");
    let v = volume(&box_drafted(20.0));
    // Each side's lean takes a wedge; the four wedges overlap at the
    // corners in pyramids: V = abh − (a + b)h²t + (4/3)h³t².
    let (h, t) = (10.0, (20.0f64.to_radians()).tan());
    let exact = a * b * h - (a + b) * h * h * t + 4.0 / 3.0 * h * h * h * t * t;
    assert!((v - exact).abs() < 1e-9 * exact, "{v} vs {exact}");
}
