//! Point classification beside a bearing race's groove (PK-431).
//!
//! A revolved groove arc is an exact torus face, full in u and trimmed
//! in v. Its chart cell does not resolve, and the ray counter poisoned
//! every ray that crossed it, so points beside the ball_bearing gallery
//! script's outer-race groove classified OnBoundary on all six
//! directions. The u-periodic chart-trim parity now decides those hits.

use popkern_bool::error::PointClassification;
use popkern_core::geometry::{Point3, PortableMath};

const PREAMBLE: &str = "let m = model();
let pitch = 0.7; let groove_r = 0.16; let lip = sqrt(groove_r * groove_r - 0.12 * 0.12);";

/// The gallery's chamfered inner race.
const INNER: &str = r#"
let race = m.profile(xz())
    .path(0.35, 0.0).line_to(0.58, 0.0).line_to(0.58, 0.2 - lip)
    .arc_to([pitch - groove_r, 0.2], [0.58, 0.2 + lip])
    .line_to(0.58, 0.4).line_to(0.35, 0.4).close()
    .sketch("s").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
race.chamfer(edges().circular().at_height(0.0).or(edges().circular().at_height(0.4)), 0.02);
m"#;

/// The gallery's chamfered outer race.
const OUTER: &str = r#"
let race = m.profile(xz())
    .path(0.82, 0.2 - lip).line_to(0.82, 0.0).line_to(1.05, 0.0)
    .line_to(1.05, 0.4).line_to(0.82, 0.4).line_to(0.82, 0.2 + lip)
    .arc_to([pitch + groove_r, 0.2], [0.82, 0.2 - lip]).close()
    .sketch("o").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
race.chamfer(edges().circular().at_height(0.0).or(edges().circular().at_height(0.4)), 0.02);
m"#;

const Z: [f64; 9] = [0.05, 0.1, 0.12, 0.15, 0.2, 0.25, 0.28, 0.3, 0.35];

/// Classify a grid of (r, z) section points at eight azimuths; returns
/// (correct, wrong, inconclusive).
fn classify_grid(
    body: &str,
    radii: &[f64],
    inside: impl Fn(f64, f64) -> bool,
) -> (usize, usize, usize) {
    let src = format!("{PREAMBLE}{body}");
    let mut run = popkern_rhai::ScriptHost::new()
        .run(&src, "race.rhai")
        .expect("race script");
    let parts = run
        .model
        .take()
        .expect("model")
        .to_parts_all()
        .expect("parts");
    let b = parts[0].brep();
    assert!(
        b.faces.values().any(|f| f
            .surface
            .as_ref()
            .is_some_and(|s| s.geometry.kind_name() == "Torus")),
        "the groove is an exact torus face"
    );
    let boxes = popkern_bool::adaptive_precision::body_face_aabbs(b);
    let (mut ok, mut wrong, mut inconclusive) = (0, 0, 0);
    for k in 0..8 {
        let phi = k as f64 * 0.785 + 0.13;
        for &r in radii {
            for z in Z {
                let p = Point3::new(r * phi.pcos(), r * phi.psin(), z);
                match popkern_bool::point_in_solid_final(&p, b, &boxes, &b.tolerance) {
                    Ok(PointClassification::OnBoundary) | Err(_) => inconclusive += 1,
                    Ok(c) if (c == PointClassification::Inside) == inside(r, z) => ok += 1,
                    Ok(_) => wrong += 1,
                }
            }
        }
    }
    (ok, wrong, inconclusive)
}

/// Inside the race's section, away from the chamfered faces.
fn outside_groove(r: f64, z: f64) -> bool {
    ((r - 0.7) * (r - 0.7) + (z - 0.2) * (z - 0.2)).sqrt() > 0.16 && z > 0.02 && z < 0.38
}

#[test]
fn points_beside_the_outer_race_groove_classify() {
    let radii = [0.7, 0.8, 0.815, 0.83, 0.845, 0.855, 0.87, 0.95, 1.0];
    let got = classify_grid(OUTER, &radii, |r, z| {
        r > 0.82 && r < 1.05 && outside_groove(r, z)
    });
    assert_eq!(got, (648, 0, 0), "(correct, wrong, inconclusive)");
}

#[test]
fn points_beside_the_inner_race_groove_classify() {
    let radii = [0.4, 0.5, 0.545, 0.556, 0.57, 0.575, 0.585, 0.6, 0.7];
    let got = classify_grid(INNER, &radii, |r, z| r < 0.58 && outside_groove(r, z));
    assert_eq!(got, (648, 0, 0), "(correct, wrong, inconclusive)");
}
