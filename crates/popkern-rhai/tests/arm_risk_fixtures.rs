//! The six-axis arm's kernel risk fixtures as standalone tests
//! (six-axis-arm §9.1, PK-717), built from the gallery script's own
//! part functions so they test exactly what the arm builds:
//!
//! - **R2**, the biarc cycloid disc: the two caps, the seven bores as
//!   cylinders, the flank as exact arcs (tangent-continuous runs share
//!   a face: 464 arcs in 146 faces), and the built flank within b/4 of
//!   the exact equidistant epitrochoid — measured on the B-rep's edges,
//!   not taken from the script's own fit report.
//! - **R4**, the 64-window encoder wheel: watertight, one planar face
//!   per polygon side, its volume the analytic one.
//!
//! R1 (tees) lives in `cylinder_tee_unions` and `fillet_tangent_chains`,
//! R3 (finned housing) in `finned_housing`; R5 (tapped holes) is timed
//! by popkern-part's `tapped_housing_taps_cost_what_the_first_did`
//! (twelve M5 taps ~13 s, 2026-10-14, against its 20 s budget). The arm
//! threads its tool flange's ISO M6 holes and keeps tap-drill bores on
//! its ~150 bolt-circle holes (PK-801, design §6).

use popkern_core::geometry::curve::CurveEval as _;
use popkern_rhai::ScriptHost;

const SCRIPT: &str = include_str!("../examples/cycloidal_joint.rhai");

/// The script's layout, flank fit and part-body functions.
fn part_functions() -> String {
    let slice = |from: &str, to: &str| -> &str {
        let a = SCRIPT.find(from).expect("section start");
        let b = SCRIPT[a..].find(to).expect("section end") + a;
        &SCRIPT[a..b]
    };
    [
        slice("// ---- Derived layout", "// ---- Part bodies"),
        slice("fn rod(t, z", "// A shaft with the magnet hub"),
        slice("// A 64-window wheel", "// Cross-roller races"),
    ]
    .concat()
}

/// The L module's size row.
const SIZE_L: &str = r#"let pi = 3.141592653589793;
let s = #{ name: "L", od: 110.0, len: 118.0, pcd: 84.0, bolts: 12, tap: 4.2, screw: "M5",
    n: 30, big_r: 40.0, rp: 2.5, e: 0.9, b: 0.02, disc_t: 6.0, cam_r: 11.0, rop: 26.0,
    op_r: 3.0, shaft_r: 8.0, d: 5.0 };
let g = layout(s);"#;

fn run(body: &str) -> (popkern_script::Model, Vec<String>) {
    let source = format!("{}\n{SIZE_L}\n{body}", part_functions());
    let mut host = ScriptHost::new();
    let outcome = host
        .run(&source, "risk.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    (outcome.model.expect("model"), outcome.printed)
}

/// The equidistant epitrochoid (the script's `cyc_pt`, ported).
fn cycloid(t: f64) -> (f64, f64) {
    let (big_r, rr, e, n) = (40.0, 2.5 + 0.02, 0.9, 30.0);
    let k = 1.0 - n;
    let psi = (k * t).sin().atan2(big_r / (e * n) - (k * t).cos());
    (
        big_r * t.cos() - rr * (t + psi).cos() - e * (n * t).cos(),
        -big_r * t.sin() + rr * (t + psi).sin() + e * (n * t).sin(),
    )
}

#[test]
fn r2_the_biarc_disc_holds_the_cycloid_to_a_quarter_backlash() {
    let (mut model, printed) = run(
        "let fit = cycloid_fit(s, pi);\nprint(fit[0].len());\nlet m = model();\n\
         disc_body(m, s, fit, pi);\nm.assert_watertight();\nm",
    );
    let arcs: usize = printed[0].trim().parse().expect("arc count");
    let part = model.to_part().expect("part");
    let brep = part.brep();
    // The seven bores are cylinders at their radii (cam bore r 11.05,
    // windows r 3 + 0.9 + 0.05); the flank is exact arcs, a run of
    // tangent-continuous arcs sharing one face (the profile compiles to
    // one exact wire and the extrude breaks faces where tangency
    // breaks): no more flank faces than arcs.
    let bores = brep
        .faces
        .values()
        .filter(|f| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(popkern_core::geometry::surface::SurfaceGeometry::Cylinder { radius, .. })
                    if (radius - 11.05).abs() < 1e-9 || (radius - 3.95).abs() < 1e-9
            )
        })
        .count();
    assert_eq!(bores, 7, "cam bore and six windows");
    let flank = brep.faces.len() - 2 - bores;
    assert!(
        flank >= 1 && flank <= arcs,
        "{flank} flank faces for {arcs} arcs"
    );
    // The exact flank, densely: 2^20 samples (chord sag ~1e-9 mm).
    let n = 1 << 20;
    let profile: Vec<(f64, f64)> = (0..n)
        .map(|i| cycloid(std::f64::consts::TAU * i as f64 / n as f64))
        .collect();
    let mut worst: f64 = 0.0;
    for (_, edge) in brep.edges.iter() {
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let (a, b) = (curve.t_start, curve.t_end);
        for k in 0..=8 {
            let p = curve.geometry.point_at(a + (b - a) * k as f64 / 8.0);
            // The flank's cap edges (the bores lie inside r = 30).
            if p.z.abs() > 1e-9 || p.x.hypot(p.y) < 33.0 {
                continue;
            }
            let d = profile
                .iter()
                .map(|q| (q.0 - p.x).hypot(q.1 - p.y))
                .fold(f64::INFINITY, f64::min);
            worst = worst.max(d);
        }
    }
    eprintln!(
        "[pk717] R2: {arcs} arcs in {flank} faces, worst flank deviation {worst:.2e} mm (b/4 = 5e-3)"
    );
    assert!(
        worst > 0.0 && worst <= 0.02 / 4.0,
        "flank deviation {worst}"
    );
}

#[test]
fn r4_the_encoder_wheel_is_watertight_with_its_analytic_volume() {
    let (mut model, printed) = run(
        "print(g.mag_in);\nlet m = model();\nencoder_body(m, g, s, pi);\n\
         m.assert_watertight();\nm",
    );
    let r1: f64 = printed[0].trim().parse().expect("mag_in");
    let part = model.to_part().expect("part");
    let brep = part.brep();
    // 256 polygon sides, top, bottom and the hub bore.
    assert_eq!(brep.faces.len(), 256 + 3, "faces");
    // The polygon (the script's points), by the shoelace formula.
    let (r0, hw) = (r1 - 3.0, 0.25 * std::f64::consts::TAU / 64.0);
    let mut pts = Vec::new();
    for k in 0..64 {
        let t0 = std::f64::consts::TAU * k as f64 / 64.0;
        for (r, t) in [(r1, t0 - hw), (r0, t0 - hw), (r0, t0 + hw), (r1, t0 + hw)] {
            pts.push((r * t.cos(), r * t.sin()));
        }
    }
    let area: f64 = (0..pts.len())
        .map(|i| {
            let (a, b) = (pts[i], pts[(i + 1) % pts.len()]);
            a.0 * b.1 - b.0 * a.1
        })
        .sum::<f64>()
        .abs()
        / 2.0;
    let bore = std::f64::consts::PI * (8.0_f64 + 0.05).powi(2);
    let want = (area - bore) * 1.5;
    let got = popkern_measure::body_volume(brep, &brep.tolerance).expect("volume");
    assert!(
        (got - want).abs() / want < 1e-4,
        "encoder volume {got} vs {want}"
    );
}
