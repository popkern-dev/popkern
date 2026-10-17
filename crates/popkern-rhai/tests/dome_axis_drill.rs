//! A revolved dome drilled along its axis, upright and turned (PK-877).
//!
//! The dome is the six-axis arm housing's shape: a cylinder r 1 to
//! z 2 under a spherical cap r 1 to z 3, drilled by an r 0.3 cylinder
//! from z 1. Before PK-877 the upright drill served only through the
//! cellular -> graph-merge fallback and the turned one declined. Three
//! defects were stacked:
//! - the ray classifier had no trim test for sphere faces, so a ray
//!   from inside the turned dome's wall crossed the sphere's untrimmed
//!   lower half and a wall fragment read Outside;
//! - the pole sample of the seam piece running from the pole to the rim
//!   kept an arbitrary azimuth (u = −π on a u = 0 meridian), so the
//!   sphere walked on the 3D-Newell fallback, which wired the closed rim
//!   backwards: the sphere zone, the hole wall and the hole bottom
//!   served with honest-looking flags on inside-out faces (OCC read the
//!   exported body as 4.02 for 7.82);
//! - the mesher's incircle conditioning (collinear nudge, cocircular
//!   jitter) stayed in the served UV, and the volume quadrature
//!   integrated the sphere zone 1.3e-8 rad off its exact rim.
//!
//! Every check here is one the inside-out serve passed or one that
//! catches it directly: no fallback, `validate_all` clean, watertight,
//! the analytic volume to 1e-9, every use pcurve imaging onto its own
//! edge in its traversal direction, latitude-exact rims on the sphere,
//! and (when the OCP venv is present) OpenCASCADE's volume of the
//! exported STEP.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::PcurveGeometry;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::BRep;
use popkern_rhai::ScriptHost;

/// The dome minus the drilled core: 8π/3 less the core's cylinder part
/// (π·0.3²·1) and its spherical part ((2π/3)(1 − 0.91^{3/2})).
fn analytic_volume() -> f64 {
    let pi = std::f64::consts::PI;
    8.0 * pi / 3.0 - pi * 0.09 - (2.0 * pi / 3.0) * (1.0 - 0.91_f64.powf(1.5))
}

/// The turn the second case applies to dome and tool alike.
const TURN: &str = ".rotate(axis_at([0.3, -0.2, 0.5], [0.6, 0.0, 0.8]), 37.0)";

fn script(turn: &str, drilled: bool) -> String {
    let cut = if drilled {
        "body = body.subtract(hole);"
    } else {
        ""
    };
    format!(
        r#"let m = model();
let body = m.profile(xz()).path(0.0, 0.0).line_to(1.0, 0.0).line_to(1.0, 2.0)
    .arc_to([0.70710678118654752, 2.70710678118654752], [0.0, 3.0]).close()
    .sketch("dome").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle()){turn};
let hole = m.tool(|t| t.profile(xy_at(1.0)).circle(0.3).sketch("hole").extrude(blind(3.0)){turn});
{cut}
m
"#
    )
}

/// The largest distance between a use pcurve's image and its edge
/// curve at matched fractions of the half-edge's traversal (a closed
/// edge's canonical half-edge runs its curve forward, the twin
/// backward).
fn worst_pcurve_image_miss(brep: &BRep) -> f64 {
    let mut worst = 0.0_f64;
    for (he, h) in &brep.half_edges {
        let (Some(pc), Some(face)) = (
            h.pcurve.as_ref(),
            brep.loops.get(h.loop_).and_then(|l| brep.faces.get(l.face)),
        ) else {
            continue;
        };
        let Some(surface) = face.surface.as_ref() else {
            continue;
        };
        let edge = &brep.edges[h.edge];
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let start = curve.geometry.point_at(curve.t_start);
        let end = curve.geometry.point_at(curve.t_end);
        let origin = brep.vertices[h.origin].position;
        let forward = if start.distance_to(&end) < 1e-9 {
            edge.half_edge == he
        } else {
            origin.distance_to(&start) <= origin.distance_to(&end)
        };
        let (ta, tb) = if forward {
            (curve.t_start, curve.t_end)
        } else {
            (curve.t_end, curve.t_start)
        };
        let (p0, p1) = pc.parameter_range();
        for k in 1..8 {
            let f = f64::from(k) / 8.0;
            let q = pc.point_at(p0 + (p1 - p0) * f);
            let image = surface.geometry.point_at(q.x, q.y);
            worst = worst.max(image.distance_to(&curve.geometry.point_at(ta + (tb - ta) * f)));
        }
    }
    worst
}

/// Every use on the sphere face is a straight chart line (latitude
/// rims and meridian seams), not a fit.
fn sphere_trims_are_chart_lines(brep: &BRep) -> bool {
    brep.half_edges.values().all(|h| {
        let on_sphere = brep
            .loops
            .get(h.loop_)
            .and_then(|l| brep.faces.get(l.face))
            .and_then(|f| f.surface.as_ref())
            .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Sphere { .. }));
        !on_sphere
            || h.pcurve.as_ref().is_some_and(|pc| match &pc.geometry {
                PcurveGeometry::Linear { .. } => true,
                PcurveGeometry::Nurbs(n) => n.degree == 1 && n.control_points.len() == 2,
                _ => false,
            })
    })
}

/// OpenCASCADE's adaptive volume of the one solid in `path`, or `None`
/// when the OCP venv is absent.
fn occ_volume(path: &std::path::Path) -> Option<f64> {
    let python = std::path::Path::new(&std::env::var("HOME").unwrap_or_default())
        .join("venvs/ocp/bin/python");
    let script = std::path::Path::new("../../tools/occ-oracle/occ_solids.py");
    if !python.exists() || !script.exists() {
        eprintln!("[dome-axis-drill] OCP oracle unavailable — OCC leg skipped");
        return None;
    }
    let output = std::process::Command::new(&python)
        .arg(script)
        .arg(path)
        .arg("--adaptive")
        .output()
        .expect("oracle invocation");
    let report = String::from_utf8_lossy(&output.stdout);
    let volumes: Vec<f64> = report
        .split("\"volume\": ")
        .skip(1)
        .filter_map(|s| s.split(['}', ',']).next()?.trim().parse().ok())
        .collect();
    assert_eq!(volumes.len(), 1, "one solid in the export: {report}");
    Some(volumes[0])
}

#[test]
fn the_drilled_dome_serves_exactly_upright_and_turned() {
    let out_dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("dome_axis_drill");
    std::fs::create_dir_all(&out_dir).expect("output dir");
    for (label, turn) in [("upright", ""), ("turned", TURN)] {
        popkern_bool::merge_fallback::take_merge_fallback_events();
        popkern_bool::chart_walk_fallback::take_chart_walk_fallback_events();
        let mut run = ScriptHost::new()
            .run(&script(turn, true), "dome.rhai")
            .unwrap_or_else(|e| panic!("{label}: {}", e.render()));
        let mut model = run.model.take().expect("model");
        let path = out_dir.join(format!("{label}.step"));
        model.to_step(&path).expect("STEP export");
        let parts = model.to_parts_all().expect("parts");
        assert_eq!(parts.len(), 1, "{label}: one body");
        let brep = parts[0].brep();

        let merge = popkern_bool::merge_fallback::take_merge_fallback_events();
        assert!(merge.is_empty(), "{label}: graph-merge fallback {merge:?}");
        let walk = popkern_bool::chart_walk_fallback::take_chart_walk_fallback_events();
        assert!(walk.is_empty(), "{label}: 3D-Newell walk fallback {walk:?}");

        assert_eq!(brep.faces.len(), 5, "{label}: faces");
        let errors = brep.validate_all();
        assert!(errors.is_empty(), "{label}: {errors:?}");
        assert_eq!(parts[0].open_edge_count().expect("mesh"), 0, "{label}");

        // A use run against its edge images a chord away (2.0 on the
        // reversed equator); honest uses sit at roundoff.
        let miss = worst_pcurve_image_miss(brep);
        let gate = brep.tolerance.linear();
        assert!(
            miss < gate,
            "{label}: a use pcurve images {miss:e} off its edge"
        );
        assert!(
            sphere_trims_are_chart_lines(brep),
            "{label}: sphere trim fit"
        );

        let volume = popkern_measure::body_volume(brep, &brep.tolerance).expect("volume");
        let exact = analytic_volume();
        assert!(
            (volume - exact).abs() < 1e-9,
            "{label}: volume {volume} vs analytic {exact}"
        );
        if let Some(occ) = occ_volume(&path) {
            assert!(
                (occ - exact).abs() < 1e-9,
                "{label}: OCC volume {occ} vs analytic {exact} (inside-out faces?)"
            );
        }
    }
}

/// `TURN` applied to a dome-frame point: 37° about (0.6, 0, 0.8)
/// through (0.3, −0.2, 0.5), right-handed (Rodrigues).
fn turned(p: [f64; 3]) -> popkern_core::geometry::Point3 {
    use popkern_core::geometry::PortableMath;
    let n = (0.6_f64 * 0.6 + 0.8 * 0.8).sqrt();
    let k = [0.6 / n, 0.0, 0.8 / n];
    let o = [0.3, -0.2, 0.5];
    let v = [p[0] - o[0], p[1] - o[1], p[2] - o[2]];
    let (s, c) = (37.0_f64.to_radians().psin(), 37.0_f64.to_radians().pcos());
    let kxv = [
        k[1] * v[2] - k[2] * v[1],
        k[2] * v[0] - k[0] * v[2],
        k[0] * v[1] - k[1] * v[0],
    ];
    let kv = k[0] * v[0] + k[1] * v[1] + k[2] * v[2];
    let r: Vec<f64> = (0..3)
        .map(|i| o[i] + v[i] * c + kxv[i] * s + k[i] * kv * (1.0 - c))
        .collect();
    popkern_core::geometry::Point3::new(r[0], r[1], r[2])
}

/// Point classification against the turned dome alone (PK-877): the
/// sphere face's solver window is its whole chart, and before the
/// sphere trim gate a ray from inside the wall below the equator also
/// crossed the sphere's untrimmed lower half next to its exit through
/// the wall (measured: (−0.21, −0.74, 1.56), at r 0.30, z 1.63 in the
/// dome's frame, voted Outside). Upright the face's box kept those rays
/// off the sphere; turned it does not.
#[test]
fn points_in_the_turned_dome_classify_through_its_sphere_trim() {
    use popkern_bool::error::PointClassification;
    for (label, turn) in [("upright", ""), ("turned", TURN)] {
        let mut run = ScriptHost::new()
            .run(&script(turn, false), "dome.rhai")
            .unwrap_or_else(|e| panic!("{label}: {}", e.render()));
        let parts = run
            .model
            .take()
            .expect("model")
            .to_parts_all()
            .expect("parts");
        let brep = parts[0].brep();
        let boxes = popkern_bool::adaptive_precision::body_face_aabbs(brep);
        let place = |p: [f64; 3]| {
            if turn.is_empty() {
                popkern_core::geometry::Point3::new(p[0], p[1], p[2])
            } else {
                turned(p)
            }
        };
        let mut wrong = Vec::new();
        for k in 0..8 {
            let phi = f64::from(k) * std::f64::consts::FRAC_PI_4 + 0.1;
            let (s, c) = (phi.sin(), phi.cos());
            for z in [0.5, 1.0, 1.5, 1.63, 1.8, 1.95] {
                for (r, expect) in [
                    (0.3, PointClassification::Inside),
                    (0.6, PointClassification::Inside),
                    (0.9, PointClassification::Inside),
                    (1.05, PointClassification::Outside),
                ] {
                    let p = place([r * c, r * s, z]);
                    let got = popkern_bool::adaptive_precision::robust_point_in_solid_final(
                        &p,
                        brep,
                        &boxes,
                        &brep.tolerance,
                    );
                    if got != expect {
                        wrong.push((r, z, phi, got));
                    }
                }
            }
        }
        assert!(wrong.is_empty(), "{label}: misclassified {wrong:?}");
    }
}
