//! A tool face kept by a boolean keeps its holes (PK-755).
//!
//! Subtracting a revolved cavity with a step leaves the step's annulus
//! (a face of the TOOL, kept and flipped) in the result. The loop
//! election filed any loop at least half the outer loop's area as a
//! separate face and skipped the winding test on the tool operand, so a
//! thin annulus came out as two overlapping discs: validate() passed,
//! the mesh filled the hole, and an interference check reported the
//! cycloid housing colliding with its own output flange.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::BRep;
use popkern_rhai::ScriptHost;

fn served(source: &str) -> BRep {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(source, "annulus.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    let mut model = outcome.model.expect("model");
    model.to_part().expect("part").brep().clone()
}

/// The planar faces lying in z = `z`, as (outer loop vertex radii,
/// inner loop count).
fn faces_at(brep: &BRep, z: f64) -> Vec<(Vec<f64>, usize)> {
    brep.faces
        .iter()
        .filter(|(_, f)| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Plane { .. })
            )
        })
        .filter_map(|(h, f)| {
            let pts = brep.face_outer_loop_positions(h);
            (!pts.is_empty() && pts.iter().all(|p| (p.z - z).abs() < 1e-6)).then(|| {
                (
                    pts.iter().map(|p| p.x.hypot(p.y)).collect(),
                    f.inner_loops.len(),
                )
            })
        })
        .collect()
}

/// A cylinder r 40 minus a revolved cavity whose step at z = 30 runs
/// from radius `outer` in to `inner`.
fn stepped(outer: f64, inner: f64) -> String {
    format!(
        r#"let m = model();
let body = m.profile(xy()).circle(40.0).sketch("solid").extrude(blind(50.0));
let cavity = m.tool(|t| {{
    let p = t.profile(xz()).path(0.0, 5.0).line_to({outer}, 5.0).line_to({outer}, 30.0);
    p = p.line_to({inner}, 30.0).line_to({inner}, 51.0).line_to(0.0, 51.0).close();
    p.sketch("cavity").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle())
}});
body.subtract(cavity, #{{ label: "cavity" }});
m"#
    )
}

#[test]
fn a_kept_tool_annulus_stays_one_face_with_its_hole() {
    // Inner/outer area ratios 0.25 (always worked), 0.71 (the cycloid
    // housing's step) and 0.98 (a hairline step).
    for (outer, inner) in [(30.0, 15.0), (32.0, 27.0), (32.0, 31.7)] {
        let brep = served(&stepped(outer, inner));
        assert!(brep.validate().is_ok(), "{outer}/{inner}: invalid");
        let step = faces_at(&brep, 30.0);
        assert_eq!(step.len(), 1, "{outer}/{inner}: step faces {step:?}");
        assert_eq!(
            step[0].1, 1,
            "{outer}/{inner}: the step lost its hole: {step:?}"
        );
        // Cylinder less the cavity (its bore ends at the top, z 50):
        // π(40²·50 − outer²·25 − inner²·20).
        let exact =
            std::f64::consts::PI * (1600.0 * 50.0 - outer * outer * 25.0 - inner * inner * 20.0);
        let volume = popkern_measure::body_volume(&brep, &brep.tolerance).expect("volume");
        assert!(
            ((volume - exact) / exact).abs() < 1e-3,
            "{outer}/{inner}: volume {volume} vs {exact}"
        );
    }
}
