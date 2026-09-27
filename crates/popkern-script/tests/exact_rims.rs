//! PK-417: a bore cut through a block leaves exact circular rims — the
//! boolean re-spells the plane × cylinder intersection curves as the
//! `Arc`s they are — so the rims select as `circular()`, carry the
//! bore's radius, take a fillet by label, and the body stays
//! watertight and STEP-clean. Union rims (a boss meeting a plate) are
//! the same class and host G2 fillets on the exact arc (PK-419).

use std::path::Path;

use popkern_core::geometry::curve::CurveGeometry;
use popkern_part::prelude::*;
use popkern_script::prelude::*;
use popkern_script::{Model, Plane};

fn bored_block(m: &mut Model) -> Body {
    let block = m
        .profile(Plane::XY)
        .rect(1.2, 0.8)
        .sketch("block_sketch")
        .expect("sketch")
        .extrude(m, End::Blind(0.5.into()))
        .expect("extrude");
    let bore = m
        .tool(|m| {
            m.profile(Plane::xy_at(-0.1))
                .circle(0.18)
                .sketch("bore_sketch")?
                .extrude(m, End::Blind(0.7.into()))
        })
        .expect("bore");
    block.subtract(m, bore).expect("subtract")
}

#[test]
fn a_bore_leaves_two_exact_circular_rims_that_select_and_fillet_by_label() {
    let mut m = Model::new();
    let mut block = bored_block(&mut m);
    let brep = m.brep_of(block).expect("brep");
    let arcs: Vec<f64> = brep
        .edges
        .values()
        .filter_map(|e| match e.curve.as_ref().map(|c| &c.geometry) {
            Some(CurveGeometry::Arc { radius, .. })
            | Some(CurveGeometry::Circle { radius, .. }) => Some(*radius),
            _ => None,
        })
        .collect();
    assert_eq!(arcs.len(), 2, "two rims, both exact: {arcs:?}");
    assert!(arcs.iter().all(|r| (r - 0.18).abs() < 1e-9), "{arcs:?}");
    let nurbs = brep
        .edges
        .values()
        .filter(|e| {
            matches!(
                e.curve.as_ref().map(|c| &c.geometry),
                Some(CurveGeometry::Nurbs(_))
            )
        })
        .count();
    assert_eq!(nurbs, 0, "no NURBS edge is left on a block with a bore");
    let issues = brep.validate_all();
    assert!(issues.is_empty(), "valid after promotion: {issues:?}");
    // The rims select by kind and by height, and take a fillet by label.
    m.label_edge("bore_rim", &edges().circular().at_height(0.5))
        .expect("the top rim is circular");
    block = block
        .fillet(&mut m, edges().labeled("bore_rim"), 0.05)
        .expect("fillet on the labeled rim");
    let _ = block;
    assert_eq!(m.report().failures().count(), 0);
    m.assert_watertight().expect("watertight");
    let step = Path::new(env!("CARGO_TARGET_TMPDIR")).join("bored_block.step");
    m.assert_step_roundtrip(&step).expect("STEP round trip");
    let oracle = Path::new(&std::env::var("HOME").unwrap_or_default()).join("venvs/ocp/bin/python");
    let check = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occ-oracle/occ_check.py");
    if oracle.exists() && check.exists() {
        let output = std::process::Command::new(&oracle)
            .arg(&check)
            .arg(&step)
            .output()
            .expect("oracle");
        let report = String::from_utf8_lossy(&output.stdout);
        assert!(report.contains("\"passed\": true"), "OCC oracle: {report}");
    }
}

/// A union's junction rim (a boss meeting a plate) is the same class
/// as a bore's: exact, selectable as `circular()`, and a fillet host
/// at G1 and at G2. Until PK-419 the union rim stayed NURBS because
/// the fillet's spine walker seeded a closed rim at its midpoint and
/// converged to the antipodal ball, tearing the G2 strip along its
/// rails (70 open mesh edges); the seed now sits at the rim vertex
/// and the walker certifies the station.
#[test]
fn a_union_junction_rim_is_exact_and_hosts_g1_and_g2_fillets() {
    for continuity in [ContinuityKind::G1, ContinuityKind::G2] {
        let mut m = Model::new();
        let plate = m
            .profile(Plane::XY)
            .rect(2.4, 1.6)
            .sketch("plate_sketch")
            .expect("sketch")
            .extrude(&mut m, End::Blind(0.4.into()))
            .expect("extrude");
        let boss = m
            .tool(|m| {
                m.profile(Plane::xy_at(0.2))
                    .circle(0.35)
                    .sketch("boss_sketch")?
                    .extrude(m, End::Blind(0.8.into()))
            })
            .expect("boss");
        let plate = plate.union(&mut m, boss).expect("union");
        let brep = m.brep_of(plate).expect("brep");
        let circular: Vec<f64> = brep
            .edges
            .values()
            .filter_map(|e| match e.curve.as_ref().map(|c| &c.geometry) {
                Some(CurveGeometry::Arc { radius, .. })
                | Some(CurveGeometry::Circle { radius, .. }) => Some(*radius),
                _ => None,
            })
            .collect();
        assert_eq!(
            circular.len(),
            2,
            "the boss top rim and the junction rim are both exact: {circular:?}"
        );
        assert!(circular.iter().all(|r| (r - 0.35).abs() < 1e-9));
        m.label_edge("junction", &edges().circular().at_height(0.4))
            .expect("the junction rim selects as circular");
        let plate = plate
            .fillet_continuity(
                &mut m,
                edges().labeled("junction"),
                0.08,
                continuity,
                "round",
            )
            .expect("fillet on the junction rim");
        let _ = plate;
        assert_eq!(m.report().failures().count(), 0, "{continuity:?}");
        m.assert_watertight()
            .unwrap_or_else(|e| panic!("{continuity:?}: {e:?}"));
        let round = m.feature("round").expect("label");
        let report = m.evidence(&round).expect("fillet evidence");
        let blend = report.blend.as_ref().expect("blend record");
        assert!(
            !blend.downgraded(),
            "{continuity:?}: {:?}",
            blend.continuity
        );
        let step =
            Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("boss_plate_{continuity:?}.step"));
        m.assert_step_roundtrip(&step).expect("STEP round trip");
    }
}

/// A bore rim (subtraction) hosts a G2 fillet on its exact arc too:
/// the walker's seed rule is spelling-agnostic.
#[test]
fn a_bore_rim_hosts_a_g2_fillet_on_its_exact_arc() {
    let mut m = Model::new();
    let block = bored_block(&mut m);
    m.label_edge("bore_rim", &edges().circular().at_height(0.5))
        .expect("the top rim is circular");
    let block = block
        .fillet_continuity(
            &mut m,
            edges().labeled("bore_rim"),
            0.05,
            ContinuityKind::G2,
            "round",
        )
        .expect("G2 fillet on the bore rim");
    let _ = block;
    assert_eq!(m.report().failures().count(), 0);
    m.assert_watertight().expect("watertight");
    let round = m.feature("round").expect("label");
    let report = m.evidence(&round).expect("fillet evidence");
    let blend = report.blend.as_ref().expect("blend record");
    assert!(!blend.downgraded(), "{:?}", blend.continuity);
    assert_eq!(
        blend.continuity.as_ref().map(|c| c.achieved),
        Some(Continuity::G2)
    );
}

/// The union rim's fillet rolls in the VOID: it adds the concave
/// root's material and its G1 surface is the torus OUTSIDE the boss
/// (major radius boss + r). Until the closed-ring embedding fix the
/// boolean shipped the boss wall and top cap inside-out (flag and
/// winding inverted together — winding-honest, still wrong), the
/// rolling ball took the boss interior, and both G1 and G2 carved a
/// groove into the boss (volume fell by ~0.006 instead of rising by
/// ~0.0032; the G1 torus sat at major 0.27 with the rim's trims on
/// the true edges, which rendered as a scalloped ring).
#[test]
fn a_union_rim_fillet_adds_material_on_the_void_side() {
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
    let (boss_r, r) = (0.35, 0.08);
    // Concave-corner fillet: cross-section r²(1 − π/4) swept around
    // its centroid radius (Pappus); centroid offset r(10 − 3π)/(12 − 3π).
    let section = r * r * (1.0 - std::f64::consts::FRAC_PI_4);
    let rho =
        boss_r + r * (10.0 - 3.0 * std::f64::consts::PI) / (12.0 - 3.0 * std::f64::consts::PI);
    let expected_gain = section * std::f64::consts::TAU * rho;
    for continuity in [ContinuityKind::G1, ContinuityKind::G2] {
        let mut m = Model::new();
        let plate = m
            .profile(Plane::XY)
            .rect(2.4, 1.6)
            .sketch("plate_sketch")
            .expect("sketch")
            .extrude(&mut m, End::Blind(0.4.into()))
            .expect("extrude");
        let boss = m
            .tool(|m| {
                m.profile(Plane::xy_at(0.2))
                    .circle(boss_r)
                    .sketch("boss_sketch")?
                    .extrude(m, End::Blind(0.8.into()))
            })
            .expect("boss");
        let plate = plate.union(&mut m, boss).expect("union");
        let united = m.brep_of(plate).expect("brep");
        // Every face of the union faces out of the solid: the boss
        // wall's flag-folded normal points away from its axis, and
        // every planar face points away from the body's interior.
        for (fh, face) in united.faces.iter() {
            let sof = face.surface.as_ref().expect("surface");
            let sign = if sof.orientation { 1.0 } else { -1.0 };
            match &sof.geometry {
                SurfaceGeometry::Cylinder { .. } => {
                    let p = sof.geometry.point_at(0.3, 0.5);
                    let n = sof.geometry.normal_at(0.3, 0.5) * sign;
                    assert!(
                        n.x * p.x + n.y * p.y > 0.0,
                        "{continuity:?}: boss wall {fh:?} faces into the boss"
                    );
                }
                SurfaceGeometry::Plane { origin, normal, .. } => {
                    let n = *normal * sign;
                    let inward_probe = *origin - n * 1.0e-3;
                    // The body occupies |x| ≤ 1.2, |y| ≤ 0.8, 0 ≤ z ≤ 1.0.
                    let inside = inward_probe.x.abs() <= 1.2
                        && inward_probe.y.abs() <= 0.8
                        && (0.0..=1.0).contains(&inward_probe.z);
                    assert!(inside, "{continuity:?}: plane {fh:?} faces into the solid");
                }
                _ => {}
            }
        }
        let tol = popkern_core::tolerance::Tolerance::default();
        let before = popkern_measure::body_volume(&united, &tol).expect("volume");
        m.label_edge("junction", &edges().circular().at_height(0.4))
            .expect("junction rim");
        let plate = plate
            .fillet_continuity(&mut m, edges().labeled("junction"), r, continuity, "round")
            .expect("fillet");
        let brep = m.brep_of(plate).expect("brep");
        let after = popkern_measure::body_volume(&brep, &tol).expect("volume");
        let gain = after - before;
        // G1 is the circular section (volumes are tessellation-
        // measured, hence the band). A G2 blend with the same contact
        // setback hugs the corner closer than the arc (measured 79 %
        // of the circular fill), so it is bounded by the arc instead.
        let (lo, hi) = match continuity {
            ContinuityKind::G1 => (0.85 * expected_gain, 1.15 * expected_gain),
            _ => (0.5 * expected_gain, 1.02 * expected_gain),
        };
        assert!(
            (lo..hi).contains(&gain),
            "{continuity:?}: the root fillet must ADD material in [{lo:.5}, {hi:.5}]; \
             measured {gain:.5} ({before:.5} -> {after:.5})"
        );
        if continuity == ContinuityKind::G1 {
            let majors: Vec<f64> = brep
                .faces
                .values()
                .filter_map(|f| match f.surface.as_ref().map(|s| &s.geometry) {
                    Some(SurfaceGeometry::Torus { major_radius, .. }) => Some(*major_radius),
                    _ => None,
                })
                .collect();
            assert_eq!(majors.len(), 1, "one root torus: {majors:?}");
            assert!(
                (majors[0] - (boss_r + r)).abs() < 1e-9,
                "the G1 root torus sits outside the boss: major {}",
                majors[0]
            );
        }
    }
}
