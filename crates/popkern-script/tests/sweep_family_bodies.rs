//! PK-440: sweep-family verbs recorded through the recorder serve
//! VALID bodies — a helical spring (cylindrical, tapered, left-hand,
//! phased), a pipe along a model edge, a swept profile that keeps its
//! exact curvature, and lofts over sketch-body sections (consumed).
//! These are the first tests the Sweep/Loft/Pipe replay arms ever had.

use popkern_script::body_sweep::HelixParams;
use popkern_script::prelude::*;
use popkern_script::{Plane, selection::edges};

fn z_axis() -> AxisSpec {
    AxisSpec::Literal {
        origin: [0.0, 0.0, 0.0],
        direction: [0.0, 0.0, 1.0],
    }
}

/// A cylindrical coil spring: circular wire section, 6.5 turns.
#[test]
fn helical_spring_serves_watertight() {
    let mut m = Model::new();
    m.profile(Plane::XY)
        .circle(0.04)
        .sketch("wire_section")
        .expect("sketch")
        .helical_sweep(&mut m, z_axis(), HelixParams::new(0.3, 0.12, 6.5))
        .expect("helical sweep");
    let part = m.to_part().expect("spring validates");
    assert!(part.face_count() > 2, "side faces plus two caps");
    m.assert_watertight().expect("watertight spring");
}

/// Taper, phase, and handedness all reach the kernel.
#[test]
fn tapered_left_hand_phased_spring_serves() {
    let mut m = Model::new();
    let helix = HelixParams {
        taper_to: Some(Scalar::from(0.15)),
        start_angle: Scalar::from(90.0),
        left_handed: true,
        ..HelixParams::new(0.3, 0.15, 4.0)
    };
    m.profile(Plane::XY)
        .circle(0.04)
        .sketch("wire_section")
        .expect("sketch")
        .helical_sweep(&mut m, z_axis(), helix)
        .expect("tapered helical sweep");
    m.assert_watertight().expect("watertight tapered spring");
    let brep = m.to_brep().expect("brep");
    // The coil tightens toward the top: max radial extent near the top
    // quarter is smaller than near the bottom quarter.
    let mut bottom_r: f64 = 0.0;
    let mut top_r: f64 = 0.0;
    let z_top = 4.0 * 0.15;
    for v in brep.vertices.values() {
        let r = (v.position.x * v.position.x + v.position.y * v.position.y).sqrt();
        if v.position.z < 0.25 * z_top {
            bottom_r = bottom_r.max(r);
        } else if v.position.z > 0.75 * z_top {
            top_r = top_r.max(r);
        }
    }
    assert!(
        top_r < bottom_r - 0.05,
        "taper must tighten the coil: bottom {bottom_r:.3} top {top_r:.3}"
    );
}

/// Two half-turn arcs phased 180° apart occupy opposite half-planes —
/// the phase parameter genuinely places the helix start.
#[test]
fn phase_places_the_helix_start() {
    let radial_extents = |phase: f64| -> (f64, f64) {
        let mut m = Model::new();
        m.profile(Plane::XY)
            .circle(0.04)
            .sketch("wire_section")
            .expect("sketch")
            .helical_sweep(
                &mut m,
                z_axis(),
                HelixParams {
                    start_angle: Scalar::from(phase),
                    ..HelixParams::new(0.3, 0.4, 0.5)
                },
            )
            .expect("helical sweep");
        let brep = m.to_brep().expect("brep");
        let xs: Vec<f64> = brep.vertices.values().map(|v| v.position.x).collect();
        let max = xs.iter().cloned().fold(f64::MIN, f64::max);
        let min = xs.iter().cloned().fold(f64::MAX, f64::min);
        (min, max)
    };
    // Phase 0: a half turn from +x through +y to -x (x spans both signs).
    // Phase 180: from -x through -y back to +x.
    let (min0, max0) = radial_extents(0.0);
    let (min180, max180) = radial_extents(180.0);
    assert!(
        (max0 - max180).abs() < 0.1 && (min0 - min180).abs() < 0.1,
        "half-turn extents should mirror: {min0:.2}..{max0:.2} vs {min180:.2}..{max180:.2}"
    );
    // The start cap sits at +x for phase 0 and -x for phase 180: check
    // via the y extents of the two bodies' first-quarter region.
    let y_signs = |phase: f64| -> f64 {
        let mut m = Model::new();
        m.profile(Plane::XY)
            .circle(0.04)
            .sketch("wire_section")
            .expect("sketch")
            .helical_sweep(
                &mut m,
                z_axis(),
                HelixParams {
                    start_angle: Scalar::from(phase),
                    ..HelixParams::new(0.3, 0.4, 0.25)
                },
            )
            .expect("helical sweep");
        let brep = m.to_brep().expect("brep");
        brep.vertices.values().map(|v| v.position.y).sum::<f64>()
            / brep.vertices.len().max(1) as f64
    };
    assert!(
        y_signs(0.0) > 0.0,
        "quarter turn from +x lives in the +y quadrant"
    );
    assert!(
        y_signs(180.0) < 0.0,
        "quarter turn from -x lives in the -y quadrant"
    );
}

/// An analytic circle profile sweeps its exact wire: the spring has a
/// handful of smooth side faces, not a polygonal tube (PK-408 applied
/// to sweep).
#[test]
fn swept_circle_keeps_exact_curvature() {
    let mut m = Model::new();
    m.profile(Plane::XY)
        .circle(0.04)
        .sketch("wire_section")
        .expect("sketch")
        .helical_sweep(&mut m, z_axis(), HelixParams::new(0.3, 0.2, 2.0))
        .expect("helical sweep");
    let part = m.to_part().expect("spring validates");
    // Exact wire: 4 quarter-arc spans → 4 side faces + 2 caps. A
    // polyline profile would shatter into dozens of facets.
    assert!(
        part.face_count() <= 8,
        "expected exact-arc side faces, got {} faces (polyline profile?)",
        part.face_count()
    );
}

/// A pipe along an edge of the model: tube around the top rim circle
/// of a cylinder.
#[test]
fn pipe_along_model_edge_serves() {
    let mut m = Model::new();
    let disc = m
        .profile(Plane::XY)
        .circle(0.5)
        .sketch("disc_sketch")
        .expect("sketch")
        .extrude(&mut m, End::Blind(Scalar::from(0.3)))
        .expect("extrude");
    let tube = disc
        .pipe(&mut m, edges().circular().at_height(0.3), 0.06)
        .expect("pipe");
    // The tube is its own body; the path body is referenced, not
    // consumed (a pipe must not replace the solid it borrowed an edge
    // from).
    assert!(disc.is_alive(&m), "path body survives the pipe");
    assert!(tube.is_alive(&m), "tube body serves");
    assert_ne!(disc.id, tube.id, "the tube is a new body");
    let part = m.to_part().expect("pipe validates");
    assert!(part.face_count() >= 1);
    m.assert_watertight().expect("watertight pipe ring");
}

/// A twisted square blade lofted through three polygon sections on
/// stacked planes; the section sketch bodies are consumed.
#[test]
fn loft_through_three_polygon_sections_serves() {
    let mut m = Model::new();
    let rotated_square = |m: &mut Model, z: f64, angle: f64, half: f64, label: &str| {
        let pts: Vec<(f64, f64)> = (0..4)
            .map(|i| {
                let a = angle + std::f64::consts::FRAC_PI_2 * i as f64;
                (half * a.cos(), half * a.sin())
            })
            .collect();
        let plane = if z == 0.0 { Plane::XY } else { Plane::xy_at(z) };
        m.profile(plane)
            .polygon(&pts)
            .sketch(label)
            .expect("section sketch")
    };
    let s0 = rotated_square(&mut m, 0.0, 0.0, 0.20, "root");
    let s1 = rotated_square(&mut m, 0.15, 0.35, 0.16, "mid");
    let s2 = rotated_square(&mut m, 0.30, 0.70, 0.10, "tip");
    let blade = m.loft(&[s0, s1, s2]).expect("loft");
    assert!(blade.is_alive(&m), "the loft serves on the first section");
    assert!(
        !s1.is_alive(&m) && !s2.is_alive(&m),
        "section sketch bodies are consumed by the loft"
    );
    m.assert_watertight().expect("watertight lofted blade");
}

/// Circle sections loft their exact wires: few smooth faces, not a
/// faceted barrel.
#[test]
fn loft_circle_sections_stay_exact() {
    let mut m = Model::new();
    let s0 = m
        .profile(Plane::XY)
        .circle(0.30)
        .sketch("base_ring")
        .expect("sketch");
    let s1 = m
        .profile(Plane::xy_at(0.25))
        .circle(0.18)
        .sketch("waist_ring")
        .expect("sketch");
    let s2 = m
        .profile(Plane::xy_at(0.5))
        .circle(0.26)
        .sketch("lip_ring")
        .expect("sketch");
    m.loft(&[s0, s1, s2]).expect("loft");
    let part = m.to_part().expect("lofted vase validates");
    assert!(
        part.face_count() <= 8,
        "expected exact-arc loft faces, got {} (polyline sections?)",
        part.face_count()
    );
    m.assert_watertight().expect("watertight loft");
}

/// Fewer than two sections is a typed argument error, transactionally
/// declined (the document is unchanged).
#[test]
fn loft_with_one_section_is_declined() {
    let mut m = Model::new();
    let s0 = m
        .profile(Plane::XY)
        .circle(0.3)
        .sketch("only_section")
        .expect("sketch");
    let before = m.bodies().len();
    let err = m.loft(&[s0]).expect_err("one section must decline");
    assert!(matches!(err, ScriptError::Argument { .. }), "{err:?}");
    assert_eq!(m.bodies().len(), before, "decline leaves the model as-is");
}

/// Sweep with twist along an open arc path. The path body is a partial
/// revolve whose outer top arc is isolated by height plus chord length
/// (a 45° arc at radius r has chord ~0.765·r, so inner and outer arcs
/// differ).
#[test]
fn sweep_along_arc_edge_with_twist_serves() {
    let build = |twist: f64| -> Result<usize, ScriptError> {
        let mut m = Model::new();
        let guide = m
            .profile(Plane::XZ)
            .rect(0.2, 0.1)
            .at(0.7, 0.05)
            .sketch("guide_section")?
            .revolve(
                &mut m,
                AxisDesc::Literal {
                    origin: [0.0, 0.0, 0.0],
                    direction: [0.0, 0.0, 1.0],
                },
                RevolveEnd::MidPlane(Scalar::from(45.0)),
                "guide_ring",
            )?;
        let sq = m
            .profile(Plane::XY)
            .rect(0.1, 0.1)
            .sketch("blade_section")?;
        sq.sweep_with(
            &mut m,
            guide,
            edges().at_height(0.1).chord_length_between(0.55, 0.7),
            twist,
            OperationMode::NewBody,
            "swept",
        )?;
        m.assert_watertight()?;
        Ok(m.to_brep()?.faces.len())
    };
    let plain = build(0.0).expect("untwisted sweep serves");
    let twisted = build(90.0).expect("twisted sweep serves");
    assert!(plain >= 1 && twisted >= 1);
}
