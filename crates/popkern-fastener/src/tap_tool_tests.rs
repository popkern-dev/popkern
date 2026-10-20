// This file is only compiled under #[cfg(test)] (see tap_tool.rs).
use crate::holes::tapped_hole;
use crate::{HoleDepth, TapOptions, ThreadSpec};
use popkern_bool::{PointClassification, face_bounding_boxes, point_in_solid_final};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// A block `t` thick whose entry face is z = t, tapped down -z at
/// (x, 0).
fn block(d: f64, t: f64) -> BRep {
    let h = 3.0 * d;
    let pts = [
        Point3::new(-h, -h, 0.0),
        Point3::new(4.0 * h, -h, 0.0),
        Point3::new(4.0 * h, h, 0.0),
        Point3::new(-h, h, 0.0),
    ];
    BRep::prism(&pts, Vec3::new(0.0, 0.0, 1.0), t, Tolerance::default()).expect("block")
}

fn tap(body: &BRep, spec: &ThreadSpec, x: f64, t: f64, depth: HoleDepth) -> BRep {
    tapped_hole(
        body,
        spec,
        Point3::new(x, 0.0, t),
        Vec3::new(0.0, 0.0, -1.0),
        depth,
        &TapOptions::default(),
        Tolerance::default(),
    )
    .expect("tapped hole")
}

/// Points around the hole at (x, 0) where the two bodies classify
/// differently (on-boundary samples excluded).
fn disagreements(a: &BRep, b: &BRep, d: f64, x: f64, t: f64) -> usize {
    let (ba, bb) = (face_bounding_boxes(a), face_bounding_boxes(b));
    let tol = Tolerance::default();
    let r = 1.2 * d;
    let mut differ = 0;
    for i in 0..14 {
        for j in 0..14 {
            for k in 0..22 {
                // Off-lattice offsets keep samples off the thread's
                // symmetric planes.
                let p = Point3::new(
                    x - r + 2.0 * r * (i as f64 + 0.37) / 14.0,
                    -r + 2.0 * r * (j as f64 + 0.61) / 14.0,
                    t + 0.3 * d - (t + 0.6 * d) * (k as f64 + 0.43) / 22.0,
                );
                let ca = point_in_solid_final(&p, a, &ba, &tol).expect("classify a");
                let cb = point_in_solid_final(&p, b, &bb, &tol).expect("classify b");
                let boundary = matches!(ca, PointClassification::OnBoundary)
                    || matches!(cb, PointClassification::OnBoundary);
                if !boundary && ca != cb {
                    differ += 1;
                }
            }
        }
    }
    differ
}

/// The one-tool cut (PK-876) removes what the bore, thread and chamfer
/// cuts remove in sequence: blind and through, three sizes, with the
/// cached tool's lid and through-reach rounded up past the material.
#[test]
fn the_one_tool_cut_matches_the_sequential_cuts() {
    for designation in ["M3", "M5", "M8"] {
        let spec = ThreadSpec::parse(designation).expect("spec");
        let d = spec.major_diameter();
        let t = 2.6 * d;
        for depth in [HoleDepth::Blind(1.7 * d), HoleDepth::Through] {
            let _ = super::take_tap_tool_fallbacks();
            let one = tap(&block(d, t), &spec, 0.0, t, depth);
            let fell_back = super::take_tap_tool_fallbacks();
            assert!(
                fell_back.is_empty(),
                "{designation} {depth:?}: {fell_back:?}"
            );
            super::SEQUENTIAL_ONLY.with(|c| c.set(true));
            let seq = tap(&block(d, t), &spec, 0.0, t, depth);
            super::SEQUENTIAL_ONLY.with(|c| c.set(false));
            assert!(
                one.validate().is_ok(),
                "{designation} {depth:?}: one-tool invalid"
            );
            assert!(
                seq.validate().is_ok(),
                "{designation} {depth:?}: sequential invalid"
            );
            let differ = disagreements(&one, &seq, d, 0.0, t);
            assert_eq!(
                differ, 0,
                "{designation} {depth:?}: {differ} samples differ"
            );
        }
    }
}

/// A second hole of the same shape reuses the first's tool: one tool
/// build, and both holes match the sequential cuts.
#[test]
fn holes_of_one_shape_share_a_tool() {
    let spec = ThreadSpec::parse("M4").expect("spec");
    let d = spec.major_diameter();
    let t = 2.6 * d;
    let x2 = 4.0 * d;
    let mut one = tap(&block(d, t), &spec, 0.0, t, HoleDepth::Blind(1.7 * d));
    let before = super::BUILDS.with(std::cell::Cell::get);
    one = tap(&one, &spec, x2, t, HoleDepth::Blind(1.7 * d));
    let after = super::BUILDS.with(std::cell::Cell::get);
    assert_eq!(after, before, "the second hole built a new tool");
    super::SEQUENTIAL_ONLY.with(|c| c.set(true));
    let mut seq = tap(&block(d, t), &spec, 0.0, t, HoleDepth::Blind(1.7 * d));
    seq = tap(&seq, &spec, x2, t, HoleDepth::Blind(1.7 * d));
    super::SEQUENTIAL_ONLY.with(|c| c.set(false));
    for x in [0.0, x2] {
        assert_eq!(disagreements(&one, &seq, d, x, t), 0, "hole at x {x}");
    }
}

/// Every coarse metric size at every common fit clearance either serves
/// a valid one-tool tapped hole or declines with a structured
/// ProfileOverlap — never a boolean decline (PK-899: an M4 at 0.15 left
/// a 0.0017 mm crest land and every cut path declined on it).
#[test]
fn every_size_and_fit_serves_or_declines_with_the_profile_reason() {
    for designation in ["M3", "M4", "M5", "M6", "M8"] {
        let spec = ThreadSpec::parse(designation).expect("spec");
        let d = spec.major_diameter();
        let t = 2.6 * d;
        for clearance in [0.0, 0.05, 0.1, 0.15, 0.2] {
            let _ = super::take_tap_tool_fallbacks();
            let result = tapped_hole(
                &block(d, t),
                &spec,
                Point3::new(0.0, 0.0, t),
                Vec3::new(0.0, 0.0, -1.0),
                HoleDepth::Blind(1.7 * d),
                &TapOptions {
                    clearance,
                    entry_chamfer: true,
                },
                Tolerance::default(),
            );
            match result {
                Ok(body) => {
                    assert!(
                        body.validate().is_ok(),
                        "{designation} {clearance}: invalid"
                    );
                    let fell_back = super::take_tap_tool_fallbacks();
                    assert!(
                        fell_back.is_empty(),
                        "{designation} {clearance}: {fell_back:?}"
                    );
                }
                Err(crate::FastenerError::ProfileOverlap { .. }) => {}
                Err(e) => panic!("{designation} {clearance}: {e}"),
            }
        }
    }
}
