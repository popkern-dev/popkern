//! PK-849: a bearing-head screw seated flush on its own tapped hole
//! clears the tap. The screw stays at full major diameter for
//! `head_side_full_major` below its head, and the tap's entry cone
//! must clear that stretch: ending at P/8 + clearance, the first tap
//! ridge met a flush M6 socket screw's top turn over a quarter turn.
//! (M3 and M4 taps into this block decline in merge, and M12 at
//! clearance 0.1 fails precision validation, with either cone: PK-895.)

use super::*;
use crate::builders::{SocketCapParams, socket_cap_screw};
use popkern_bool::{PointClassification, face_bounding_boxes, point_in_solid_final};
use popkern_core::geometry::{Point3, Vec3};

/// `point`'s ray-vote classification against `body`.
fn classify(point: &Point3, body: &BRep) -> PointClassification {
    point_in_solid_final(
        point,
        body,
        &face_bounding_boxes(body),
        &Tolerance::default(),
    )
    .expect("classification")
}

/// A block whose entry face is z = 0, material in z > 0, tapped along
/// +z (the screws' insertion frame: bearing plane z = 0, shank +z).
fn tapped_block(spec: &ThreadSpec, clearance: f64) -> BRep {
    let d = spec.major_diameter();
    let half = 2.0 * d;
    let pts = [
        Point3::new(-half, -half, 0.0),
        Point3::new(half, -half, 0.0),
        Point3::new(half, half, 0.0),
        Point3::new(-half, half, 0.0),
    ];
    let block = BRep::prism(
        &pts,
        Vec3::new(0.0, 0.0, 1.0),
        3.0 * d,
        Tolerance::default(),
    )
    .expect("block");
    let options = TapOptions {
        clearance,
        ..TapOptions::default()
    };
    tapped_hole(
        &block,
        spec,
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        HoleDepth::Blind(2.0 * d),
        &options,
        Tolerance::default(),
    )
    .unwrap_or_else(|e| panic!("{spec:?} clearance {clearance}: {e:?}"))
}

#[test]
fn a_flush_socket_screw_clears_its_tapped_hole() {
    for (designation, clearances) in [("M6", &[0.1, 0.2][..]), ("M8", &[0.2]), ("M12", &[0.2])] {
        let spec = ThreadSpec::parse(designation).expect("spec");
        let full_major = thread_form::head_side_full_major(&spec).expect("runout");
        let screw = socket_cap_screw(
            &SocketCapParams {
                spec: spec.clone(),
                length: 1.5 * spec.major_diameter(),
                thread_length: None,
                tip_chamfer: true,
                cosmetic_thread: false,
            },
            Tolerance::default(),
        )
        .expect("screw");
        for &clearance in clearances {
            let block = tapped_block(&spec, clearance);
            // At the screw's crest radius through the entry region
            // (three leads), no point is inside both screw and tap. A
            // positive clearance: at zero the mating flanks coincide and
            // read inside both wherever the threads engage.
            let r = spec.major_diameter() / 2.0 - 0.02 * spec.pitch();
            let reach = 3.0 * spec.lead();
            let (mut screw_full, mut ridge) = (0, 0);
            for k in 1..=36 {
                let z = reach * k as f64 / 36.0;
                for j in 0..24 {
                    let a = std::f64::consts::TAU * j as f64 / 24.0;
                    let p = Point3::new(r * a.cos(), r * a.sin(), z);
                    let in_screw = classify(&p, &screw) == PointClassification::Inside;
                    let in_block = classify(&p, &block) == PointClassification::Inside;
                    assert!(
                        !(in_screw && in_block),
                        "{designation} clearance {clearance}: screw and tap overlap at \
                         z {z:.3} (full-major runout {full_major:.3}), angle {a:.2}"
                    );
                    screw_full += usize::from(in_screw && z < full_major);
                    ridge += usize::from(in_block);
                }
            }
            // Controls: the probe sees the screw's full-major runout and
            // the tap's ridges it must clear.
            assert!(
                screw_full > 0,
                "{designation}: no full-major runout sampled"
            );
            assert!(
                ridge > 0,
                "{designation}: no tap ridge reaches the crest radius"
            );
        }
    }
}
