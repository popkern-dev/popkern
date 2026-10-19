//! PK-849: a bearing-head screw seated flush on its own tapped hole
//! clears the tap. The screw stays at full major diameter for
//! `head_side_full_major` below its head, and the tap's entry cone
//! must clear that stretch: ending at P/8 + clearance, the first tap
//! ridge met a flush M6 socket screw's top turn over a quarter turn.

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
    tapped_block_at(spec, clearance, 2.0 * spec.major_diameter())
}

/// [`tapped_block`] with a blind `depth`.
fn tapped_block_at(spec: &ThreadSpec, clearance: f64, depth: f64) -> BRep {
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
        HoleDepth::Blind(depth),
        &options,
        Tolerance::default(),
    )
    .unwrap_or_else(|e| panic!("{spec:?} clearance {clearance}: {e:?}"))
}

#[test]
fn a_flush_socket_screw_clears_its_tapped_hole() {
    for (designation, clearances) in [
        ("M4", &[0.1][..]),
        ("M6", &[0.1, 0.2]),
        ("M8", &[0.2]),
        ("M12", &[0.1, 0.2]),
    ] {
        let spec = ThreadSpec::parse(designation).expect("spec");
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
            assert_seated_clear(&spec, clearance, &screw, designation);
        }
    }
}

/// A hex bolt bears on its flange face like a socket head (head clear
/// z = 0, the same runout), so it seats flush the same way.
#[test]
fn a_flush_hex_bolt_clears_its_tapped_hole() {
    use crate::builders::{HexBoltParams, hex_bolt};
    for designation in ["M6", "M8"] {
        let spec = ThreadSpec::parse(designation).expect("spec");
        let bolt = hex_bolt(
            &HexBoltParams {
                spec: spec.clone(),
                length: 1.5 * spec.major_diameter(),
                thread_length: None,
                across_flats: None,
                head_height: None,
                head_chamfer: true,
                tip_chamfer: true,
                cosmetic_thread: false,
            },
            Tolerance::default(),
        )
        .expect("bolt");
        assert_seated_clear(&spec, 0.2, &bolt, designation);
    }
}

/// At the screw's crest radius through the entry region (three leads),
/// no point is inside both `screw` (seated flush, insertion frame) and
/// its tapped block. A positive `clearance`: at zero the mating flanks
/// coincide and read inside both wherever the threads engage.
fn assert_seated_clear(spec: &ThreadSpec, clearance: f64, screw: &BRep, label: &str) {
    let block = tapped_block(spec, clearance);
    let full_major = thread_form::head_side_full_major(spec).expect("runout");
    let r = spec.major_diameter() / 2.0 - 0.02 * spec.pitch();
    let reach = 3.0 * spec.lead();
    let (mut screw_full, mut ridge) = (0, 0);
    for k in 1..=36 {
        let z = reach * k as f64 / 36.0;
        for j in 0..24 {
            let a = std::f64::consts::TAU * j as f64 / 24.0;
            let p = Point3::new(r * a.cos(), r * a.sin(), z);
            let in_screw = classify(&p, screw) == PointClassification::Inside;
            let in_block = classify(&p, &block) == PointClassification::Inside;
            assert!(
                !(in_screw && in_block),
                "{label} clearance {clearance}: screw and tap overlap at z {z:.3} \
                 (full-major runout {full_major:.3}), angle {a:.2}"
            );
            screw_full += usize::from(in_screw && z < full_major);
            ridge += usize::from(in_block);
        }
    }
    // Controls: the probe sees the screw's full-major runout and the
    // tap's ridges it must clear.
    assert!(screw_full > 0, "{label}: no full-major runout sampled");
    assert!(ridge > 0, "{label}: no tap ridge reaches the crest radius");
}

/// PK-895: a blind tap's thread tool stays clear of the hole's floor at
/// every depth phase. Its helix runs half a lead past the nominal end,
/// rounds up to a quarter turn, and the ridge reaches its half width
/// beyond: ending ¾ lead above the floor, the tip met the floor at every
/// whole-pitch M3 depth, on the floor rim's seam vertex, and the cut
/// declined. Twelve depth phases per size and clearance serve.
#[test]
fn blind_taps_serve_at_every_depth_phase() {
    for designation in ["M3", "M4", "M5"] {
        let spec = ThreadSpec::parse(designation).expect("spec");
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
        for clearance in [0.0, 0.1] {
            for k in 0..12 {
                let depth = 1.7 * d + spec.pitch() * k as f64 / 6.0;
                let options = TapOptions {
                    clearance,
                    entry_chamfer: true,
                };
                tapped_hole(
                    &block,
                    &spec,
                    Point3::new(0.0, 0.0, 0.0),
                    Vec3::new(0.0, 0.0, 1.0),
                    HoleDepth::Blind(depth),
                    &options,
                    Tolerance::default(),
                )
                .unwrap_or_else(|e| {
                    panic!("{designation} clearance {clearance} depth {depth}: {e:?}")
                });
            }
        }
    }
}

/// PK-898: an M12 blind tap serves at clearance 0 and 0.1. The entry
/// chamfer's cut trims each thread flank x bore record at the cone, and
/// the assembly heals the trimmed piece onto its nodes by interpolating
/// samples of the record; 17 fixed samples left a 0.7-turn piece, about
/// 22 long, 2e-4 off both surfaces, and the hole failed the 1e-4 edge
/// precision bound at every depth phase.
#[test]
fn an_m12_blind_tap_serves_at_tight_clearances() {
    let spec = ThreadSpec::parse("M12").expect("spec");
    for clearance in [0.0, 0.1] {
        for k in 0..4 {
            let depth = 1.7 * spec.major_diameter() + spec.pitch() * k as f64 / 4.0;
            tapped_block_at(&spec, clearance, depth);
        }
    }
}
