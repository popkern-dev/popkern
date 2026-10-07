//! PK-701: a hinge as an assembly — leaves and pin as parts mated about
//! the pin axis, the open angle an assembly parameter that re-places
//! leaf B without rebuilding any part, limits, the swing sweep against
//! an analytic first contact, and `.pkn` round trips of mate drives.

use popkern_assembly::{Assembly, Document};
use popkern_core::geometry::Point3;
use popkern_history::feature_hinge::HingeKindDesc;
use popkern_script::body_hinge::HingeOptions;
use popkern_script::scalar::Scalar;

fn hinge(kind: HingeKindDesc) -> Assembly {
    let mut asm = Assembly::new("hinged");
    asm.add_hinge(
        "door",
        kind,
        Scalar::Literal(60.0),
        &HingeOptions::default(),
    )
    .expect("hinge");
    asm
}

/// Where leaf B's part-local point lands in the world.
fn leaf_b_world(asm: &Assembly, p: Point3) -> Point3 {
    let instance = asm
        .instances(popkern_assembly::Parent::Root)
        .expect("root")
        .iter()
        .find(|i| i.name == "door_leaf_b")
        .expect("leaf b");
    instance.placement.apply_point(&p)
}

#[test]
fn a_butt_hinge_assembly_has_three_mated_parts_and_an_angle_slider() {
    let mut asm = hinge(HingeKindDesc::Butt);
    assert_eq!(asm.definitions().len(), 3);
    assert!(asm.params().iter().any(|p| p.name == "door_angle"));
    // Flat: leaf B's local +X tip stays on +X.
    let tip = Point3::new(15.0, 0.0, 30.0);
    let flat = leaf_b_world(&asm, tip);
    assert!(
        (flat.x - 15.0).abs() < 1e-6 && flat.y.abs() < 1e-6,
        "{flat:?}"
    );
    // Open to 90°: leaf B turned 90° about +Z — its tip on +Y.
    let rebuilt = asm.set_assembly_param("door_angle", 90.0).expect("90");
    assert!(
        rebuilt.is_empty(),
        "an angle edit replays no part: {rebuilt:?}"
    );
    let open = leaf_b_world(&asm, tip);
    assert!(
        open.x.abs() < 1e-6 && (open.y - 15.0).abs() < 1e-6,
        "{open:?}"
    );
    assert!((open.z - 30.0).abs() < 1e-9, "the pin axis holds z");
    // Folded shut: tip on −X.
    asm.set_assembly_param("door_angle", 0.0).expect("0");
    let shut = leaf_b_world(&asm, tip);
    assert!(
        (shut.x + 15.0).abs() < 1e-6 && shut.y.abs() < 1e-6,
        "{shut:?}"
    );
}

#[test]
fn the_angle_is_limited() {
    let mut asm = hinge(HingeKindDesc::PrintInPlace);
    assert_eq!(asm.definitions().len(), 2, "print-in-place has no pin part");
    for bad in [-5.0, 400.0] {
        let err = asm
            .set_assembly_param("door_angle", bad)
            .expect_err("limited");
        assert!(err.to_string().contains("outside its limits"), "{err}");
    }
    // A refused value leaves the angle where it was.
    let p = asm
        .params()
        .iter()
        .find(|p| p.name == "door_angle")
        .expect("p");
    assert!((p.display_value() - 180.0).abs() < 1e-9);
}

#[test]
fn mate_drives_and_limits_survive_a_pkn_round_trip() {
    let mut asm = hinge(HingeKindDesc::Butt);
    asm.set_assembly_param("door_angle", 120.0).expect("120");
    let bytes = asm.to_pkn_bytes().expect("save");
    let Document::Assembly(mut loaded) = Document::from_pkn_bytes(&bytes).expect("load") else {
        panic!("an assembly document");
    };
    assert!(loaded.drives_mates("door_angle"));
    let tip = Point3::new(15.0, 0.0, 30.0);
    let before = leaf_b_world(&asm, tip);
    assert!((leaf_b_world(&loaded, tip) - before).length() < 1e-9);
    loaded
        .set_assembly_param("door_angle", 90.0)
        .expect("slider after load");
    assert!(loaded.set_assembly_param("door_angle", 361.0).is_err());
}

/// Whether two convex polygons overlap (separating-axis theorem).
fn overlap(a: &[[f64; 2]], b: &[[f64; 2]]) -> bool {
    for poly in [a, b] {
        for i in 0..poly.len() {
            let (p, q) = (poly[i], poly[(i + 1) % poly.len()]);
            let n = [q[1] - p[1], p[0] - q[0]];
            let proj = |s: &[[f64; 2]]| {
                s.iter()
                    .map(|v| v[0] * n[0] + v[1] * n[1])
                    .fold((f64::MAX, f64::MIN), |(lo, hi), x| (lo.min(x), hi.max(x)))
            };
            let ((alo, ahi), (blo, bhi)) = (proj(a), proj(b));
            if ahi <= blo || bhi <= alo {
                return false;
            }
        }
    }
    true
}

fn rect(x0: f64, x1: f64, y0: f64, y1: f64) -> Vec<[f64; 2]> {
    vec![[x0, y0], [x1, y0], [x1, y1], [x0, y1]]
}

/// The first opening angle past flat where leaf B's plate meets leaf
/// A's, from the cross-sections alone: knuckles are rotation-invariant
/// and the notches clear them, so contact is plate against plate —
/// in A's knuckle zones (A's plate reaches the axis, B's stops at the
/// notch) and in B's (the reverse).
fn analytic_first_contact(r: f64, t: f64, w: f64, notch: f64) -> f64 {
    let zones = [
        (rect(-w, 0.0, -r, -r + t), rect(notch, w, -r, -r + t)),
        (rect(-w, -notch, -r, -r + t), rect(0.0, w, -r, -r + t)),
    ];
    let hits = |alpha: f64| {
        let (s, c) = alpha.to_radians().sin_cos();
        zones.iter().any(|(a, b)| {
            // Opening past flat turns leaf B clockwise by alpha.
            let turned: Vec<[f64; 2]> = b
                .iter()
                .map(|p| [c * p[0] + s * p[1], -s * p[0] + c * p[1]])
                .collect();
            overlap(a, &turned)
        })
    };
    let mut alpha = 0.0;
    while !hits(alpha) {
        alpha += 0.01;
        assert!(alpha < 180.0, "no contact");
    }
    let (mut lo, mut hi) = (alpha - 0.01, alpha);
    while hi - lo > 1e-6 {
        let mid = 0.5 * (lo + hi);
        if hits(mid) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    180.0 + hi
}

#[test]
fn the_swing_finds_the_analytic_first_contact() {
    let mut asm = hinge(HingeKindDesc::PrintInPlace);
    // Default 60 print-in-place hinge: t = 2, clearance 0.3,
    // R = pin/2 + c + t = 2 + 0.3 + 2, leaf width 15.
    let (t, c) = (2.0, 0.3);
    let r = 2.0 + c + t;
    let expected = analytic_first_contact(r, t, 15.0, r + c);
    let report = asm
        .swing("door_angle", 180.0, 359.0, 5.0, 0.1)
        .expect("swing");
    let found = report
        .first_contact
        .expect("the leaves collide opening past flat");
    eprintln!(
        "swing: first contact {found:.2}° (analytic {expected:.2}°), {} samples",
        report.samples.len()
    );
    assert!(
        (found - expected).abs() < 0.5,
        "found {found}, analytic {expected}"
    );
    // Closing from flat to shut never collides (the leaves fold onto the
    // pin's diameter apart).
    let closing = asm
        .swing("door_angle", 180.0, 0.0, 10.0, 0.1)
        .expect("close");
    assert_eq!(closing.first_contact, None, "{:?}", closing.samples);
    // The sweep restores the angle.
    let p = asm
        .params()
        .iter()
        .find(|p| p.name == "door_angle")
        .expect("p");
    assert!((p.display_value() - 180.0).abs() < 1e-9);
}

/// PK-696: a print-in-place hinge's leaves are apart by exactly the
/// printed clearance — the pin in its bore, the knuckles across their
/// gaps, each leaf's notch round the other's knuckles — measured by the
/// interference check's certified bounds, at several clearances and
/// knuckle counts, and at an open angle.
#[test]
fn print_in_place_leaves_are_apart_by_the_printed_clearance() {
    use popkern_assembly::{InterferenceOptions, PairStatus};
    for (knuckles, c) in [(3.0, 0.15), (5.0, 0.3), (9.0, 0.6)] {
        let mut asm = Assembly::new("pip");
        asm.add_hinge(
            "h",
            HingeKindDesc::PrintInPlace,
            Scalar::Literal(60.0),
            &HingeOptions {
                knuckles: Some(Scalar::Literal(knuckles)),
                clearance: Some(Scalar::Literal(c)),
                gap: Some(Scalar::Literal(c)),
                ..HingeOptions::default()
            },
        )
        .expect("hinge");
        for angle in [180.0, 95.0] {
            asm.set_assembly_param("h_angle", angle).expect("angle");
            let report = asm
                .interference(&InterferenceOptions {
                    clearance: 2.0 * c,
                    ..InterferenceOptions::default()
                })
                .expect("interference");
            assert_eq!(
                report.pairs.len(),
                1,
                "{knuckles} knuckles, c {c}: {:?}",
                report.pairs
            );
            match &report.pairs[0].status {
                PairStatus::Clear {
                    distance,
                    lower_bound,
                    refined,
                } => {
                    assert!(*refined, "the distance is refined onto the true faces");
                    assert!(
                        (distance - c).abs() < 1e-6,
                        "{knuckles} knuckles at {angle}°: distance {distance} vs clearance {c}"
                    );
                    assert!(
                        *lower_bound > 0.0,
                        "certified apart: lower bound {lower_bound}"
                    );
                }
                other => panic!("{knuckles} knuckles, c {c}, {angle}°: {other:?}"),
            }
        }
    }
}
