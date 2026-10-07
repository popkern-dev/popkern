// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use std::f64::consts::PI;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fastener::ThreadSpec;

use crate::leaf::hinge_component;
use crate::params::{
    HingeComponent, HingeHoles, HingeKind, HingeLayout, HingeParams, hinge_layout,
};

/// Tessellated volume against the exact one: the measure crate's fine
/// tessellation (chord height 1e-3) on cylinders of radius ≥ 1;
/// the 54-hinge sweep measured 4.0e-4 at worst (2026-10-07).
const VOLUME_REL_TOL: f64 = 1.0e-3;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn volume(body: &BRep) -> f64 {
    popkern_measure::body_volume(body, &tol()).expect("volume")
}

fn open_edges(body: &BRep) -> usize {
    let mesh = popkern_mesh::tessellate(body, &popkern_mesh::TessellationOptions::display())
        .expect("mesh");
    popkern_mesh::welded_open_edges(&mesh)
}

/// Area of the disc of radius `r` below the line `y = c`.
fn disc_below(r: f64, c: f64) -> f64 {
    let c = c.clamp(-r, r);
    c * (r * r - c * c).sqrt() + r * r * (c / r).asin() + r * r * PI / 2.0
}

/// The exact volume of a flat, hole-free leaf: the profile area along
/// the length, the other leaf's zones notched back to the notch line,
/// the bore (or, print-in-place leaf A, the integral pin added through
/// the zones).
fn leaf_volume(
    kind: HingeKind,
    component: HingeComponent,
    p: &HingeParams,
    l: &HingeLayout,
) -> f64 {
    let (r, t, w) = (l.knuckle_radius, p.thickness, p.leaf_width);
    let lower_left_segment = disc_below(r, -r + t) / 2.0;
    let area = w * t + PI * r * r - lower_left_segment;
    let near = l.notch * t + PI * r * r - lower_left_segment;
    let other = if component == HingeComponent::LeafA {
        HingeComponent::LeafB
    } else {
        HingeComponent::LeafA
    };
    let zones: f64 = l
        .owned(other)
        .iter()
        .map(|&(z0, z1)| (z1 + p.gap).min(p.length) - (z0 - p.gap).max(0.0))
        .sum();
    let kept = p.length - zones;
    if kind == HingeKind::PrintInPlace && component == HingeComponent::LeafA {
        let pin = PI * (p.pin_diameter / 2.0).powi(2);
        kept * area + zones * (area - near) + zones * pin
    } else {
        let bore = PI * l.bore_radius.powi(2);
        kept * (area - bore) + zones * (area - near)
    }
}

fn check(kind: HingeKind, p: &HingeParams) -> f64 {
    let layout = hinge_layout(kind, p).expect("layout");
    let mut worst: f64 = 0.0;
    let mut components = vec![HingeComponent::LeafA, HingeComponent::LeafB];
    if kind == HingeKind::Butt {
        components.push(HingeComponent::Pin);
    }
    for c in components {
        let body = hinge_component(kind, c, p, tol())
            .unwrap_or_else(|e| panic!("{kind:?} {c:?} {p:?}: {e}"));
        assert_eq!(open_edges(&body), 0, "{kind:?} {c:?} watertight {p:?}");
        let expected = match c {
            HingeComponent::Pin => PI * (p.pin_diameter / 2.0).powi(2) * p.length,
            leaf => leaf_volume(kind, leaf, p, &layout),
        };
        let got = volume(&body);
        let rel = (got - expected).abs() / expected;
        assert!(
            rel < VOLUME_REL_TOL,
            "{kind:?} {c:?}: volume {got} vs exact {expected} (rel {rel:.2e}) {p:?}"
        );
        worst = worst.max(rel);
    }
    worst
}

/// The design review's consumer example, as the first test.
#[test]
fn the_default_butt_hinge_builds_three_exact_bodies() {
    let worst = check(HingeKind::Butt, &HingeParams::new(60.0));
    assert!(worst < VOLUME_REL_TOL);
}

#[test]
fn a_parameter_sweep_builds_every_combination() {
    // 54 deterministic combinations over knuckles 2..9, pin 2..8,
    // clearance 0..0.5, thickness 1..4 on a 60 hinge.
    let knuckles = [2, 3, 4, 5, 7, 9];
    let pins = [2.0, 4.0, 6.0];
    let clearances = [0.0, 0.25, 0.5];
    let thicknesses = [1.0, 2.5, 4.0];
    let mut worst: f64 = 0.0;
    let mut built = 0;
    for (i, &n) in knuckles.iter().enumerate() {
        for (j, &pin) in pins.iter().enumerate() {
            for (k, &c) in clearances.iter().enumerate() {
                // One thickness per (n, pin, c), rotating, keeps the
                // sweep at 54 combinations while covering all three.
                let t = thicknesses[(i + j + k) % 3];
                let p = HingeParams {
                    knuckles: n,
                    pin_diameter: pin,
                    clearance: c,
                    thickness: t,
                    ..HingeParams::new(60.0)
                };
                worst = worst.max(check(HingeKind::Butt, &p));
                built += 1;
            }
        }
    }
    assert_eq!(built, 54);
    eprintln!("butt sweep: {built} hinges, worst volume deviation {worst:.2e}");
}

#[test]
fn builds_are_deterministic() {
    let p = HingeParams::new(60.0);
    for c in [HingeComponent::LeafA, HingeComponent::LeafB] {
        let a = hinge_component(HingeKind::Butt, c, &p, tol()).expect("a");
        let b = hinge_component(HingeKind::Butt, c, &p, tol()).expect("b");
        assert_eq!(a.faces.len(), b.faces.len());
        let pa: Vec<_> = a.vertices.values().map(|v| v.position).collect();
        let pb: Vec<_> = b.vertices.values().map(|v| v.position).collect();
        assert_eq!(pa, pb, "{c:?} vertex positions bit-identical");
    }
}

#[test]
fn countersunk_screw_holes_cut_both_leaves() {
    let p = HingeParams {
        holes: Some(HingeHoles {
            per_leaf: 2,
            thread: ThreadSpec::parse("M3").expect("M3"),
            countersink: true,
        }),
        ..HingeParams::new(60.0)
    };
    let plain = HingeParams::new(60.0);
    for c in [HingeComponent::LeafA, HingeComponent::LeafB] {
        let holed = hinge_component(HingeKind::Butt, c, &p, tol()).expect("holed");
        let bare = hinge_component(HingeKind::Butt, c, &plain, tol()).expect("bare");
        assert_eq!(open_edges(&holed), 0);
        assert!(
            holed.faces.len() > bare.faces.len(),
            "{c:?} gained hole faces"
        );
        let removed = volume(&bare) - volume(&holed);
        // Each hole removes at least its 3.4 bore through 2.0 of plate.
        let bore = PI * 1.7 * 1.7 * 2.0;
        assert!(removed > 2.0 * bore, "{c:?} removed {removed}");
    }
}

#[test]
fn leaf_b_turns_to_the_open_angle() {
    let p = HingeParams {
        angle_deg: 90.0,
        ..HingeParams::new(60.0)
    };
    let b = hinge_component(HingeKind::Butt, HingeComponent::LeafB, &p, tol()).expect("b");
    // Flat leaf B reaches x = +15; turned 90° about +Z it reaches
    // y = +15 instead, and stays within the knuckle radius in +x.
    let max_x = b
        .vertices
        .values()
        .map(|v| v.position.x)
        .fold(f64::MIN, f64::max);
    let max_y = b
        .vertices
        .values()
        .map(|v| v.position.y)
        .fold(f64::MIN, f64::max);
    assert!((max_y - 15.0).abs() < 1e-9, "max y {max_y}");
    assert!(max_x <= 4.0 + 1e-9, "max x {max_x}");
}
