// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use std::f64::consts::PI;
use std::time::Instant;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::HingeError;
use crate::living::{LivingHingeParams, LivingHingeStyle, living_hinge};

/// Fine-tessellation volume against the exact one.
const VOLUME_REL_TOL: f64 = 1.0e-3;

fn build(p: &LivingHingeParams) -> BRep {
    let body = living_hinge(p, Tolerance::default()).unwrap_or_else(|e| panic!("{p:?}: {e}"));
    let mesh = popkern_mesh::tessellate(&body, &popkern_mesh::TessellationOptions::display())
        .expect("mesh");
    assert_eq!(
        popkern_mesh::welded_open_edges(&mesh),
        0,
        "watertight {p:?}"
    );
    body
}

fn volume(body: &BRep) -> f64 {
    popkern_measure::body_volume(body, &Tolerance::default()).expect("volume")
}

/// Profile area × length: panels and the zone at full thickness, minus
/// the gap above the web, plus the material the root fillets keep.
fn straight_volume(p: &LivingHingeParams) -> f64 {
    let r = p.root_radius;
    let area = (2.0 * p.panel_width + p.span) * p.thickness - p.span * (p.thickness - p.web)
        + 2.0 * (r * r - PI * r * r / 4.0);
    area * p.length
}

/// The plate minus every slot (a rectangle with one semicircular end).
fn slotted_volume(p: &LivingHingeParams) -> f64 {
    let w = p.slot_width;
    let slot = w * (p.length - p.bridge - w / 2.0) + PI * w * w / 8.0;
    let plate = (2.0 * p.panel_width + p.span) * p.length;
    (plate - slot * p.slot_centres().len() as f64) * p.thickness
}

fn assert_volume(body: &BRep, expected: f64, what: &str) -> f64 {
    let got = volume(body);
    let rel = (got - expected).abs() / expected;
    assert!(
        rel < VOLUME_REL_TOL,
        "{what}: {got} vs exact {expected} ({rel:.2e})"
    );
    rel
}

#[test]
fn straight_living_hinges_sweep_with_exact_volumes() {
    let mut worst: f64 = 0.0;
    for (span, web_frac, r_frac) in [
        (10.0, 0.2, 0.5),
        (10.0, 0.1, 0.9),
        (6.0, 0.3, 0.0),
        (20.0, 0.15, 0.4),
        (4.0, 0.25, 0.7),
    ] {
        let mut p = LivingHingeParams::straight(span, 40.0);
        p.web = p.thickness * web_frac;
        p.root_radius = ((p.thickness - p.web) * r_frac).min(span / 2.0 * 0.9);
        let body = build(&p);
        worst = worst.max(assert_volume(&body, straight_volume(&p), "straight"));
    }
    eprintln!("straight sweep: worst volume deviation {worst:.2e}");
}

#[test]
fn slotted_living_hinges_sweep_with_exact_volumes() {
    let mut worst: f64 = 0.0;
    for (span, pitch, width, bridge) in [
        (10.0, 1.25, 0.4, 4.0),
        (12.0, 1.0, 0.5, 2.0),
        (8.0, 2.0, 1.0, 6.0),
        (20.0, 0.8, 0.3, 3.0),
    ] {
        let p = LivingHingeParams {
            slot_pitch: pitch,
            slot_width: width,
            bridge,
            ..LivingHingeParams::slotted(span, 40.0)
        };
        let body = build(&p);
        worst = worst.max(assert_volume(&body, slotted_volume(&p), "slotted"));
    }
    eprintln!("slotted sweep: worst volume deviation {worst:.2e}");
}

#[test]
fn two_hundred_slots_build_in_one_extrusion() {
    let params = |span: f64| LivingHingeParams {
        slot_pitch: 0.5,
        slot_width: 0.2,
        bridge: 2.0,
        ..LivingHingeParams::slotted(span, 40.0)
    };
    let (p50, p200) = (params(25.0), params(100.0));
    assert_eq!(p50.slot_centres().len(), 50);
    assert_eq!(p200.slot_centres().len(), 200);
    // Recorded budget (design §8): one extrusion whose cost grows
    // linearly in the slot count (0.29 s for 200 slots measured
    // 2026-10-07; 58.7 s before PK-709 made the Euler checks linear).
    // Gated as an interleaved ratio, min of three runs, so suite
    // contention scales both sides (PK-810: an absolute 1.0 s budget
    // tripped at 1.11 s under full-suite load). 50 -> 200 slots reads
    // 10.8x today (growth ~n^1.7); a quadratic reads 16x and the
    // PK-709 regression 200x, which the 5 s backstop also catches.
    // Load can still land on one side of a single interleaved
    // measurement (nextest's parallel crate run read 14.5x once), so
    // up to three measurements run and the best ratio counts: a real
    // growth regression reads high on every one.
    let mut body = None;
    let (mut ratio, mut small, mut large) = (f64::MAX, f64::MAX, f64::MAX);
    for _ in 0..3 {
        let (mut s, mut l) = (f64::MAX, f64::MAX);
        for _ in 0..3 {
            let start = Instant::now();
            living_hinge(&p50, Tolerance::default()).expect("50 slots");
            s = s.min(start.elapsed().as_secs_f64());
            let start = Instant::now();
            body = Some(living_hinge(&p200, Tolerance::default()).expect("200 slots"));
            l = l.min(start.elapsed().as_secs_f64());
        }
        if l / s.max(1e-9) < ratio {
            (ratio, small, large) = (l / s.max(1e-9), s, l);
        }
        if ratio < 14.0 {
            break;
        }
    }
    let body = body.expect("200 slots built");
    eprintln!(
        "200-slot living hinge: {large:.3} s, 50 slots {small:.3} s, ratio {ratio:.1}, {} faces",
        body.faces.len()
    );
    assert!(
        ratio < 14.0,
        "200 slots cost {ratio:.1}x 50 slots ({large:.3} s vs {small:.3} s): growth is quadratic"
    );
    assert!(large < 5.0, "200 slots took {large:.3} s");
    assert_volume(&body, slotted_volume(&p200), "200 slots");
}

#[test]
fn builds_are_deterministic() {
    for p in [
        LivingHingeParams::straight(10.0, 40.0),
        LivingHingeParams::slotted(10.0, 40.0),
    ] {
        let a = living_hinge(&p, Tolerance::default()).expect("a");
        let b = living_hinge(&p, Tolerance::default()).expect("b");
        let pa: Vec<_> = a.vertices.values().map(|v| v.position).collect();
        let pb: Vec<_> = b.vertices.values().map(|v| v.position).collect();
        assert_eq!(pa, pb);
    }
}

#[test]
fn declines_name_the_field() {
    let decline = |p: LivingHingeParams| living_hinge(&p, Tolerance::default()).expect_err("no");
    let s = LivingHingeParams::straight(10.0, 40.0);
    assert!(matches!(
        decline(LivingHingeParams {
            root_radius: s.thickness,
            ..s.clone()
        }),
        HingeError::InvalidSpec {
            field: "root_radius",
            ..
        }
    ));
    assert!(matches!(
        decline(LivingHingeParams {
            web: 0.0,
            ..s.clone()
        }),
        HingeError::InvalidSpec { field: "web", .. }
    ));
    let t = LivingHingeParams::slotted(10.0, 40.0);
    assert!(matches!(
        decline(LivingHingeParams {
            slot_width: t.slot_pitch,
            ..t.clone()
        }),
        HingeError::InvalidSpec {
            field: "slot_width",
            ..
        }
    ));
    assert!(matches!(
        decline(LivingHingeParams {
            bridge: 40.0,
            ..t.clone()
        }),
        HingeError::InvalidSpec {
            field: "bridge",
            ..
        }
    ));
    for style in [LivingHingeStyle::Straight, LivingHingeStyle::Slotted] {
        assert_eq!(LivingHingeStyle::parse(style.spelling()), Some(style));
    }
}
