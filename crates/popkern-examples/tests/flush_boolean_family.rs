//! PK-184 — the surviving coincident/flush boolean sub-class, pinned
//! at the boolean substrate (no dress-up stage in the way).
//!
//! Two configurations the earlier coincident campaigns never reached:
//!
//! * **Triple-flush partial union** (mounting-bracket substrate): the
//!   wall is flush with the base on THREE planes at once, one operand
//!   face fully contained in the other's coplanar face on two of
//!   them. The chart substrate must resolve boundary-substitute
//!   record directions by original-loop flow (the balance search was
//!   ambiguous), and the weld's collinear valence-2 residue vertices
//!   must dissolve so the physical corner line is ONE edge again.
//!
//! * **Disconnecting flush subtract** (slotted-block substrate): a
//!   through-slot flush on four planes splits the block in two. The
//!   tool's side walls lose their ENTIRE boundary registration to
//!   welded record pairs (the anchorless case — resolved by
//!   direction agreement with the original loop, no anchor needed),
//!   and the assembly must file one shell/solid per component or the
//!   correct population fails Euler and falls to graph-merge.

use popkern_part::prelude::*;

#[test]
fn triple_flush_l_union_is_watertight_with_merged_corners() {
    let base = Profile::on(Plane::XY)
        .rect(60.0, 40.0)
        .extrude(6.0)
        .expect("base");
    let wall = Profile::on(Plane::XY)
        .at(-27.0, 0.0)
        .rect(6.0, 40.0)
        .extrude(35.0)
        .expect("wall");
    let body = base.union(&wall).expect("flush L union");
    assert_eq!(
        body.open_edge_count().expect("opens"),
        0,
        "flush L union must be watertight"
    );
    assert_eq!(
        body.brep().faces.len(),
        8,
        "L-shape: 6 box faces + wall top/x=-24, coplanar skins merged"
    );

    // The weld residue class: the tall corner at (-30, -20) must be
    // ONE 35-length edge, not a 6 + 29 collinear split. A leftover
    // valence-2 vertex there leaks segment semantics into chord
    // selectors (the fillet stage filleted HALF the corner and the
    // continuation sat r(√2−1) off the trimmed support).
    let brep = body.brep();
    let mut corner_spans: Vec<f64> = Vec::new();
    for (_, edge) in brep.edges.iter() {
        let he = &brep.half_edges[edge.half_edge];
        let a = brep.vertices[he.origin].position;
        let b = brep.vertices[brep.half_edges[he.twin].origin].position;
        if (a.x + 30.0).abs() < 1e-9
            && (b.x + 30.0).abs() < 1e-9
            && (a.y + 20.0).abs() < 1e-9
            && (b.y + 20.0).abs() < 1e-9
        {
            corner_spans.push((a.z - b.z).abs());
        }
    }
    assert_eq!(
        corner_spans.len(),
        1,
        "the (-30,-20) corner line must be a single edge, got spans {corner_spans:?}"
    );
    assert!(
        (corner_spans[0] - 35.0).abs() < 1e-9,
        "the merged corner spans the full wall height, got {}",
        corner_spans[0]
    );
}

#[test]
fn flush_l_union_fillet_selects_only_full_corners() {
    // The end-to-end chain the corpus part runs: with the corner
    // lines merged, the chord∈[5,7] selector reaches ONLY the two
    // base-plate corners and the fillet serves watertight.
    let base = Profile::on(Plane::XY)
        .rect(60.0, 40.0)
        .extrude(6.0)
        .expect("base");
    let wall = Profile::on(Plane::XY)
        .at(-27.0, 0.0)
        .rect(6.0, 40.0)
        .extrude(35.0)
        .expect("wall");
    let body = base.union(&wall).expect("union");
    let filleted = body
        .fillet(3.0, |e| {
            e.straight()
                .direction(Dir::Z)
                .chord_length_between(5.0, 7.0)
        })
        .expect("fillet on the two base corners");
    assert_eq!(
        filleted.open_edge_count().expect("opens"),
        0,
        "filleted L must be watertight"
    );
}

/// Moat cut (PK-217 floating-piece leg): subtracting a full-height
/// square tube (walls only) leaves an outer ring AND a floating core
/// whose AABB nests inside the ring's. The shell partition must
/// classify the core by parity against the ring's face subset —
/// even crossings ⇒ not enclosed ⇒ its OWN solid, not a cavity.
#[test]
fn moat_cut_serves_floating_core_as_own_solid() {
    let block = Profile::on(Plane::XY)
        .rect(40.0, 40.0)
        .extrude(10.0)
        .expect("block");
    let tube_outer = Profile::on(Plane::XY)
        .rect(20.0, 20.0)
        .extrude(10.0)
        .expect("tube outer");
    let tube_inner = Profile::on(Plane::XY)
        .rect(12.0, 12.0)
        .extrude(10.0)
        .expect("tube inner");
    let tube = tube_outer.subtract(&tube_inner).expect("tube");
    let result = block.subtract(&tube).expect("moat cut");
    assert_eq!(result.open_edge_count().expect("opens"), 0);
    assert_eq!(result.brep().shells.len(), 2, "ring + core shells");
    assert_eq!(result.brep().solids.len(), 2, "core is a separate solid");
    assert!(
        result.brep().validate().is_ok(),
        "{:?}",
        result.brep().validate()
    );
    let volume = popkern_measure::body_volume(
        result.brep(),
        &popkern_core::tolerance::Tolerance::default(),
    )
    .expect("volume");
    // Ring (40² − 20²)·10 + core 12²·10 = 12000 + 1440.
    assert!(
        (volume - 13440.0).abs() < 1e-6,
        "moat volume {volume} != 13440"
    );
}

/// Cavity union (PK-217 cavity leg): two blocks with facing blind
/// pockets, unioned so the pockets close into a fully-enclosed void.
/// The void's shell reads Inside the outer component's parity ⇒ a
/// CAVITY shell of the SAME solid, and the measured volume excludes
/// the void.
#[test]
fn cavity_union_serves_void_as_inner_shell_of_one_solid() {
    let lower = Profile::on(Plane::XY)
        .rect(20.0, 20.0)
        .extrude(10.0)
        .expect("lower block")
        .subtract(
            &Profile::on(Plane::xy_at(4.0))
                .rect(12.0, 12.0)
                .extrude(6.0)
                .expect("lower pocket"),
        )
        .expect("lower U");
    let upper = Profile::on(Plane::xy_at(6.0))
        .rect(20.0, 20.0)
        .extrude(10.0)
        .expect("upper block")
        .subtract(
            &Profile::on(Plane::xy_at(6.0))
                .rect(12.0, 12.0)
                .extrude(6.0)
                .expect("upper pocket"),
        )
        .expect("upper U");
    let result = lower.union(&upper).expect("cavity union");
    assert_eq!(result.open_edge_count().expect("opens"), 0);
    assert_eq!(result.brep().shells.len(), 2, "outer + cavity shells");
    assert_eq!(result.brep().solids.len(), 1, "one solid with a void");
    assert!(
        result.brep().validate().is_ok(),
        "{:?}",
        result.brep().validate()
    );
    let volume = popkern_measure::body_volume(
        result.brep(),
        &popkern_core::tolerance::Tolerance::default(),
    )
    .expect("volume");
    // Outer 20·20·16 minus the enclosed void 12·12·8.
    assert!(
        (volume - 5248.0).abs() < 1e-6,
        "cavity volume {volume} != 5248"
    );
}

#[test]
fn disconnecting_flush_slot_serves_two_closed_shells() {
    let block = Profile::on(Plane::XY)
        .rect(70.0, 30.0)
        .extrude(12.0)
        .expect("block");
    let slot = Profile::on(Plane::XY)
        .rect(12.0, 30.0)
        .extrude(12.0)
        .expect("slot");
    let result = block.subtract(&slot).expect("flush through-slot subtract");
    assert_eq!(
        result.open_edge_count().expect("opens"),
        0,
        "both pieces must be watertight"
    );
    assert_eq!(
        result.brep().faces.len(),
        12,
        "two disjoint boxes: 6 + 6 faces"
    );
    assert_eq!(
        result.brep().shells.len(),
        2,
        "one shell per disconnected component"
    );
    assert_eq!(
        result.brep().solids.len(),
        2,
        "mutually exterior components are separate solids"
    );
    assert!(
        result.brep().validate().is_ok(),
        "multi-shell result passes validation: {:?}",
        result.brep().validate()
    );
}
