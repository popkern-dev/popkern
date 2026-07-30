//! Strip tessellation measure-fidelity pins (PK-218).
//!
//! A blend strip's chart mixes units (u = radians, v = millimetres).
//! Before PK-218 the boundary-first protection filters ran in raw UV,
//! so on long strips the clearance band exceeded the chart's metric
//! width, evacuated the interior lattice, and the boundary-only
//! Delaunay chorded rail-to-rail across the arc: every mesh vertex
//! lay ON the cylinder while the mesh volume read 0.5% low (OCC read
//! the exact smooth value from the same body). The fix is twofold:
//! metric-scaled protection filters, and the structured tensor mesh
//! for clean four-sided ruled charts whose error is the adaptive
//! edge sampling's own chord tolerance by construction.

use popkern_part::prelude::*;

/// The PK-218 reproducer: 80×50×5 plate, both 80-rims filleted at
/// r = 1.5. Analytic smooth volume 19922.7433; the rail-chord defect
/// read 19823.69 (−99 = two arc-chord lunes). The structured strip
/// mesh must land within the arc-chord budget at chord_height 0.001
/// (measured 0.11 total), and stay watertight.
#[test]
fn plate_rim_strip_mesh_volume_matches_analytic() {
    let plate = Profile::on(Plane::XY)
        .rect(80.0, 50.0)
        .extrude(5.0)
        .unwrap();
    let part = plate
        .fillet(1.5, |e| {
            e.straight()
                .direction(Dir::X)
                .chord_length_between(76.0, 84.0)
                .at_height(5.0)
        })
        .unwrap();
    let brep = part.brep();
    assert!(brep.validate_all().is_empty(), "plate must validate");
    let opens = part.open_edge_count().expect("open edge count");
    assert_eq!(opens, 0, "structured strip mesh must stay watertight");

    let expected = 80.0 * 50.0 * 5.0 - 2.0 * (1.0 - std::f64::consts::PI / 4.0) * 1.5 * 1.5 * 80.0;
    let volume = popkern_measure::body_volume(brep, &popkern_core::tolerance::Tolerance::default())
        .expect("volume");
    assert!(
        (volume - expected).abs() < 0.2,
        "plate rim-strip volume {volume} != analytic {expected} within the arc-chord budget"
    );
}
