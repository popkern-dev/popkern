//! Machinery checks for the offset classification oracle: it must
//! bless a correct offset and flag a wrong-shape result. The full
//! convex-body Steiner calibration battery lives with the fixture
//! corpus (popkern-examples).

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::{
    OffsetOracleOptions, OffsetSemantics, classify_offset_agreement, grid_offset_volume,
    steiner_offset_volume,
};
use crate::TessellationOptions;
use crate::point_query::PointMeshIndex;

fn sphere_index(radius: f64) -> PointMeshIndex {
    let brep = BRep::sphere(radius, Tolerance::default()).expect("sphere");
    let mesh = crate::tessellate(&brep, &TessellationOptions::default()).expect("tessellate");
    PointMeshIndex::build(&mesh)
}

#[test]
fn oracle_blesses_a_correct_sphere_dilation() {
    let source = sphere_index(1.0);
    let grown = sphere_index(1.25);
    let report = classify_offset_agreement(
        &source,
        &grown,
        OffsetSemantics::Grow { distance: 0.25 },
        &OffsetOracleOptions::default(),
        |_| false,
    );
    assert!(
        report.is_clean(),
        "disagreements: {:?}",
        report.disagreements
    );
    assert!(
        report.classified_fraction() > 0.5,
        "sweep barely classified anything: {report:?}"
    );
}

#[test]
fn oracle_flags_a_wrong_shape_served_as_the_offset() {
    // A sphere grown by 0.25 served as if it were the 0.10 offset:
    // every point in the spherical shell between r=1.10 and r=1.25
    // is predicted OUTSIDE but actually inside.
    let source = sphere_index(1.0);
    let wrong = sphere_index(1.25);
    let report = classify_offset_agreement(
        &source,
        &wrong,
        OffsetSemantics::Grow { distance: 0.10 },
        &OffsetOracleOptions::default(),
        |_| false,
    );
    assert!(
        report.disagreement_count > 100,
        "wrong shape must produce mass disagreement: {report:?}"
    );
    assert!(!report.disagreements.is_empty());
}

#[test]
fn oracle_flags_a_wrong_erosion() {
    let source = sphere_index(1.0);
    let wrong_inner = sphere_index(0.95);
    let report = classify_offset_agreement(
        &source,
        &wrong_inner,
        OffsetSemantics::Shrink { distance: 0.25 },
        &OffsetOracleOptions::default(),
        |_| false,
    );
    assert!(
        report.disagreement_count > 100,
        "under-eroded body must fail the shrink predicate: {report:?}"
    );
}

#[test]
fn grid_volume_matches_steiner_on_the_sphere() {
    let source = sphere_index(1.0);
    let d = 0.3;
    let grid = grid_offset_volume(&source, d, 64);
    let pi = std::f64::consts::PI;
    let steiner = steiner_offset_volume(4.0 * pi / 3.0, 4.0 * pi, 4.0 * pi, 2.0, d);
    let exact = 4.0 * pi / 3.0 * (1.0 + d) * (1.0 + d) * (1.0 + d);
    assert!(
        (steiner - exact).abs() < 1e-12,
        "Steiner constants for the sphere are exact"
    );
    let rel = (grid - steiner).abs() / steiner;
    assert!(
        rel < 0.01,
        "grid-integrated dilation volume {grid} vs Steiner {steiner} (rel {rel})"
    );
}
