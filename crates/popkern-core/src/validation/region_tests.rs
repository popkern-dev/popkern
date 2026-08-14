//! Tests for region-scoped fast validation (PK-273).

use crate::geometry::Point3;
use crate::geometry::gregory::{GregoryPatch, RectangularGregory};
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle};

/// A folded rectangular Gregory (opposite corners swapped) — the
/// same construction the Gregory self-intersection tests use.
fn folded_gregory() -> GregoryPatch {
    let lin = |u: f64, v: f64| Point3::new(u, v, 0.0);
    let mut boundary_cps = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for (i, row) in boundary_cps.iter_mut().enumerate() {
        for (j, cp) in row.iter_mut().enumerate() {
            *cp = lin(i as f64 / 3.0, j as f64 / 3.0);
        }
    }
    let mut inner_cps = [
        (lin(1.0 / 3.0, 1.0 / 3.0), lin(1.0 / 3.0, 1.0 / 3.0)),
        (lin(1.0 / 3.0, 2.0 / 3.0), lin(1.0 / 3.0, 2.0 / 3.0)),
        (lin(2.0 / 3.0, 2.0 / 3.0), lin(2.0 / 3.0, 2.0 / 3.0)),
        (lin(2.0 / 3.0, 1.0 / 3.0), lin(2.0 / 3.0, 1.0 / 3.0)),
    ];
    boundary_cps[0][0] = Point3::new(1.0, 1.0, 0.0);
    boundary_cps[3][3] = Point3::new(0.0, 0.0, 0.0);
    inner_cps[0].0 = Point3::new(2.0 / 3.0, 2.0 / 3.0, 0.0);
    inner_cps[0].1 = Point3::new(2.0 / 3.0, 2.0 / 3.0, 0.0);
    inner_cps[2].0 = Point3::new(1.0 / 3.0, 1.0 / 3.0, 0.0);
    inner_cps[2].1 = Point3::new(1.0 / 3.0, 1.0 / 3.0, 0.0);
    GregoryPatch::Rectangular(Box::new(RectangularGregory {
        boundary_cps,
        inner_cps,
    }))
}

fn put_folded_gregory_on(brep: &mut BRep, face: FaceHandle) {
    brep.faces[face].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Gregory(folded_gregory()),
        orientation: true,
    });
}

/// On a cube, a face's one-ring by shared vertices is every face
/// except the opposite one.
#[test]
fn validation_region_dilates_by_shared_vertices() {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let dirty = cube.faces.keys().next().unwrap();
    let region = cube.validation_region(&[dirty]);
    assert_eq!(
        region.len(),
        5,
        "cube face + 4 vertex-sharing neighbours; opposite face excluded"
    );
    assert!(region.contains(&dirty));
}

#[test]
fn validation_region_ignores_dead_handles_and_empty_input() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let dirty = cube.faces.keys().next().unwrap();
    cube.faces.remove(dirty);
    assert!(cube.validation_region(&[dirty]).is_empty());
    assert!(cube.validation_region(&[]).is_empty());
}

/// Row 404: a folded Gregory on a NEIGHBOUR of the dirty face must
/// be caught — the one-ring dilation is what finds it.
#[test]
fn region_scoped_sweep_catches_boundary_straddling_defect() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let dirty = cube.faces.keys().next().unwrap();
    let region = cube.validation_region(&[dirty]);
    let neighbour = *region
        .iter()
        .find(|f| **f != dirty)
        .expect("cube face has neighbours");
    put_folded_gregory_on(&mut cube, neighbour);

    let dilated = cube.validation_region(&[dirty]);
    assert!(
        cube.validate_gregory_self_intersection_among(&dilated)
            .is_err(),
        "one-ring dilation must reach the neighbour's folded patch"
    );
    assert!(
        cube.validate_gregory_self_intersection_among(&[dirty])
            .is_ok(),
        "without dilation the defect is invisible — the dilation is load-bearing"
    );
}

/// Documented containment: a defect wholly OUTSIDE the region is not
/// reported by the scoped sweep but IS reported by the full one.
#[test]
fn region_scoped_sweep_documented_containment_outside_region() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let dirty = cube.faces.keys().next().unwrap();
    let region = cube.validation_region(&[dirty]);
    let opposite = cube
        .faces
        .keys()
        .find(|f| !region.contains(f))
        .expect("cube has an opposite face outside the one-ring");
    put_folded_gregory_on(&mut cube, opposite);

    let dilated = cube.validation_region(&[dirty]);
    assert!(
        cube.validate_gregory_self_intersection_among(&dilated)
            .is_ok(),
        "scoped sweep must not see the opposite face"
    );
    assert!(
        cube.validate_gregory_self_intersection().is_err(),
        "the full sweep remains the backstop"
    );
}

/// Parity: over the full face set, the region variant reports
/// exactly what `validate_collect_fast` reports.
#[test]
fn region_collect_fast_over_all_faces_matches_whole_model() {
    let mut cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let some_face = cube.faces.keys().nth(2).unwrap();
    put_folded_gregory_on(&mut cube, some_face);

    let all_faces: Vec<FaceHandle> = cube.faces.keys().collect();
    let scoped = cube.validate_region_collect_fast(&all_faces);
    let whole = cube.validate_collect_fast();
    assert_eq!(format!("{scoped:?}"), format!("{whole:?}"));
    assert!(
        whole
            .iter()
            .any(|e| matches!(e, crate::error::ValidationError::TwistedSurface { face } if *face == some_face)),
        "the folded patch must be reported by both"
    );
}

/// The verdict cache returns identical reports on repeat calls and
/// records hits (PK-273 structured observability).
#[test]
fn gregory_verdict_cache_hits_on_repeat_and_preserves_verdicts() {
    GregoryPatch::self_intersection_cache_reset();
    let folded = folded_gregory();
    let first = folded.self_intersection_report(1e-6);
    let (hits_0, misses_0) = GregoryPatch::self_intersection_cache_stats();
    assert_eq!(hits_0, 0);
    assert!(misses_0 >= 1);
    let second = folded.self_intersection_report(1e-6);
    assert_eq!(first, second);
    let (hits_1, _) = GregoryPatch::self_intersection_cache_stats();
    assert!(hits_1 >= 1, "repeat call with identical content must hit");
    assert!(folded.self_intersects_with_tolerance(1e-6));

    // A different tolerance is a different key — no stale verdict.
    let (_, misses_before) = GregoryPatch::self_intersection_cache_stats();
    let _ = folded.self_intersection_report(1e-7);
    let (_, misses_after) = GregoryPatch::self_intersection_cache_stats();
    assert!(misses_after > misses_before);
    GregoryPatch::self_intersection_cache_reset();
}
