//! §5.3 item #3 extended — Karciauskas-Peters biquintic-full
//! Gregory dispatch integration test.
//!
//! Verifies that opting into `FilletOptions::use_biquintic_full_gregory`
//! on a mixed-radius trihedral vertex routes the three-quad
//! sub-patches through `RectangularGregoryBiquinticFull` instead of
//! the default `RectangularGregoryDeg5`, without regressing the
//! downstream face / mesh contract.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn fillet_mixed_corner(use_biquintic_full: bool) -> BRep {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    assert_eq!(incident_edges.len(), 3);

    let radii = [0.08, 0.10, 0.12];
    let specs: Vec<FilletSpec> = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let opts = FilletOptions {
        use_rational_gregory: true,
        use_deg5_gregory: true,
        use_biquintic_full_gregory: use_biquintic_full,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&cube, &specs, &opts, Tolerance::default(), None)
        .expect("fillet_edges should succeed on mixed-radius cube corner");
    result.brep
}

/// Flag off → baseline deg5 Gregory sub-patches. Flag on → biquintic-full
/// variant appears on at least one face. Face/topology counts agree
/// (the solver does not change the surface topology, only the interior
/// control-point geometry).
#[test]
fn biquintic_full_flag_routes_to_biquintic_full_variant() {
    let brep_deg5 = fillet_mixed_corner(false);
    let brep_full = fillet_mixed_corner(true);

    assert_eq!(
        brep_deg5.faces.len(),
        brep_full.faces.len(),
        "biquintic-full opt-in must not change face count"
    );

    let mut deg5_count = 0usize;
    let mut biquintic_full_count = 0usize;
    for (_, face) in brep_deg5.faces.iter() {
        if let Some(sof) = face.surface.as_ref()
            && let SurfaceGeometry::Gregory(gp) = &sof.geometry
            && matches!(
                gp,
                popkern_core::geometry::gregory::GregoryPatch::RectangularDeg5(_)
            )
        {
            deg5_count += 1;
        }
    }
    for (_, face) in brep_full.faces.iter() {
        if let Some(sof) = face.surface.as_ref()
            && let SurfaceGeometry::Gregory(gp) = &sof.geometry
            && matches!(
                gp,
                popkern_core::geometry::gregory::GregoryPatch::RectangularBiquinticFull(_)
            )
        {
            biquintic_full_count += 1;
        }
    }
    assert!(
        deg5_count >= 1,
        "baseline (flag off) must produce at least one RectangularDeg5 sub-patch"
    );
    assert_eq!(
        biquintic_full_count, deg5_count,
        "flag on must replace every deg5 sub-patch with the biquintic-full variant"
    );
}

/// The biquintic-full opt-in must not break the BRep's validation
/// invariants (Euler characteristic, manifold edges, etc). Uses the
/// cheap `validate()` path (structural + geometric, no self-intersection).
#[test]
fn biquintic_full_dispatch_produces_valid_brep() {
    let brep = fillet_mixed_corner(true);
    brep.validate()
        .expect("biquintic-full vertex blend must yield a valid BRep");
}
