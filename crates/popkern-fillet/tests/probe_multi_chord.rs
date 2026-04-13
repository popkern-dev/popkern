//! Probe test to verify the multi-chord Phase B decomposition of
//! a 12-edge cube fillet produces the expected face decomposition.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn count_face_types(brep: &BRep) -> (usize, usize, usize, usize) {
    let mut planes = 0;
    let mut tori = 0;
    let mut spheres = 0;
    let mut others = 0;
    for fh in brep.faces.keys() {
        let face = &brep.faces[fh];
        match face.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Plane { .. }) => planes += 1,
            Some(SurfaceGeometry::Torus { .. }) => tori += 1,
            Some(SurfaceGeometry::Sphere { .. }) => spheres += 1,
            _ => others += 1,
        }
    }
    (planes, tori, spheres, others)
}

/// Expected: 6 plane centres + 12 torus strips + 8 sphere vertex
/// blends at cube corners = 26 faces.
#[test]
fn twelve_edge_cube_fillet_face_type_count() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).unwrap();
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let result = fillet_edges(&cube, &specs, &FilletOptions::default(), tol, None).unwrap();

    let (planes, tori, spheres, others) = count_face_types(&result.brep);
    assert_eq!(planes, 6, "expected 6 plane centres");
    assert_eq!(tori, 12, "expected 12 torus strip faces");
    assert_eq!(spheres, 8, "expected 8 sphere vertex-blend faces");
    assert_eq!(others, 0, "expected no other surface types");
}
