//! PK-787: sphere faces keep their material side through a popkern STEP
//! round trip.
//!
//! The kernel's orientation flag says a face's outer loop runs
//! counter-clockwise in its surface chart; STEP's ADVANCED_FACE
//! `same_sense` says it runs counter-clockwise about the written
//! surface's normal. The sphere chart (colatitude) has an inward
//! `∂u × ∂v` while SPHERICAL_SURFACE's normal points outward, so the flag
//! and `same_sense` must be negations. The importer negated them; the
//! exporter did not, so every exported sphere face came back reversed
//! and its tessellation filled the complement of the trim (boss_plate's
//! corner patches meshed 0.0703 mm² for a 0.01005 mm² face).

use std::sync::Arc;

use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::geometry::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step, import_step,
};

fn document_with_body(brep: BRep) -> StepDocument {
    StepDocument {
        products: vec![StepProductNode {
            name: "body".into(),
            id: "body".into(),
            children: Vec::new(),
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::Other("SHAPE_REPRESENTATION".into()),
                solids: vec![brep],
                context: RepresentationContext::default(),
                non_brep_items: Vec::new(),
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    }
}

/// Mesh area of every sphere face, and the whole mesh's volume.
fn sphere_area_and_volume(brep: &BRep) -> (f64, f64) {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .expect("tessellates");
    let mut sphere_area = 0.0;
    for fm in &mesh.face_meshes {
        let is_sphere = brep.faces[fm.face]
            .surface
            .as_ref()
            .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Sphere { .. }));
        if !is_sphere {
            continue;
        }
        for t in &fm.triangles {
            let p = |i: usize| fm.vertices[t[i]].position;
            sphere_area += 0.5 * (p(1) - p(0)).cross(&(p(2) - p(0))).length();
        }
    }
    let volume: f64 = mesh
        .triangles
        .iter()
        .map(|t| {
            let p = |i: usize| mesh.vertices[t[i]].position;
            let (a, b, c) = (p(0), p(1), p(2));
            (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                + a.z * (b.x * c.y - b.y * c.x))
                / 6.0
        })
        .sum();
    (sphere_area, volume.abs())
}

fn round_trip(brep: &BRep) -> BRep {
    let bytes = export_step(
        &document_with_body(brep.clone()),
        &StepExportOptions::default(),
    )
    .expect("exports");
    let imported = import_step(&bytes, &StepImportOptions::default()).expect("imports");
    imported.document.products[0].shapes[0].solids[0].clone()
}

fn assert_sphere_faces_survive(name: &str, native: &BRep) {
    let (native_area, native_volume) = sphere_area_and_volume(native);
    assert!(native_area > 0.0, "{name}: fixture has sphere faces");
    let imported = round_trip(native);
    let (area, volume) = sphere_area_and_volume(&imported);
    assert!(
        (area - native_area).abs() <= 0.02 * native_area,
        "{name}: imported sphere faces mesh {area:.5} vs native {native_area:.5}"
    );
    assert!(
        (volume - native_volume).abs() <= 1e-3 * native_volume,
        "{name}: imported volume {volume:.6} vs native {native_volume:.6}"
    );
}

#[test]
fn filleted_cube_corner_patches_round_trip_with_their_material_side() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .map(|edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let filleted = fillet_edges(&cube, &specs, &FilletOptions::default(), tol, None)
        .expect("all-edge fillet")
        .brep;
    assert_sphere_faces_survive("filleted cube", &filleted);
}
