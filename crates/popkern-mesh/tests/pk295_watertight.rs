//! PK-295 regression pin: the served rim-filleted-cylinder minus
//! strip-crossing-cuboid subtraction must TESSELLATE WATERTIGHT.
//! Promoted from the PK-295 diagnostic probe — the historical result
//! shipped 5 faces and 20 open mesh edges with a clean validate_all,
//! and mid-campaign states shipped 57 open edges along the seam-side
//! bulge until the slit-seam dissolution pass landed.

use popkern_core::geometry::{Point3, Transform3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use std::collections::HashMap;

fn rim_filleted_cylinder(tol: Tolerance) -> BRep {
    let brep = BRep::cylinder(1.0, 2.0, tol).expect("cylinder");
    let mut rim = Vec::new();
    for (e, d) in brep.edges.iter() {
        let he = &brep.half_edges[d.half_edge];
        let a = brep.vertices[he.origin].position;
        let b = brep.vertices[brep.half_edges[he.twin].origin].position;
        let f1 = brep.loops[he.loop_].face;
        let f2 = brep.loops[brep.half_edges[he.twin].loop_].face;
        if f1 != f2 && (a.z - 2.0).abs() < 1e-9 && (b.z - 2.0).abs() < 1e-9 {
            rim.push(e);
        }
    }
    let specs: Vec<FilletSpec> = rim
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.2),
        })
        .collect();
    fillet_edges(&brep, &specs, &FilletOptions::default(), tol, None)
        .expect("rim fillet must complete")
        .brep
}

#[test]
fn pk295_served_subtraction_watertightness() {
    let tol = Tolerance::default();
    let part = rim_filleted_cylinder(tol);
    let cutter = BRep::unit_cube(tol)
        .expect("cube")
        .transform(&Transform3::scale_uniform(1.8))
        .expect("scale")
        .transform(&Transform3::translation(-0.9, -0.9, 0.6))
        .expect("translate");
    let result = popkern_bool::subtraction(part, cutter, tol, None).expect("difference");

    let mesh = popkern_mesh::tessellate(&result, &popkern_mesh::TessellationOptions::default())
        .expect("tessellation");
    let mut key: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut canon: Vec<Point3> = Vec::new();
    let id = |p: Point3, key: &mut HashMap<(i64, i64, i64), usize>, canon: &mut Vec<Point3>| {
        let k = (
            (p.x * 1.0e6).round() as i64,
            (p.y * 1.0e6).round() as i64,
            (p.z * 1.0e6).round() as i64,
        );
        *key.entry(k).or_insert_with(|| {
            canon.push(p);
            canon.len() - 1
        })
    };
    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        let w = [
            id(mesh.vertices[t[0]].position, &mut key, &mut canon),
            id(mesh.vertices[t[1]].position, &mut key, &mut canon),
            id(mesh.vertices[t[2]].position, &mut key, &mut canon),
        ];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    let open = counts.values().filter(|&&c| c == 1).count();
    assert_eq!(
        open, 0,
        "served subtraction must tessellate watertight (PK-295), {open} open mesh edges"
    );
}
