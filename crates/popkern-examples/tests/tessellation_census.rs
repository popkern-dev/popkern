//! PK-259 acceptance pin: per-face tessellation census on a served
//! fillet-corner body, through the same serde round-trip the viz
//! loader uses. Every surfaced production face must receive
//! triangles. (The bare-cube renders that motivated PK-259 were the
//! viz screenshot CLI stub rendering a hardcoded unit cube —
//! tessellation was measured innocent by this census, and the pin
//! keeps it that way.)

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec};

#[test]
fn fillet_corner_face_triangle_census() {
    let tol = Tolerance::default();
    let brep = BRep::unit_cube(tol).unwrap();
    let corner = brep
        .vertices
        .iter()
        .find(|(_, v)| {
            v.position.x.abs() < 1e-9 && v.position.y.abs() < 1e-9 && v.position.z.abs() < 1e-9
        })
        .map(|(vh, _)| vh)
        .expect("corner");
    let mut edges = Vec::new();
    for (e, d) in brep.edges.iter() {
        let he = &brep.half_edges[d.half_edge];
        let a = he.origin;
        let b = brep.half_edges[he.twin].origin;
        if a == corner || b == corner {
            edges.push(e);
        }
    }
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip([0.15, 0.1, 0.1])
        .map(|(&edge, r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let result = popkern_fillet::fillet_edges(&brep, &specs, &FilletOptions::default(), tol, None)
        .expect("cube fillet");
    // The viz screenshot path loads the body from BRep JSON — run the
    // census on the round-tripped body to cover the same surface.
    let json = serde_json::to_string(&result.brep).expect("serialize");
    let body: BRep = serde_json::from_str(&json).expect("deserialize");
    let body = &body;
    let mesh = popkern_mesh::tessellate(body, &popkern_mesh::TessellationOptions::default())
        .expect("tessellate");
    let mut lines = Vec::new();
    let mut missing = Vec::new();
    for (face, data) in body.faces.iter() {
        let Some(surface) = data.surface.as_ref() else {
            continue;
        };
        let kind = format!("{:?}", std::mem::discriminant(&surface.geometry));
        let triangles = mesh
            .face_meshes
            .iter()
            .filter(|fm| fm.face == face)
            .map(|fm| fm.triangles.len())
            .sum::<usize>();
        lines.push(format!(
            "  {face:?} role={:?} surface={kind} triangles={triangles}",
            data.role
        ));
        if triangles == 0 {
            missing.push(face);
        }
    }
    eprintln!("[pk259-census]\n{}", lines.join("\n"));
    assert!(
        missing.is_empty(),
        "every surfaced face must receive triangles; empty: {missing:?}"
    );
}
