use super::*;
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_mesh::MeshVertex;

/// PK-1138: the mesher closes a curved rim with triangles whose three
/// vertices lie on one chart row: collinear in (u, v), they cover no
/// surface, while their facets span the arc's sagitta off it. Taken
/// flat, two spherical sectors' common read 1.6e-3 low with its area
/// exact. A cone's rim chord closed by such a triangle leaves the
/// face's volume and area terms as they were.
#[test]
fn a_triangle_collinear_in_uv_adds_only_its_slivers() {
    let tol = Tolerance::default();
    let brep = BRep::cone(4.0, 4.0, tol)
        .and_then(|c| c.transform(&Transform3::translation(1.0, 2.0, 3.0)))
        .expect("cone");
    let mesh = crate::body::tessellate_fine(&brep).expect("mesh");
    let (fm, surface) = mesh
        .face_meshes
        .iter()
        .find_map(|fm| {
            let s = &brep.faces[fm.face].surface.as_ref()?.geometry;
            matches!(s, SurfaceGeometry::Cone { .. }).then_some((fm, s))
        })
        .expect("cone face mesh");
    let gate = 100.0 * tol.linear();
    // Three points of the rim row, its arc's sagitta between them.
    let rim = fm
        .vertices
        .iter()
        .map(|v| v.uv.1)
        .fold(f64::NEG_INFINITY, f64::max);
    let mut closed = fm.clone();
    let first = closed.vertices.len();
    for u in [1.0, 1.05, 1.1] {
        closed.vertices.push(MeshVertex {
            position: surface.point_at(u, rim),
            uv: (u, rim),
            ..fm.vertices[0].clone()
        });
    }
    closed.triangles.push([first, first + 1, first + 2]);
    for kind in [Integrand::Volume, Integrand::Area] {
        let before = face_quadrature(&brep, fm, Some(surface), gate, kind);
        let after = face_quadrature(&brep, &closed, Some(surface), gate, kind);
        assert!(
            (after - before).abs() < 1e-12 * before.abs().max(1.0),
            "{kind:?}: {before} -> {after}"
        );
    }
}
