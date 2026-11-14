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

/// Sphere, cylinder and cone primitives far from the origin, their
/// meshes reading their faces, measure their closed forms to 1e-10
/// relative (PK-1186: the exact term stands in only for a misread face).
#[test]
fn analytic_faces_far_from_the_origin_read_their_exact_volume() {
    let tol = Tolerance::default();
    let pi = std::f64::consts::PI;
    let place = |b: BRep| {
        b.transform(&Transform3::rotation_x(0.7))
            .and_then(|b| b.transform(&Transform3::translation(100.0, -60.0, 35.0)))
            .expect("placed")
    };
    let cases = [
        (BRep::sphere(3.0, tol).expect("ball"), 4.0 / 3.0 * pi * 27.0),
        (BRep::cylinder(1.5, 4.0, tol).expect("rod"), pi * 2.25 * 4.0),
        (
            BRep::cone(2.0, 5.0, tol).expect("cone"),
            pi * 4.0 * 5.0 / 3.0,
        ),
    ];
    for (body, want) in cases {
        let v = crate::body_volume(&place(body), &tol).expect("volume");
        assert!((v - want).abs() < 1e-10 * want, "volume {v} for {want}");
    }
}

/// A measuring mesh that covers part of a sphere face twice (as the
/// fallback grid did beside a pole) moves the mesh's quadrature, and the
/// face takes its exact term, which reads only the mesh's orientation
/// (PK-1186).
#[test]
fn a_double_covering_mesh_leaves_the_exact_term() {
    let tol = Tolerance::default();
    let body = BRep::sphere(3.0, tol)
        .and_then(|b| b.transform(&Transform3::translation(5.0, 5.0, 13.0)))
        .expect("ball");
    let mut mesh = crate::body::tessellate_fine(&body).expect("mesh");
    let fm = &mut mesh.face_meshes[0];
    let extra: Vec<[usize; 3]> = fm
        .triangles
        .iter()
        .take(fm.triangles.len() / 10)
        .copied()
        .collect();
    fm.triangles.extend(extra);
    let fm = &mesh.face_meshes[0];
    let sof = body
        .faces
        .get(fm.face)
        .and_then(|f| f.surface.as_ref())
        .expect("surface");
    let mesh_area = || face_quadrature(&body, fm, Some(&sof.geometry), 1e-4, Integrand::Area);
    let exact = analytic::misread_face_term(&body, fm, sof, &mesh_area).expect("exact term");
    let quad = face_quadrature(&body, fm, Some(&sof.geometry), 1e-4, Integrand::Volume);
    let want = 4.0 / 3.0 * std::f64::consts::PI * 27.0;
    assert!(
        (exact.abs() - want).abs() < 1e-10 * want,
        "exact {exact} for {want}"
    );
    assert!(
        (quad.abs() - want).abs() > 1e-3 * want,
        "quadrature {quad} unmoved"
    );
}
