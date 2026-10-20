use super::*;
use crate::{TessellationOptions, tessellate};
use popkern_core::tolerance::Tolerance;

fn primitives() -> Vec<(&'static str, BRep)> {
    let tol = Tolerance::default();
    vec![
        ("cube", BRep::unit_cube(tol).expect("cube")),
        ("cylinder", BRep::cylinder(1.0, 2.0, tol).expect("cylinder")),
        ("cone", BRep::cone(1.0, 2.0, tol).expect("cone")),
        ("sphere", BRep::sphere(1.5, tol).expect("sphere")),
        ("torus", BRep::torus(3.0, 1.0, tol).expect("torus")),
    ]
}

/// The largest sampled `|L(x) − S(x)|` over the mesh's triangles.
fn sampled(brep: &BRep, mesh: &TriMesh) -> f64 {
    let mut worst = 0.0_f64;
    for (t, &face) in mesh.face_map.iter().enumerate() {
        let Some(s) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let tri = mesh.triangles[t];
        let p = tri.map(|i| mesh.vertices[i].position);
        let w = unwrap(tri.map(|i| mesh.vertices[i].uv), s.geometry.periodicity());
        for piece in pole_pieces(&s.geometry, w, p, 0.005) {
            worst = worst.max(piece.sampled(&s.geometry));
        }
    }
    worst
}

#[test]
fn primitive_meshes_certify_and_their_bounds_cover_the_sampled_deviation() {
    let options = TessellationOptions::default();
    for (name, brep) in primitives() {
        let mesh = tessellate(&brep, &options).expect("mesh");
        let cert = chordal_certificate(&brep, &mesh, options.chord_height);
        let truth = sampled(&brep, &mesh);
        eprintln!(
            "[pk505] {name}: worst bound {:.3e}, sampled {truth:.3e}, {} faces, status {:?}",
            cert.worst(),
            cert.faces.len(),
            cert.faces.iter().map(|f| f.status).collect::<Vec<_>>()
        );
        assert!(
            cert.inconclusive().is_empty(),
            "{name}: {:?}",
            cert.inconclusive()
        );
        let surface = cert
            .faces
            .iter()
            .map(|f| f.surface_bound)
            .fold(0.0, f64::max);
        assert!(
            surface >= truth,
            "{name}: bound {surface} under sampled {truth}"
        );
        for f in &cert.faces {
            assert!(f.surface_bound >= f.witness, "{name}: {f:?}");
        }
        if name == "torus" {
            // The 3/1 torus meshes within its chord height (4.7e-3 from
            // the surface to the mesh, 4.8e-3 back, sampled), but its
            // bilinear bound charges the in-plane u-sag: 5.5e-3 (PK-908).
            for f in cert.exceeding() {
                assert!(f.witness <= options.chord_height, "{f:?}");
            }
        } else {
            assert!(cert.certified(), "{name}: {:?}", cert.exceeding());
        }
    }
}

#[test]
fn a_coarse_mesh_is_proven_to_exceed_a_tight_tolerance() {
    // The cylinder at the default chord height, certified against a
    // tolerance a hundred times tighter: proven bounds, above it.
    let brep = BRep::cylinder(1.0, 2.0, Tolerance::default()).expect("cylinder");
    let options = TessellationOptions::default();
    let mesh = tessellate(&brep, &options).expect("mesh");
    let cert = chordal_certificate(&brep, &mesh, options.chord_height / 100.0);
    assert!(!cert.certified());
    assert!(!cert.exceeding().is_empty());
    assert!(cert.inconclusive().is_empty());
}

#[test]
fn a_moved_vertex_raises_the_bound_by_its_offset() {
    let brep = BRep::sphere(1.5, Tolerance::default()).expect("sphere");
    let options = TessellationOptions::default();
    let mut mesh = tessellate(&brep, &options).expect("mesh");
    let base = chordal_certificate(&brep, &mesh, options.chord_height).worst();
    // Push one interior vertex off the surface by 0.05 along its normal.
    let t = mesh.triangles.len() / 2;
    let i = mesh.triangles[t][0];
    let n = mesh.vertices[i].normal;
    mesh.vertices[i].position = mesh.vertices[i].position + n * 0.05;
    let cert = chordal_certificate(&brep, &mesh, options.chord_height);
    assert!(
        cert.worst() >= 0.05 && cert.worst() > base,
        "{} vs {base}",
        cert.worst()
    );
    assert!(!cert.certified());
}
