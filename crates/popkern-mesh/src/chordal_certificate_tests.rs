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

/// A sampled two-sided distance between each face and its triangles:
/// surface samples to the face's nearest triangle, and triangle samples
/// to their projection on the surface. Bent pairings (PK-908) bound
/// this, not the affine `|L(x) − S(x)|`. Surface samples come from every
/// fourth triangle (a maximum over fewer samples is still a lower
/// bound; each is measured against every triangle).
fn sampled(brep: &BRep, mesh: &TriMesh) -> f64 {
    use popkern_core::geometry::SurfaceEval;
    let mut worst = 0.0_f64;
    for f in &mesh.face_meshes {
        let Some(s) = brep.faces.get(f.face).and_then(|x| x.surface.as_ref()) else {
            continue;
        };
        let tris: Vec<usize> = (0..mesh.triangles.len())
            .filter(|&t| mesh.face_map[t] == f.face)
            .collect();
        let all: Vec<[Point3; 3]> = tris
            .iter()
            .map(|&t| mesh.triangles[t].map(|i| mesh.vertices[i].position))
            .collect();
        for (k, (&t, p)) in tris.iter().zip(&all).enumerate() {
            let tri = mesh.triangles[t];
            let w = unwrap(tri.map(|i| mesh.vertices[i].uv), s.geometry.periodicity());
            for piece in pole_pieces(&s.geometry, w, *p, 0.005)
                .iter()
                .filter(|_| k % 4 == 0)
            {
                worst = worst.max(piece.sampled_against(&s.geometry, &all));
            }
            let n = 4;
            for i in 0..=n {
                for j in 0..=n - i {
                    let (b1, b2) = (i as f64 / n as f64, j as f64 / n as f64);
                    let b0 = 1.0 - b1 - b2;
                    let q = Point3::new(
                        b0 * p[0].x + b1 * p[1].x + b2 * p[2].x,
                        b0 * p[0].y + b1 * p[1].y + b2 * p[2].y,
                        b0 * p[0].z + b1 * p[1].z + b2 * p[2].z,
                    );
                    if let Some((u, v)) = s.geometry.closest_point(&q) {
                        worst = worst.max(q.distance_to(&s.geometry.point_at(u, v)));
                    }
                }
            }
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
        // The 3/1 torus meshes within its chord height (4.7e-3 from the
        // surface to the mesh, 4.8e-3 back, sampled); its affine and
        // bilinear pairings charged the in-plane u-sag (5.5e-3), the bent
        // one does not (PK-908).
        assert!(cert.certified(), "{name}: {:?}", cert.exceeding());
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

/// PK-905: faces are certified in parallel and collected in order, so
/// the certificate is the same at any thread count.
#[cfg(feature = "parallel")]
#[test]
fn thread_count_does_not_change_the_certificate() {
    let brep = BRep::torus(3.0, 1.0, Tolerance::default()).expect("torus");
    let options = TessellationOptions::default();
    let mesh = tessellate(&brep, &options).expect("mesh");
    let at = |threads: usize| {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("pool")
            .install(|| chordal_certificate(&brep, &mesh, options.chord_height))
    };
    let (single, many) = (at(1), at(8));
    assert_eq!(single, many);
    assert!(single.faces.len() > 1);
}

#[test]
fn closed_form_feet_are_the_nearest_points_in_the_window() {
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
    let arc = CurveGeometry::Arc {
        center: Point3::new(1.0, -2.0, 0.5),
        normal: Vec3::new(0.0, 0.6, 0.8),
        radius: 1.5,
        start_angle: 0.0,
        end_angle: 1.0,
    };
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: Vec3::new(2.0, 0.0, 0.0),
    };
    // A window on another turn than the arc's span, and a line window.
    for (curve, window) in [
        (
            &arc,
            (std::f64::consts::TAU + 0.2, std::f64::consts::TAU + 0.9),
        ),
        (&arc, (-3.0, -2.5)),
        (&line, (-1.0, 2.0)),
    ] {
        for k in 0..40 {
            let s = k as f64 * 0.37;
            let p = Point3::new(s.sin() * 3.0, s.cos() * 2.0 - 2.0, (2.0 * s).sin());
            let t = chordal_boundary::window_foot(curve, p, window);
            assert!((window.0..=window.1).contains(&t), "{t} outside {window:?}");
            let sampled = (0..=2000)
                .map(|j| window.0 + (window.1 - window.0) * j as f64 / 2000.0)
                .map(|x| curve.point_at(x).distance_to(&p))
                .fold(f64::INFINITY, f64::min);
            assert!(
                curve.point_at(t).distance_to(&p) <= sampled + 1e-12,
                "{p:?}"
            );
        }
    }
}
