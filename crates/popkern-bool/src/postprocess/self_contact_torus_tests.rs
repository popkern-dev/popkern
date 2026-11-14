use super::*;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::EdgeHandle;

fn edges_of(brep: &BRep, face: FaceHandle) -> Vec<EdgeHandle> {
    let mut out = Vec::new();
    let f = &brep.faces[face];
    for l in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let start = brep.loops[l].half_edge;
        let mut h = start;
        loop {
            out.push(brep.half_edges[h].edge);
            h = brep.half_edges[h].next;
            if h == start {
                break;
            }
        }
    }
    out
}

/// The cosine's range over an interval reaches 1 and -1 exactly where a
/// multiple of the period or a half period lies inside.
#[test]
fn cosine_ranges() {
    let pi = std::f64::consts::PI;
    let close =
        |(a, b): (f64, f64), (c, d): (f64, f64)| (a - c).abs() < 1e-12 && (b - d).abs() < 1e-12;
    assert!(close(cos_range(0.1, 0.2), (0.2f64.cos(), 0.1f64.cos())));
    assert!(close(cos_range(-0.1, 0.2), (0.2f64.cos(), 1.0)));
    assert!(close(
        cos_range(3.0, 3.5),
        (-1.0, 3.0f64.cos().max(3.5f64.cos()))
    ));
    assert!(close(
        cos_range(2.0 * pi - 0.1, 2.0 * pi + 0.1),
        ((0.1f64).cos(), 1.0)
    ));
    assert!(close(
        cos_range(-3.5, -3.0),
        (-1.0, 3.0f64.cos().max(3.5f64.cos()))
    ));
    assert!(close(cos_range(0.0, 7.0), (-1.0, 1.0)));
}

/// Every face of the primitive torus, turned, lies in its closed-form
/// box (dense uv samples over its boundary's uv box), and each box is
/// smaller than the whole torus's.
#[test]
fn a_torus_faces_box_holds_it() {
    let tol = Tolerance::default();
    let turn = Transform3::rotation_x(0.7).compose(&Transform3::rotation_z(-1.3));
    let torus = BRep::torus(1.0, 0.3, tol)
        .and_then(|b| b.transform(&turn))
        .expect("torus");
    for face in torus.faces.keys() {
        let surface = &torus.faces[face]
            .surface
            .as_ref()
            .expect("surface")
            .geometry;
        let bx = torus_face_box(&torus, face, surface).expect("box");
        let tau = std::f64::consts::TAU;
        let mut held = 0;
        for i in 0..=120 {
            for j in 0..=120 {
                let (u, v) = (tau * f64::from(i) / 120.0, tau * f64::from(j) / 120.0);
                if torus.face_uv_trim_containment(face, surface, u, v)
                    != Some(popkern_core::validation::FaceUvContainment::Interior)
                {
                    continue;
                }
                held += 1;
                let p = surface.point_at(u, v);
                assert!(
                    p.x >= bx.min.x
                        && p.y >= bx.min.y
                        && p.z >= bx.min.z
                        && p.x <= bx.max.x
                        && p.y <= bx.max.y
                        && p.z <= bx.max.z,
                    "{p:?} outside {bx:?}"
                );
            }
        }
        assert!(held > 100, "{held} samples on the face");
        let size = bx.min.distance_to(&bx.max);
        assert!(size < 0.9 * 2.0 * 3f64.sqrt() * 1.3, "box {size}");
    }
}

/// A box over the torus clears every plane x torus pair when its bottom
/// lies a tenth above the torus's top, and not when it rests on the top
/// (the plane is tangent to the tube's upper circle, where the normals
/// are parallel).
#[test]
fn a_plane_over_a_torus_is_cleared_unless_it_rests_on_it() {
    let tol = Tolerance::default();
    let linear = tol.linear();
    for (gap, cleared) in [(0.1, true), (0.0, false)] {
        let torus = BRep::torus(1.0, 0.3, tol).expect("torus");
        let slab = BRep::prism(
            &[
                Point3::new(-2.0, -2.0, 0.3 + gap),
                Point3::new(2.0, -2.0, 0.3 + gap),
                Point3::new(2.0, 2.0, 0.3 + gap),
                Point3::new(-2.0, 2.0, 0.3 + gap),
            ],
            Vec3::new(0.0, 0.0, 1.0),
            1.0,
            tol,
        )
        .expect("slab");
        let both = crate::operation::helpers::combine_breps(torus, slab).expect("both");
        let mut decided = Vec::new();
        for f in both.faces.keys() {
            let sf = &both.faces[f].surface.as_ref().expect("surface").geometry;
            if !matches!(sf, SurfaceGeometry::Plane { .. }) {
                continue;
            }
            let (box_f, edges_f) = (
                super::super::containing_box(&both, f).expect("box"),
                edges_of(&both, f),
            );
            for g in both.faces.keys() {
                let sg = &both.faces[g].surface.as_ref().expect("surface").geometry;
                if !matches!(sg, SurfaceGeometry::Torus { .. }) {
                    continue;
                }
                let (box_g, edges_g) = (
                    super::super::containing_box(&both, g).expect("box"),
                    edges_of(&both, g),
                );
                decided.push(plane_torus_clear(
                    &both,
                    (f, sf, &box_f, &edges_f),
                    (g, sg, &box_g, &edges_g),
                    linear,
                ));
            }
        }
        assert_eq!(decided.len(), 6 * 16);
        if cleared {
            assert!(decided.iter().all(|&c| c), "gap {gap}");
        } else {
            assert!(!decided.iter().all(|&c| c), "gap {gap}");
        }
    }
}

/// The least covering arc: intervals either side of the seam join
/// across it, a gap anywhere else is kept out, and intervals covering
/// the circle give the whole period.
#[test]
fn covering_arcs() {
    let tau = std::f64::consts::TAU;
    let close =
        |(a, b): (f64, f64), (c, d): (f64, f64)| (a - c).abs() < 1e-12 && (b - d).abs() < 1e-12;
    // One edge at 0, another at the period, the rest beside them.
    let beside_seam = [(5.0, tau), (0.0, 0.5), (tau, tau), (0.0, 0.0)];
    assert!(close(
        covering_arc(beside_seam.into_iter(), tau).expect("arc"),
        (5.0, tau + 0.5)
    ));
    assert!(close(
        covering_arc([(1.0, 2.0), (1.5, 3.0)].into_iter(), tau).expect("arc"),
        (1.0, 3.0)
    ));
    assert!(close(
        covering_arc([(0.0, 3.5), (3.0, tau)].into_iter(), tau).expect("arc"),
        (0.0, tau)
    ));
    assert!(close(
        covering_arc([(-1.3, 0.27)].into_iter(), tau).expect("arc"),
        (tau - 1.3, tau + 0.27)
    ));
}
