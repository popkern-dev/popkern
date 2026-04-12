//! Tests that trigger fillet, offset, and swept detection and repair paths,
//! exercising the repair bodies that cube-based tests cannot reach.

use super::*;
use crate::constructors::builder::BRepBuilder;
use crate::geometry::Point3;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::tolerance::Tolerance;

// -- Geometry builders -------------------------------------------------------

/// Three coplanar quads on z=0 sharing edges. Middle face has G1 with
/// both neighbors (all normals +Z), so it is detected as a fillet.
fn make_three_coplanar_faces() -> BRep {
    let mut b = BRepBuilder::new(Tolerance::default());

    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v4 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v8 = b.add_vertex(Point3::new(2.0, 0.0, 0.0));
    let v12 = b.add_vertex(Point3::new(3.0, 0.0, 0.0));
    let v3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));
    let v7 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let v11 = b.add_vertex(Point3::new(2.0, 1.0, 0.0));
    let v15 = b.add_vertex(Point3::new(3.0, 1.0, 0.0));

    let solid = b.add_solid();
    let shell = b.add_shell(solid);

    let (_, he_04, _) = b.add_edge(v0, v4);
    let (_, he_47, he_74) = b.add_edge(v4, v7);
    let (_, he_73, _) = b.add_edge(v7, v3);
    let (_, he_30, _) = b.add_edge(v3, v0);
    let (_, he_48, _) = b.add_edge(v4, v8);
    let (_, he_8_11, he_11_8) = b.add_edge(v8, v11);
    let (_, he_11_7, _) = b.add_edge(v11, v7);
    let (_, he_8_12, _) = b.add_edge(v8, v12);
    let (_, he_12_15, _) = b.add_edge(v12, v15);
    let (_, he_15_11, _) = b.add_edge(v15, v11);

    let face_a = b.add_face(shell);
    b.add_loop(&[he_04, he_47, he_73, he_30], face_a)
        .expect("face A loop");

    let face_b = b.add_face(shell);
    b.add_loop(&[he_48, he_8_11, he_11_7, he_74], face_b)
        .expect("face B loop");

    let face_c = b.add_face(shell);
    b.add_loop(&[he_8_12, he_12_15, he_15_11, he_11_8], face_c)
        .expect("face C loop");

    for (fh, x_off) in [(face_a, 0.0), (face_b, 1.0), (face_c, 2.0)] {
        b.set_face_surface(
            fh,
            SurfaceOnFace {
                geometry: SurfaceGeometry::Nurbs(flat_nurbs(x_off, 0.0, 1.0, 1.0)),
                orientation: true,
            },
        );
    }
    b.build_unchecked()
}

/// Two faces sharing an edge, one at z=0 and one at z=`d` (constant offset).
fn make_offset_pair(d: f64) -> BRep {
    let mut b = BRepBuilder::new(Tolerance::default());

    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let v3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));
    let v4 = b.add_vertex(Point3::new(0.0, 0.0, d));
    let v5 = b.add_vertex(Point3::new(1.0, 0.0, d));

    let solid = b.add_solid();
    let shell = b.add_shell(solid);

    let (_, he_01, he_10) = b.add_edge(v0, v1);
    let (_, he_12, _) = b.add_edge(v1, v2);
    let (_, he_23, _) = b.add_edge(v2, v3);
    let (_, he_30, _) = b.add_edge(v3, v0);
    let (_, he_04, _) = b.add_edge(v0, v4);
    let (_, he_45, _) = b.add_edge(v4, v5);
    let (_, he_51, _) = b.add_edge(v5, v1);

    let face_lo = b.add_face(shell);
    b.add_loop(&[he_01, he_12, he_23, he_30], face_lo)
        .expect("lower loop");

    let face_hi = b.add_face(shell);
    b.add_loop(&[he_10, he_04, he_45, he_51], face_hi)
        .expect("upper loop");

    b.set_face_surface(
        face_lo,
        SurfaceOnFace {
            geometry: SurfaceGeometry::Nurbs(flat_nurbs(0.0, 0.0, 1.0, 1.0)),
            orientation: true,
        },
    );
    b.set_face_surface(
        face_hi,
        SurfaceOnFace {
            geometry: SurfaceGeometry::Nurbs(flat_nurbs_z(0.0, 0.0, d, 1.0, 1.0)),
            orientation: true,
        },
    );
    b.build_unchecked()
}

/// Single face with a 5x5 swept-like NURBS (extrusion along z, arc profile).
fn make_swept_face() -> BRep {
    let mut b = BRepBuilder::new(Tolerance::default());
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let v3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));
    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let (_, he01, _) = b.add_edge(v0, v1);
    let (_, he12, _) = b.add_edge(v1, v2);
    let (_, he23, _) = b.add_edge(v2, v3);
    let (_, he30, _) = b.add_edge(v3, v0);
    let face = b.add_face(shell);
    b.add_loop(&[he01, he12, he23, he30], face)
        .expect("swept loop");

    // Each u-row identical profile in XY, offset in z by small amount.
    let (nu, nv) = (5, 5);
    let u_span = 0.08; // small => shape deviation < 0.1
    let mut cps = Vec::with_capacity(nu * nv);
    let mut wts = Vec::with_capacity(nu * nv);
    for i in 0..nu {
        let z = u_span * i as f64 / (nu - 1) as f64;
        for j in 0..nv {
            let t = j as f64 / (nv - 1) as f64;
            cps.push(Point3::new(t, 0.2 * (std::f64::consts::PI * t).sin(), z));
            wts.push(1.0);
        }
    }
    b.set_face_surface(
        face,
        SurfaceOnFace {
            geometry: SurfaceGeometry::Nurbs(NurbsSurface {
                degree_u: 2,
                degree_v: 2,
                knots_u: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
                knots_v: vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0],
                control_points: cps,
                num_u: nu,
                num_v: nv,
                weights: wts,
            }),
            orientation: true,
        },
    );
    b.build_unchecked()
}

/// Single face with a 2x2 NURBS (too small for swept detection).
fn make_tiny_nurbs_face() -> (BRep, FaceHandle) {
    let mut b = BRepBuilder::new(Tolerance::default());
    let v0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let v1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let v2 = b.add_vertex(Point3::new(1.0, 1.0, 0.0));
    let v3 = b.add_vertex(Point3::new(0.0, 1.0, 0.0));
    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let (_, h01, _) = b.add_edge(v0, v1);
    let (_, h12, _) = b.add_edge(v1, v2);
    let (_, h23, _) = b.add_edge(v2, v3);
    let (_, h30, _) = b.add_edge(v3, v0);
    let face = b.add_face(shell);
    b.add_loop(&[h01, h12, h23, h30], face).expect("loop");
    b.set_face_surface(
        face,
        SurfaceOnFace {
            geometry: SurfaceGeometry::Nurbs(NurbsSurface {
                degree_u: 1,
                degree_v: 1,
                knots_u: vec![0.0, 0.0, 1.0, 1.0],
                knots_v: vec![0.0, 0.0, 1.0, 1.0],
                control_points: vec![
                    Point3::new(0.0, 0.0, 0.0),
                    Point3::new(0.0, 1.0, 0.0),
                    Point3::new(1.0, 0.0, 0.0),
                    Point3::new(1.0, 1.0, 0.0),
                ],
                num_u: 2,
                num_v: 2,
                weights: vec![1.0; 4],
            }),
            orientation: true,
        },
    );
    (b.build_unchecked(), face)
}

/// Flat degree-3 NURBS patch (4x4) at z=0.
fn flat_nurbs(x_off: f64, y_off: f64, w: f64, h: f64) -> NurbsSurface {
    flat_nurbs_z(x_off, y_off, 0.0, w, h)
}

/// Flat degree-3 NURBS patch (4x4) at given z height.
fn flat_nurbs_z(x_off: f64, y_off: f64, z: f64, w: f64, h: f64) -> NurbsSurface {
    let (nu, nv) = (4, 4);
    let mut cps = Vec::with_capacity(nu * nv);
    let mut wts = Vec::with_capacity(nu * nv);
    for i in 0..nu {
        for j in 0..nv {
            cps.push(Point3::new(
                x_off + w * i as f64 / 3.0,
                y_off + h * j as f64 / 3.0,
                z,
            ));
            wts.push(1.0);
        }
    }
    NurbsSurface {
        degree_u: 3,
        degree_v: 3,
        knots_u: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: cps,
        num_u: nu,
        num_v: nv,
        weights: wts,
    }
}

// -- Fillet tests ------------------------------------------------------------

#[test]
fn detect_fillet_on_coplanar_middle_face() {
    let brep = make_three_coplanar_faces();
    let tol = Tolerance::default();
    let mut found = false;
    for (fh, _) in &brep.faces {
        let r = detect_surface_type(&brep, fh, &tol);
        if r.surface_type == SurfaceType::Fillet {
            found = true;
            assert!(r.confidence > 0.0);
            assert!(r.description.contains("G1-continuous"));
        }
    }
    assert!(found, "expected at least one face detected as fillet");
}

#[test]
fn repair_fillet_exercises_enforcement_body() {
    let mut brep = make_three_coplanar_faces();
    let tol = Tolerance::default();
    let fh = brep
        .faces
        .keys()
        .find(|fh| detect_surface_type(&brep, *fh, &tol).surface_type == SurfaceType::Fillet)
        .expect("no fillet face detected");
    let r = repair_fillet(&mut brep, fh, &tol, 1.0);
    assert_eq!(r.surface_type, SurfaceType::Fillet);
    assert!(
        r.description.contains("fillet repair"),
        "desc: {}",
        r.description
    );
    assert!(r.max_displacement.is_finite());
}

#[test]
fn repair_all_finds_fillet_faces() {
    let mut brep = make_three_coplanar_faces();
    let tol = Tolerance::default();
    let results = repair_all_type_specific(&mut brep, &tol, 1.0);
    assert!(
        results
            .iter()
            .any(|r| r.surface_type == SurfaceType::Fillet),
        "repair_all should find fillet faces"
    );
}

// -- Offset tests ------------------------------------------------------------

#[test]
fn detect_offset_on_parallel_faces() {
    let brep = make_offset_pair(0.5);
    let tol = Tolerance::default();
    let mut found = false;
    for (fh, _) in &brep.faces {
        let r = detect_surface_type(&brep, fh, &tol);
        if r.surface_type == SurfaceType::Offset {
            found = true;
            assert!(r.confidence > 0.0);
            assert!(r.description.contains("constant offset"));
        }
    }
    assert!(found, "expected at least one face detected as offset");
}

#[test]
fn repair_offset_exercises_adjustment_body() {
    let mut brep = make_offset_pair(0.5);
    let tol = Tolerance::default();
    let fh = brep
        .faces
        .keys()
        .find(|fh| detect_surface_type(&brep, *fh, &tol).surface_type == SurfaceType::Offset)
        .expect("no offset face detected");
    let r = repair_offset(&mut brep, fh, &tol, 1.0);
    assert_eq!(r.surface_type, SurfaceType::Offset);
    assert!(
        r.description.contains("offset repair"),
        "desc: {}",
        r.description
    );
    assert!(r.max_displacement.is_finite());
}

#[test]
fn repair_all_finds_offset_faces() {
    let mut brep = make_offset_pair(0.5);
    let tol = Tolerance::default();
    let results = repair_all_type_specific(&mut brep, &tol, 1.0);
    assert!(
        results
            .iter()
            .any(|r| r.surface_type == SurfaceType::Offset),
        "repair_all should find offset faces"
    );
}

// -- Swept tests -------------------------------------------------------------

#[test]
fn detect_swept_on_uniform_cross_section_face() {
    let brep = make_swept_face();
    let tol = Tolerance::default();
    let face = brep.faces.keys().next().expect("test has one face");
    let r = detect_surface_type(&brep, face, &tol);
    assert_eq!(
        r.surface_type,
        SurfaceType::Swept,
        "expected Swept, got {:?}: {}",
        r.surface_type,
        r.description
    );
    assert!(r.confidence > 0.0);
}

#[test]
fn repair_swept_exercises_cross_section_body() {
    let mut brep = make_swept_face();
    let tol = Tolerance::default();
    let face = brep.faces.keys().next().expect("test has one face");
    let r = repair_swept(&mut brep, face, &tol, 1.0);
    assert_eq!(r.surface_type, SurfaceType::Swept);
    assert!(
        r.description.contains("swept repair"),
        "desc: {}",
        r.description
    );
    assert!(r.max_displacement.is_finite());
}

#[test]
fn repair_all_finds_swept_faces() {
    let mut brep = make_swept_face();
    let tol = Tolerance::default();
    let results = repair_all_type_specific(&mut brep, &tol, 1.0);
    assert!(
        results.iter().any(|r| r.surface_type == SurfaceType::Swept),
        "repair_all should find swept faces"
    );
}

// -- Edge cases --------------------------------------------------------------

#[test]
fn detect_all_on_three_face_strip() {
    let brep = make_three_coplanar_faces();
    let tol = Tolerance::default();
    let results = detect_all_surface_types(&brep, &tol);
    assert_eq!(results.len(), 3);
    let n = results
        .iter()
        .filter(|r| r.surface_type == SurfaceType::Fillet)
        .count();
    assert!(n >= 1, "expected >= 1 fillet, got {n}");
}

#[test]
fn detect_all_on_offset_pair() {
    let brep = make_offset_pair(0.5);
    let tol = Tolerance::default();
    let results = detect_all_surface_types(&brep, &tol);
    assert_eq!(results.len(), 2);
    let n = results
        .iter()
        .filter(|r| r.surface_type == SurfaceType::Offset)
        .count();
    assert!(n >= 1, "expected >= 1 offset, got {n}");
}

#[test]
fn swept_detection_needs_min_grid_size() {
    let (brep, face) = make_tiny_nurbs_face();
    let r = detect_surface_type(&brep, face, &Tolerance::default());
    assert_ne!(
        r.surface_type,
        SurfaceType::Swept,
        "2x2 should not be swept"
    );
}

#[test]
fn offset_not_detected_at_tiny_distance() {
    let brep = make_offset_pair(1e-7); // below linear()*5 threshold
    let tol = Tolerance::default();
    for (fh, _) in &brep.faces {
        let r = detect_surface_type(&brep, fh, &tol);
        assert_ne!(
            r.surface_type,
            SurfaceType::Offset,
            "tiny offset: no detect"
        );
    }
}
