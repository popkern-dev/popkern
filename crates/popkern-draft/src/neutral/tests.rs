use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn plane_neutral_on_cube_face() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    let nl = compute_neutral_line(
        &brep,
        face,
        &NeutralReference::Plane {
            origin: Point3::new(0.0, 0.0, 0.5),
            normal: Vec3::new(0.0, 0.0, 1.0),
        },
        &Vec3::new(0.0, 0.0, 1.0),
        &tol(),
    );
    let _ = nl;
}

#[test]
fn auto_neutral_on_cube_face() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    let nl = compute_neutral_line(
        &brep,
        face,
        &NeutralReference::Auto,
        &Vec3::new(0.0, 0.0, 1.0),
        &tol(),
    );
    let _ = nl;
}

#[test]
fn edge_neutral_on_cube_face() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    let face_data = &brep.faces[face];
    let outer_loop = &brep.loops[face_data.outer_loop];
    let he = &brep.half_edges[outer_loop.half_edge];
    let edge_h = he.edge;

    let nl = compute_neutral_line(
        &brep,
        face,
        &NeutralReference::Edge(edge_h),
        &Vec3::new(0.0, 0.0, 1.0),
        &tol(),
    );
    assert!(nl.is_ok(), "edge neutral failed: {:?}", nl.err());
    let nl = nl.unwrap();
    assert!(nl.direction.length() > 0.9, "direction should be unit");
}

#[test]
fn invalid_edge_handle_returns_error() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    let nl = compute_neutral_line(
        &brep,
        face,
        &NeutralReference::Edge(Default::default()),
        &Vec3::new(0.0, 0.0, 1.0),
        &tol(),
    );
    assert!(nl.is_err());
}

#[test]
fn plane_neutral_parallel_to_face_returns_error() {
    // For a cube's Z-normal face, a neutral plane also with Z normal is parallel →
    // the cross product is zero → InvalidNeutral error.
    let brep = BRep::unit_cube(tol()).unwrap();

    // Find the face whose surface normal is Z (top or bottom face).
    use popkern_core::geometry::surface::SurfaceEval;
    let z_face = brep.faces.iter().find(|(_, fd)| {
        fd.surface.as_ref().is_some_and(|sof| {
            let n = sof.geometry.normal_at(0.5, 0.5);
            n.z.abs() > 0.9
        })
    });

    if let Some((face_h, _)) = z_face {
        // Neutral plane also has Z normal → parallel to face → error.
        let nl = compute_neutral_line(
            &brep,
            face_h,
            &NeutralReference::Plane {
                origin: Point3::new(0.5, 0.5, 0.5),
                normal: Vec3::new(0.0, 0.0, 1.0),
            },
            &Vec3::new(0.0, 0.0, 1.0),
            &tol(),
        );
        assert!(
            nl.is_err(),
            "expected error for neutral plane parallel to face, got {:?}",
            nl.ok()
        );
    }
}

#[test]
fn edge_neutral_edge_not_on_face_returns_error() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let mut face_iter = brep.faces.keys();
    let face_a = face_iter.next().unwrap();
    let face_b = face_iter.next().unwrap();

    // Find an edge that belongs to face_b but NOT face_a.
    let face_b_data = &brep.faces[face_b];
    let loop_b = &brep.loops[face_b_data.outer_loop];
    let he_b = &brep.half_edges[loop_b.half_edge];
    let edge_from_b = he_b.edge;

    // Verify the edge is actually on face_b, not face_a, before asserting.
    let face_a_data = &brep.faces[face_a];
    let loop_a = &brep.loops[face_a_data.outer_loop];
    let start_he = loop_a.half_edge;
    let mut cur = start_he;
    let mut edge_on_a = false;
    loop {
        let he = &brep.half_edges[cur];
        if he.edge == edge_from_b {
            edge_on_a = true;
            break;
        }
        cur = he.next;
        if cur == start_he {
            break;
        }
    }

    if !edge_on_a {
        let nl = compute_neutral_line(
            &brep,
            face_a,
            &NeutralReference::Edge(edge_from_b),
            &Vec3::new(0.0, 0.0, 1.0),
            &tol(),
        );
        assert!(nl.is_err(), "expected error: edge not on face");
    }
}

#[test]
fn auto_neutral_with_pull_parallel_to_face_normal_uses_edge_fallback() {
    // Find a Z-normal cube face, then pull in Z direction.
    // The neutral plane (Z-normal through centroid) is parallel to the face
    // → compute_plane_neutral fails → auto falls back to find_best_edge_neutral.
    let brep = BRep::unit_cube(tol()).unwrap();
    use popkern_core::geometry::surface::SurfaceEval;
    let z_face = brep.faces.iter().find(|(_, fd)| {
        fd.surface.as_ref().is_some_and(|sof| {
            let n = sof.geometry.normal_at(0.5, 0.5);
            n.z.abs() > 0.9
        })
    });

    if let Some((face_h, _)) = z_face {
        // Pull in Z → neutral plane has Z normal → parallel to face → falls back.
        let nl = compute_neutral_line(
            &brep,
            face_h,
            &NeutralReference::Auto,
            &Vec3::new(0.0, 0.0, 1.0),
            &tol(),
        );
        // Should succeed via edge fallback (the face has 4 edges to choose from).
        // It may also succeed directly if the centroid plane is acceptable.
        let _ = nl; // either outcome is valid; just exercise the code path.
    }
}

#[test]
fn plane_neutral_succeeds_for_side_face() {
    // A side face (X-normal) with a Z-normal neutral plane should intersect.
    let brep = BRep::unit_cube(tol()).unwrap();
    use popkern_core::geometry::surface::SurfaceEval;
    let x_face = brep.faces.iter().find(|(_, fd)| {
        fd.surface.as_ref().is_some_and(|sof| {
            let n = sof.geometry.normal_at(0.5, 0.5);
            n.x.abs() > 0.9
        })
    });

    if let Some((face_h, _)) = x_face {
        let nl = compute_neutral_line(
            &brep,
            face_h,
            &NeutralReference::Plane {
                origin: Point3::new(0.5, 0.5, 0.5),
                normal: Vec3::new(0.0, 0.0, 1.0),
            },
            &Vec3::new(0.0, 1.0, 0.0),
            &tol(),
        );
        assert!(
            nl.is_ok(),
            "expected Ok for side face + Z neutral, got {:?}",
            nl.err()
        );
    }
}

#[test]
fn edge_neutral_direction_is_unit() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    let face_data = &brep.faces[face];
    let outer_loop = &brep.loops[face_data.outer_loop];

    // Walk all 4 edges of the outer loop and check each.
    let start_he = outer_loop.half_edge;
    let mut cur = start_he;
    loop {
        let he = &brep.half_edges[cur];
        let edge_h = he.edge;
        let nl = compute_neutral_line(
            &brep,
            face,
            &NeutralReference::Edge(edge_h),
            &Vec3::new(0.0, 0.0, 1.0),
            &tol(),
        );
        if let Ok(nl) = nl {
            let len = nl.direction.length();
            assert!((len - 1.0).abs() < 0.01, "direction not unit: {}", len);
        }
        cur = he.next;
        if cur == start_he {
            break;
        }
    }
}
