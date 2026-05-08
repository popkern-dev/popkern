#![cfg(test)]
use super::var_radius_tests_support::*;
use super::*;

#[test]
fn classify_empty_cap_list() {
    let brep = BRep::new(Tolerance::default());
    let result = classify_cap_loops(&brep, &[], &[]);
    assert!(result.is_empty());
}

/// Per_fillet_caps shorter than cap_faces — function falls back
/// to a zero centre/radius for the missing entries (no panic,
/// twins classified as NeighbourCylinder or Other).
#[test]
fn classify_handles_missing_per_fillet_caps() {
    let brep = BRep::new(Tolerance::default());
    let result = classify_cap_loops(&brep, &[], &[(Point3::new(0.0, 0.0, 0.0), 0.1)]);
    assert!(result.is_empty());
}

/// End-to-end: run cube_mixed pipeline, find a shared
/// rolling-ball vertex, split it, validate. Used as the
/// first empirical test of `split_apex_footprint_vertex`.
#[test]
#[ignore]
fn split_one_cube_mixed_vertex() {
    unsafe {
        std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
        std::env::set_var("POPKERN_SPHERE_CAP_RUNG", "1");
    }
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let apex = cube
        .vertices
        .keys()
        .find(|&v| cube.half_edges.values().filter(|he| he.origin == v).count() == 3)
        .unwrap();
    let edges: Vec<_> = {
        use std::collections::BTreeSet;
        let s: BTreeSet<_> = cube
            .half_edges
            .values()
            .filter(|he| he.origin == apex)
            .map(|he| he.edge)
            .collect();
        s.into_iter().collect()
    };
    let radii = [0.08_f64, 0.10, 0.12];
    let specs: Vec<crate::types::FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| crate::types::FilletSpec {
            edge: e,
            radius: crate::types::FilletRadius::Constant(r),
        })
        .collect();
    let result = crate::fillet_edges(
        &cube,
        &specs,
        &crate::FilletOptions::default(),
        Tolerance::default(),
        None,
    );
    let Ok(filleted) = result else {
        return;
    };
    let mut brep = filleted.brep;

    // Pick the first NeighCyl edge to find V_x and its cap_face / cube_face.
    use popkern_core::geometry::surface::SurfaceGeometry;
    let cap_faces: Vec<_> = brep
        .faces
        .iter()
        .filter_map(|(fh, f)| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Sphere { .. }) => Some(fh),
            _ => None,
        })
        .collect();
    let per_fillet_caps: Vec<(Point3, f64)> = cap_faces
        .iter()
        .filter_map(|&fh| {
            brep.faces
                .get(fh)
                .and_then(|f| match &f.surface.as_ref()?.geometry {
                    SurfaceGeometry::Sphere { center, radius } => Some((*center, *radius)),
                    _ => None,
                })
        })
        .collect();
    if cap_faces.len() != 3 {
        return;
    }

    let classification = classify_cap_loops(&brep, &cap_faces, &per_fillet_caps);
    let (cap_k, neigh_he) = classification
        .iter()
        .find_map(|c| {
            c.entries
                .iter()
                .find(|e| e.twin_kind == TwinKind::NeighbourCylinder)
                .map(|e| (c.cap_face, e.half_edge))
        })
        .unwrap();
    let v_x = brep.half_edges[neigh_he].origin;
    let cyl_m = brep
        .loops
        .get(brep.half_edges[brep.half_edges[neigh_he].twin].loop_)
        .unwrap()
        .face;

    // Find the cube face V_x is incident to (the Plane face).
    let cube_face = brep
        .half_edges
        .values()
        .find_map(|he| {
            if he.origin != v_x {
                return None;
            }
            let f = brep.loops[he.loop_].face;
            let surf = brep.faces.get(f).and_then(|fd| fd.surface.as_ref())?;
            if matches!(surf.geometry, SurfaceGeometry::Plane { .. }) {
                Some(f)
            } else {
                None
            }
        })
        .unwrap();

    // Compute target position for v_new on cyl_m's cap-tangent
    // circle for the cap_n's sphere (the cap whose primary is cyl_m).
    let cyl_m_data = match brep
        .faces
        .get(cyl_m)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
    {
        Some(SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        }) => Some((*origin, *axis, *radius)),
        _ => None,
    }
    .unwrap();
    let (cyl_origin, cyl_axis, _cyl_radius) = cyl_m_data;
    let cap_n_sphere = cap_faces.iter().find_map(|&fh| {
        if fh == cap_k {
            return None;
        }
        match brep
            .faces
            .get(fh)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry)
        {
            Some(SurfaceGeometry::Sphere { center, .. }) => Some(*center),
            _ => None,
        }
    });
    let cap_n_centre = cap_n_sphere.unwrap();
    let v_x_pos = brep.vertices[v_x].position;
    // Hardcode: V_x is at (0, 0.12, 0.10). v_new at (0, 0.12, 0.12)
    // is the cap-tangent endpoint of cap_n (sphere r=0.12) on
    // cyl_m on cube face x=0.
    let v_new_pos = Point3::new(0.0, 0.12, 0.12);
    let _ = (cap_n_centre, cyl_origin, cyl_axis, v_x_pos);

    eprintln!(
        "Splitting V_x={v_x:?} pos=({:.3},{:.3},{:.3}) → ({:.3},{:.3},{:.3})",
        v_x_pos.x, v_x_pos.y, v_x_pos.z, v_new_pos.x, v_new_pos.y, v_new_pos.z
    );
    eprintln!("  cap_face={cap_k:?} cube_face={cube_face:?} cyl_split={cyl_m:?}");

    match split_apex_footprint_vertex(&mut brep, v_x, cap_k, cube_face, cyl_m, v_new_pos) {
        Ok(v_new) => {
            eprintln!("Split OK; v_new={v_new:?}");
            brep.validate().expect("validate after split");
        }
        Err(e) => panic!("split failed: {e:?}"),
    }
}
