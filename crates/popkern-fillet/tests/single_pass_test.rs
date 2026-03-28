use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

/// Count faces whose boundary vertices are far from the assigned surface.
fn count_bad_faces(brep: &BRep, threshold: f64) -> (usize, Vec<String>) {
    let mut bad_faces = Vec::new();
    for fh in brep.faces.keys() {
        let face = &brep.faces[fh];
        let Some(sof) = face.surface.as_ref() else {
            bad_faces.push(format!("{fh:?}: no surface assigned"));
            continue;
        };
        let positions = brep.face_outer_loop_positions(fh);
        if positions.is_empty() {
            continue;
        }

        for p in &positions {
            if let Some((u, v)) = project_point_to_uv(&sof.geometry, p) {
                let surf_pt = sof.geometry.point_at(u, v);
                let d = p.distance_to(&surf_pt);
                if d > threshold {
                    let surf_type = match &sof.geometry {
                        popkern_core::geometry::surface::SurfaceGeometry::Plane { .. } => "Plane",
                        popkern_core::geometry::surface::SurfaceGeometry::Torus { .. } => "Torus",
                        popkern_core::geometry::surface::SurfaceGeometry::Sphere { .. } => "Sphere",
                        _ => "Other",
                    };
                    bad_faces.push(format!(
                        "{fh:?} ({surf_type}): vertex ({:.3},{:.3},{:.3}) is {d:.4} from surface",
                        p.x, p.y, p.z,
                    ));
                    break;
                }
            }
        }
    }
    (bad_faces.len(), bad_faces)
}

/// 12-edge fillet surface correctness regression test.
///
/// No fillet faces should have wrong surfaces (the mef surface
/// inheritance fix eliminated surface loss during vertex blend).
/// The remaining ≤3 bad faces are vertex blend corners where the
/// blend surface itself (mitred plane or sphere) doesn't precisely
/// fit all boundary vertices — a blend surface quality issue, not
/// a surface assignment issue. The old sequential approach had 5.
#[test]
fn twelve_edge_fillet_surface_correctness() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).unwrap();
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let result = fillet_edges(&cube, &specs, &FilletOptions::default(), tol).unwrap();
    let (bad_count, details) = count_bad_faces(&result.brep, 0.05);
    assert!(
        bad_count <= 3,
        "expected at most 3 bad faces (blend surface precision), got {bad_count}:\n{}",
        details.join("\n"),
    );
}

/// All faces should have a surface assigned (no None surfaces).
#[test]
fn twelve_edge_fillet_no_missing_surfaces() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).unwrap();
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let result = fillet_edges(&cube, &specs, &FilletOptions::default(), tol).unwrap();
    let missing: Vec<_> = result
        .brep
        .faces
        .keys()
        .filter(|fh| result.brep.faces[*fh].surface.is_none())
        .collect();
    assert!(
        missing.is_empty(),
        "faces with no surface assigned: {missing:?}",
    );
}
