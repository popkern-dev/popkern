use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

/// Every face's boundary vertices must lie on its assigned surface.
///
/// This catches the misassignment bug: when face_remap errors cause a
/// face to end up with a surface from a different original face, the
/// vertices won't project correctly onto that surface.
#[test]
fn all_12_edge_fillet_faces_have_correct_surfaces() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).unwrap();
    let edges: Vec<_> = cube.edges.keys().collect();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let result = fillet_edges(&cube, &specs, &FilletOptions::default(), tol).unwrap();
    let brep = &result.brep;

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
                if d > 0.05 {
                    bad_faces.push(format!(
                        "{fh:?}: vertex ({:.3},{:.3},{:.3}) is {d:.4} from surface",
                        p.x, p.y, p.z,
                    ));
                    break; // One bad vertex per face is enough.
                }
            }
        }
    }

    assert!(
        bad_faces.is_empty(),
        "faces with incorrect surface assignments:\n{}",
        bad_faces.join("\n"),
    );
}
