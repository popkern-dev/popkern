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
                    bad_faces.push(format!(
                        "{fh:?}: vertex ({:.3},{:.3},{:.3}) is {d:.4} from surface",
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
/// The three-phase approach (split all, then merge all) should produce
/// at most 2 misassigned faces. The old sequential approach produced 5.
/// The remaining ≤2 are vertex blend corner cases where the blend
/// surface restoration assigns a plane instead of the correct surface.
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
        bad_count <= 2,
        "expected at most 2 bad faces (vertex blend corners), got {bad_count}:\n{}",
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
