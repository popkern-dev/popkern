#![cfg(test)]
use super::var_radius_tests_support::*;
use super::*;

#[test]
#[ignore]
fn dump_neighcyl_connectivity_cube_mixed() {
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
    let brep = filleted.brep;
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
    dump_neighcyl_connectivity(&brep, &cap_faces, &per_fillet_caps);
}

/// Phase 6c-D investigation: walk each cap's outer loop in both
/// rung-on and rung-off modes side-by-side.  Reveals exactly
/// where the topology branches between the two paths.  Skips
/// reroute (we want the *pre-surgery* topology each path
/// produces).
#[test]
#[ignore]
fn dump_cap_loops_rung_compare() {
    use popkern_core::geometry::surface::SurfaceGeometry;

    for (label, rung) in [("RUNG_OFF", false), ("RUNG_ON", true)] {
        unsafe {
            std::env::set_var("POPKERN_SPHERE_FIT_BOUNDARY", "1");
            if rung {
                std::env::set_var("POPKERN_SPHERE_CAP_RUNG", "1");
            } else {
                std::env::remove_var("POPKERN_SPHERE_CAP_RUNG");
            }
            std::env::remove_var("POPKERN_VAR_RADIUS_REROUTE");
            // Clean other env vars that could affect topology.
            std::env::remove_var("POPKERN_APEX_TOPOLOGY_REBUILD");
            std::env::remove_var("POPKERN_OUTER_FILLET_G1");
            std::env::remove_var("POPKERN_EDGE_TRIPLE");
        }
        eprintln!("\n========== {label} ==========");

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
        let Ok(filleted) = crate::fillet_edges(
            &cube,
            &specs,
            &crate::FilletOptions::default(),
            Tolerance::default(),
            None,
        ) else {
            eprintln!("  fillet_edges failed");
            continue;
        };
        let brep = filleted.brep;
        eprintln!(
            "BRep totals: {} faces, {} edges, {} half_edges, {} vertices",
            brep.faces.len(),
            brep.edges.len(),
            brep.half_edges.len(),
            brep.vertices.len()
        );
        eprintln!("corner_patches: {}", filleted.corner_patches.len());
        for cp in &filleted.corner_patches {
            eprintln!(
                "  {:?} method={:?} face={:?}",
                cp.classification, cp.method, cp.face
            );
        }
        // Also list ALL faces in the BRep with their surface type.
        eprintln!("All faces:");
        for (fh, f) in brep.faces.iter() {
            let kind = match f.surface.as_ref().map(|s| &s.geometry) {
                Some(SurfaceGeometry::Plane { .. }) => "Plane",
                Some(SurfaceGeometry::Cylinder { .. }) => "Cyl",
                Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
                Some(SurfaceGeometry::Gregory(_)) => "Gregory",
                Some(SurfaceGeometry::PolyhedralSpline(_)) => "PolyhedralSpline",
                Some(SurfaceGeometry::Torus { .. }) => "Torus",
                Some(SurfaceGeometry::Cone { .. }) => "Cone",
                Some(SurfaceGeometry::Nurbs(_)) => "Nurbs",
                None => "?",
            };
            eprintln!("  {fh:?} {kind}");
        }

        // Identify cap faces by their cap_sphere data on
        // CornerPatchInfo (works regardless of surface type).
        let mut cap_info: Vec<(FaceHandle, Option<(Point3, f64)>)> = filleted
            .corner_patches
            .iter()
            .map(|p| (p.face, p.cap_sphere))
            .collect();
        cap_info.sort_by_key(|(fh, _)| {
            // Stable order: by sphere centre x then radius
            let (cx, r) = (*fh, 0.0);
            let _ = (cx, r);
            *fh
        });
        for (k, (fh, sphere)) in cap_info.iter().enumerate() {
            let surface_kind = match brep
                .faces
                .get(*fh)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry)
            {
                Some(SurfaceGeometry::Sphere { center, radius }) => format!(
                    "Sphere centre=({:.3},{:.3},{:.3}) r={:.3}",
                    center.x, center.y, center.z, radius
                ),
                Some(SurfaceGeometry::Gregory(_)) => "Gregory".to_string(),
                Some(SurfaceGeometry::PolyhedralSpline(_)) => "PolyhedralSpline".to_string(),
                Some(other) => format!("{other:?}"),
                None => "<no surface>".to_string(),
            };
            eprintln!(
                "\ncap[{k}] {fh:?} cap_sphere={:?} surface={surface_kind}",
                sphere
                    .map(|(c, r)| format!("centre=({:.3},{:.3},{:.3}) r={:.3}", c.x, c.y, c.z, r))
            );
            let Some(face_data) = brep.faces.get(*fh) else {
                continue;
            };
            let outer = face_data.outer_loop;
            let Some(loop_data) = brep.loops.get(outer) else {
                continue;
            };
            let start = loop_data.half_edge;
            let mut cursor = start;
            let mut step = 0;
            loop {
                if step > 16 {
                    break;
                }
                let Some(he) = brep.half_edges.get(cursor) else {
                    break;
                };
                let edge = he.edge;
                let v_origin = he.origin;
                let pos = brep.vertices[v_origin].position;
                let twin = he.twin;
                let twin_face = brep
                    .half_edges
                    .get(twin)
                    .map(|t| brep.loops[t.loop_].face)
                    .unwrap_or(*fh);
                let twin_kind = match brep
                    .faces
                    .get(twin_face)
                    .and_then(|f| f.surface.as_ref())
                    .map(|s| &s.geometry)
                {
                    Some(SurfaceGeometry::Plane { .. }) => "Plane",
                    Some(SurfaceGeometry::Cylinder { .. }) => "Cyl",
                    Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
                    Some(SurfaceGeometry::Gregory(_)) => "Greg",
                    Some(SurfaceGeometry::PolyhedralSpline(_)) => "Greg",
                    _ => "?",
                };
                eprintln!(
                    "  step {step}: he={cursor:?} edge={edge:?} V={v_origin:?}=({:.3},{:.3},{:.3}) → twin face {twin_face:?} ({twin_kind})",
                    pos.x, pos.y, pos.z
                );
                cursor = he.next;
                step += 1;
                if cursor == start {
                    break;
                }
            }
        }
    }
}
