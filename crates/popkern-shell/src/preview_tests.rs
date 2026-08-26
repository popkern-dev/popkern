use super::*;

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn cube_outward_preview_carries_faces_and_strips() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let mesh = shell_preview(&brep, 0.15, ShellDirection::Outward, tol()).unwrap();
    assert!(!mesh.triangles.is_empty());
    let face_tris = mesh.face_map.iter().filter(|f| f.is_some()).count();
    let strip_tris = mesh.face_map.iter().filter(|f| f.is_none()).count();
    assert!(face_tris > 0, "face grids present");
    assert!(strip_tris > 0, "gap strips approximated on a convex grow");
    // Face samples realize the offset distance exactly (the preview
    // samples the exact offset map, never a fit).
    for (t, fh) in mesh.triangles.iter().zip(&mesh.face_map) {
        if fh.is_none() {
            continue;
        }
        for &vi in t {
            let p = mesh.positions[vi];
            // Unit cube: outward face samples sit 0.15 outside the box.
            let dx = (-p.x).max(p.x - 1.0).max(0.0);
            let dy = (-p.y).max(p.y - 1.0).max(0.0);
            let dz = (-p.z).max(p.z - 1.0).max(0.0);
            let d = (dx * dx + dy * dy + dz * dz).sqrt();
            assert!(
                (d - 0.15).abs() < 1e-9,
                "face sample must realize the offset: {d}"
            );
        }
    }
}

#[test]
fn cube_inward_preview_has_no_strips() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let mesh = shell_preview(&brep, 0.1, ShellDirection::Inward, tol()).unwrap();
    assert!(mesh.face_map.iter().any(|f| f.is_some()));
    assert!(
        mesh.face_map.iter().all(|f| f.is_some()),
        "inward cube offsets intersect: no gap strips"
    );
}

#[test]
fn midplane_preview_combines_both_skins() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let inward = shell_preview(&brep, 0.2, ShellDirection::Inward, tol()).unwrap();
    let mid = shell_preview(&brep, 0.2, ShellDirection::Midplane, tol()).unwrap();
    assert!(mid.triangles.len() > inward.triangles.len());
}
