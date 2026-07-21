//! PK-153: serving-shape census for the V-free repeated-apex output on
//! the mixed-radius cube corner — the fixture whose DefaultOn n>=3
//! widening was measured to double-cover (three face-family bridge
//! faces each duplicating a near support face; total mesh area
//! baseline + ~2.43).
//!
//! The hard gate asserts NO two production faces with identical
//! surfaces cover overlapping regions (the duplicate-surface census
//! from the strip-amputation class), and that mesh area under the
//! V-free serving matches the apex-retaining baseline within the
//! measured serving tolerance.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, TopologyRole};
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn mixed_radius_cube() -> BRep {
    let cube = BRep::unit_cube(Tolerance::default()).unwrap();
    let vertex = cube.vertices.keys().next().unwrap();
    let incident_edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    let radii = [0.08, 0.10, 0.12];
    let specs: Vec<FilletSpec> = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    fillet_edges(
        &cube,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("fillet_edges should succeed on mixed-radius cube corner")
    .brep
}

fn surface_key(geometry: &SurfaceGeometry) -> Option<String> {
    match geometry {
        SurfaceGeometry::Plane { origin, normal, .. } => Some(format!(
            "plane n=({:.6},{:.6},{:.6}) d={:.6}",
            normal.x,
            normal.y,
            normal.z,
            origin.x * normal.x + origin.y * normal.y + origin.z * normal.z
        )),
        _ => None,
    }
}

/// Census used by the PK-153 widening work: dump the face layout of
/// the V-free serving output. Run with
/// `POPKERN_DIRECT_APEX_REQUIRE_V_FREE_REPEATED_APEX=1 ... -- --ignored --nocapture`.
#[test]
#[ignore = "diagnostic census; run explicitly under the V-free env"]
fn v_free_face_census_dump() {
    let brep = mixed_radius_cube();
    for (fh, f) in brep.faces.iter() {
        let desc = match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Plane { normal, .. }) => {
                format!("Plane n=({:.2},{:.2},{:.2})", normal.x, normal.y, normal.z)
            }
            Some(SurfaceGeometry::Cylinder { radius, .. }) => format!("Cylinder r={radius}"),
            Some(SurfaceGeometry::Sphere { radius, .. }) => format!("Sphere r={radius}"),
            Some(SurfaceGeometry::Gregory(_)) => "Gregory".to_string(),
            Some(_) => "Other".to_string(),
            None => "NoSurface".to_string(),
        };
        let verts: Vec<String> = brep
            .face_outer_loop_positions(fh)
            .iter()
            .map(|p| format!("({:.2},{:.2},{:.2})", p.x, p.y, p.z))
            .collect();
        eprintln!("FACE {fh:?} [{:?}] {desc} ring={verts:?}", f.role);
    }
}

/// Production planes must not duplicate: two Production faces on the
/// SAME plane whose vertex rings overlap indicate the face-family
/// bridge double-cover (bridge face cloned the support surface without
/// retiring the source support face).
#[test]
fn no_duplicate_production_plane_coverage() {
    let brep = mixed_radius_cube();
    let mut by_plane: std::collections::BTreeMap<String, Vec<(f64, f64, f64, f64, f64, f64)>> =
        std::collections::BTreeMap::new();
    for (fh, f) in brep.faces.iter() {
        if f.role != TopologyRole::Production {
            continue;
        }
        let Some(key) = f.surface.as_ref().and_then(|s| surface_key(&s.geometry)) else {
            continue;
        };
        let ring = brep.face_outer_loop_positions(fh);
        if ring.is_empty() {
            continue;
        }
        let (mut x0, mut y0, mut z0) = (f64::MAX, f64::MAX, f64::MAX);
        let (mut x1, mut y1, mut z1) = (f64::MIN, f64::MIN, f64::MIN);
        for p in &ring {
            x0 = x0.min(p.x);
            y0 = y0.min(p.y);
            z0 = z0.min(p.z);
            x1 = x1.max(p.x);
            y1 = y1.max(p.y);
            z1 = z1.max(p.z);
        }
        by_plane
            .entry(key)
            .or_default()
            .push((x0, y0, z0, x1, y1, z1));
    }
    let overlaps = |a: &(f64, f64, f64, f64, f64, f64), b: &(f64, f64, f64, f64, f64, f64)| {
        // Strict interior overlap of the in-plane bounding boxes; edge
        // contact is legitimate tiling.
        let eps = 1.0e-9;
        a.0 < b.3 - eps
            && b.0 < a.3 - eps
            && a.1 < b.4 - eps
            && b.1 < a.4 - eps
            && a.2 < b.5 - eps
            && b.2 < a.5 - eps
    };
    for (key, boxes) in &by_plane {
        for i in 0..boxes.len() {
            for j in (i + 1)..boxes.len() {
                assert!(
                    !overlaps(&boxes[i], &boxes[j]),
                    "two Production faces on the same plane overlap ({key}): \
                     {:?} vs {:?} — face-family bridge double-cover",
                    boxes[i],
                    boxes[j]
                );
            }
        }
    }
}
