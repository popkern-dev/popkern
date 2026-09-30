//! PK-598: the fitted cap's scoped self-intersection postcondition
//! still declines a cap that interpenetrates its region.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::fit_cap_region_self_intersection;
use crate::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};

/// The octagonal pyramid apex (h 0.8, mixed radii) under a G2 request,
/// served by the tangent fit (the curvature-weighted fit declines) with the KP20 cap forced off: the served body and its
/// cap sectors.
fn fitted_octagon_apex() -> (BRep, Vec<FaceHandle>) {
    let tol = Tolerance::default();
    let body = BRep::n_gonal_pyramid(8, 1.0, 0.8, tol).expect("pyramid");
    let apex = body
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(v, _)| v)
        .expect("apex");
    let mut edges: Vec<_> = body
        .half_edges
        .values()
        .filter(|he| he.origin == apex)
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip([0.04, 0.06, 0.05, 0.07].iter().cycle())
        .map(|(&edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..Default::default()
    };
    let result = crate::dispatch_overrides::with_kp20_cap_forced_to_decline(|| {
        fillet_edges(&body, &specs, &options, tol, None)
    })
    .expect("fillet");
    let caps: Vec<FaceHandle> = result
        .corner_patches
        .iter()
        .filter(|p| p.setback_cap)
        .map(|p| p.face)
        .collect();
    assert_eq!(caps.len(), 8, "the fitted cap serves 8 sectors");
    (result.brep, caps)
}

#[test]
fn a_served_fit_cap_passes_and_a_displaced_sector_declines() {
    let (mut brep, caps) = fitted_octagon_apex();
    assert!(
        fit_cap_region_self_intersection(&brep, &caps, &[]).is_ok(),
        "the served cap is self-intersection clean"
    );
    // Negative control: carry sector 0 across the cap, twice its offset
    // from the pyramid axis to the far side. Its rim stays topologically
    // where it was, but its surface now runs through the sectors and
    // strips opposite (not its edge neighbours, whose contact every
    // scan exempts: sinking it along the axis crossed only those).
    let Some(SurfaceGeometry::Nurbs(nurbs)) = brep.faces[caps[0]]
        .surface
        .as_mut()
        .map(|s| &mut s.geometry)
    else {
        panic!("fitted sectors are NURBS");
    };
    let count = nurbs.control_points.len() as f64;
    let (cx, cy) = nurbs
        .control_points
        .iter()
        .fold((0.0, 0.0), |(x, y), p| (x + p.x / count, y + p.y / count));
    for p in &mut nurbs.control_points {
        *p = *p - Vec3::new(2.0 * cx, 2.0 * cy, 0.0);
    }
    let verdict = fit_cap_region_self_intersection(&brep, &caps, &[]);
    assert!(
        verdict.is_err(),
        "a cap interpenetrating its region must decline"
    );
}
