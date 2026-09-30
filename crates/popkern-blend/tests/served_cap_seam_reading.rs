//! PK-607: the served whole-seam cap reading sees a crease the
//! mid-seam number cannot.
//!
//! The cube corner 0.05/0.2/0.3 under a G1 request, observed at
//! construction (serving retry off), commits a three-sector
//! direct-apex cap whose sector seam certifies a normal-sine lower
//! bound of 0.9998 at fraction 0.875-0.94 (an independent projection
//! oracle agrees). The mid-seam number (fractions 0.4-0.6) reads 0.245
//! there; the whole-seam `served_cap_seam_g1` must read the crease.

use popkern_blend::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
fn the_whole_seam_reading_sees_a_crease_the_mid_seam_misses() {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    let corner = cube
        .vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&Point3::new(1.0, 1.0, 1.0)) < 1e-9)
        .map(|(h, _)| h)
        .expect("corner");
    let mut edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == corner)
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip([0.05, 0.2, 0.3])
        .map(|(&edge, r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G1,
        serving_self_intersection_retry: false,
        ..Default::default()
    };
    let result = fillet_edges(&cube, &specs, &options, tol, None).expect("fillet");
    let details: Vec<String> = result
        .repair_report
        .apex_routes
        .iter()
        .flat_map(|r| r.exact_boundary_report.diagnostics.iter())
        .filter(|d| d.code == "exact_boundary_cap.mid_seam_normal")
        .filter_map(|d| d.detail_code.clone())
        .collect();
    let field = |detail: &str, key: &str| -> f64 {
        detail
            .split(';')
            .find_map(|part| part.strip_prefix(key)?.strip_prefix('=')?.parse().ok())
            .unwrap_or(f64::NAN)
    };
    let [detail] = details.as_slice() else {
        panic!("one committed direct-apex cap: {details:?}");
    };
    let mid = field(detail, "served_mid_seam_normal");
    let whole = field(detail, "served_cap_seam_g1");
    assert!(mid < 0.3, "the mid-seam number reads 0.245: {detail}");
    assert!(
        whole > 0.9,
        "the whole-seam reading sees the crease: {detail}"
    );
}
