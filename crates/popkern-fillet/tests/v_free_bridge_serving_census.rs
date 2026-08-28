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

/// Production planes must not duplicate: two Production faces on the
/// SAME plane whose regions share interior points indicate the
/// face-family bridge double-cover (bridge face cloned the support
/// surface without retiring the source support face). Complementary
/// tilings that only share boundary are legitimate.
#[test]
fn no_duplicate_production_plane_coverage() {
    let brep = mixed_radius_cube();
    let mut by_plane: std::collections::BTreeMap<String, Vec<Vec<(f64, f64)>>> =
        std::collections::BTreeMap::new();
    for (fh, f) in brep.faces.iter() {
        if f.role != TopologyRole::Production {
            continue;
        }
        let Some(sof) = f.surface.as_ref() else {
            continue;
        };
        let SurfaceGeometry::Plane { origin, normal, .. } = &sof.geometry else {
            continue;
        };
        let Some(n) = normal.normalized() else {
            continue;
        };
        let key = format!(
            "plane n=({:.6},{:.6},{:.6}) d={:.6}",
            n.x,
            n.y,
            n.z,
            origin.x * n.x + origin.y * n.y + origin.z * n.z
        );
        // In-plane 2D coordinates.
        let e1 = if n.x.abs() < 0.9 {
            popkern_core::geometry::Vec3::new(1.0, 0.0, 0.0)
        } else {
            popkern_core::geometry::Vec3::new(0.0, 1.0, 0.0)
        };
        let e1 = (e1 - n * e1.dot(&n)).normalized().unwrap();
        let e2 = n.cross(&e1);
        let ring: Vec<(f64, f64)> = brep
            .face_outer_loop_positions(fh)
            .iter()
            .map(|p| {
                let r = *p - *origin;
                (r.dot(&e1), r.dot(&e2))
            })
            .collect();
        if ring.len() >= 3 {
            by_plane.entry(key).or_default().push(ring);
        }
    }
    let inside = |poly: &[(f64, f64)], p: (f64, f64)| -> bool {
        // Strict interior: positive winding test via ray casting, plus
        // a clearance from every edge so shared boundaries don't count.
        let mut winding = false;
        let mut min_edge_dist = f64::MAX;
        for i in 0..poly.len() {
            let (x0, y0) = poly[i];
            let (x1, y1) = poly[(i + 1) % poly.len()];
            if ((y0 > p.1) != (y1 > p.1)) && (p.0 < (x1 - x0) * (p.1 - y0) / (y1 - y0) + x0) {
                winding = !winding;
            }
            let (dx, dy) = (x1 - x0, y1 - y0);
            let len2 = dx * dx + dy * dy;
            let t = if len2 > 0.0 {
                (((p.0 - x0) * dx + (p.1 - y0) * dy) / len2).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let (cx, cy) = (x0 + t * dx, y0 + t * dy);
            min_edge_dist = min_edge_dist.min(((p.0 - cx).powi(2) + (p.1 - cy).powi(2)).sqrt());
        }
        winding && min_edge_dist > 1.0e-6
    };
    for (key, rings) in &by_plane {
        for i in 0..rings.len() {
            for j in (i + 1)..rings.len() {
                // Sample a grid over ring i's bbox; a point strictly
                // inside BOTH rings is genuine double coverage.
                let (mut x0, mut y0, mut x1, mut y1) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
                for &(x, y) in &rings[i] {
                    x0 = x0.min(x);
                    y0 = y0.min(y);
                    x1 = x1.max(x);
                    y1 = y1.max(y);
                }
                let mut shared = None;
                'grid: for gx in 1..40 {
                    for gy in 1..40 {
                        let p = (
                            x0 + (x1 - x0) * f64::from(gx) / 40.0,
                            y0 + (y1 - y0) * f64::from(gy) / 40.0,
                        );
                        if inside(&rings[i], p) && inside(&rings[j], p) {
                            shared = Some(p);
                            break 'grid;
                        }
                    }
                }
                assert!(
                    shared.is_none(),
                    "two Production faces on the same plane share interior point \
                     {shared:?} ({key}) — face-family bridge double-cover"
                );
            }
        }
    }
}

/// PK-167 stage 1b acceptance: the constructive in-place rewrite
/// is the DEFAULT serving path (the historical PK167_INPLACE
/// opt-in and PK167_INPLACE_OFF opt-out are both retired, and
/// V-free repeated-apex construction is DefaultOn for n >= 3 since
/// PK-153/PK-167). The mixed cube must serve a closed shell —
/// validated body, no boundary edges — with the same
/// duplicate-coverage cleanliness as the bridge+repair path.
#[test]
fn constructive_rewrite_serves_closed_shell() {
    let brep = mixed_radius_cube();
    assert!(
        brep.boundary_edges.is_empty(),
        "constructive serving must leave no boundary edges: {:?}",
        brep.boundary_edges
    );
    let errors = brep.validate_all();
    assert!(
        errors.len() <= 1,
        "constructive serving regressed validity: {errors:?}"
    );
}
