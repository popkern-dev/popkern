//! PK-141 hole-corner cap geometry: the classic rolling-ball corner
//! solved on curved supports.
//!
//! At a true corner junction of two fillets sharing one support (the
//! boss-on-wall class: both roll on the wall, one against the top
//! plane and one against the side plane), the fillets stop at
//! chain-corner level and the territory between their end stations
//! ships as a mitre gap. The honest closure is:
//!
//! 1. the CORNER BALL: center `c` equidistant `r` from all three
//!    support surfaces (solved by alternating projection onto the
//!    three offset surfaces — exact for the analytic supports);
//! 2. two STRIP EXTENSIONS: each fillet's envelope continued from
//!    its current end station to the `c` station (centers on the
//!    shared-support offset intersected with the fillet's own plane
//!    offset — sampled and projected, then swept with the standard
//!    fillet-surface machinery);
//! 3. the SPHERE WEDGE at `c` between the two extended end arcs.
//!
//! This module is the GEOMETRY solve only, gated by its own G0
//! acceptance (extension boundaries against the strips and supports)
//! — the topology carve consumes it separately.

// The carve consumer lands in the follow-up stage; until it does,
// the solved-geometry API is exercised by its own gate only.
#![allow(dead_code)]

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, SurfaceEval};

/// Solved corner-cap geometry for one hole corner.
pub(crate) struct CornerCapGeometry {
    /// The corner-ball center.
    pub center: Point3,
    /// Contact of the corner ball on the shared support.
    pub shared_contact: Point3,
    /// Contact on each fillet's own (non-shared) support.
    pub own_contacts: [Point3; 2],
    /// Sampled center paths for the two strip extensions, each from
    /// the fillet's end station to the corner-ball center.
    pub extension_centers: [Vec<Point3>; 2],
    /// The shared support (the surface both fillets roll on — the
    /// wall of the corner). The carve projects wall trims onto it.
    pub shared_surface: SurfaceGeometry,
}

/// Distance from `p` to the offset surface of `s` at distance `r`
/// (unsigned), together with the projection of `p` onto that offset.
fn project_to_offset(s: &SurfaceGeometry, p: &Point3, r: f64) -> Option<(f64, Point3)> {
    let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(s, p)?;
    let q = s.point_at(u, v);
    let d = p.distance_to(&q);
    if d < 1e-12 {
        return None;
    }
    let dir = (*p - q) * (1.0 / d);
    let off = q + dir * r;
    Some(((d - r).abs(), off))
}

/// Solve the corner-ball center: the point at distance `r` from all
/// three supports, seeded near the junction. Alternating projection
/// onto the three offset surfaces, residual-terminated.
pub(crate) fn solve_corner_ball(
    supports: [&SurfaceGeometry; 3],
    seed: Point3,
    r: f64,
    tol: f64,
) -> Option<Point3> {
    let mut p = seed;
    for _ in 0..500 {
        let mut worst = 0.0_f64;
        for s in supports {
            let (_, off) = project_to_offset(s, &p, r)?;
            p = off;
        }
        for s in supports {
            let (d, _) = project_to_offset(s, &p, r)?;
            worst = worst.max(d);
        }
        if worst < tol {
            return Some(p);
        }
    }
    None
}

/// Contact of the ball centered at `c` on support `s`.
pub(crate) fn ball_contact(s: &SurfaceGeometry, c: &Point3) -> Option<Point3> {
    let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(s, c)?;
    Some(s.point_at(u, v))
}

/// Sample the extension center path from `end_center` to `c`: linear
/// blend refined onto the intersection of the two offset surfaces
/// (shared support + the fillet's own support), residual-terminated.
pub(crate) fn extension_center_path(
    shared: &SurfaceGeometry,
    own: &SurfaceGeometry,
    end_center: Point3,
    c: Point3,
    r: f64,
    samples: usize,
    tol: f64,
) -> Option<Vec<Point3>> {
    let n = samples.max(3);
    let mut path = Vec::with_capacity(n);
    for i in 0..n {
        let t = i as f64 / (n - 1) as f64;
        let mut p = Point3 {
            x: end_center.x + t * (c.x - end_center.x),
            y: end_center.y + t * (c.y - end_center.y),
            z: end_center.z + t * (c.z - end_center.z),
        };
        let mut ok = false;
        for _ in 0..500 {
            let (_, off_a) = project_to_offset(shared, &p, r)?;
            p = off_a;
            let (_, off_b) = project_to_offset(own, &p, r)?;
            p = off_b;
            let (da, _) = project_to_offset(shared, &p, r)?;
            let (db, _) = project_to_offset(own, &p, r)?;
            if da.max(db) < tol {
                ok = true;
                break;
            }
        }
        if !ok {
            return None;
        }
        path.push(p);
    }
    Some(path)
}

/// Solve the full corner-cap geometry. `end_centers` are the two
/// fillets' cross-section centers at their corner-end stations;
/// `own_supports` their non-shared supports in the same order.
pub(crate) fn solve_corner_cap(
    shared: &SurfaceGeometry,
    own_supports: [&SurfaceGeometry; 2],
    end_centers: [Point3; 2],
    r: f64,
    tol: f64,
) -> Option<CornerCapGeometry> {
    let seed = Point3 {
        x: 0.5 * (end_centers[0].x + end_centers[1].x),
        y: 0.5 * (end_centers[0].y + end_centers[1].y),
        z: 0.5 * (end_centers[0].z + end_centers[1].z),
    };
    let center = solve_corner_ball([shared, own_supports[0], own_supports[1]], seed, r, tol)?;
    let shared_contact = ball_contact(shared, &center)?;
    let own_contacts = [
        ball_contact(own_supports[0], &center)?,
        ball_contact(own_supports[1], &center)?,
    ];
    let extension_centers = [
        extension_center_path(shared, own_supports[0], end_centers[0], center, r, 9, tol)?,
        extension_center_path(shared, own_supports[1], end_centers[1], center, r, 9, tol)?,
    ];
    Some(CornerCapGeometry {
        center,
        shared_contact,
        own_contacts,
        extension_centers,
        shared_surface: shared.clone(),
    })
}

/// The end cross-section center of a fillet at the junction: the
/// point at distance `r` from both of the fillet's supports nearest
/// the two given trim endpoints (reconstructed by offset projection
/// from the trim-end midpoint pushed off the supports).
pub(crate) fn end_station_center(
    support_1: &SurfaceGeometry,
    support_2: &SurfaceGeometry,
    trim_end_1: Point3,
    trim_end_2: Point3,
    r: f64,
    tol: f64,
) -> Option<Point3> {
    let mid = Point3 {
        x: 0.5 * (trim_end_1.x + trim_end_2.x),
        y: 0.5 * (trim_end_1.y + trim_end_2.y),
        z: 0.5 * (trim_end_1.z + trim_end_2.z),
    };
    // Seed off the supports: push the midpoint away from support 1
    // along the direction that increases distance to it.
    let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(support_1, &mid)?;
    let q = support_1.point_at(u, v);
    let d = mid.distance_to(&q);
    let seed = if d > 1e-9 {
        q + (mid - q) * (r / d)
    } else {
        mid
    };
    let mut p = seed;
    for _ in 0..500 {
        let (_, off_a) = project_to_offset(support_1, &p, r)?;
        p = off_a;
        let (_, off_b) = project_to_offset(support_2, &p, r)?;
        p = off_b;
        let (da, _) = project_to_offset(support_1, &p, r)?;
        let (db, _) = project_to_offset(support_2, &p, r)?;
        if da.max(db) < tol {
            // Anchor near the trim ends: the two-offset intersection
            // is a curve; pull the iterate's nearest point toward the
            // station by re-seeding is unnecessary — the seed already
            // sits at the station and the crawl stays local for the
            // corner-scale spans this solves.
            return Some(p);
        }
    }
    None
}

/// Convenience for tests and the carve: distance from a point to a
/// surface via projection.
pub(crate) fn off_surface(s: &SurfaceGeometry, p: &Point3) -> Option<f64> {
    let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(s, p)?;
    Some(p.distance_to(&s.point_at(u, v)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Vec3;

    fn wall() -> SurfaceGeometry {
        SurfaceGeometry::Cylinder {
            origin: Point3 {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            axis: Vec3 {
                x: 0.0,
                y: 0.0,
                z: 1.0,
            },
            radius: 1.0,
        }
    }

    fn plane(origin: Point3, normal: Vec3) -> SurfaceGeometry {
        let u_axis = if normal.z.abs() > 0.9 {
            Vec3 {
                x: 1.0,
                y: 0.0,
                z: 0.0,
            }
        } else {
            Vec3 {
                x: 0.0,
                y: 0.0,
                z: 1.0,
            }
        };
        let v_axis = normal.cross(&u_axis);
        SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        }
    }

    /// The boss -y top corner (the measured defect corner): wall
    /// r=1, top plane z=1.25, side plane y=-0.25, fillet r=0.05.
    /// The solved ball must land at the hand-derived center
    /// (1.0063, -0.30, 1.30) with contacts on all three supports.
    #[test]
    fn boss_defect_corner_ball_and_extensions_solve() {
        let w = wall();
        let top = plane(
            Point3 {
                x: 1.0,
                y: -0.2,
                z: 1.25,
            },
            Vec3 {
                x: 0.0,
                y: 0.0,
                z: 1.0,
            },
        );
        let side = plane(
            Point3 {
                x: 1.0,
                y: -0.25,
                z: 1.0,
            },
            Vec3 {
                x: 0.0,
                y: -1.0,
                z: 0.0,
            },
        );
        let r = 0.05;
        // End centers from the measured trim endpoints.
        let top_end = end_station_center(
            &w,
            &top,
            Point3 {
                x: 0.9682,
                y: -0.25,
                z: 1.30,
            },
            Point3 {
                x: 1.0167,
                y: -0.2625,
                z: 1.25,
            },
            r,
            1e-9,
        )
        .expect("top end station");
        let vert_end = end_station_center(
            &w,
            &side,
            Point3 {
                x: 0.9583,
                y: -0.2857,
                z: 1.25,
            },
            Point3 {
                x: 1.0062,
                y: -0.25,
                z: 1.25,
            },
            r,
            1e-9,
        )
        .expect("vertical end station");
        let cap = solve_corner_cap(&w, [&top, &side], [top_end, vert_end], r, 1e-9)
            .expect("corner cap solves");
        // Hand-derived ball center for this corner.
        let expect = Point3 {
            x: (1.05_f64 * 1.05 - 0.09).sqrt(),
            y: -0.30,
            z: 1.30,
        };
        assert!(
            cap.center.distance_to(&expect) < 1e-6,
            "ball center {:?} vs {:?}",
            cap.center,
            expect
        );
        // Every extension center is exactly r from the shared wall
        // and its own plane.
        for (path, own) in cap.extension_centers.iter().zip([&top, &side]) {
            for p in path {
                let dw = off_surface(&w, p).unwrap();
                let dp = off_surface(own, p).unwrap();
                assert!(
                    (dw - r).abs() < 1e-8 && (dp - r).abs() < 1e-8,
                    "extension center off offsets: wall {dw:.2e} own {dp:.2e}"
                );
            }
        }
        // Contacts lie on their supports.
        assert!(off_surface(&w, &cap.shared_contact).unwrap() < 1e-8);
        assert!(off_surface(&top, &cap.own_contacts[0]).unwrap() < 1e-8);
        assert!(off_surface(&side, &cap.own_contacts[1]).unwrap() < 1e-8);
    }
}
