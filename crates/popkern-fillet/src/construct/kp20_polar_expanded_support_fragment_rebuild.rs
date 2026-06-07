//! Boundary-and-tangent constrained macro-fragment rebuild.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::Kp20Bi6Patch;

const FREE_INTERIOR_REBUILD_STRENGTH: f64 = 0.85;

pub(super) fn rebuild_macro_fragments_from_fixed_boundary(
    fragments: [Kp20Bi6Patch; 3],
) -> [Kp20Bi6Patch; 3] {
    fragments.map(rebuild_fragment_from_boundary)
}

fn rebuild_fragment_from_boundary(fragment: Kp20Bi6Patch) -> Kp20Bi6Patch {
    let mut rebuilt = fragment;
    for u in 2..5 {
        for v in 2..5 {
            let boundary_point = tangent_ribbon_point(&fragment, u, v);
            rebuilt.cps[u][v] = fragment.cps[u][v]
                + (boundary_point - fragment.cps[u][v]) * FREE_INTERIOR_REBUILD_STRENGTH;
        }
    }
    rebuilt
}

fn tangent_ribbon_point(fragment: &Kp20Bi6Patch, u: usize, v: usize) -> Point3 {
    let tu = u as f64 / 6.0;
    let tv = v as f64 / 6.0;
    let row_point = hermite_point(
        fragment.cps[0][v],
        (fragment.cps[1][v] - fragment.cps[0][v]) * 6.0,
        fragment.cps[6][v],
        (fragment.cps[6][v] - fragment.cps[5][v]) * 6.0,
        tu,
    );
    let column_point = hermite_point(
        fragment.cps[u][0],
        (fragment.cps[u][1] - fragment.cps[u][0]) * 6.0,
        fragment.cps[u][6],
        (fragment.cps[u][6] - fragment.cps[u][5]) * 6.0,
        tv,
    );
    lerp_point(row_point, column_point, 0.5)
}

fn hermite_point(
    p0: Point3,
    d0: popkern_core::geometry::Vec3,
    p1: Point3,
    d1: popkern_core::geometry::Vec3,
    t: f64,
) -> Point3 {
    let t2 = t * t;
    let t3 = t2 * t;
    let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
    let h10 = t3 - 2.0 * t2 + t;
    let h01 = -2.0 * t3 + 3.0 * t2;
    let h11 = t3 - t2;
    Point3::new(
        h00 * p0.x + h10 * d0.x + h01 * p1.x + h11 * d1.x,
        h00 * p0.y + h10 * d0.y + h01 * p1.y + h11 * d1.y,
        h00 * p0.z + h10 * d0.z + h01 * p1.z + h11 * d1.z,
    )
}

fn lerp_point(a: Point3, b: Point3, t: f64) -> Point3 {
    a + (b - a) * t
}
