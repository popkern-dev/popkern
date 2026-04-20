//! Tests for [`super::build_guide_from_hermite`].
//!
//! The **reconstruction invariant** is the load-bearing contract:
//! the built patch must exactly reproduce the Hermite data at each
//! of the 4 corners (position, u-tangent, v-tangent, mixed twist).
//! If any of these drift, the higher-level K-P orchestration
//! (Shirman-Séquin enclosure at the central vertex + K-P 2004
//! §4.2 outer-corner formula) won't compose correctly — adjacent
//! guides computed from compatible Hermite data wouldn't actually
//! produce G1-matching surfaces.

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::gregory::KpGuidePatch;

const EPS: f64 = 1e-12;

fn vec_max(a: Vec3, b: Vec3) -> f64 {
    ((a.x - b.x).abs())
        .max((a.y - b.y).abs())
        .max((a.z - b.z).abs())
}

fn point_max(a: Point3, b: Point3) -> f64 {
    ((a.x - b.x).abs())
        .max((a.y - b.y).abs())
        .max((a.z - b.z).abs())
}

/// A set of 4 Hermite corners describing an asymmetric bicubic
/// patch with non-zero twist at every corner. Used as the
/// reference fixture for reconstruction-invariant tests.
fn asymmetric_hermite() -> [GuideCornerData; 4] {
    [
        // c0 at (0, 0)
        GuideCornerData {
            position: Point3::new(0.0, 0.0, 0.0),
            tangent_u: Vec3::new(1.0, 0.2, 0.1),
            tangent_v: Vec3::new(0.1, 1.2, -0.2),
            twist: Vec3::new(0.3, -0.15, 0.4),
        },
        // c1 at (1, 0)
        GuideCornerData {
            position: Point3::new(1.0, 0.0, 0.2),
            tangent_u: Vec3::new(1.1, -0.3, 0.4),
            tangent_v: Vec3::new(0.2, 1.3, 0.1),
            twist: Vec3::new(-0.2, 0.25, -0.3),
        },
        // c2 at (1, 1)
        GuideCornerData {
            position: Point3::new(1.2, 1.3, -0.1),
            tangent_u: Vec3::new(0.9, 0.4, -0.2),
            tangent_v: Vec3::new(-0.3, 0.8, 0.3),
            twist: Vec3::new(0.15, 0.2, -0.35),
        },
        // c3 at (0, 1)
        GuideCornerData {
            position: Point3::new(-0.1, 1.2, 0.15),
            tangent_u: Vec3::new(1.2, -0.1, 0.05),
            tangent_v: Vec3::new(0.0, 1.1, -0.15),
            twist: Vec3::new(-0.1, -0.25, 0.3),
        },
    ]
}

// ---- Reconstruction invariant: Hermite data → CPs → Hermite data ----

#[test]
fn position_recovered_at_each_corner() {
    let corners = asymmetric_hermite();
    let g = build_guide_from_hermite(corners);
    assert!(
        point_max(g.point_at(0.0, 0.0), corners[0].position) < EPS,
        "c0 position drift"
    );
    assert!(
        point_max(g.point_at(1.0, 0.0), corners[1].position) < EPS,
        "c1 position drift"
    );
    assert!(
        point_max(g.point_at(1.0, 1.0), corners[2].position) < EPS,
        "c2 position drift"
    );
    assert!(
        point_max(g.point_at(0.0, 1.0), corners[3].position) < EPS,
        "c3 position drift"
    );
}

#[test]
fn tangent_u_recovered_at_each_corner() {
    let corners = asymmetric_hermite();
    let g = build_guide_from_hermite(corners);
    assert!(
        vec_max(g.partial_u(0.0, 0.0), corners[0].tangent_u) < EPS,
        "c0 ∂u drift"
    );
    assert!(
        vec_max(g.partial_u(1.0, 0.0), corners[1].tangent_u) < EPS,
        "c1 ∂u drift"
    );
    assert!(
        vec_max(g.partial_u(1.0, 1.0), corners[2].tangent_u) < EPS,
        "c2 ∂u drift"
    );
    assert!(
        vec_max(g.partial_u(0.0, 1.0), corners[3].tangent_u) < EPS,
        "c3 ∂u drift"
    );
}

#[test]
fn tangent_v_recovered_at_each_corner() {
    let corners = asymmetric_hermite();
    let g = build_guide_from_hermite(corners);
    assert!(
        vec_max(g.partial_v(0.0, 0.0), corners[0].tangent_v) < EPS,
        "c0 ∂v drift"
    );
    assert!(
        vec_max(g.partial_v(1.0, 0.0), corners[1].tangent_v) < EPS,
        "c1 ∂v drift"
    );
    assert!(
        vec_max(g.partial_v(1.0, 1.0), corners[2].tangent_v) < EPS,
        "c2 ∂v drift"
    );
    assert!(
        vec_max(g.partial_v(0.0, 1.0), corners[3].tangent_v) < EPS,
        "c3 ∂v drift"
    );
}

#[test]
fn twist_recovered_at_each_corner() {
    let corners = asymmetric_hermite();
    let g = build_guide_from_hermite(corners);
    assert!(
        vec_max(g.second_partial_uv(0.0, 0.0), corners[0].twist) < EPS,
        "c0 ∂²uv drift"
    );
    assert!(
        vec_max(g.second_partial_uv(1.0, 0.0), corners[1].twist) < EPS,
        "c1 ∂²uv drift"
    );
    assert!(
        vec_max(g.second_partial_uv(1.0, 1.0), corners[2].twist) < EPS,
        "c2 ∂²uv drift"
    );
    assert!(
        vec_max(g.second_partial_uv(0.0, 1.0), corners[3].twist) < EPS,
        "c3 ∂²uv drift"
    );
}

// ---- Special cases ----

#[test]
fn flat_unit_square_from_hermite_data() {
    // Corners on the z = 0 plane, consistent tangents in the +u
    // and +v directions, zero twist → the built guide should be
    // the flat unit square (matching `flat_guide()` in
    // kp_guide_patch_tests).
    let corners = [
        GuideCornerData::no_twist(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ),
        GuideCornerData::no_twist(
            Point3::new(1.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ),
        GuideCornerData::no_twist(
            Point3::new(1.0, 1.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ),
        GuideCornerData::no_twist(
            Point3::new(0.0, 1.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ),
    ];
    let g = build_guide_from_hermite(corners);
    for &u in &[0.0, 0.25, 0.5, 0.75, 1.0] {
        for &v in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            let pt = g.point_at(u, v);
            assert!(pt.z.abs() < EPS, "z drift at ({u}, {v}): {pt:?}");
            assert!((pt.x - u).abs() < EPS);
            assert!((pt.y - v).abs() < EPS);
        }
    }
}

#[test]
fn no_twist_helper_zeroes_twist_field() {
    let c = GuideCornerData::no_twist(
        Point3::new(1.0, 2.0, 3.0),
        Vec3::new(0.5, 0.0, 0.0),
        Vec3::new(0.0, 0.5, 0.0),
    );
    assert_eq!(c.twist.x, 0.0);
    assert_eq!(c.twist.y, 0.0);
    assert_eq!(c.twist.z, 0.0);
}

#[test]
fn degenerate_zero_tangents_collapses_to_bilinear_corner_patch() {
    // Corners on a square with zero tangents + zero twist. The
    // inner CPs all collapse to the corner positions. This is
    // the "all zero Hermite data except positions" degenerate
    // case — the patch is still valid (all 16 CPs defined) but
    // highly pinched at the corners.
    let corners = [
        GuideCornerData::no_twist(Point3::new(0.0, 0.0, 0.0), Vec3::zero(), Vec3::zero()),
        GuideCornerData::no_twist(Point3::new(3.0, 0.0, 0.0), Vec3::zero(), Vec3::zero()),
        GuideCornerData::no_twist(Point3::new(3.0, 3.0, 0.0), Vec3::zero(), Vec3::zero()),
        GuideCornerData::no_twist(Point3::new(0.0, 3.0, 0.0), Vec3::zero(), Vec3::zero()),
    ];
    let g = build_guide_from_hermite(corners);
    // All 8 edge CPs collapse to the adjacent corner.
    assert!(point_max(g.grid[1][0], Point3::new(0.0, 0.0, 0.0)) < EPS);
    assert!(point_max(g.grid[2][0], Point3::new(3.0, 0.0, 0.0)) < EPS);
    assert!(point_max(g.grid[3][1], Point3::new(3.0, 0.0, 0.0)) < EPS);
    assert!(point_max(g.grid[3][2], Point3::new(3.0, 3.0, 0.0)) < EPS);
    assert!(point_max(g.grid[2][3], Point3::new(3.0, 3.0, 0.0)) < EPS);
    assert!(point_max(g.grid[1][3], Point3::new(0.0, 3.0, 0.0)) < EPS);
    assert!(point_max(g.grid[0][2], Point3::new(0.0, 3.0, 0.0)) < EPS);
    assert!(point_max(g.grid[0][1], Point3::new(0.0, 0.0, 0.0)) < EPS);
    // All 4 inner CPs collapse to the adjacent corner.
    assert!(point_max(g.grid[1][1], Point3::new(0.0, 0.0, 0.0)) < EPS);
    assert!(point_max(g.grid[2][1], Point3::new(3.0, 0.0, 0.0)) < EPS);
    assert!(point_max(g.grid[2][2], Point3::new(3.0, 3.0, 0.0)) < EPS);
    assert!(point_max(g.grid[1][2], Point3::new(0.0, 3.0, 0.0)) < EPS);
}

// ---- Round-trip: extract Hermite from known patch, rebuild, compare ----

#[test]
fn round_trip_known_bicubic_patch() {
    // Start from a hand-picked asymmetric bicubic patch. Extract
    // the 4 corners of Hermite data by evaluating the patch's
    // partials, rebuild via `build_guide_from_hermite`, and
    // check the grids agree exactly (bit-level, modulo floating-
    // point rounding).
    let mut source = KpGuidePatch::zeros();
    // Arbitrary CPs with z-variation.
    for i in 0..4 {
        for j in 0..4 {
            let u = i as f64 / 3.0;
            let v = j as f64 / 3.0;
            source.grid[i][j] = Point3::new(
                u + 0.1 * v,
                v - 0.05 * u,
                0.2 * u * (1.0 - u) + 0.1 * v * (1.0 - v) + 0.05 * u * v,
            );
        }
    }
    let hermite = [
        GuideCornerData {
            position: source.point_at(0.0, 0.0),
            tangent_u: source.partial_u(0.0, 0.0),
            tangent_v: source.partial_v(0.0, 0.0),
            twist: source.second_partial_uv(0.0, 0.0),
        },
        GuideCornerData {
            position: source.point_at(1.0, 0.0),
            tangent_u: source.partial_u(1.0, 0.0),
            tangent_v: source.partial_v(1.0, 0.0),
            twist: source.second_partial_uv(1.0, 0.0),
        },
        GuideCornerData {
            position: source.point_at(1.0, 1.0),
            tangent_u: source.partial_u(1.0, 1.0),
            tangent_v: source.partial_v(1.0, 1.0),
            twist: source.second_partial_uv(1.0, 1.0),
        },
        GuideCornerData {
            position: source.point_at(0.0, 1.0),
            tangent_u: source.partial_u(0.0, 1.0),
            tangent_v: source.partial_v(0.0, 1.0),
            twist: source.second_partial_uv(0.0, 1.0),
        },
    ];
    let rebuilt = build_guide_from_hermite(hermite);
    // A bicubic is fully determined by its 16 CPs, and a bicubic
    // is equivalently determined by its 4 corners × 4 Hermite
    // values each — same dimensionality. Round-trip equals the
    // source exactly.
    for i in 0..4 {
        for j in 0..4 {
            assert!(
                point_max(rebuilt.grid[i][j], source.grid[i][j]) < EPS,
                "grid[{i}][{j}] drift: source {:?}, rebuilt {:?}",
                source.grid[i][j],
                rebuilt.grid[i][j]
            );
        }
    }
}
