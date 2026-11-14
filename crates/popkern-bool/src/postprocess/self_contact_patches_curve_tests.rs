use super::super::window_patches;
use super::*;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};

/// A bicubic sheet over [0, 4]², waved in z by up to 0.5, unit weights
/// but one of 1.6.
fn sheet() -> NurbsSurface {
    let (p, n) = (3, 6);
    let knots = vec![0.0, 0.0, 0.0, 0.0, 0.3, 0.6, 1.0, 1.0, 1.0, 1.0];
    let mut control_points = Vec::new();
    let mut weights = Vec::new();
    for i in 0..n {
        for j in 0..n {
            let z = 0.25 * ((i * 3 + j * 5) % 5) as f64 - 0.5;
            control_points.push(Point3::new(0.8 * i as f64, 0.8 * j as f64, z));
            weights.push(if (i, j) == (2, 3) { 1.6 } else { 1.0 });
        }
    }
    NurbsSurface {
        degree_u: p,
        degree_v: p,
        knots_u: knots.clone(),
        knots_v: knots,
        control_points,
        num_u: n,
        num_v: n,
        weights,
    }
}

/// A small deterministic generator in [0, 1).
fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    (*seed >> 11) as f64 / (1u64 << 53) as f64
}

/// Segments and arcs about the sheet: whatever `curve_apart` clears
/// lies farther than the pad from the sheet at a 120² lattice of it and
/// 200 points of the curve (sampled distances bound the true one from
/// above, so a cleared curve sampled within the pad would be a false
/// clearance), and both outcomes occur.
#[test]
fn curves_cleared_of_the_patches_lie_farther_than_the_pad() {
    let s = sheet();
    let window = ((0.0, 1.0), (0.0, 1.0));
    let patches = window_patches(&s, &window).expect("patches");
    let surface = SurfaceGeometry::Nurbs(s);
    let n = 120;
    let lattice: Vec<Point3> = (0..=n)
        .flat_map(|i| (0..=n).map(move |j| (i, j)))
        .map(|(i, j)| surface.point_at(i as f64 / n as f64, j as f64 / n as f64))
        .collect();
    let pad = 0.1;
    let (mut cleared, mut kept) = (0, 0);
    let mut seed = 5u64;
    for k in 0..200 {
        let at = Point3::new(
            4.0 * lcg(&mut seed),
            4.0 * lcg(&mut seed),
            1.6 * lcg(&mut seed) - 0.8,
        );
        let d = Vec3::new(
            lcg(&mut seed) - 0.5,
            lcg(&mut seed) - 0.5,
            0.3 * (lcg(&mut seed) - 0.5),
        );
        let (curve, range) = if k % 2 == 0 {
            (
                CurveGeometry::Line {
                    origin: at,
                    direction: d * 3.0,
                },
                (0.0, 1.0),
            )
        } else {
            (
                CurveGeometry::Circle {
                    center: at,
                    normal: d.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0)),
                    radius: 0.2 + 0.8 * lcg(&mut seed),
                },
                (0.0, 4.0),
            )
        };
        let nurbs = curve.to_nurbs_in_range(range);
        if curve_apart(&nurbs, &patches, pad, 8) {
            cleared += 1;
            for i in 0..=200 {
                let p = curve.point_at(range.0 + (range.1 - range.0) * f64::from(i) / 200.0);
                let near = lattice
                    .iter()
                    .map(|q| q.distance_to(&p))
                    .fold(f64::INFINITY, f64::min);
                assert!(near > pad, "curve {k} cleared {near} from the sheet");
            }
        } else {
            kept += 1;
        }
    }
    assert!(cleared > 30 && kept > 30, "cleared {cleared}, kept {kept}");
}
