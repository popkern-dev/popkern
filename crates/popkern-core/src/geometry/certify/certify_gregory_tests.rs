//! Rational Gregory exact forms: enclosures match the patch away from
//! blended corners and decline (never guess) at them.

use super::super::*;
use super::{Rng, binonic};
use crate::geometry::gregory::GregoryPatch;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};

fn rational_gregories() -> Vec<SurfaceGeometry> {
    use crate::geometry::gregory::{
        DualJetCorner, RectangularGregory, RectangularGregoryBiNonicDualJet,
    };
    let SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonic(plain)) = binonic() else {
        return Vec::new();
    };
    let g = plain.grid;
    let nudge = |p: Point3, dz: f64| Point3::new(p.x, p.y, p.z + dz);
    let claims = |i: usize, j: usize, di: isize, dj: isize| {
        let at =
            |a: isize, b: isize| g[(i as isize + a * di) as usize][(j as isize + b * dj) as usize];
        DualJetCorner {
            cell_11: (nudge(at(1, 1), 0.05), nudge(at(1, 1), -0.04)),
            cell_12: (nudge(at(1, 2), 0.03), at(1, 2)),
            cell_21: (at(2, 1), nudge(at(2, 1), -0.02)),
            cell_22: (nudge(at(2, 2), 0.01), nudge(at(2, 2), 0.02)),
        }
    };
    let dual_jet = RectangularGregoryBiNonicDualJet {
        grid: g,
        corners: [
            Some(claims(0, 0, 1, 1)),
            None,
            Some(claims(9, 9, -1, -1)),
            None,
        ],
    };
    let mut boundary_cps = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
    for (i, row) in boundary_cps.iter_mut().enumerate() {
        for (j, p) in row.iter_mut().enumerate() {
            *p = Point3::new(i as f64, j as f64, ((i * 3 + j) % 4) as f64 * 0.2);
        }
    }
    let inner = |i: usize, j: usize| {
        (
            nudge(boundary_cps[i][j], 0.1),
            nudge(boundary_cps[i][j], -0.1),
        )
    };
    let cubic = RectangularGregory {
        boundary_cps,
        inner_cps: [inner(1, 1), inner(1, 2), inner(2, 2), inner(2, 1)],
    };
    vec![
        SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonicDualJet(Box::new(dual_jet))),
        SurfaceGeometry::Gregory(GregoryPatch::Rectangular(Box::new(cubic))),
    ]
}

#[test]
fn rational_gregory_forms_enclose_the_patch_and_decline_at_blended_corners() {
    let mut rng = Rng(17);
    for surface in rational_gregories() {
        let mut source = SurfaceJetSource::new(&surface).unwrap_or_else(|e| panic!("{e}"));
        for _ in 0..60 {
            let (u, v) = (rng.range(0.03, 0.97), rng.range(0.03, 0.97));
            let jet = source
                .jet(Iv::point(u), Iv::point(v))
                .unwrap_or_else(|e| panic!("point ({u}, {v}): {e}"));
            let p = surface.point_at(u, v);
            let err = (jet.partial(0, 0).mid() - Vec3::new(p.x, p.y, p.z)).length();
            assert!(err < 1e-11, "Gregory value at ({u}, {v}) off by {err:e}");
            let (u0, v0) = (rng.range(0.03, 0.8), rng.range(0.03, 0.8));
            let boxed = source
                .jet(Iv::new(u0, u0 + 0.15), Iv::new(v0, v0 + 0.15))
                .unwrap_or_else(|e| panic!("box at ({u0}, {v0}): {e}"));
            for k in 0..=6 {
                for l in 0..=6 {
                    let (su, sv) = (u0 + 0.025 * k as f64, v0 + 0.025 * l as f64);
                    let q = surface.point_at(su, sv);
                    let pos = boxed.partial(0, 0);
                    let slack = 1e-12;
                    assert!(
                        pos.x.lo - slack <= q.x
                            && q.x <= pos.x.hi + slack
                            && pos.z.lo - slack <= q.z
                            && q.z <= pos.z.hi + slack,
                        "Gregory box misses ({su}, {sv})"
                    );
                }
            }
        }
        // A box touching a blended corner cannot bound the quotient.
        assert_eq!(
            source.jet(Iv::new(0.0, 0.1), Iv::new(0.0, 0.1)),
            Err(CertifyDecline::DenominatorContainsZero)
        );
    }
}

/// A pcurve on a patch edge (`u = 0` exactly) yields a subnormal-width
/// box straddling 0. Restricting to it by the `alpha/beta` ratio
/// divided by a subnormal interval and put zero in a polynomial's
/// unit weight hull (measured on pentagon_tower cap seams); the
/// narrow-box path must enclose it tightly instead.
#[test]
fn a_patch_edge_subnormal_box_encloses_tightly() {
    let surface = binonic();
    let mut source = SurfaceJetSource::new(&surface).unwrap_or_else(|e| panic!("{e}"));
    let edge = Iv::new(-1.5e-323, 1.5e-323);
    let near_half = Iv::new(0.5f64.next_down().next_down(), 0.5f64.next_up().next_up());
    for (u, v) in [(edge, near_half), (near_half, edge), (edge, edge)] {
        let jet = source
            .jet(u, v)
            .unwrap_or_else(|e| panic!("patch-edge box declined: {e}"));
        let p = jet.partial(0, 0);
        assert!(
            p.x.width() < 1e-13 && p.y.width() < 1e-13 && p.z.width() < 1e-13,
            "{p:?}"
        );
    }
}
