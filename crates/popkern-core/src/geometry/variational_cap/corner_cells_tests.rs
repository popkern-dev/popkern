use super::*;
use crate::geometry::Vec3;
use crate::geometry::surface::SurfaceEval;
use crate::geometry::variational_cap::space::{SectorSpace, apply3};

/// A deterministic coefficient set over a space.
fn coefficients(space: &SectorSpace) -> [Vec<f64>; 3] {
    let total = space.total();
    [
        (0..total).map(|k| (k as f64 * 0.37).sin()).collect(),
        (0..total).map(|k| (k as f64 * 0.11).cos()).collect(),
        (0..total).map(|k| (k as f64 * 0.23).sin() * 0.5).collect(),
    ]
}

/// Copy each dual cell's u-claim into its v-claim.
fn equal_claims(space: &SectorSpace, coords: &mut [Vec<f64>; 3]) {
    for s in 0..space.n {
        for slot in 0..CELLS {
            let (i, j) = cell_control(space.m, slot);
            for c in coords.iter_mut() {
                c[space.claim_index(s, slot)] = c[space.index(s, i, j)];
            }
        }
    }
}

const SAMPLES: [(f64, f64); 6] = [
    (0.1, 0.7),
    (0.5, 0.5),
    (0.93, 0.97),
    (0.999, 0.3),
    (0.6, 0.999),
    (0.99, 0.99),
];

#[test]
fn equal_claims_reproduce_the_polynomial_sector() {
    let plain = SectorSpace::bezier(4, 9);
    let cells = SectorSpace::bezier(4, 9).with_corner_cells();
    assert!(cells.corner_cells);
    let mut x = coefficients(&cells);
    equal_claims(&cells, &mut x);
    // The plain space's coefficients: each sector's m² block.
    let y: [Vec<f64>; 3] = std::array::from_fn(|c| {
        (0..plain.n)
            .flat_map(|s| (0..plain.m * plain.m).map(move |l| (s, l)))
            .map(|(s, l)| x[c][s * cells.per() + l])
            .collect()
    });
    for s in 0..cells.n {
        for &(u, v) in &SAMPLES {
            for (du, dv) in [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)] {
                let a = apply3(&cells.functional(s, u, v, du, dv), &x);
                let b = apply3(&plain.functional(s, u, v, du, dv), &y);
                assert!(
                    (a - b).length() < 1e-9 * (1.0 + b.length()),
                    "s{s} ({u}, {v}) d({du},{dv}): {a:?} vs {b:?}"
                );
            }
        }
    }
}

#[test]
fn each_rim_reads_only_its_own_claim() {
    let space = SectorSpace::bezier(5, 9).with_corner_cells();
    let m = space.m;
    let claims_v: Vec<usize> = (0..CELLS).map(|k| space.claim_index(2, k)).collect();
    let claims_u: Vec<usize> = (0..CELLS)
        .map(|k| {
            let (i, j) = cell_control(m, k);
            space.index(2, i, j)
        })
        .collect();
    let derivatives = [(0, 0), (1, 0), (0, 1), (2, 0), (1, 1), (0, 2)];
    for t in [0.2, 0.7, 0.95, 0.999, 1.0] {
        for &(du, dv) in &derivatives {
            // The rim u = 1, the corner included under its owner.
            let on_u = space.functional_at(2, 1.0, t, du, dv, CornerOwner::URim);
            assert!(
                on_u.iter().all(|(i, w)| !claims_v.contains(i) || *w == 0.0),
                "u-rim t {t} d({du},{dv}) reads a v-claim: {on_u:?}"
            );
            let on_v = space.functional_at(2, t, 1.0, du, dv, CornerOwner::VRim);
            assert!(
                on_v.iter().all(|(i, w)| !claims_u.contains(i) || *w == 0.0),
                "v-rim t {t} d({du},{dv}) reads a u-claim: {on_v:?}"
            );
        }
    }
}

#[test]
fn the_face_surface_is_the_dual_jet_patch_of_the_blend() {
    let space = SectorSpace::bezier(3, 9).with_corner_cells();
    let x = coefficients(&space);
    for s in 0..space.n {
        let face = space.face_surface(&x, s);
        assert!(matches!(
            face,
            crate::geometry::SurfaceGeometry::Gregory(
                crate::geometry::gregory::GregoryPatch::RectangularBiNonicDualJet(_)
            )
        ));
        for &(u, v) in &SAMPLES {
            for (du, dv) in [(0, 0), (1, 0), (0, 1)] {
                let a = apply3(&space.functional(s, u, v, du, dv), &x);
                let b = match (du, dv) {
                    (0, 0) => {
                        let p = face.point_at(u, v);
                        Vec3::new(p.x, p.y, p.z)
                    }
                    (1, 0) => face.partial_u(u, v),
                    _ => face.partial_v(u, v),
                };
                assert!(
                    (a - b).length() < 1e-9 * (1.0 + a.length()),
                    "s{s} ({u}, {v}) d({du},{dv}): blend {a:?} patch {b:?}"
                );
            }
        }
    }
}
