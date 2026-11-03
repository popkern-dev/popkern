//! Partitions where many bodies crowd one place, and the N-ary entry's
//! recorded pairwise fallback.

use popkern_core::geometry::{PortableMath, Transform3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::partition;
use crate::error::BoolError;

fn cube_at(x: f64, y: f64, z: f64) -> BRep {
    BRep::unit_cube(Tolerance::default())
        .and_then(|c| c.transform(&Transform3::translation(x, y, z)))
        .expect("cube")
}

/// A planar body's volume: the divergence theorem over its faces' fan
/// triangles (exact on convex planar faces).
fn volume(b: &BRep) -> f64 {
    use popkern_core::geometry::Vec3;
    let mut v = 0.0;
    for (fh, _) in &b.faces {
        let ring = b.face_outer_loop_positions(fh);
        let Some(p0) = ring.first() else { continue };
        let o = Vec3::new(p0.x, p0.y, p0.z);
        for w in ring[1..].windows(2) {
            let p1 = Vec3::new(w[0].x, w[0].y, w[0].z);
            let p2 = Vec3::new(w[1].x, w[1].y, w[1].z);
            v += o.dot(&p1.cross(&p2)) / 6.0;
        }
    }
    v.abs()
}

fn star_rods(n: usize, tol: Tolerance) -> Vec<BRep> {
    (0..n)
        .map(|k| {
            BRep::cylinder(0.4 + 0.03 * k as f64, 3.0, tol)
                .and_then(|c| c.transform(&Transform3::translation(0.0, 0.017 * k as f64, -1.5)))
                .and_then(|c| c.transform(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2)))
                .and_then(|c| c.transform(&Transform3::rotation_z(0.3 * k as f64 + 0.05)))
                .expect("rod")
        })
        .collect()
}

/// Whether `p` lies inside rod `k` of [`star_rods`], from the rod's
/// own placement (independent of any boolean), with the distance to
/// its boundary.
fn in_star_rod(k: usize, p: &popkern_core::geometry::Point3) -> (bool, f64) {
    let (s, c) = (0.3 * k as f64 + 0.05).psin_cos();
    // Undo the z rotation, then the y quarter turn: the rod's axis is
    // the world x axis after them, offset along y.
    let (x, y) = (c * p.x + s * p.y, -s * p.x + c * p.y);
    let (axial, rho) = (x, (y - 0.017 * k as f64).phypot(p.z));
    let r = 0.4 + 0.03 * k as f64;
    let margin = (r - rho).abs().min((1.5 - axial.abs()).abs());
    (rho < r && axial.abs() < 1.5, margin)
}

/// Six rods crossing near one point crowd triple points of different
/// pairs within tolerances of each other: their crossings share one
/// vertex (else a micro edge the radial walk ties on). The union
/// agrees with the rods themselves at sample points.
#[test]
fn six_rods_crossing_near_one_point_partition() {
    let tol = Tolerance::default();
    let rods = star_rods(6, tol);
    let refs: Vec<&BRep> = rods.iter().collect();
    let p = partition(&refs, tol).expect("partition");
    let union = p
        .select(&super::SetExpr::Any, &crate::BoolOptions::default())
        .expect("union");
    assert!(
        union.validate_all().is_empty(),
        "{:?}",
        union.validate_all()
    );
    let boxes = crate::intersect::compute_face_aabbs(&union);
    let mut checked = 0;
    for i in 0..9 {
        for j in 0..9 {
            for l in 0..5 {
                let at = popkern_core::geometry::Point3::new(
                    -1.6 + 0.4 * i as f64 + 0.013,
                    -1.6 + 0.4 * j as f64 + 0.007,
                    -0.5 + 0.25 * l as f64 + 0.003,
                );
                let rods_at: Vec<(bool, f64)> = (0..6).map(|k| in_star_rod(k, &at)).collect();
                if rods_at.iter().any(|(_, m)| *m < 1e-3) {
                    continue;
                }
                let expect = rods_at.iter().any(|(inside, _)| *inside);
                let got = crate::classify::point_in_solid_final(&at, &union, &boxes, &tol)
                    .expect("classified");
                assert_eq!(
                    got == crate::error::PointClassification::Inside,
                    expect,
                    "at {at:?}"
                );
                checked += 1;
            }
        }
    }
    assert!(checked > 300, "{checked}");
}

/// A declined partition falls back to two-operand booleans over the
/// expression, recorded as a `PartitionFallback` shadow line; an answer
/// (a non-manifold result) passes through unrescued.
#[test]
fn a_declined_partition_falls_back_to_pairwise_booleans() {
    let tol = Tolerance::default();
    let (a, b, c) = (
        cube_at(0.0, 0.0, 0.0),
        cube_at(0.5, 0.0, 0.0),
        cube_at(0.25, 0.5, 0.0),
    );
    let refs = [&a, &b, &c];
    let decline = || {
        Err(BoolError::MergeFailed {
            stage: "complex".into(),
            reason: "a declined stage".into(),
        })
    };
    let options = crate::BoolOptions::default();
    let a_not_c = super::SetExpr::Minus(
        Box::new(super::SetExpr::Body(0)),
        Box::new(super::SetExpr::Body(2)),
    );
    for expr in [super::SetExpr::Any, super::SetExpr::All, a_not_c] {
        let fell_back =
            super::nary::or_fold(decline(), &refs, &expr, tol, &options).expect("fallback");
        let served = partition(&refs, tol)
            .and_then(|p| p.select(&expr, &options))
            .expect("partition");
        assert!(fell_back.validate_all().is_empty(), "{expr:?}");
        assert_eq!(fell_back.faces.len(), served.faces.len(), "{expr:?}");
        assert!(
            (volume(&fell_back) - volume(&served)).abs() < 1e-9,
            "{expr:?}"
        );
    }
    let answer = super::nary::or_fold(
        Err(BoolError::Cancelled),
        &refs,
        &super::SetExpr::Any,
        tol,
        &options,
    );
    assert!(matches!(answer, Err(BoolError::Cancelled)));
}
