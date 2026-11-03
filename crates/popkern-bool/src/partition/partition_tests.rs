//! Partitions of overlapping solids: every region's membership matches
//! its interior point classified against each body alone.

use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::{BodySet, partition};
use crate::PointClassification;

fn cube_at(x: f64, y: f64, z: f64) -> BRep {
    BRep::unit_cube(Tolerance::default())
        .and_then(|c| c.transform(&Transform3::translation(x, y, z)))
        .expect("cube")
}

/// Membership by point classification against each body alone: the
/// oracle the partition's own membership must equal.
fn classified(bodies: &[&BRep], p: &popkern_core::geometry::Point3) -> BodySet {
    let tol = Tolerance::default();
    let inside: Vec<usize> = bodies
        .iter()
        .enumerate()
        .filter(|(_, b)| {
            let boxes = crate::intersect::compute_face_aabbs(b);
            matches!(
                crate::classify::point_in_solid_final(p, b, &boxes, &tol),
                Ok(PointClassification::Inside)
            )
        })
        .map(|(i, _)| i)
        .collect();
    BodySet::of(&inside)
}

fn check(bodies: &[&BRep], want: &[&[usize]]) {
    let p = partition(bodies, Tolerance::default()).unwrap_or_else(|e| panic!("{e}"));
    let mut got: Vec<BodySet> = p.regions().iter().map(|r| r.bodies).collect();
    for r in p.regions() {
        assert_eq!(
            classified(bodies, &r.interior_point),
            r.bodies,
            "region {:?} at {:?}",
            r.bodies,
            r.interior_point
        );
    }
    got.sort();
    let mut want: Vec<BodySet> = want.iter().map(|w| BodySet::of(w)).collect();
    want.sort();
    assert_eq!(got, want);
}

#[test]
fn two_overlapping_cubes_cut_three_regions() {
    let (a, b) = (cube_at(0.0, 0.0, 0.0), cube_at(0.5, 0.3, 0.2));
    check(&[&a, &b], &[&[0], &[1], &[0, 1]]);
}

/// Three cubes each offset along a diagonal: seven regions, including
/// the one all three hold, and records of different pairs crossing on
/// shared faces.
#[test]
fn three_overlapping_cubes_cut_seven_regions() {
    let a = cube_at(0.0, 0.0, 0.0);
    let b = cube_at(0.5, 0.3, 0.2);
    let c = cube_at(0.25, 0.6, 0.45);
    check(
        &[&a, &b, &c],
        &[&[0], &[1], &[2], &[0, 1], &[0, 2], &[1, 2], &[0, 1, 2]],
    );
}

#[test]
fn a_sheet_is_refused_by_kind() {
    let mut open = cube_at(0.0, 0.0, 0.0);
    let face = open.faces.keys().next().expect("face");
    open.retire_faces_to_boundary(&[face]).expect("retire");
    let solid = cube_at(0.5, 0.5, 0.5);
    match partition(&[&solid, &open], Tolerance::default()) {
        Err(crate::BoolError::UnsupportedBodyKind { operand: 1, .. }) => {}
        other => panic!("{other:?}"),
    }
}

/// Three orthogonal cylinders through one point (the Steinmetz
/// tricylinder): their pairwise records are curves crossing at the
/// points where all three walls meet. Regions are connected components:
/// each rod alone keeps its two ends; each pair without the third is
/// four pieces, one per sign quadrant of the two axes it spans (z ∩ x
/// outside y needs |x| and |z| both large); the three meet in one.
#[test]
fn three_orthogonal_cylinders_cut_nineteen_regions() {
    let tol = Tolerance::default();
    let rod = |turn: Transform3| {
        BRep::cylinder(0.5, 3.0, tol)
            .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, -1.5)))
            .and_then(|c| c.transform(&turn))
            .expect("rod")
    };
    let z = rod(Transform3::identity());
    let x = rod(Transform3::rotation_y(std::f64::consts::FRAC_PI_2));
    let y = rod(Transform3::rotation_x(std::f64::consts::FRAC_PI_2));
    let pairs: [&[usize]; 4] = [&[0, 1]; 4];
    let mut want: Vec<&[usize]> = vec![&[0], &[0], &[1], &[1], &[2], &[2], &[0, 1, 2]];
    want.extend(pairs);
    want.extend([&[0usize, 2][..]; 4]);
    want.extend([&[1usize, 2][..]; 4]);
    check(&[&z, &x, &y], &want);
}

/// Membership on a 9³ lattice over the bodies' boxes: `body` holds each
/// probe exactly when `want` holds for the bodies containing it (probes
/// either classification leaves undecided are skipped). `want` is
/// written independently of `SetExpr::holds`, so the expression's own
/// semantics are under test too.
fn agrees_on_lattice(bodies: &[&BRep], want: &dyn Fn(&BodySet) -> bool, body: &BRep) {
    let tol = Tolerance::default();
    let (mut lo, mut hi) = ([f64::MAX; 3], [f64::MIN; 3]);
    for b in bodies {
        for v in b.vertices.values() {
            for (k, c) in [v.position.x, v.position.y, v.position.z]
                .into_iter()
                .enumerate()
            {
                lo[k] = lo[k].min(c);
                hi[k] = hi[k].max(c);
            }
        }
    }
    let boxes = crate::intersect::compute_face_aabbs(body);
    let n = 9;
    let mut checked = 0;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let f = |a: usize, s: f64, d: usize| {
                    lo[d] + (hi[d] - lo[d]) * (a as f64 + s) / n as f64
                };
                let p = popkern_core::geometry::Point3::new(
                    f(i, 0.381_966, 0),
                    f(j, 0.414_214, 1),
                    f(k, 0.732_051, 2),
                );
                let want = want(&classified(bodies, &p));
                let got = if body.faces.is_empty() {
                    Some(false)
                } else {
                    match crate::classify::point_in_solid_final(&p, body, &boxes, &tol) {
                        Ok(PointClassification::Inside) => Some(true),
                        Ok(PointClassification::Outside) => Some(false),
                        _ => None,
                    }
                };
                if let Some(got) = got {
                    assert_eq!(got, want, "at {p:?}");
                    checked += 1;
                }
            }
        }
    }
    assert!(checked > n * n * n / 2, "{checked} probes decided");
}

/// Every set expression over three overlapping cubes selects a valid
/// body whose membership matches the expression on the lattice.
#[test]
fn three_cubes_select_every_expression() {
    use super::SetExpr;
    let a = cube_at(0.0, 0.0, 0.0);
    let b = cube_at(0.5, 0.3, 0.2);
    let c = cube_at(0.25, 0.6, 0.45);
    let bodies = [&a, &b, &c];
    let p = partition(&bodies, Tolerance::default()).expect("partition");
    let split = crate::BoolOptions {
        non_manifold: crate::NonManifoldPolicy::SplitToManifold,
    };
    // Odd keeps regions that meet only along edges (the pieces of a ⊕
    // touch where the bodies' boundaries cross): it declines by
    // default and serves as manifold pieces under the split policy.
    assert!(matches!(
        p.select(&SetExpr::Odd, &crate::BoolOptions::default()),
        Err(crate::BoolError::NonManifoldResult { .. })
    ));
    let plain = crate::BoolOptions::default();
    type Want = fn(&BodySet) -> bool;
    let cases: [(SetExpr, crate::BoolOptions, Want); 6] = [
        (SetExpr::Any, plain, |m| !m.is_empty()),
        (SetExpr::All, plain, |m| m.len() == 3),
        (SetExpr::AtLeast(2), plain, |m| m.len() >= 2),
        (SetExpr::Odd, split, |m| m.len() % 2 == 1),
        (
            SetExpr::body(0)
                .and(SetExpr::body(1))
                .minus(SetExpr::body(2)),
            plain,
            |m| m.contains(0) && m.contains(1) && !m.contains(2),
        ),
        (SetExpr::body(1).minus(SetExpr::body(0)), plain, |m| {
            m.contains(1) && !m.contains(0)
        }),
    ];
    for (expr, options, want) in cases {
        let body = p
            .select(&expr, &options)
            .unwrap_or_else(|e| panic!("{expr:?}: {e}"));
        assert!(
            body.validate_all().is_empty(),
            "{expr:?}: {:?}",
            body.validate_all()
        );
        agrees_on_lattice(&bodies, &want, &body);
    }
}

/// Every region of three overlapping cubes as its own valid solid, and
/// every shared face an interface: the two faces coincide (an interior
/// point of one lies on the other) and face opposite ways.
#[test]
fn three_cubes_give_region_solids_and_their_interfaces() {
    use popkern_core::geometry::surface::SurfaceEval;
    let a = cube_at(0.0, 0.0, 0.0);
    let b = cube_at(0.5, 0.3, 0.2);
    let c = cube_at(0.25, 0.6, 0.45);
    let bodies = [&a, &b, &c];
    let p = partition(&bodies, Tolerance::default()).expect("partition");
    let cells = p
        .region_solids(&super::SetExpr::Any)
        .expect("region solids");
    assert_eq!(cells.solids.len(), 7);
    for s in &cells.solids {
        assert!(
            s.body.validate_all().is_empty(),
            "{:?}: {:?}",
            s.bodies,
            s.body.validate_all()
        );
        let inside = p.regions()[s.region].interior_point;
        assert_eq!(classified(&bodies, &inside), s.bodies);
    }
    assert!(!cells.interfaces.is_empty());
    for (k, _) in cells.solids.iter().enumerate() {
        assert!(
            cells
                .interfaces
                .iter()
                .any(|f| f.solids.0 == k || f.solids.1 == k),
            "region solid {k} shares no face"
        );
    }
    let normal_at = |body: &BRep, face, p: &popkern_core::geometry::Point3| {
        let f = &body.faces[face];
        let s = f.surface.as_ref().expect("surface");
        let (u, v) = s.geometry.closest_point(p).expect("uv");
        let n = s.geometry.normal_at(u, v);
        if s.orientation { n } else { n * -1.0 }
    };
    for f in &cells.interfaces {
        let (si, sj) = (
            &cells.solids[f.solids.0].body,
            &cells.solids[f.solids.1].body,
        );
        let q = crate::operation::contact::face_interior_point(si, f.faces.0).expect("interior");
        let (u, v) = sj.faces[f.faces.1]
            .surface
            .as_ref()
            .expect("surface")
            .geometry
            .closest_point(&q)
            .expect("uv");
        assert_eq!(
            crate::classify_uv_loops::point_in_face_region_exact(sj, f.faces.1, u, v),
            Some(true),
            "{f:?}"
        );
        assert!(
            normal_at(si, f.faces.0, &q).dot(&normal_at(sj, f.faces.1, &q)) < -0.99,
            "{f:?}"
        );
    }
}

fn block(lo: [f64; 3], hi: [f64; 3]) -> BRep {
    use popkern_core::geometry::{Point3, Vec3};
    BRep::prism(
        &[
            Point3::new(lo[0], lo[1], lo[2]),
            Point3::new(hi[0], lo[1], lo[2]),
            Point3::new(hi[0], hi[1], lo[2]),
            Point3::new(lo[0], hi[1], lo[2]),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        hi[2] - lo[2],
        Tolerance::default(),
    )
    .expect("block")
}

/// A slab through the middle of a block cuts it in three: the two
/// outer pieces outside the slab, the middle one inside it, sharing a
/// face with each.
#[test]
fn a_slab_splits_a_block_into_three_pieces() {
    let target = block([0.0; 3], [2.0; 3]);
    let slab = block([0.8, -1.0, -1.0], [1.2, 3.0, 3.0]);
    let cut = super::split(&target, &[&slab], Tolerance::default()).expect("split");
    assert_eq!(cut.pieces.len(), 3);
    let inside: Vec<bool> = cut
        .pieces
        .iter()
        .map(|p| p.inside_tools.contains(0))
        .collect();
    assert_eq!(inside.iter().filter(|&&i| i).count(), 1);
    for p in &cut.pieces {
        assert!(
            p.body.validate_all().is_empty(),
            "{:?}",
            p.body.validate_all()
        );
    }
    assert_eq!(cut.interfaces.len(), 2);
    let middle = inside.iter().position(|&i| i).expect("middle");
    assert!(
        cut.interfaces
            .iter()
            .all(|f| f.solids.0 == middle || f.solids.1 == middle)
    );
}

#[test]
fn a_tool_that_misses_leaves_one_piece() {
    let target = block([0.0; 3], [1.0; 3]);
    let far = block([5.0; 3], [6.0; 3]);
    let cut = super::split(&target, &[&far], Tolerance::default()).expect("split");
    assert_eq!(cut.pieces.len(), 1);
    assert!(cut.pieces[0].inside_tools.is_empty());
    assert_eq!(cut.pieces[0].body.faces.len(), 6);
    assert!(cut.interfaces.is_empty());
}

/// A block resting on a plate shares part of the plate's top: two
/// regions, their union one valid body, and as region solids one
/// interface, the contact patch.
#[test]
fn a_block_resting_on_a_plate_partitions_with_its_contact() {
    let plate = block([0.0, 0.0, 0.0], [4.0, 4.0, 1.0]);
    let resting = block([1.0, 1.0, 1.0], [2.0, 2.0, 2.0]);
    let bodies = [&plate, &resting];
    check(&bodies, &[&[0], &[1]]);
    let p = partition(&bodies, Tolerance::default()).expect("partition");
    let union = p
        .select(&super::SetExpr::Any, &crate::BoolOptions::default())
        .expect("union");
    assert!(
        union.validate_all().is_empty(),
        "{:?}",
        union.validate_all()
    );
    let cells = p
        .region_solids(&super::SetExpr::Any)
        .expect("region solids");
    assert_eq!(cells.solids.len(), 2);
    assert_eq!(cells.interfaces.len(), 1);
}

/// Identical cubes are one region both hold.
#[test]
fn identical_cubes_are_one_region() {
    let (a, b) = (cube_at(0.0, 0.0, 0.0), cube_at(0.0, 0.0, 0.0));
    check(&[&a, &b], &[&[0, 1]]);
}
