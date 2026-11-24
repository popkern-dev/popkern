//! A rod lying along a block's top, touching exactly and overlapping by
//! 5e-7 (within the tolerance): the operands meet in one tangent record,
//! a contact curve (W3c). The union declines on the line under the
//! default policy and serves one general body under the general one.
//! Before, the record entered the arrangement as a crossing: the split
//! pieces kept their copies of it and the general glue declined (a wedge
//! between the fins in two regions; the corpus's PY-K cells).

use popkern_core::geometry::{Point3, PortableMath, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::{BoolError, BoolOperation, BoolOptions, NonManifoldPolicy, boolean_with_options};

fn block_and_rod(gap: f64) -> (BRep, BRep) {
    let tol = Tolerance::default();
    let rect: Vec<Point3> = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]
        .iter()
        .map(|&(x, y)| Point3::new(x, y, 0.0))
        .collect();
    let block = BRep::prism(&rect, Vec3::new(0.0, 0.0, 1.0), 2.0, tol).expect("block");
    let rod = BRep::cylinder(0.5, 3.0, tol)
        .and_then(|c| c.transform(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2)))
        .and_then(|c| c.transform(&Transform3::translation(0.5, 2.0, 2.5 + gap)))
        .expect("rod");
    (block, rod)
}

#[test]
fn a_rod_lying_along_a_top_touches_along_a_contact_curve() {
    let tol = Tolerance::default();
    for gap in [0.0, -5e-7] {
        let (block, rod) = block_and_rod(gap);
        let decline = BoolOptions {
            non_manifold: NonManifoldPolicy::Decline,
        };
        match boolean_with_options(
            BoolOperation::Union,
            block.clone(),
            rod.clone(),
            tol,
            None,
            &decline,
        ) {
            Err(BoolError::NonManifoldResult { contacts }) => {
                assert!(
                    contacts
                        .iter()
                        .all(|c| (c.location.z - 2.0).abs() < 1e-5
                            && (c.location.y - 2.0).abs() < 1e-3),
                    "gap {gap}: {contacts:?}"
                );
            }
            other => panic!("gap {gap}: {:?}", other.map(|b| b.faces.len())),
        }
        let general = BoolOptions {
            non_manifold: NonManifoldPolicy::General,
        };
        let body = boolean_with_options(BoolOperation::Union, block, rod, tol, None, &general)
            .unwrap_or_else(|e| panic!("gap {gap}: {e}"));
        assert!(
            body.validate_all().is_empty(),
            "gap {gap}: {:?}",
            body.validate_all().first()
        );
        assert_eq!(body.solid_regions().len(), 2, "gap {gap}");
    }
}

/// The rod turned about its axis, its seam off the contact line: the
/// analytic arm is tangency-silent and no edge lies on the line, so the
/// operands meet in no record. The face-pair pass traces the
/// tangent-curve contact into one exact line along the top
/// (`ff_parallel_contact_curves`), and the union declines on it, read from
/// the topology. Before, construction found no curve: the general glue
/// had only the rims' touches.
#[test]
fn a_turned_rod_on_a_top_is_traced_into_a_contact_line() {
    let tol = Tolerance::default();
    let rect: Vec<Point3> = [(0.0, 0.0), (4.0, 0.0), (4.0, 4.0), (0.0, 4.0)]
        .iter()
        .map(|&(x, y)| Point3::new(x, y, 0.0))
        .collect();
    let block = BRep::prism(&rect, Vec3::new(0.0, 0.0, 1.0), 2.0, tol).expect("block");
    let rod = BRep::cylinder(0.5, 3.0, tol)
        .and_then(|c| c.transform(&Transform3::rotation_z(0.7)))
        .and_then(|c| c.transform(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2)))
        .and_then(|c| c.transform(&Transform3::translation(0.5, 2.0, 2.5)))
        .expect("rod");
    let ix = crate::intersect::find_intersections(&block, &rod, &tol).expect("intersect");
    assert!(ix.records.is_empty(), "{} records", ix.records.len());
    assert_eq!(ix.contact_curves.len(), 1, "{:?}", ix.contact_curves);
    let c = &ix.contact_curves[0];
    use popkern_core::geometry::curve::CurveEval;
    for t in [c.t_start, 0.5 * (c.t_start + c.t_end), c.t_end] {
        let p = c.geometry.point_at(t);
        assert!(
            (p.y - 2.0).abs() < 1e-6 && (p.z - 2.0).abs() < 1e-9,
            "{p:?}"
        );
    }
    match crate::boolean(BoolOperation::Union, block, rod, tol, None) {
        Err(BoolError::NonManifoldResult { contacts }) => {
            assert!(
                contacts.iter().any(|k| (k.location.z - 2.0).abs() < 1e-6),
                "{contacts:?}"
            );
        }
        other => panic!("{:?}", other.map(|b| b.faces.len())),
    }
}

/// Two balls of radius 1000 overlapping by two tolerances cross at
/// 9e-5 rad, under the angular tolerance, round a lens two tolerances
/// thick: a crossing, never a contact curve, and the cell complex's
/// radial order alternates the balls' faces at its rim. Read by the
/// bare angle, the record went to the contact reading and then the
/// order tied the faces by curvature, and the common was served empty
/// (SS-E 695). Off both balls' seams and poles.
#[test]
fn balls_a_thousand_across_two_tolerances_deep_meet_in_a_lens() {
    let tol = Tolerance::default();
    let r = 1000.0;
    let d = Vec3::new(1.0_f64.pcos(), 1.0_f64.psin(), 0.3)
        .normalized()
        .expect("dir");
    let reach = 2.0 * r - 2.0 * tol.linear();
    let a = BRep::sphere(r, tol).expect("a");
    let b = BRep::sphere(r, tol)
        .and_then(|s| {
            s.transform(&Transform3::translation(
                d.x * reach,
                d.y * reach,
                d.z * reach,
            ))
        })
        .expect("b");
    let common = boolean_with_options(
        BoolOperation::Intersection,
        a,
        b,
        tol,
        None,
        &BoolOptions::default(),
    )
    .expect("common");
    assert_eq!(common.faces.len(), 2, "the lens: {}", common.faces.len());
    common.validate().expect("valid");
}
