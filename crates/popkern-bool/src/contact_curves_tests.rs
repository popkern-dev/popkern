//! A rod lying along a block's top, touching exactly and overlapping by
//! 5e-7 (within the tolerance): the operands meet in one tangent record,
//! a contact curve (W3c). The union declines on the line under the
//! default policy and serves one general body under the general one.
//! Before, the record entered the arrangement as a crossing: the split
//! pieces kept their copies of it and the general glue declined (a wedge
//! between the fins in two regions; the corpus's PY-K cells).

use popkern_core::geometry::{Point3, Transform3, Vec3};
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
