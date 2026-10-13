use crate::error::ValidationError;
use crate::geometry::Point3;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

fn pcurve_errors(brep: &BRep) -> Vec<ValidationError> {
    brep.validate_pcurve_consistency_all()
}

/// Primitives on every chart kind carry consistent pcurves.
#[test]
fn primitives_have_consistent_pcurves() {
    let tol = Tolerance::default();
    for (name, body) in [
        ("cube", BRep::unit_cube(tol)),
        ("cylinder", BRep::cylinder(0.5, 2.0, tol)),
        ("cone", BRep::cone(1.0, 2.0, tol)),
        ("sphere", BRep::sphere(1.0, tol)),
        ("torus", BRep::torus(2.0, 0.5, tol)),
    ] {
        let body = body.expect(name);
        assert!(
            pcurve_errors(&body).is_empty(),
            "{name}: {:?}",
            pcurve_errors(&body)
        );
    }
}

/// Negative control: an open edge's pcurve reversed against its
/// half-edge (its trace still coincides) is reported.
#[test]
fn a_reversed_open_pcurve_is_reported() {
    let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube");
    let heh = brep.half_edges.keys().next().expect("half-edge");
    let pcurve = brep.half_edges[heh].pcurve.as_mut().expect("pcurve");
    let nurbs = pcurve.as_nurbs_mut().expect("cube pcurves are NURBS");
    nurbs.control_points.reverse();
    nurbs.weights.reverse();
    let errors = pcurve_errors(&brep);
    assert!(
        matches!(errors.as_slice(), [ValidationError::PcurveDeviation { half_edge, .. }] if *half_edge == heh),
        "{errors:?}"
    );
}

/// Negative control: a closed edge's pcurve moved off its circle (an
/// offset that kept the source pcurve) is reported.
#[test]
fn a_closed_pcurve_off_its_edge_is_reported() {
    let mut brep = BRep::cylinder(0.5, 2.0, Tolerance::default()).expect("cylinder");
    let (heh, shifted) = brep
        .half_edges
        .iter()
        .find_map(|(h, he)| {
            let pc = he.pcurve.as_ref()?;
            let (a, b) = (pc.point_at(pc.t_start), pc.point_at(pc.t_end));
            ((b.x - a.x).abs() > 3.0).then(|| {
                // A cylinder rim: slide it 0.1 along the axis.
                let n = NurbsCurve {
                    degree: 1,
                    knots: vec![0.0, 0.0, 1.0, 1.0],
                    control_points: vec![
                        Point3::new(a.x, a.y + 0.1, 0.0),
                        Point3::new(b.x, b.y + 0.1, 0.0),
                    ],
                    weights: vec![1.0, 1.0],
                };
                (h, n)
            })
        })
        .expect("a cylinder rim pcurve");
    assert!(pcurve_errors(&brep).is_empty());
    brep.half_edges[heh].pcurve = Some(shifted.into());
    let errors = pcurve_errors(&brep);
    assert!(
        errors.iter().any(
            |e| matches!(e, ValidationError::PcurveDeviation { half_edge, .. } if *half_edge == heh)
        ),
        "{errors:?}"
    );
}
