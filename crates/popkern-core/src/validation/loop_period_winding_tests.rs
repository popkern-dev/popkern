use super::*;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::{Point3, Transform3};
use crate::tolerance::Tolerance;

fn winding_errors(brep: &BRep) -> Vec<ValidationError> {
    brep.validate_loop_period_winding_all()
}

/// Valid bodies on every periodic chart kind read clean: cylinder,
/// sphere (pole caps), cone (apex), torus (doubly periodic).
#[test]
fn periodic_primitives_bound_regions() {
    let tol = Tolerance::default();
    let bodies = [
        ("cylinder", BRep::cylinder(0.5, 2.0, tol)),
        ("sphere", BRep::sphere(1.0, tol)),
        ("cone", BRep::cone(1.0, 2.0, tol)),
        ("torus", BRep::torus(2.0, 0.5, tol)),
        (
            "moved cylinder",
            BRep::cylinder(0.5, 2.0, tol)
                .and_then(|b| b.transform(&Transform3::translation(3.0, -1.0, 0.5))),
        ),
    ];
    for (name, body) in bodies {
        let body = body.expect(name);
        assert!(
            winding_errors(&body).is_empty(),
            "{name}: {:?}",
            winding_errors(&body)
        );
    }
}

/// Negative control: a cylinder side's rim pcurve stretched the long
/// way round (twice the period) leaves the face a net winding.
#[test]
fn a_rim_run_the_long_way_round_is_reported() {
    let tol = Tolerance::default();
    let mut body = BRep::cylinder(0.5, 2.0, tol).expect("cylinder");
    let side = body
        .faces
        .iter()
        .find(|(_, f)| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Cylinder { .. })
            )
        })
        .map(|(h, _)| h)
        .expect("side face");
    assert!(winding_errors(&body).is_empty());
    // The side's rim half-edge: the pcurve whose u sweeps a full turn.
    let start = body.loops[body.faces[side].outer_loop].half_edge;
    let mut current = start;
    let rim = loop {
        let he = &body.half_edges[current];
        let pc = he.pcurve.as_ref().expect("pcurve");
        let (a, b) = (pc.point_at(pc.t_start), pc.point_at(pc.t_end));
        if (b.x - a.x).abs() > 3.0 {
            break current;
        }
        current = he.next;
        assert_ne!(current, start, "a rim pcurve");
    };
    let pc = body.half_edges[rim].pcurve.as_ref().expect("pcurve");
    let (a, b) = (pc.point_at(pc.t_start), pc.point_at(pc.t_end));
    let long = Point3::new(a.x + 2.0 * (b.x - a.x), b.y, 0.0);
    body.half_edges[rim].pcurve = Some(
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(a.x, a.y, 0.0), long],
            weights: vec![1.0, 1.0],
        }
        .into(),
    );
    let errors = winding_errors(&body);
    assert!(
        matches!(
            errors.as_slice(),
            [ValidationError::LoopPeriodWinding { face, winding_u, winding_v: 0 }]
                if *face == side && winding_u.abs() == 1
        ),
        "{errors:?}"
    );
}
