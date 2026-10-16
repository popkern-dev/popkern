//! The served-body loop winding gate (PK-862).

use super::gate_loop_period_winding;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::{Point3, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn all_winding(body: &BRep) -> usize {
    let faces: Vec<_> = body.faces.keys().collect();
    body.validate_loop_period_winding_among(&faces).len()
}

/// A cylinder side whose rim pcurve runs the long way round the period
/// (a net winding, the cooriented-coincident-face signature) is
/// re-derived by the gate and the body reads clean; a clean body is
/// left as it is.
#[test]
fn a_rim_run_the_long_way_round_is_rederived() {
    let tol = Tolerance::default();
    let clean = BRep::cylinder(0.5, 2.0, tol).expect("cylinder");
    let mut untouched = clean.clone();
    assert!(gate_loop_period_winding(&mut untouched).is_ok());
    assert_eq!(
        format!("{:?}", untouched.half_edges),
        format!("{:?}", clean.half_edges)
    );

    let mut body = clean;
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
    body.half_edges[rim].pcurve = Some(
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(a.x, a.y, 0.0),
                Point3::new(a.x + 2.0 * (b.x - a.x), b.y, 0.0),
            ],
            weights: vec![1.0, 1.0],
        }
        .into(),
    );
    assert_eq!(all_winding(&body), 1, "the corrupted rim winds");
    assert!(gate_loop_period_winding(&mut body).is_ok());
    assert_eq!(all_winding(&body), 0, "re-derived");
    assert!(body.validate().is_ok());
}
