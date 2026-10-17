//! The split point refuses a chord whose copies miss a face use of its
//! edge (PK-873): landing it on fewer faces opens the weld.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::{CurveSamples, edge_curve_midpoint, face_uses};

#[test]
fn a_chord_splits_only_with_a_copy_in_every_face_of_its_edge() {
    let brep = BRep::cylinder(1.0, 2.0, Tolerance::default()).expect("cylinder");
    // A rim: an edge with a curve whose two half-edges lie in different
    // faces.
    let (edge, faces) = brep
        .edges
        .iter()
        .find_map(|(eh, e)| {
            let h = &brep.half_edges[e.half_edge];
            let t = &brep.half_edges[h.twin];
            let (f0, f1) = (brep.loops[h.loop_].face, brep.loops[t.loop_].face);
            (e.curve.is_some() && f0 != f1).then_some((eh, [f0, f1]))
        })
        .expect("a rim edge");
    assert_eq!(face_uses(&brep, edge), 2);
    let curve = brep.edges[edge].curve.as_ref().expect("curve");
    let at = |f: f64| {
        curve
            .geometry
            .point_at(curve.t_start + (curve.t_end - curve.t_start) * f)
    };
    let (a, b) = (at(0.1), at(0.2));
    let gate = 1e-7;
    let mut cache = CurveSamples::new();
    for one in faces {
        assert!(
            edge_curve_midpoint(&brep, &[one], a, b, gate, &mut cache, None).is_none(),
            "a split with one face's copy landed"
        );
    }
    let m = edge_curve_midpoint(&brep, &faces, a, b, gate, &mut cache, None)
        .expect("a split with both copies");
    assert!(m.distance_to(&at(0.15)) < 1e-9, "midpoint {m:?}");
}
