use popkern_core::geometry::Point3;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// A cone and the same cone turned about its axis (its seam elsewhere):
/// the containment probe's outward side came from the outer loop's vertex
/// polygon, which on the cone's lateral face is its rim vertex and its
/// apex alone, so the shortcut never applied and the full pipeline
/// declined the coincident apex (loop successor; 32 scenes of the
/// generated corpus's coincident cone cell).
#[test]
fn a_cone_and_itself_turned_about_its_axis_unite_and_meet_in_the_cone() {
    for size in [1.0, 1e-3] {
        cone_and_turned_copy(size);
    }
}

/// The cone of base radius `0.8 size` and its turned copy (a thousandth
/// size met in no record, and the no-records path declined).
fn cone_and_turned_copy(size: f64) {
    let tol = Tolerance::default();
    let a = BRep::cone(0.8 * size, 1.6 * size, tol).expect("cone");
    let turn = Transform3::rotation_around_axis(
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
        2.94,
    )
    .expect("turn");
    let b = a.transform(&turn).expect("turned");
    for (name, result) in [
        ("union", crate::union(a.clone(), b.clone(), tol, None)),
        (
            "common",
            crate::intersection(a.clone(), b.clone(), tol, None),
        ),
    ] {
        let body = result.unwrap_or_else(|e| panic!("{size} {name}: {e}"));
        let (v, e, f) = (body.vertices.len(), body.edges.len(), body.faces.len());
        assert_eq!(
            (v, e, f),
            (2, 2, 2),
            "{size} {name}: the cone's own topology"
        );
    }
    let rest = crate::subtraction(a, b, tol, None);
    assert!(
        rest.as_ref().map_or(true, |body| body.faces.is_empty()),
        "difference: nothing left, got {:?}",
        rest.map(|body| body.faces.len())
    );
}
