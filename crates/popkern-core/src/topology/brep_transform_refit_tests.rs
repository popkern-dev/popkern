use crate::geometry::transform::Transform3;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

/// A cube whose one vertex sits `miss` off its edges' curve ends.
fn cube_with_vertex_off(miss: f64) -> (BRep, crate::topology::VertexHandle) {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let v = brep.vertices.keys().next().unwrap();
    let p = brep.vertices[v].position;
    brep.vertices[v].position = crate::geometry::Point3::new(p.x, p.y, p.z + miss);
    (brep, v)
}

/// A curve end 1.2e-9 off its vertex at unit size is 1.2e-6 off at a
/// thousand, past the bare tolerance: the scale records it on the vertex
/// and its edges (generated corpus W at scale 1000, a box split by a
/// NURBS dome, failed closure before any boolean ran).
#[test]
fn a_scaled_curve_end_off_its_vertex_widens_the_vertex() {
    let (brep, v) = cube_with_vertex_off(1.2e-9);
    assert!(brep.validate_all().is_empty(), "valid at unit size");
    let scaled = brep.transform(&Transform3::scale_uniform(1000.0)).unwrap();
    let issues = scaled.validate_all();
    assert!(issues.is_empty(), "scaled body invalid: {issues:?}");
    let tol = scaled.vertex_tolerance(v);
    assert!(tol >= 1.2e-6 * (1.0 - 1e-6), "vertex tolerance {tol}");
    // The edges at the vertex hold it; the vertices away from it keep
    // the bare tolerance.
    let far = scaled
        .vertices
        .keys()
        .filter(|&w| w != v)
        .all(|w| scaled.vertex_tolerance(w) <= Tolerance::default().linear());
    assert!(far, "only the vertex off its curves widens");
}

/// Negative control: the scaled body's curve end really lies past the
/// bare tolerance, so a body left at it fails closure.
#[test]
fn the_scaled_miss_is_past_the_bare_tolerance() {
    let (brep, v) = cube_with_vertex_off(1.2e-9);
    let mut scaled = brep.transform(&Transform3::scale_uniform(1000.0)).unwrap();
    let bare = Tolerance::default().linear();
    scaled.vertices[v].tolerance = Some(bare);
    for e in scaled.edges.keys().collect::<Vec<_>>() {
        scaled.edges[e].tolerance = None;
    }
    let issues = scaled.validate_all();
    assert!(
        issues.iter().any(|i| i.to_string().contains("Shell gap")),
        "expected a shell gap at the bare tolerance: {issues:?}"
    );
}
