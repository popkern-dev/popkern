use super::*;
use crate::tolerance::Tolerance;

#[test]
fn returns_none_for_stale_half_edge() {
    let brep = BRep::new(Tolerance::default());
    let stale = HalfEdgeHandle::default();
    assert!(surface_derived_split_position(&brep, stale, 0.5).is_none());
}

#[test]
fn returns_none_when_half_edge_has_no_pcurve() {
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    // mev edges don't have pcurves attached, so helper returns None.
    assert!(surface_derived_split_position(&brep, mev.half_edge, 0.5).is_none());
}

#[test]
fn returns_some_on_planar_face_equal_to_chord_midpoint() {
    // Cube faces are planar — the surface-derived split position
    // coincides with the chord midpoint within fp precision.
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let half_edge = brep.half_edges.keys().next().unwrap();
    let result = surface_derived_split_position(&brep, half_edge, 0.5);
    let Some(point) = result else {
        panic!("expected Some on a cube half-edge (Plane surface + pcurve)");
    };
    // Compute the chord midpoint for comparison.
    let he = &brep.half_edges[half_edge];
    let next = &brep.half_edges[he.next];
    let pos_a = brep.vertices[he.origin].position;
    let pos_b = brep.vertices[next.origin].position;
    let chord_mid = Point3::new(
        0.5 * (pos_a.x + pos_b.x),
        0.5 * (pos_a.y + pos_b.y),
        0.5 * (pos_a.z + pos_b.z),
    );
    let drift = (point - chord_mid).length();
    assert!(
        drift < 1.0e-12,
        "surface-derived point {point:?} should equal chord midpoint {chord_mid:?} on a planar face; drift={drift:.3e}"
    );
}

#[test]
fn returns_some_at_endpoint_t_equals_endpoint_position() {
    // t = 0 should return the half-edge's origin vertex position
    // (within fp precision); a basic sanity check on the pcurve
    // parameter mapping.
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let half_edge = brep.half_edges.keys().next().unwrap();
    let Some(point) = surface_derived_split_position(&brep, half_edge, 0.0) else {
        panic!("expected Some at t=0");
    };
    let origin = brep.vertices[brep.half_edges[half_edge].origin].position;
    let drift = (point - origin).length();
    assert!(
        drift < 1.0e-12,
        "surface-derived point at t=0 should equal origin vertex; drift={drift:.3e}"
    );
}

#[test]
fn returns_some_on_curved_surface_equals_surface_point_at_pcurve_uv() {
    // On a cylinder face the chord midpoint lies INSIDE the
    // cylindrical surface (below the curve), so chord and
    // surface-derived diverge. We verify the surface-derived value
    // is on the cylinder.
    let brep = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    // Find a half-edge on the lateral (cylindrical) face.
    let half_edge = brep
        .half_edges
        .iter()
        .find_map(|(handle, data)| {
            let face_handle = brep.loops.get(data.loop_)?.face;
            let face = brep.faces.get(face_handle)?;
            let surface = face.surface.as_ref()?;
            data.pcurve.as_ref()?;
            matches!(
                surface.geometry,
                crate::geometry::surface::SurfaceGeometry::Cylinder { .. }
            )
            .then_some(handle)
        })
        .expect("cylinder should expose at least one lateral half-edge with pcurve");
    let Some(point) = surface_derived_split_position(&brep, half_edge, 0.5) else {
        panic!("expected Some on cylinder lateral half-edge");
    };
    // Confirm the result is on the cylinder (distance from axis ≈ radius).
    let axis_dist = (point.x.powi(2) + point.y.powi(2)).sqrt();
    assert!(
        (axis_dist - 1.0).abs() < 1.0e-6,
        "surface-derived point {point:?} should be on cylinder of radius 1.0; axis_dist={axis_dist:.6e}"
    );
}
