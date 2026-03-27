use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceEval;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

#[test]
fn recompute_pcurves_on_cube_is_idempotent() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // All half-edges already have pcurves from unit_cube.
    brep.recompute_all_pcurves().unwrap();
    for (_, he) in &brep.half_edges {
        assert!(he.pcurve.is_some(), "pcurve should still be present");
    }
}

#[test]
fn recompute_fills_missing_pcurves() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // Strip all pcurves.
    let he_keys: Vec<_> = brep.half_edges.keys().collect();
    for heh in &he_keys {
        brep.half_edges[*heh].pcurve = None;
    }
    // Recompute should fill them back.
    brep.recompute_all_pcurves().unwrap();
    for (_, he) in &brep.half_edges {
        assert!(he.pcurve.is_some(), "pcurve should be recomputed");
    }
}

#[test]
fn pcurve_endpoints_match_vertices() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // Strip and recompute to exercise the fix path.
    let he_keys: Vec<_> = brep.half_edges.keys().collect();
    for heh in &he_keys {
        brep.half_edges[*heh].pcurve = None;
    }
    brep.recompute_all_pcurves().unwrap();

    for fh in brep.faces.keys().collect::<Vec<_>>() {
        let face = &brep.faces[fh];
        let surface = match &face.surface {
            Some(s) => &s.geometry,
            None => continue,
        };
        let loop_ = &brep.loops[face.outer_loop];
        let mut cursor = loop_.half_edge;
        loop {
            let he = &brep.half_edges[cursor];
            if let Some(ref pcurve) = he.pcurve {
                let (t0, t1) = pcurve.parameter_range();
                // Check start
                let uv_start = pcurve.point_at(t0);
                let pos_from_uv = surface.point_at(uv_start.x, uv_start.y);
                let actual_start = brep.vertices[he.origin].position;
                let err_start = pos_from_uv.distance_to(&actual_start);
                assert!(err_start < 1e-6, "pcurve start error: {err_start}");
                // Check end
                let uv_end = pcurve.point_at(t1);
                let pos_from_uv_end = surface.point_at(uv_end.x, uv_end.y);
                let next_he = &brep.half_edges[he.next];
                let actual_end = brep.vertices[next_he.origin].position;
                let err_end = pos_from_uv_end.distance_to(&actual_end);
                assert!(err_end < 1e-6, "pcurve end error: {err_end}");
            }
            cursor = brep.half_edges[cursor].next;
            if cursor == loop_.half_edge {
                break;
            }
        }
    }
}

#[test]
fn cylinder_seam_pcurve_is_degree_1() {
    use crate::geometry::surface::SurfaceGeometry;

    let mut brep = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    // Clear all pcurves and recompute.
    for heh in brep.half_edges.keys().collect::<Vec<_>>() {
        brep.half_edges[heh].pcurve = None;
    }
    brep.recompute_all_pcurves().unwrap();

    // The seam edges on the lateral (cylindrical) face should get degree-1
    // pcurves since their UV path is a straight vertical line.
    let lateral_face = brep.faces.keys().find(|&fh| {
        matches!(
            brep.faces[fh].surface.as_ref().map(|s| &s.geometry),
            Some(SurfaceGeometry::Cylinder { .. })
        )
    });
    if let Some(fh) = lateral_face {
        let loop_ = &brep.loops[brep.faces[fh].outer_loop];
        let mut cursor = loop_.half_edge;
        let mut found_degree_1 = false;
        loop {
            let he = &brep.half_edges[cursor];
            if let Some(ref pcurve) = he.pcurve {
                if pcurve.degree == 1 && pcurve.control_points.len() == 2 {
                    found_degree_1 = true;
                }
            }
            cursor = he.next;
            if cursor == loop_.half_edge {
                break;
            }
        }
        assert!(
            found_degree_1,
            "cylinder lateral face should have at least one degree-1 pcurve (seam)"
        );
    }
}

#[test]
fn cylinder_recomputed_pcurves_are_accurate() {
    use crate::geometry::surface::SurfaceGeometry;

    let mut brep = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    // Clear all pcurves and recompute.
    for heh in brep.half_edges.keys().collect::<Vec<_>>() {
        brep.half_edges[heh].pcurve = None;
    }
    brep.recompute_all_pcurves().unwrap();

    // Check that pcurve-evaluated 3D positions match vertex positions on
    // the cylindrical face.
    let lateral_face = brep.faces.keys().find(|&fh| {
        matches!(
            brep.faces[fh].surface.as_ref().map(|s| &s.geometry),
            Some(SurfaceGeometry::Cylinder { .. })
        )
    });
    if let Some(fh) = lateral_face {
        let surface = &brep.faces[fh].surface.as_ref().unwrap().geometry;
        let loop_ = &brep.loops[brep.faces[fh].outer_loop];
        let mut cursor = loop_.half_edge;
        loop {
            let he = &brep.half_edges[cursor];
            if let Some(ref pcurve) = he.pcurve {
                let (t0, t1) = pcurve.parameter_range();
                let uv_start = pcurve.point_at(t0);
                let pos_start = surface.point_at(uv_start.x, uv_start.y);
                let actual_start = brep.vertices[he.origin].position;
                let err = pos_start.distance_to(&actual_start);
                assert!(err < 1e-6, "cylinder pcurve start error: {err}");

                let uv_end = pcurve.point_at(t1);
                let pos_end = surface.point_at(uv_end.x, uv_end.y);
                let next_he = &brep.half_edges[he.next];
                let actual_end = brep.vertices[next_he.origin].position;
                let err_end = pos_end.distance_to(&actual_end);
                assert!(err_end < 1e-6, "cylinder pcurve end error: {err_end}");
            }
            cursor = brep.half_edges[cursor].next;
            if cursor == loop_.half_edge {
                break;
            }
        }
    }
}
