//! Tests for the variable-radius in-place path (Round 5h + 5i).
//! Split out so neither `modify_inplace_tests.rs` nor the
//! production module exceeds the 500-line structural cap.

use crate::types::{FilletOptions, FilletRadius, FilletSpec};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube() -> BRep {
    BRep::unit_cube(tol()).unwrap()
}

/// Round 5i: variable in-place reassigns the fillet face's surface
/// to a true NURBS sweep (not a Torus with averaged minor radius),
/// confirming geometric exactness.
#[test]
fn modify_open_in_place_variable_assigns_nurbs_sweep_surface() {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let brep = cube();
    let edge = brep.edges.keys().next().unwrap();
    let opts = FilletOptions {
        propagate_tangent: false,
        vertex_mode: crate::types::VertexMode::Mitred,
    };
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let r1 = crate::fillet_edges(&brep, &specs, &opts, tol(), None).unwrap();
    let blend = r1
        .brep
        .faces
        .keys()
        .find(|fh| {
            matches!(
                r1.brep.faces[*fh].surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Torus { .. })
            )
        })
        .expect("post-fillet brep should contain a Torus face");

    let modified =
        crate::modify_inplace::modify_open_in_place_variable(&r1.brep, blend, 0.05, 0.25)
            .expect("variable in-place should succeed");

    let surf = modified.faces[blend]
        .surface
        .as_ref()
        .expect("modified blend should have a surface");
    match &surf.geometry {
        SurfaceGeometry::Nurbs(_) => { /* expected */ }
        other => panic!("variable in-place should assign a NURBS sweep, got {other:?}"),
    }
}

/// Round 5h: `modify_open_in_place_variable` moves corners using
/// per-end radii. The start corners must be at distance `r_start`
/// from the start endpoint and the end corners at distance `r_end`
/// from the end endpoint.
#[test]
fn modify_open_in_place_variable_uses_per_end_radii() {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let brep = cube();
    let edge = brep.edges.keys().next().unwrap();
    let opts = FilletOptions {
        propagate_tangent: false,
        vertex_mode: crate::types::VertexMode::Mitred,
    };
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let r1 = crate::fillet_edges(&brep, &specs, &opts, tol(), None).unwrap();
    let blend = r1
        .brep
        .faces
        .keys()
        .find(|fh| {
            matches!(
                r1.brep.faces[*fh].surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Torus { .. })
            )
        })
        .expect("post-fillet brep should contain a Torus face");

    let r_start = 0.05;
    let r_end = 0.25;
    let modified =
        crate::modify_inplace::modify_open_in_place_variable(&r1.brep, blend, r_start, r_end)
            .expect("variable in-place should succeed for plane-plane fillet");

    let plan = crate::unfillet::compute_unfillet_plan(&r1.brep, blend).unwrap();
    let p_start = plan.restored_edge_start;
    let p_end = plan.restored_edge_end;

    let mut start_distances: Vec<f64> = Vec::new();
    let mut end_distances: Vec<f64> = Vec::new();
    for (_, v) in &modified.vertices {
        let pos = v.position;
        let d_start = ((pos.x - p_start.x).powi(2)
            + (pos.y - p_start.y).powi(2)
            + (pos.z - p_start.z).powi(2))
        .sqrt();
        let d_end =
            ((pos.x - p_end.x).powi(2) + (pos.y - p_end.y).powi(2) + (pos.z - p_end.z).powi(2))
                .sqrt();
        if d_start < d_end && d_start > 1e-9 {
            start_distances.push(d_start);
        } else if d_end < d_start && d_end > 1e-9 {
            end_distances.push(d_end);
        }
    }

    let near_start_at_radius = start_distances
        .iter()
        .filter(|d| (**d - r_start).abs() < 1e-6)
        .count();
    let near_end_at_radius = end_distances
        .iter()
        .filter(|d| (**d - r_end).abs() < 1e-6)
        .count();
    assert!(
        near_start_at_radius >= 2,
        "should have ≥2 corners at distance {r_start} from start; got {near_start_at_radius} (distances: {start_distances:?})"
    );
    assert!(
        near_end_at_radius >= 2,
        "should have ≥2 corners at distance {r_end} from end; got {near_end_at_radius} (distances: {end_distances:?})"
    );
}
