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

/// Round 5m: `modify_open_in_place_bypoints` accepts a `ByPoints`
/// radius profile. Corners are positioned by the radii at
/// parameters 0 and 1, and the fillet surface is a NURBS sweep
/// shaped by the full piecewise-linear profile.
#[test]
fn modify_open_in_place_bypoints_assigns_nurbs_sweep() {
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

    // 3-point profile: 0.05 → 0.10 → 0.15 along the edge.
    let pts = vec![(0.0, 0.05), (0.5, 0.10), (1.0, 0.15)];
    let modified = crate::modify_inplace::modify_open_in_place_bypoints(&r1.brep, blend, &pts)
        .expect("ByPoints in-place should succeed");

    let surf = modified.faces[blend]
        .surface
        .as_ref()
        .expect("modified blend should have a surface");
    match &surf.geometry {
        SurfaceGeometry::Nurbs(_) => { /* expected */ }
        other => panic!("ByPoints in-place should assign a NURBS sweep, got {other:?}"),
    }

    // Corners at the start should be at distance r_start (0.05),
    // corners at the end at distance r_end (0.15).
    let plan = crate::unfillet::compute_unfillet_plan(&r1.brep, blend).unwrap();
    let p_start = plan.restored_edge_start;
    let p_end = plan.restored_edge_end;
    let mut start_at_radius = 0;
    let mut end_at_radius = 0;
    for (_, v) in &modified.vertices {
        let pos = v.position;
        let d_start = ((pos.x - p_start.x).powi(2)
            + (pos.y - p_start.y).powi(2)
            + (pos.z - p_start.z).powi(2))
        .sqrt();
        let d_end =
            ((pos.x - p_end.x).powi(2) + (pos.y - p_end.y).powi(2) + (pos.z - p_end.z).powi(2))
                .sqrt();
        if (d_start - 0.05).abs() < 1e-6 {
            start_at_radius += 1;
        }
        if (d_end - 0.15).abs() < 1e-6 {
            end_at_radius += 1;
        }
    }
    assert!(
        start_at_radius >= 2,
        "should have ≥2 corners at distance 0.05 from start; got {start_at_radius}"
    );
    assert!(
        end_at_radius >= 2,
        "should have ≥2 corners at distance 0.15 from end; got {end_at_radius}"
    );
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
