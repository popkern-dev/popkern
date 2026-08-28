//! PK-324 acceptance gates: curved edge recompute staleness (row 470)
//! and vertex recompute determinism at valence>3 (row 471).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::classify::classify_faces;
use popkern_shell::edge_recompute::recompute_edges;
use popkern_shell::inner_skin::build_inner_skin;
use popkern_shell::surface_offset::compute_offset_surfaces;
use popkern_shell::vertex_recompute::recompute_vertices;
use popkern_shell::{ShellDirection, ShellOptions, thicken::thicken};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn curved_bodies(t: Tolerance) -> Vec<(&'static str, BRep)> {
    vec![
        ("cylinder", BRep::cylinder(0.5, 1.0, t).unwrap()),
        ("cone", BRep::cone(0.5, 1.0, t).unwrap()),
        ("sphere", BRep::sphere(0.5, t).unwrap()),
        ("torus", BRep::torus(0.6, 0.25, t).unwrap()),
    ]
}

/// Row 470: no clone-curve survivals. After the inward recompute
/// ladder, every inner-skin edge adjacent to a non-planar offset
/// surface must carry a curve measurably DIFFERENT from the outer
/// clone it started as — a stale clone sits a full offset distance off
/// the inner skin and previously escaped whenever validation had no
/// reachable check for it.
#[test]
fn no_clone_curve_survives_curved_recompute() {
    let thickness = 0.04;
    for (name, brep) in &curved_bodies(tol()) {
        let classification =
            classify_faces(brep, thickness, &[], ShellDirection::Inward, &[]).unwrap();
        let offset_result = compute_offset_surfaces(brep, &classification, true, &tol()).unwrap();
        let (mut inner, skin_map) = build_inner_skin(brep, &[]).unwrap();
        let clones: Vec<_> = inner
            .edges
            .iter()
            .map(|(eh, e)| (eh, e.curve.clone()))
            .collect();
        recompute_vertices(&mut inner, &offset_result.surfaces, &skin_map, brep, &tol()).unwrap();
        recompute_edges(&mut inner, &offset_result.surfaces, &skin_map, brep, &tol()).unwrap();

        let reverse_face = skin_map.reverse_face_map();
        let mut checked = 0usize;
        for (eh, clone) in clones {
            let Some(clone) = clone else { continue };
            let edge = &inner.edges[eh];
            let Some(new_curve) = &edge.curve else {
                panic!("{name}: edge {eh:?} lost its curve in recompute");
            };
            // Only edges with a non-planar adjacent offset surface are
            // in this pin's scope (plane-plane lines are covered by the
            // exact arm and endpoint validation).
            let he = &inner.half_edges[edge.half_edge];
            let twin = &inner.half_edges[he.twin];
            let non_planar_adjacent = [he.loop_, twin.loop_].iter().any(|&lp| {
                inner
                    .loops
                    .get(lp)
                    .and_then(|l| reverse_face.get(&l.face))
                    .and_then(|fh| offset_result.surfaces.get(fh))
                    .is_some_and(|s| !matches!(s, SurfaceGeometry::Plane { .. }))
            });
            if !non_planar_adjacent {
                continue;
            }
            let mid_old = clone
                .geometry
                .point_at(clone.t_start + (clone.t_end - clone.t_start) * 0.5);
            let mid_new = new_curve
                .geometry
                .point_at(new_curve.t_start + (new_curve.t_end - new_curve.t_start) * 0.5);
            assert!(
                mid_old.distance_to(&mid_new) > thickness * 0.25,
                "{name}: edge {eh:?} kept its stale clone curve (moved {})",
                mid_old.distance_to(&mid_new)
            );
            checked += 1;
        }
        assert!(checked > 0, "{name}: pin checked no edges — vacuous");
    }
}

/// Row 471: vertex recompute is deterministic at valence>3 — repeated
/// thickens of the 6-gonal pyramid (apex valence 6, least-squares over
/// all incident plane constraints in sorted-handle order) produce
/// bit-identical vertex positions.
#[test]
fn valence_gt3_vertex_recompute_is_bit_identical_across_runs() {
    let build = || {
        let pyramid = BRep::n_gonal_pyramid(6, 1.0, 1.2, tol()).unwrap();
        let result = thicken(&pyramid, 0.04, &ShellOptions::default(), tol()).unwrap();
        let mut keys: Vec<_> = result.brep.vertices.keys().collect();
        keys.sort();
        keys.iter()
            .map(|&vh| {
                let p = result.brep.vertices[vh].position;
                (p.x.to_bits(), p.y.to_bits(), p.z.to_bits())
            })
            .collect::<Vec<_>>()
    };
    let first = build();
    assert!(!first.is_empty());
    for run in 0..3 {
        assert_eq!(first, build(), "vertex bits drifted on repeat run {run}");
    }
}

/// Row 471 support: the pyramid apex actually exercises the
/// least-squares arm (6 planes meet there) and the solved apex lies
/// within tolerance of every incident offset plane — the constraints
/// are consistent for a pyramid apex, so least-squares must be exact.
#[test]
fn pyramid_apex_least_squares_hits_all_planes() {
    let pyramid = BRep::n_gonal_pyramid(6, 1.0, 1.2, tol()).unwrap();
    let thickness = 0.04;
    let classification =
        classify_faces(&pyramid, thickness, &[], ShellDirection::Inward, &[]).unwrap();
    let offset_result = compute_offset_surfaces(&pyramid, &classification, true, &tol()).unwrap();
    let (mut inner, skin_map) = build_inner_skin(&pyramid, &[]).unwrap();
    recompute_vertices(
        &mut inner,
        &offset_result.surfaces,
        &skin_map,
        &pyramid,
        &tol(),
    )
    .unwrap();
    // The apex is the vertex with max z.
    let apex = inner
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(_, v)| v.position)
        .unwrap();
    let mut lateral_planes = 0usize;
    for surf in offset_result.surfaces.values() {
        let SurfaceGeometry::Plane { origin, normal, .. } = surf else {
            continue;
        };
        // Lateral planes only: the base plane is not incident at the apex.
        if normal.z.abs() > 0.99 {
            continue;
        }
        lateral_planes += 1;
        let d = (apex - *origin).dot(normal).abs();
        assert!(d < 1e-9, "apex off an incident lateral offset plane by {d}");
    }
    assert_eq!(lateral_planes, 6, "expected 6 lateral planes at the apex");
}
