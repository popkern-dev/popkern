//! PK-244: winding-flag honesty audit across cap-author families.
//!
//! Every consumer of the stored orientation flag (edge convexity
//! classification, mesher triangle winding, STEP `same_sense`,
//! measured volume sign) rests on the interior-left invariant the
//! core audit measures: outer loop CCW in UV ⇔ flag true. The strip
//! authors were made honest in PK-238/PK-241; this file measures the
//! CAP authors, one fixture per family. A failure lists the
//! dishonest faces with their surface class — that IS the defect
//! record for the author behind them.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges, vertex_blend_3d};

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{fillet_corner_full, incident_edges, pick_vertex_with_valence};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Format every dishonest face with its surface class for the
/// assertion message; empty = the fixture's authors are honest.
fn dishonest_report(brep: &BRep) -> Vec<String> {
    brep.audit_winding_honesty()
        .into_iter()
        .filter(|v| v.honest == Some(false))
        .map(|v| {
            let surface = brep
                .faces
                .get(v.face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| format!("{:?}", std::mem::discriminant(&s.geometry)))
                .unwrap_or_else(|| "<no-surface>".into());
            format!(
                "{:?} [{surface}] area {:.3e} flag {}",
                v.face, v.signed_uv_area, v.orientation
            )
        })
        .collect()
}

fn assert_honest(brep: &BRep, fixture: &str) {
    let report = dishonest_report(brep);
    assert!(
        report.is_empty(),
        "{fixture}: {} dishonest face(s):\n  {}",
        report.len(),
        report.join("\n  ")
    );
}

/// Cube corner, uniform radii: setback/sphere-promotion cap family.
#[test]
fn cube_corner_uniform_caps_are_winding_honest() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let vertex = pick_vertex_with_valence(&brep, 3).unwrap();
    let edges = incident_edges(&brep, vertex);
    let (out, _) = fillet_corner_full(&brep, &edges, &[0.2, 0.2, 0.2], true, false);
    assert_honest(&out, "cube_corner_uniform");
}

/// Cube corner, mixed radii: variable-radius setback caps.
#[test]
fn cube_corner_mixed_caps_are_winding_honest() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let vertex = pick_vertex_with_valence(&brep, 3).unwrap();
    let edges = incident_edges(&brep, vertex);
    let (out, _) = fillet_corner_full(&brep, &edges, &[0.2, 0.3, 0.25], true, false);
    assert_honest(&out, "cube_corner_mixed");
}

/// Square pyramid apex (N=4): direct-apex partition cap family.
#[test]
fn pyramid4_apex_caps_are_winding_honest() {
    let brep = BRep::square_pyramid(1.0, 0.8, tol()).unwrap();
    let vertex = pick_vertex_with_valence(&brep, 4).unwrap();
    let edges = incident_edges(&brep, vertex);
    let (out, _) = fillet_corner_full(&brep, &edges, &[0.15; 4], true, false);
    assert_honest(&out, "pyramid4_apex");
}

/// Pentagonal pyramid apex (N=5): KP20 polar cap family.
#[test]
fn pentagon5_apex_caps_are_winding_honest() {
    let brep = BRep::n_gonal_pyramid(5, 0.6, 0.8, tol()).unwrap();
    let vertex = pick_vertex_with_valence(&brep, 5).unwrap();
    let edges = incident_edges(&brep, vertex);
    let (out, _) = fillet_corner_full(&brep, &edges, &[0.1; 5], true, true);
    assert_honest(&out, "pentagon5_apex");
}

/// Single open-edge fillet: terminal strip-end treatment authors.
#[test]
fn single_edge_open_strip_is_winding_honest() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let edge = brep.edges.keys().next().unwrap();
    let specs = [FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.2),
    }];
    let result = fillet_edges(&brep, &specs, &FilletOptions::default(), tol(), None).unwrap();
    assert_honest(&result.brep, "single_edge_open_strip");
}

/// Closed-rim fillet on a cylinder top edge: closed-annulus author.
#[test]
fn cylinder_rim_closed_annulus_is_winding_honest() {
    let brep = BRep::cylinder(0.5, 1.0, tol()).unwrap();
    // The top rim: the non-seam (non-Line) edge whose vertices sit on
    // the top cap (the constructor authors rim circles as Nurbs).
    let mut rim = None;
    for (eh, edge) in &brep.edges {
        if let Some(curve) = edge.curve.as_ref()
            && !matches!(
                curve.geometry,
                popkern_core::geometry::CurveGeometry::Line { .. }
            )
        {
            let he = &brep.half_edges[edge.half_edge];
            let z = brep.vertices[he.origin].position.z;
            if z > 0.5 {
                rim = Some(eh);
                break;
            }
        }
    }
    let rim = rim.expect("cylinder has a top rim edge");
    let specs = [FilletSpec {
        edge: rim,
        radius: FilletRadius::Constant(0.1),
    }];
    let result = fillet_edges(&brep, &specs, &FilletOptions::default(), tol(), None).unwrap();
    assert_honest(&result.brep, "cylinder_rim_closed_annulus");
}

/// Blend-on-blend chain: vertex blend then a secondary fillet on a
/// strip-support edge — torus promotion + section caps (the honest
/// template author) plus the sphere/cylinder vertex-blend authors.
#[test]
fn blend_on_blend_chain_is_winding_honest() {
    let cube = BRep::unit_cube(tol()).unwrap();
    let vertex = pick_vertex_with_valence(&cube, 3).unwrap();
    let blended = vertex_blend_3d(&cube, vertex, 0.3, &FilletOptions::default(), tol(), None)
        .unwrap()
        .brep;
    assert_honest(&blended, "vertex_blend_3d");

    // Shortest cylinder-plane edge = the far-corner blend-on-blend
    // candidate (same selection as the §12.3 pin).
    let mut best: Option<(popkern_core::topology::EdgeHandle, f64)> = None;
    for (eh, edge) in &blended.edges {
        let he = &blended.half_edges[edge.half_edge];
        let twin = &blended.half_edges[he.twin];
        let f1 = blended.loops[he.loop_].face;
        let f2 = blended.loops[twin.loop_].face;
        let s1 = blended.faces[f1].surface.as_ref().map(|s| &s.geometry);
        let s2 = blended.faces[f2].surface.as_ref().map(|s| &s.geometry);
        use popkern_core::geometry::surface::SurfaceGeometry as Sg;
        let cyl_plane = matches!(
            (s1, s2),
            (Some(Sg::Cylinder { .. }), Some(Sg::Plane { .. }))
                | (Some(Sg::Plane { .. }), Some(Sg::Cylinder { .. }))
        );
        if cyl_plane {
            let a = blended.vertices[he.origin].position;
            let b = blended.vertices[twin.origin].position;
            let len = a.distance_to(&b);
            if best.is_none_or(|(_, l)| len < l) {
                best = Some((eh, len));
            }
        }
    }
    let (edge, _) = best.expect("blend output has a cylinder-plane edge");
    let specs = [FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.02),
    }];
    let result = fillet_edges(&blended, &specs, &FilletOptions::default(), tol(), None).unwrap();
    assert_honest(&result.brep, "blend_on_blend_chain");
}
