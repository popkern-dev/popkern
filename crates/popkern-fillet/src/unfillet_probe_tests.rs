//! Diagnostic probes for the closed-fillet construct path. All
//! tests in this file are `#[ignore]` — they document specific
//! known limitations and produce inspection output rather than
//! asserting behaviour. Run with
//! `cargo test --lib unfillet -- --ignored --nocapture` to
//! observe their output.

use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

// ---------------------------------------------------------------------------
// Round 5c-2 probe: cylinder-cylinder unfillet
// ---------------------------------------------------------------------------
//
// Both probes are kept ignored after Round 5c-2 (pivoted) because
// they exercise pre-existing popkern-fillet construct-phase
// limitations: cylinder-cylinder edges from a boolean union are
// classified as seams, and plane-cylinder fillets require the
// closed-fillet construction path which is documented as a future
// gap in the diagnostic returned by helpers::split_perpendicular_edges_by_radius.

#[test]
#[ignore = "popkern-fillet does not yet support cylinder-cylinder seam edges; tracked as Round 5c construct-phase work"]
fn probe_5c2_cylinder_cylinder_fillet_then_unfillet() {
    use popkern_bool::{BoolOperation, boolean};
    use popkern_core::geometry::Transform3;

    let cyl_a = BRep::cylinder(0.5, 1.0, tol()).unwrap();
    let cyl_b = BRep::cylinder(0.3, 1.5, tol()).unwrap();
    let cyl_b = cyl_b
        .transform(&Transform3::translation(0.4, 0.0, 0.0))
        .unwrap();

    let united = match boolean(BoolOperation::Union, cyl_a, cyl_b, tol(), None) {
        Ok(b) => {
            eprintln!("probe2: union OK, {} faces", b.faces.len());
            b
        }
        Err(e) => {
            eprintln!("probe2: union FAILED: {e:?}");
            return;
        }
    };

    use popkern_core::geometry::surface::SurfaceGeometry;
    let cyl_cyl_edge = united.edges.keys().find(|e| {
        let he = &united.half_edges[united.edges[*e].half_edge];
        let face_a = united.loops[he.loop_].face;
        let face_b = united.loops[united.half_edges[he.twin].loop_].face;
        let surf_a = united.faces[face_a].surface.as_ref().map(|s| &s.geometry);
        let surf_b = united.faces[face_b].surface.as_ref().map(|s| &s.geometry);
        matches!(
            (surf_a, surf_b),
            (
                Some(SurfaceGeometry::Cylinder { .. }),
                Some(SurfaceGeometry::Cylinder { .. })
            )
        )
    });
    let edge = match cyl_cyl_edge {
        Some(e) => {
            eprintln!("probe2: found cyl-cyl edge {e:?}");
            e
        }
        None => {
            eprintln!("probe2: NO cyl-cyl edge found");
            return;
        }
    };

    let filleted = match crate::fillet_edge(&united, edge, 0.05, tol(), None) {
        Ok(b) => {
            eprintln!(
                "probe2: fillet OK, result has {} faces, {} edges, {} vertices",
                b.brep.faces.len(),
                b.brep.edges.len(),
                b.brep.vertices.len()
            );
            b.brep
        }
        Err(e) => {
            eprintln!(
                "probe2: fillet FAILED: {}",
                match &e {
                    BlendError::PartialFailure { failures, .. } => format!(
                        "{:?}",
                        failures
                            .iter()
                            .map(|(_, f)| format!("{f:?}"))
                            .collect::<Vec<_>>()
                    ),
                    other => format!("{other:?}"),
                }
            );
            return;
        }
    };
    let _ = filleted;
}

// ---------------------------------------------------------------------------
// Round 5c-2 probe: plane-cylinder unfillet
// ---------------------------------------------------------------------------

/// Inspect the loop structure of the closed fillet edge so the
/// Round 5c-3 surgery code knows what topology shape to handle.
/// Prints the half-edges of each side of the fillet edge plus the
/// loop chain on each side.
#[test]
#[ignore = "diagnostic probe — produces inspection output only"]
fn probe_5c3_loop_structure_of_closed_edge() {
    use popkern_bool::{BoolOperation, boolean};
    use popkern_core::geometry::Transform3;

    let cube = BRep::unit_cube(tol()).unwrap();
    let cyl = BRep::cylinder(0.3, 2.0, tol()).unwrap();
    let translated = cyl
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();
    let subtracted = match boolean(BoolOperation::Subtraction, cube, translated, tol(), None) {
        Ok(b) => b,
        Err(_) => return,
    };

    use popkern_core::geometry::surface::SurfaceGeometry;
    let edge_h = subtracted.edges.keys().find(|e| {
        let he = &subtracted.half_edges[subtracted.edges[*e].half_edge];
        let fa = subtracted.loops[he.loop_].face;
        let fb = subtracted.loops[subtracted.half_edges[he.twin].loop_].face;
        let sa = subtracted.faces[fa].surface.as_ref().map(|s| &s.geometry);
        let sb = subtracted.faces[fb].surface.as_ref().map(|s| &s.geometry);
        matches!(
            (sa, sb),
            (
                Some(SurfaceGeometry::Plane { .. }),
                Some(SurfaceGeometry::Cylinder { .. })
            ) | (
                Some(SurfaceGeometry::Cylinder { .. }),
                Some(SurfaceGeometry::Plane { .. })
            )
        )
    });
    let edge_h = match edge_h {
        Some(e) => e,
        None => return,
    };

    let edata = &subtracted.edges[edge_h];
    let he_a = edata.half_edge;
    let he_a_data = &subtracted.half_edges[he_a];
    let he_b = he_a_data.twin;
    let he_b_data = &subtracted.half_edges[he_b];

    let loop_a = he_a_data.loop_;
    let loop_b = he_b_data.loop_;
    let face_a = subtracted.loops[loop_a].face;
    let face_b = subtracted.loops[loop_b].face;

    eprintln!(
        "5c3 loop probe: edge {edge_h:?}, he_a={he_a:?} (loop={loop_a:?}, face={face_a:?}), he_b={he_b:?} (loop={loop_b:?}, face={face_b:?})"
    );
    eprintln!(
        "   he_a: origin={:?}, next={:?}, prev={:?}",
        he_a_data.origin, he_a_data.next, he_a_data.prev
    );
    eprintln!(
        "   he_b: origin={:?}, next={:?}, prev={:?}",
        he_b_data.origin, he_b_data.next, he_b_data.prev
    );

    eprintln!("   loop_a chain:");
    let mut cursor = subtracted.loops[loop_a].half_edge;
    let start = cursor;
    for _ in 0..20 {
        let h = &subtracted.half_edges[cursor];
        let v_origin = subtracted.vertices[h.origin].position;
        eprintln!(
            "     {cursor:?}: edge={:?}, origin={:?}@({:.3},{:.3},{:.3}), next={:?}",
            h.edge, h.origin, v_origin.x, v_origin.y, v_origin.z, h.next
        );
        cursor = h.next;
        if cursor == start {
            break;
        }
    }

    eprintln!("   loop_b chain:");
    let mut cursor = subtracted.loops[loop_b].half_edge;
    let start = cursor;
    for _ in 0..20 {
        let h = &subtracted.half_edges[cursor];
        let v_origin = subtracted.vertices[h.origin].position;
        eprintln!(
            "     {cursor:?}: edge={:?}, origin={:?}@({:.3},{:.3},{:.3}), next={:?}",
            h.edge, h.origin, v_origin.x, v_origin.y, v_origin.z, h.next
        );
        cursor = h.next;
        if cursor == start {
            break;
        }
    }
}

/// Probe whether `compute_fillet_geometry` itself works on a
/// closed plane-cylinder fillet edge. If the geometry pipeline
/// produces a valid `EdgeGeometry` but only `apply_blend_geometry`
/// fails, the Round 5c-3 scope shrinks to "topology surgery only".
#[test]
#[ignore = "diagnostic probe — produces inspection output only"]
fn probe_5c3_geometry_works_on_closed_fillet() {
    use crate::analysis::analyze_edges;
    use crate::geometry_compute::compute_fillet_geometry;
    use crate::types::FilletRadius;
    use popkern_bool::{BoolOperation, boolean};
    use popkern_core::geometry::Transform3;

    let cube = BRep::unit_cube(tol()).unwrap();
    let cyl = BRep::cylinder(0.3, 2.0, tol()).unwrap();
    let translated = cyl
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();
    let subtracted = match boolean(BoolOperation::Subtraction, cube, translated, tol(), None) {
        Ok(b) => b,
        Err(_) => return,
    };

    use popkern_core::geometry::surface::SurfaceGeometry;
    let edge = subtracted.edges.keys().find(|e| {
        let he = &subtracted.half_edges[subtracted.edges[*e].half_edge];
        let fa = subtracted.loops[he.loop_].face;
        let fb = subtracted.loops[subtracted.half_edges[he.twin].loop_].face;
        let sa = subtracted.faces[fa].surface.as_ref().map(|s| &s.geometry);
        let sb = subtracted.faces[fb].surface.as_ref().map(|s| &s.geometry);
        matches!(
            (sa, sb),
            (
                Some(SurfaceGeometry::Plane { .. }),
                Some(SurfaceGeometry::Cylinder { .. })
            ) | (
                Some(SurfaceGeometry::Cylinder { .. }),
                Some(SurfaceGeometry::Plane { .. })
            )
        )
    });
    let edge = match edge {
        Some(e) => e,
        None => return,
    };

    // Confirm closed: v_start == v_end.
    let edata = &subtracted.edges[edge];
    let he_data = &subtracted.half_edges[edata.half_edge];
    let twin = &subtracted.half_edges[he_data.twin];
    eprintln!(
        "5c3 probe: edge {edge:?}, v_start={:?}, v_end={:?}, closed={}",
        he_data.origin,
        twin.origin,
        he_data.origin == twin.origin
    );

    // Run the analysis + geometry pipeline.
    let analyses = match analyze_edges(&subtracted, &[edge], false) {
        Ok(a) => {
            eprintln!("5c3 probe: analyze_edges OK ({} entries)", a.len());
            a
        }
        Err(e) => {
            eprintln!("5c3 probe: analyze_edges FAILED: {e}");
            return;
        }
    };

    if analyses.is_empty() {
        eprintln!("5c3 probe: no analyses produced");
        return;
    }

    let geom = match compute_fillet_geometry(
        &analyses[0],
        &FilletRadius::Constant(0.05),
        &subtracted,
        &tol(),
        crate::types::ContinuityLevel::G1,
        crate::types::CapKind::Open,
    ) {
        Ok(g) => g,
        Err(e) => {
            eprintln!("5c3 probe: compute_fillet_geometry FAILED: {e:?}");
            return;
        }
    };
    eprintln!(
        "5c3 probe: compute_fillet_geometry OK — surface kind={:?}, effective_radius={}",
        match &geom.surface {
            SurfaceGeometry::Plane { .. } => "Plane",
            SurfaceGeometry::Cylinder { .. } => "Cylinder",
            SurfaceGeometry::Sphere { .. } => "Sphere",
            SurfaceGeometry::Torus { .. } => "Torus",
            SurfaceGeometry::Cone { .. } => "Cone",
            SurfaceGeometry::Nurbs(_) => "Nurbs",
            SurfaceGeometry::Gregory(_) => "Gregory",
            SurfaceGeometry::PolyhedralSpline(_) => "PolyhedralSpline",
        },
        geom.effective_radius
    );
    use popkern_core::geometry::curve::CurveEval;
    let (t01, t11) = geom.trim_3d_1.parameter_range();
    let p01 = geom.trim_3d_1.point_at(t01);
    let p11 = geom.trim_3d_1.point_at(t11);
    let (t02, t12) = geom.trim_3d_2.parameter_range();
    let p02 = geom.trim_3d_2.point_at(t02);
    let p12 = geom.trim_3d_2.point_at(t12);
    eprintln!(
        "5c3 probe: trim_3d_1 from ({:.3},{:.3},{:.3}) to ({:.3},{:.3},{:.3})",
        p01.x, p01.y, p01.z, p11.x, p11.y, p11.z
    );
    eprintln!(
        "5c3 probe: trim_3d_2 from ({:.3},{:.3},{:.3}) to ({:.3},{:.3},{:.3})",
        p02.x, p02.y, p02.z, p12.x, p12.y, p12.z
    );
}

/// Probe how far the current unfillet pipeline gets on a
/// plane-cylinder fillet. After Round 5c-2 (pivoted), the fillet
/// engine returns a clean "closed fillet edges not yet supported"
/// diagnostic instead of the historical "no perpendicular edge
/// found" error. Tracked as Round 5c construct-phase work.
#[test]
#[ignore = "popkern-fillet does not yet support closed fillet edges (cylindrical hole rims); diagnostic improved in Round 5c-2"]
fn probe_5c2_plane_cylinder_unfillet() {
    use popkern_bool::{BoolOperation, boolean};
    use popkern_core::geometry::Transform3;

    let cube = BRep::unit_cube(tol()).unwrap();
    let cyl = BRep::cylinder(0.3, 2.0, tol()).unwrap();
    let translated = cyl
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();

    let subtracted = match boolean(BoolOperation::Subtraction, cube, translated, tol(), None) {
        Ok(b) => {
            eprintln!("probe: subtraction OK, {} faces", b.faces.len());
            b
        }
        Err(e) => {
            eprintln!("probe: subtraction FAILED: {e:?}");
            return;
        }
    };

    // Pick a plane-cylinder edge.
    use popkern_core::geometry::surface::SurfaceGeometry;
    let edge = subtracted.edges.keys().find(|e| {
        let he = &subtracted.half_edges[subtracted.edges[*e].half_edge];
        let face_a = subtracted.loops[he.loop_].face;
        let face_b = subtracted.loops[subtracted.half_edges[he.twin].loop_].face;
        let surf_a = subtracted.faces[face_a]
            .surface
            .as_ref()
            .map(|s| &s.geometry);
        let surf_b = subtracted.faces[face_b]
            .surface
            .as_ref()
            .map(|s| &s.geometry);
        matches!(
            (surf_a, surf_b),
            (
                Some(SurfaceGeometry::Plane { .. }),
                Some(SurfaceGeometry::Cylinder { .. })
            ) | (
                Some(SurfaceGeometry::Cylinder { .. }),
                Some(SurfaceGeometry::Plane { .. })
            )
        )
    });
    let edge = match edge {
        Some(e) => {
            eprintln!("probe: found plane-cylinder edge {e:?}");
            e
        }
        None => {
            eprintln!("probe: NO plane-cylinder edge found");
            return;
        }
    };

    let filleted = match crate::fillet_edge(&subtracted, edge, 0.05, tol(), None) {
        Ok(b) => {
            eprintln!("probe: fillet OK");
            b.brep
        }
        Err(e) => {
            eprintln!("probe: fillet FAILED: {e:?}");
            return;
        }
    };

    eprintln!(
        "probe: plane-cylinder fillet result has {} faces, {} edges, {} vertices",
        filleted.faces.len(),
        filleted.edges.len(),
        filleted.vertices.len()
    );

    // Find any new face (one of them should be the blend).
    let pre_face_count = subtracted.faces.len();
    let new_faces: Vec<_> = filleted
        .faces
        .keys()
        .filter(|fh| {
            // Heuristic: face index > pre_face_count would be new.
            // Actually we'd need pre/post sets. Just pick faces that
            // identify as blends.
            crate::identify_blend(&filleted, *fh).is_some()
        })
        .collect();

    eprintln!(
        "probe: identified {} blend face(s) out of {} new faces (pre count {})",
        new_faces.len(),
        filleted.faces.len() - pre_face_count,
        pre_face_count
    );

    for fh in &new_faces {
        eprintln!("probe: trying compute_unfillet_plan on {fh:?}");
        match compute_unfillet_plan(&filleted, *fh) {
            Ok(plan) => {
                eprintln!(
                    "  plan OK: support_a={:?}, support_b={:?}, edge from {:?} to {:?}",
                    plan.support_a,
                    plan.support_b,
                    plan.restored_edge_start,
                    plan.restored_edge_end
                );
                eprintln!(
                    "  trying unfillet_face on {fh:?}: {:?}",
                    unfillet_face(&filleted, *fh)
                        .map(|b| b.faces.len())
                        .map_err(|e| format!("{e:?}"))
                );
            }
            Err(e) => {
                eprintln!("  plan FAILED: {e:?}");
            }
        }
    }
}
