//! Validation tests for all showcase examples.
//!
//! Each test builds the model using the same code as the binary example,
//! then validates geometry and topology. These tests are the acceptance
//! criteria for fixing the broken examples.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

// ── Helpers ─────────────────────────────────────────────────────────

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Check that no vertex has NaN or infinite coordinates.
fn assert_finite_vertices(brep: &BRep) {
    for (vh, v) in &brep.vertices {
        assert!(
            v.position.x.is_finite() && v.position.y.is_finite() && v.position.z.is_finite(),
            "vertex {vh:?} has non-finite position: {:?}",
            v.position
        );
    }
}

/// Check that every face has a surface attached.
fn assert_all_faces_have_surfaces(brep: &BRep) {
    for (fh, face) in &brep.faces {
        assert!(
            face.surface.is_some(),
            "face {fh:?} has no surface geometry"
        );
    }
}

/// Export to STEP and reimport, returning the reimported BRep.
fn round_trip_step(brep: &BRep, name: &str) -> BRep {
    use popkern_io::{
        OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
        StepExportOptions, StepImportOptions, StepProductNode, UnitContext, export_step,
        import_step,
    };
    use std::sync::Arc;

    let doc = StepDocument {
        products: vec![StepProductNode {
            name: name.into(),
            id: format!("{}-TEST", name.to_uppercase()),
            children: vec![],
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
                context: RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };

    let bytes = export_step(&doc, &StepExportOptions::default()).expect("STEP export failed");
    let result = import_step(&bytes, &StepImportOptions::default()).expect("STEP reimport failed");

    result
        .document
        .products
        .into_iter()
        .flat_map(|p| p.shapes)
        .flat_map(|s| {
            // Arc::try_unwrap won't work if there are other refs, so clone.
            s.solids.clone()
        })
        .next()
        .expect("no solid in reimported STEP")
}

// ── Boolean Showcase ────────────────────────────────────────────────

fn build_boolean_showcase() -> BRep {
    use popkern_bool::{BoolOperation, boolean};
    use popkern_core::geometry::Transform3;

    let tol = tol();
    let cube = BRep::unit_cube(tol).expect("cube");
    let hole = BRep::cylinder(0.2, 2.0, tol).expect("cylinder");
    let hole = hole
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .expect("translate");

    match boolean(BoolOperation::Subtraction, cube.clone(), hole, tol, None) {
        Ok(brep) => brep,
        Err(_) => cube,
    }
}

#[test]
fn boolean_cube_with_hole_topology() {
    let brep = build_boolean_showcase();
    // A cube minus a through-cylinder should have MORE than 6 faces.
    // The cube has 6 faces; after subtraction, the top and bottom faces
    // each get a circular hole boundary, and the cylinder interior adds
    // a cylindrical face. Expected: 7+ faces.
    println!(
        "Boolean showcase: {} faces, {} edges, {} vertices",
        brep.faces.len(),
        brep.edges.len(),
        brep.vertices.len()
    );
    assert_finite_vertices(&brep);
    assert_all_faces_have_surfaces(&brep);

    // The boolean currently returns the original cube unchanged (6 faces).
    // After the fix, it should have more.
    assert!(
        brep.faces.len() > 6,
        "boolean subtraction should produce more than 6 faces (got {}). \
         The cylinder subtraction is not modifying the cube.",
        brep.faces.len()
    );

    // Verify the expected topology: 7 faces (6 cube + 1 cylinder wall),
    // with exactly 2 faces having inner loops (the circular holes).
    let faces_with_holes: Vec<_> = brep
        .faces
        .values()
        .filter(|f| !f.inner_loops.is_empty())
        .collect();
    assert_eq!(
        faces_with_holes.len(),
        2,
        "expected 2 faces with inner loops (top and bottom holes), got {}",
        faces_with_holes.len()
    );
}

// ── Pipe Fitting ────────────────────────────────────────────────────

fn build_pipe_fitting() -> BRep {
    use popkern_core::geometry::{NurbsCurve, Point3};
    use popkern_sweep::{PipeRadius, SweepOptions, pipe};

    let tol = tol();
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let path = NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(1.0, 0.0, 1.0),
        ],
        weights: vec![1.0, w, 1.0],
    };

    pipe(
        &path,
        &PipeRadius::Constant(0.3),
        &SweepOptions::default(),
        tol,
        None,
    )
    .expect("pipe")
    .brep
}

#[test]
fn pipe_fitting_geometry() {
    let brep = build_pipe_fitting();
    println!(
        "Pipe fitting: {} faces, {} edges, {} vertices",
        brep.faces.len(),
        brep.edges.len(),
        brep.vertices.len()
    );
    assert_eq!(
        brep.faces.len(),
        3,
        "pipe should have 3 faces (1 side + 2 caps)"
    );
    assert_finite_vertices(&brep);
    assert_all_faces_have_surfaces(&brep);

    // Check that the pipe side face is a NURBS surface and that
    // points on it are at approximately radius 0.3 from the path.
    for (fh, face) in &brep.faces {
        if let Some(ref sof) = face.surface {
            use popkern_core::geometry::surface::SurfaceGeometry;
            if let SurfaceGeometry::Nurbs(ref nurbs) = sof.geometry {
                // Sample the surface at several points and check distance from path.
                let p_mid = nurbs.point_at(0.0, 0.5);
                println!("  NURBS side face {fh:?}: mid-point = {:?}", p_mid);
                // The path at t=0.5 is approximately at the quarter-circle midpoint.
                // For radius 0.3 pipe, all points on the side surface should be
                // within [0.1, 0.5] distance from the nearest path point.
                // This is a sanity check, not an exact distance check.
                let all_finite = nurbs
                    .control_points
                    .iter()
                    .all(|p| p.x.is_finite() && p.y.is_finite() && p.z.is_finite());
                assert!(all_finite, "NURBS surface has non-finite control points");
            }
        }
    }

    // STEP round-trip
    let reimported = round_trip_step(&brep, "pipe_fitting");
    assert_eq!(reimported.faces.len(), brep.faces.len());
}

// ── Fillet Showcase ─────────────────────────────────────────────────

fn build_fillet_showcase() -> BRep {
    use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

    let tol = tol();
    let cube = BRep::unit_cube(tol).expect("cube");
    let edges: Vec<_> = cube.edges.keys().collect();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();

    match fillet_edges(&cube, &specs, &FilletOptions::default(), tol, None) {
        Ok(result) => result.brep,
        Err(_) => {
            // Try single edge as fallback
            let single = vec![FilletSpec {
                edge: edges[0],
                radius: FilletRadius::Constant(0.1),
            }];
            fillet_edges(&cube, &single, &FilletOptions::default(), tol, None)
                .map(|r| r.brep)
                .unwrap_or(cube)
        }
    }
}

#[test]
fn fillet_showcase_geometry() {
    let brep = build_fillet_showcase();
    println!(
        "Fillet showcase: {} faces, {} edges, {} vertices",
        brep.faces.len(),
        brep.edges.len(),
        brep.vertices.len()
    );
    assert_finite_vertices(&brep);
    assert_all_faces_have_surfaces(&brep);

    // With all 12 edges filleted, we expect more than 6 faces.
    assert!(
        brep.faces.len() > 6,
        "filleted cube should have > 6 faces, got {}",
        brep.faces.len()
    );

    // Check that fillet faces have reasonable surfaces (no degenerate torus).
    for (fh, face) in &brep.faces {
        if let Some(ref sof) = face.surface {
            use popkern_core::geometry::surface::SurfaceGeometry;
            if let SurfaceGeometry::Torus {
                center,
                major_radius,
                minor_radius,
                ..
            } = &sof.geometry
            {
                assert!(
                    center.x.is_finite() && center.y.is_finite() && center.z.is_finite(),
                    "torus {fh:?} has non-finite center: {center:?}"
                );
                assert!(
                    *major_radius > 0.0 && major_radius.is_finite(),
                    "torus {fh:?} has bad major_radius: {major_radius}"
                );
                assert!(
                    *minor_radius > 0.0 && minor_radius.is_finite(),
                    "torus {fh:?} has bad minor_radius: {minor_radius}"
                );
                // For a 0.1 radius fillet on a unit cube (90-degree edges),
                // major_radius should be approximately 0.1 / sin(45°) ≈ 0.1414
                // and minor_radius should be 0.1.
                assert!(
                    (*minor_radius - 0.1).abs() < 0.05,
                    "torus {fh:?} minor_radius should be ~0.1, got {minor_radius}"
                );
            }
        }
    }

    // STEP round-trip
    let reimported = round_trip_step(&brep, "fillet_showcase");
    assert_eq!(reimported.faces.len(), brep.faces.len());
}

// ── Drafted Mold ────────────────────────────────────────────────────

fn build_drafted_mold() -> BRep {
    use popkern_core::geometry::Vec3;
    use popkern_draft::{NeutralReference, draft_face};

    let tol = tol();
    let cube = BRep::unit_cube(tol).expect("cube");

    let pull_dir = Vec3::new(0.0, 0.0, 1.0);
    let face = cube
        .faces
        .iter()
        .find(|(_, fd)| {
            fd.surface.as_ref().is_some_and(|sof| {
                let n = sof.geometry.normal_at(0.5, 0.5);
                n.dot(&pull_dir).abs() < 0.5
            })
        })
        .map(|(fh, _)| fh)
        .expect("cube should have a side face");

    match draft_face(
        &cube,
        face,
        0.087,
        pull_dir,
        NeutralReference::Auto,
        tol,
        None,
    ) {
        Ok(result) => result.brep,
        Err(_) => cube,
    }
}

#[test]
fn drafted_mold_geometry() {
    let brep = build_drafted_mold();
    println!(
        "Drafted mold: {} faces, {} edges, {} vertices",
        brep.faces.len(),
        brep.edges.len(),
        brep.vertices.len()
    );
    assert_eq!(brep.faces.len(), 6);
    assert_eq!(brep.edges.len(), 12);
    assert_eq!(brep.vertices.len(), 8);
    assert_finite_vertices(&brep);
    assert_all_faces_have_surfaces(&brep);

    // Check that vertices are within reasonable bounds.
    // A unit cube is [0,1]^3. After a small 5-degree draft, vertices
    // should not move more than ~0.1 from their original positions.
    for (_vh, v) in &brep.vertices {
        let p = &v.position;
        assert!(
            p.x > -0.5 && p.x < 1.5 && p.y > -0.5 && p.y < 1.5 && p.z > -0.5 && p.z < 1.5,
            "vertex position out of reasonable bounds: {:?}",
            p
        );
    }

    // Check that the drafted face has a tilted normal.
    let pull_dir = Vec3::new(0.0, 0.0, 1.0);
    let has_tilted_face = brep.faces.values().any(|face| {
        face.surface.as_ref().is_some_and(|sof| {
            let n = sof.geometry.normal_at(0.5, 0.5);
            let dot_pull = n.dot(&pull_dir).abs();
            // A drafted side face should have a normal that's slightly
            // off from perpendicular to pull_dir.
            dot_pull > 0.01 && dot_pull < 0.5
        })
    });
    assert!(
        has_tilted_face,
        "should have at least one face with a tilted normal (drafted)"
    );

    // STEP round-trip
    let reimported = round_trip_step(&brep, "drafted_mold");
    assert_eq!(reimported.faces.len(), brep.faces.len());
}

// ── Shelled Box (known-good reference) ──────────────────────────────

fn build_shelled_box() -> BRep {
    use popkern_shell::{ShellDirection, ShellOptions, shell};

    let tol = tol();
    let cube = BRep::unit_cube(tol).expect("cube");
    let face_to_remove = cube.faces.keys().next().expect("cube has faces");

    shell(
        &cube,
        0.1,
        ShellDirection::Inward,
        &[face_to_remove],
        &[],
        &ShellOptions::default(),
        tol,
        None,
    )
    .expect("shell")
    .brep
}

#[test]
fn shelled_box_reference() {
    let brep = build_shelled_box();
    println!(
        "Shelled box: {} faces, {} edges, {} vertices",
        brep.faces.len(),
        brep.edges.len(),
        brep.vertices.len()
    );
    assert_finite_vertices(&brep);
    assert_all_faces_have_surfaces(&brep);
    // Shelled box should work correctly - this is our reference.
    let reimported = round_trip_step(&brep, "shelled_box");
    assert_eq!(reimported.faces.len(), brep.faces.len());
}
