//! Unit tests: profile compilation, part features, selectors, gates.

use crate::prelude::*;

fn block() -> Part {
    Profile::on(Plane::XY)
        .rect(2.0, 1.0)
        .extrude(0.5)
        .expect("block extrude")
}

#[test]
fn extrude_rect_makes_a_box() {
    let part = block();
    assert_eq!(part.face_count(), 6);
    assert_eq!(part.edge_count(), 12);
    assert_eq!(part.open_edge_count().expect("mesh"), 0);
}

#[test]
fn extrude_circle_makes_a_cylinder_with_true_cylinder_wall() {
    let part = Profile::on(Plane::XY)
        .circle(0.5)
        .extrude(1.0)
        .expect("cylinder extrude");
    // Canonical recognition must produce an analytic cylinder wall:
    // the top rim borders (Cylinder, Plane). (Blending a bare swept
    // cylinder rim is a separate fillet-kernel frontier — the
    // selector resolution is what this test pins.)
    let err = part
        .fillet(9999.0, |e| {
            e.between(SurfaceKind::Cylinder, SurfaceKind::Plane)
                .at_height(1.0)
        })
        .unwrap_err();
    assert!(
        !matches!(err, PartError::EmptySelection { .. }),
        "rim edge must match the Cylinder×Plane selector, got {err:?}"
    );
}

#[test]
fn empty_sketch_is_a_typed_error() {
    let err = Profile::on(Plane::XY).extrude(1.0).unwrap_err();
    assert!(matches!(err, PartError::NotApplicable(_)), "{err:?}");
}

#[test]
fn union_and_features_log() {
    let a = block();
    // Penetrating configuration (the flagship pattern): the tool
    // body is embedded in the base footprint and pierces its top.
    // Side-overlap and coplanar-face booleans on swept bodies are
    // known kernel frontiers, tracked in the ledger — not API scope.
    let b = Profile::on(Plane::xy_at(0.2))
        .rect(0.5, 0.4)
        .extrude(0.7)
        .expect("second block");
    let joined = a.union(&b).expect("union");
    assert!(joined.open_edge_count().expect("mesh") == 0);
    let names: Vec<_> = joined.features().entries().iter().map(|f| f.name).collect();
    assert!(names.contains(&"union"));
    assert_eq!(names.iter().filter(|&&n| n == "extrude").count(), 2);
}

#[test]
fn empty_selection_is_an_error_with_context() {
    let err = block().fillet(0.1, |e| e.at_height(99.0)).unwrap_err();
    match err {
        PartError::EmptySelection {
            selection,
            entity,
            available,
        } => {
            assert!(selection.contains("at_height(99)"));
            assert_eq!(entity, "edges");
            assert_eq!(available, 12);
        }
        other => panic!("expected EmptySelection, got {other:?}"),
    }
}

#[test]
fn allow_empty_turns_no_match_into_a_recorded_noop() {
    let part = block()
        .fillet(0.1, |e| e.at_height(99.0).allow_empty())
        .expect("no-op fillet");
    assert_eq!(part.face_count(), 6);
    assert!(
        part.features()
            .entries()
            .iter()
            .any(|f| f.summary.contains("no-op"))
    );
}

#[test]
fn direction_selector_finds_the_vertical_edges() {
    let part = block();
    // A box has 4 vertical edges; fillet them and expect 4 new faces.
    let filleted = part
        .fillet(0.1, |e| e.straight().direction(Dir::Z))
        .expect("vertical fillets");
    assert_eq!(filleted.face_count(), 6 + 4);
    assert_eq!(filleted.open_edge_count().expect("mesh"), 0);
}

#[test]
fn or_unions_selector_branches() {
    let part = block();
    let top_and_vertical = part
        .fillet(0.1, |e| {
            e.straight()
                .direction(Dir::Z)
                .or(|e| e.straight().at_height(0.5))
        })
        .expect("chain fillet");
    // 4 verticals + 4 top edges blended, with corner patches.
    assert!(top_and_vertical.face_count() > 6 + 8);
    assert_eq!(top_and_vertical.open_edge_count().expect("mesh"), 0);
}

#[test]
fn seam_edges_are_skipped_not_fatal() {
    let cyl = Profile::on(Plane::XY)
        .circle(0.5)
        .extrude(1.0)
        .expect("cylinder");
    // direction(Z) matches only the cylinder's seam line, which can
    // never blend: with allow_empty the operation records a no-op
    // rather than failing.
    let out = cyl.fillet(0.1, |e| e.straight().direction(Dir::Z).allow_empty());
    assert!(out.is_ok(), "{out:?}");
}

#[test]
fn polygon_profile_extrudes() {
    let tri = Profile::on(Plane::XY)
        .polygon(&[(0.0, 0.0), (2.0, 0.0), (0.0, 1.5)])
        .extrude(0.4)
        .expect("triangle prism");
    assert_eq!(tri.face_count(), 5);
    assert_eq!(tri.open_edge_count().expect("mesh"), 0);
}

/// Regression test for the sweep longitudinal-edge fix: every edge
/// of a RAW sweep output (no boolean, no fillet) must carry a curve
/// whose trimmed endpoints coincide with the edge's vertices. Bare
/// sweep solids used to leave longitudinal edges curveless (invisible
/// until a consumer needed canonical samples), and closed cap rims
/// carried arcs whose start angle was measured in the wrong frame.
#[test]
fn raw_sweep_edges_carry_matching_curves() {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let r = 0.5;
    let wire = popkern_core::geometry::NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: vec![
            popkern_core::geometry::Point3::new(r, 0.0, 0.0),
            popkern_core::geometry::Point3::new(r, r, 0.0),
            popkern_core::geometry::Point3::new(0.0, r, 0.0),
            popkern_core::geometry::Point3::new(-r, r, 0.0),
            popkern_core::geometry::Point3::new(-r, 0.0, 0.0),
            popkern_core::geometry::Point3::new(-r, -r, 0.0),
            popkern_core::geometry::Point3::new(0.0, -r, 0.0),
            popkern_core::geometry::Point3::new(r, -r, 0.0),
            popkern_core::geometry::Point3::new(r, 0.0, 0.0),
        ],
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    };
    let path = popkern_core::geometry::NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
            popkern_core::geometry::Point3::new(0.0, 0.0, 0.5),
        ],
        weights: vec![1.0, 1.0],
    };
    let result = popkern_sweep::sweep(
        &popkern_sweep::SweepProfile::Wire(wire),
        &path,
        &popkern_sweep::SweepOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("raw sweep");

    use popkern_core::geometry::curve::CurveEval;
    let brep = &result.brep;
    for (eh, e) in brep.edges.iter() {
        let he = e.half_edge;
        let v0 = brep.vertices[brep.half_edges[he].origin].position;
        let v1 = brep.vertices[brep.half_edges[brep.half_edges[he].twin].origin].position;
        let c = e
            .curve
            .as_ref()
            .unwrap_or_else(|| panic!("edge {eh:?} has no curve on raw sweep output"));
        let ceval: &dyn CurveEval = &c.geometry;
        let p0 = ceval.point_at(c.t_start);
        let p1 = ceval.point_at(c.t_end);
        let closed = v0.distance_to(&v1) < 1e-9;
        if closed {
            // Closed rim: both trim ends must land on the seam vertex.
            assert!(
                p0.distance_to(&v0) < 1e-9 && p1.distance_to(&v0) < 1e-9,
                "edge {eh:?}: closed-rim trim off its seam vertex (p0 {p0:?}, p1 {p1:?}, v {v0:?})"
            );
        } else {
            let forward = p0.distance_to(&v0) < 1e-9 && p1.distance_to(&v1) < 1e-9;
            let reversed = p0.distance_to(&v1) < 1e-9 && p1.distance_to(&v0) < 1e-9;
            assert!(
                forward || reversed,
                "edge {eh:?}: curve endpoints do not match vertices \
                 (p0 {p0:?}, p1 {p1:?}, v0 {v0:?}, v1 {v1:?})"
            );
        }
    }
}

#[test]
fn chamfer_verticals() {
    let part = block()
        .chamfer(0.08, |e| e.straight().direction(Dir::Z))
        .expect("chamfer");
    // Face count and validation gate the API contract; chamfer
    // output tessellation is not yet watertight on swept boxes (a
    // chamfer-kernel frontier recorded in the PK-100 ledger notes,
    // parallel to the rim-fillet gap).
    assert_eq!(part.face_count(), 6 + 4);
}

#[test]
fn shell_surfaces_the_kernel_validity_gap_as_a_typed_error() {
    // The shell kernel's open-shell (cup) output does not currently
    // satisfy strict Euler–Poincaré validation (the shelled_box
    // example works only because it never validates: 16 − 40 + 14 ≠
    // 2). The Part invariant — every Ok(Part) wraps a validated
    // body — makes that a typed Validation error instead of a
    // silently-invalid part. Recorded as a kernel frontier in the
    // PK-100 ledger notes; this test pins BOTH the selector
    // resolution (exactly one +Z face) and the honest failure mode,
    // and will fail deliberately when the kernel gap closes so the
    // expectation gets flipped to success.
    let cube = popkern_core::topology::BRep::unit_cube(Tolerance::default()).expect("unit cube");
    let part = Part::from_brep(cube, Tolerance::default(), "unit cube").expect("wrap");
    let err = part.shell(0.08, |f| f.facing(Dir::Z)).unwrap_err();
    match err {
        PartError::Validation { feature, .. } => {
            assert!(
                feature.contains("1 openings"),
                "selector must match exactly one face: {feature}"
            );
        }
        other => panic!("expected Validation (kernel gap), got {other:?}"),
    }
}

#[test]
fn translate_moves_the_features_log_along() {
    let part = block().translate(5.0, 0.0, 0.0).expect("translate");
    assert_eq!(part.face_count(), 6);
    let names: Vec<_> = part.features().entries().iter().map(|f| f.name).collect();
    assert_eq!(names, vec!["extrude", "transform"]);
}

#[test]
fn torus_passes_genus_aware_validation() {
    // validate() rejected the kernel's own torus constructor until
    // the Euler check became genus-aware; pin both forms.
    let torus = popkern_core::topology::BRep::torus(2.0, 0.5, Tolerance::default()).expect("torus");
    assert!(
        torus.validate().is_ok(),
        "genus-aware validate must accept a torus"
    );
    assert!(
        torus.validate_euler_poincare().is_err(),
        "strict genus-0 form still rejects it (direct callers rely on that)"
    );
}

#[test]
fn profile_placement_is_exact() {
    // Two defects hid here for the kernel's whole life: the profile
    // centroid was biased by a double-counted seam sample (bodies
    // landed r/65 off), and the sweep frame rotated profiles 90°
    // about the path. Pin exact placement and orientation.
    let drill = Profile::on(Plane::XY)
        .circle(0.25)
        .at(1.0, 0.0)
        .extrude(0.5)
        .expect("drill");
    let mut found_cylinder = false;
    for (_, f) in drill.brep().faces.iter() {
        if let Some(popkern_core::geometry::surface::SurfaceGeometry::Cylinder {
            origin,
            radius,
            ..
        }) = f.surface.as_ref().map(|s| &s.geometry)
        {
            found_cylinder = true;
            assert!(
                (origin.x - 1.0).abs() < 1e-9
                    && origin.y.abs() < 1e-9
                    && (radius - 0.25).abs() < 1e-9,
                "cylinder must sit exactly at (1, 0) r=0.25; got ({}, {}) r={}",
                origin.x,
                origin.y,
                radius
            );
        }
    }
    assert!(
        found_cylinder,
        "canonical recognition must produce a Cylinder wall"
    );

    let block = Profile::on(Plane::XY)
        .rect(2.0, 1.0)
        .extrude(0.5)
        .expect("block");
    let (mut max_x, mut max_y) = (0.0_f64, 0.0_f64);
    for (_, v) in block.brep().vertices.iter() {
        max_x = max_x.max(v.position.x.abs());
        max_y = max_y.max(v.position.y.abs());
    }
    assert!(
        (max_x - 1.0).abs() < 1e-9 && (max_y - 0.5).abs() < 1e-9,
        "rect(2, 1) must span x ±1, y ±0.5 (frame rotation); got x ±{max_x}, y ±{max_y}"
    );
}

#[test]
fn chamfered_box_is_watertight() {
    // Chamfer bevels used to be authored perpendicular to the true
    // bevel (offset row walked OUT of the face), then wired around
    // retained corner spikes with no pcurves — three fixes pinned
    // by one gate.
    let part = block()
        .chamfer(0.08, |e| e.straight().direction(Dir::Z))
        .expect("chamfer");
    assert_eq!(part.face_count(), 6 + 4);
    assert_eq!(part.open_edge_count().expect("mesh"), 0);
}

#[test]
fn pocketed_and_chamfered_plate_is_watertight() {
    // The chamfer-bracket example's shape: penetrating rectangular
    // pocket + outer vertical chamfers.
    let plate = Profile::on(Plane::XY)
        .rect(3.0, 1.8)
        .extrude(0.5)
        .expect("plate");
    let tool = Profile::on(Plane::xy_at(0.2))
        .rect(2.0, 1.0)
        .extrude(0.5)
        .expect("tool");
    let pocketed = plate.subtract(&tool).expect("pocket");
    assert_eq!(pocketed.open_edge_count().expect("mesh"), 0);
    let chamfered = pocketed
        .chamfer(0.18, |e| {
            e.straight()
                .direction(Dir::Z)
                .chord_length_between(0.4, 1.0)
        })
        .expect("chamfer");
    assert_eq!(chamfered.face_count(), 15);
    assert_eq!(chamfered.open_edge_count().expect("mesh"), 0);
}

/// FRONTIER PIN (flips to a watertight assert when the kernel gap
/// closes): interior through-hole walls have never tessellated —
/// the periodic wall ring's u-winding comes out ±2 periods instead
/// of 0 and boundary-first declines. Constructor bodies, no sweep
/// involvement. Ledger PK-100.
#[test]
fn through_hole_tessellation_gap_is_pinned() {
    use popkern_core::topology::BRep;
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let t = popkern_core::geometry::transform::Transform3::translation(0.5, 0.5, -0.5);
    let drill = BRep::cylinder(0.2, 2.0, Tolerance::default())
        .expect("cyl")
        .transform(&t)
        .expect("move");
    let out = popkern_bool::boolean(
        popkern_bool::BoolOperation::Subtraction,
        cube,
        drill,
        Tolerance::default(),
        None,
    )
    .expect("subtract");
    let part = Part::from_brep(out, Tolerance::default(), "cube with through-hole")
        .expect("hole solid validates (genus-aware)");
    let opens = part.open_edge_count().expect("mesh");
    assert!(
        opens > 0,
        "through-hole tessellation gap has closed ({opens} opens) — \
         flip this pin to assert_eq!(opens, 0) and update the ledger"
    );
}
