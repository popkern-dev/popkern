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

/// 12-vertex hexagram (6-point star) for the polygon-profile pins:
/// outer radius `r`, inner radius `r/√3`, starting at angle 0.
fn hexagram_points(r: f64) -> Vec<(f64, f64)> {
    (0..12)
        .map(|k| {
            let angle = std::f64::consts::TAU * k as f64 / 12.0;
            let radius = if k % 2 == 0 { r } else { r / 3.0_f64.sqrt() };
            (radius * angle.cos(), radius * angle.sin())
        })
        .collect()
}

/// Shoelace area of a simple polygon.
fn polygon_area(points: &[(f64, f64)]) -> f64 {
    let n = points.len();
    (0..n)
        .map(|i| {
            let (x0, y0) = points[i];
            let (x1, y1) = points[(i + 1) % n];
            x0 * y1 - x1 * y0
        })
        .sum::<f64>()
        .abs()
        / 2.0
}

/// Mesh volume via signed tetrahedra (watertight mesh required).
fn mesh_volume(mesh: &popkern_mesh::TriMesh) -> f64 {
    mesh.triangles
        .iter()
        .map(|t| {
            let a = mesh.vertices[t[0]].position;
            let b = mesh.vertices[t[1]].position;
            let c = mesh.vertices[t[2]].position;
            (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                + a.z * (b.x * c.y - b.y * c.x))
                / 6.0
        })
        .sum::<f64>()
        .abs()
}

/// W0 pin: NON-CONVEX (reflex-vertex) polygon profiles must extrude
/// watertight with correct volume. The hexagram has six 240° reflex
/// vertices; the laptop-stand outlines are reflex too.
#[test]
fn nonconvex_polygon_extrude_is_watertight_with_exact_volume() {
    let star = hexagram_points(12.0);
    let part = Profile::on(Plane::XY)
        .polygon(&star)
        .extrude(30.0)
        .expect("hexagram extrude");
    assert_eq!(part.open_edge_count().expect("mesh"), 0, "star prism opens");
    let mesh = part
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let expected = polygon_area(&star) * 30.0;
    let got = mesh_volume(&mesh);
    assert!(
        (got - expected).abs() < expected * 1e-3,
        "star prism volume {got:.3} vs analytic {expected:.3}"
    );

    // Reflex L-outline as the second regime.
    let ell = [
        (0.0, 0.0),
        (40.0, 0.0),
        (40.0, 12.0),
        (14.0, 12.0),
        (14.0, 30.0),
        (0.0, 30.0),
    ];
    let part = Profile::on(Plane::XY)
        .polygon(&ell)
        .extrude(8.0)
        .expect("L extrude");
    assert_eq!(part.open_edge_count().expect("mesh"), 0, "L prism opens");
    let mesh = part
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let expected = polygon_area(&ell) * 8.0;
    let got = mesh_volume(&mesh);
    assert!(
        (got - expected).abs() < expected * 1e-3,
        "L prism volume {got:.3} vs analytic {expected:.3}"
    );
}

/// W1 pin: polygon-profile BLIND pockets (the hexagram sockets).
/// Star pocket in a cylinder cap, star pocket in a prism — each cut
/// watertight. (Every socket here is the FIRST boolean on its body;
/// the double-socket column is pinned separately as a frontier.)
#[test]
fn star_blind_pockets_are_watertight() {
    let star = hexagram_points(12.0);

    // Socket in a cylinder cap (the column-end class).
    let cylinder = Profile::on(Plane::XY)
        .circle(35.0)
        .extrude(60.0)
        .expect("column blank");
    let tool = Profile::on(Plane::xy_at(-1.0))
        .polygon(&star)
        .extrude(17.0)
        .expect("socket tool");
    let socketed = cylinder.subtract(&tool).expect("bottom socket");
    assert_eq!(socketed.open_edge_count().expect("mesh"), 0);
    assert!(
        socketed.face_count() > cylinder.face_count(),
        "socket must add faces (silent no-op guard): {} vs {}",
        socketed.face_count(),
        cylinder.face_count()
    );

    // Socket in a prism (the base-boss regime).
    let prism = Profile::on(Plane::XY)
        .rect(60.0, 60.0)
        .extrude(22.0)
        .expect("boss blank");
    let tool = Profile::on(Plane::xy_at(5.0))
        .polygon(&star)
        .extrude(18.0)
        .expect("boss socket tool");
    let socketed = prism.subtract(&tool).expect("boss socket");
    assert_eq!(socketed.open_edge_count().expect("mesh"), 0);
}

/// W1 pin: the COLUMN configuration — sockets cut into BOTH ends of
/// one body (the second boolean runs on the first's output). Closed
/// by the closed-ring guard in fix_pcurve_endpoints: the fixer was
/// dragging a transplanted closed rim pcurve's end control points to
/// the seam-vertex UV (a closed ring's start is an arbitrary anchor,
/// not a constraint), warping the ring so the next boolean's clip saw
/// a phantom crossing and dropped socket segments.
#[test]
fn double_socket_column_is_watertight() {
    let star = hexagram_points(12.0);
    let cylinder = Profile::on(Plane::XY)
        .circle(35.0)
        .extrude(60.0)
        .expect("blank");
    let bot = Profile::on(Plane::xy_at(-1.0))
        .polygon(&star)
        .extrude(17.0)
        .expect("bottom tool");
    let top = Profile::on(Plane::xy_at(44.0))
        .polygon(&star)
        .extrude(17.0)
        .expect("top tool");
    let column = cylinder
        .subtract(&bot)
        .expect("bottom socket")
        .subtract(&top)
        .expect("top socket");
    assert_eq!(column.open_edge_count().expect("mesh"), 0);
    // Both sockets must actually exist (silent no-op guard): the
    // 28mm web between the two 16mm-deep star pockets.
    let mesh = column
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let volume = mesh_volume(&mesh);
    let blank = std::f64::consts::PI * 35.0 * 35.0 * 60.0;
    let pocket = polygon_area(&star) * 16.0;
    let expected = blank - 2.0 * pocket;
    assert!(
        (volume - expected).abs() / expected < 0.02,
        "column volume {volume:.0} vs expected {expected:.0} (blank {blank:.0}, 2 pockets)"
    );
}

/// FRONTIER PIN (W3 — flips to a watertight assert when closed):
/// union of two COAXIAL star prisms (the peg's stepped pilot tips —
/// a star boss piercing a star end face). Cellular assembly fails
/// with a NON-MANIFOLD EDGE (a fragment of the inner prism's wall
/// that lies strictly inside the outer body appears to stay wired),
/// and the graph-merge fallback produces an OpenLoop whose gap is
/// exactly the offset between the two star outlines (0.6), i.e. a
/// loop mixing the two rings. Until closed, the laptop-stand peg
/// ships as a plain prism (containment recorded in
/// laptop_stand_api.rs; the stepped lead-in is a fit nicety, not
/// structural).
#[test]
fn star_boss_union_gap_is_pinned() {
    let body_star = hexagram_points(11.8);
    let pilot_star = hexagram_points(11.2);
    let body = Profile::on(Plane::xy_at(2.0))
        .polygon(&body_star)
        .extrude(26.0)
        .expect("body");
    let pilot = Profile::on(Plane::XY)
        .polygon(&pilot_star)
        .extrude(6.0)
        .expect("pilot");
    match body.union(&pilot) {
        Ok(peg) => {
            let opens = peg.open_edge_count().expect("mesh");
            assert!(
                opens > 0,
                "star-boss union gap has closed ({opens} opens) — flip this \
                 pin to a watertight assert and restore the peg's pilot tips"
            );
        }
        Err(_) => {
            // Current behavior: validation failure (OpenLoop from the
            // mixed-ring loop). Either failure mode pins the gap.
        }
    }
}

/// FRONTIER PIN (flips to a watertight assert when closed): union of
/// a cylinder boss onto a REFLEX polygon plate. The same union onto a
/// RECT plate is watertight, but on the reflex outline the cellular
/// walker merges the plate-top's boss ring and its 9-vertex outline
/// into ONE 13-half-edge loop (no inner loop) — Euler-Poincaré
/// violation (expected 2, got 4) — and the silent graph-merge
/// fallback emits garbage pcurves (mesh vertices off-outline).
/// Until closed, the laptop-stand base ships without its socket boss
/// (thicker plate, socket cut directly — containment recorded in
/// laptop_stand_api.rs).
#[test]
fn cylinder_boss_on_reflex_plate_union_gap_is_pinned() {
    let outline: Vec<(f64, f64)> = vec![
        (-45.0, 0.0),
        (45.0, 0.0),
        (45.0, 8.0),
        (100.0, 212.0),
        (58.0, 212.0),
        (0.0, 70.0),
        (-58.0, 212.0),
        (-100.0, 212.0),
        (-45.0, 8.0),
    ];
    let plate = Profile::on(Plane::XY)
        .polygon(&outline)
        .extrude(10.0)
        .expect("plate");
    let rect = Profile::on(Plane::XY)
        .rect(200.0, 200.0)
        .extrude(10.0)
        .expect("rect");
    let boss = Profile::on(Plane::xy_at(9.0))
        .circle(30.0)
        .extrude(23.0)
        .expect("boss")
        .translate(0.0, 35.0, 0.0)
        .expect("boss placement");

    // Control: the same boss on a rect plate is closed today.
    let rect_union = rect.union(&boss).expect("rect union");
    assert_eq!(rect_union.open_edge_count().expect("mesh"), 0);

    match plate.union(&boss) {
        Ok(based) => {
            let opens = based.open_edge_count().expect("mesh");
            assert!(
                opens > 0,
                "reflex-plate boss union gap has closed ({opens} opens) — flip \
                 this pin to a watertight assert and restore the base's boss"
            );
        }
        Err(_) => {}
    }
}

/// STEP export parity for plain extruded cylinders: the cap↔wall rim
/// edges must export as SHARED EDGE_CURVEs (a Ø70×60 cylinder has 3
/// edges — two rims + seam — so exactly 3 EDGE_CURVE entities).
/// Closed 2026-07-07: the exporter's face-local edge fork fired for
/// cap rims because (a) pcurve synthesis collapsed closed rims on
/// planes to a degenerate point-line (Plane fast-path keyed on the
/// surface kind instead of the sampled points), and (b) the fork
/// gate let a phase-shifted stored pcurve shadow a synthesized one
/// that validates. OpenCASCADE read each rim as two once-used edges
/// — a bare cylinder and the laptop-stand column imported into
/// FreeCAD as an open pipe shell plus detached socket prisms.
#[test]
fn plain_cylinder_step_exports_shared_rim_edges() {
    let cyl = Profile::on(Plane::XY)
        .circle(35.0)
        .extrude(60.0)
        .expect("cyl");
    let dir = std::env::temp_dir().join("popkern_rim_fork_pin");
    std::fs::create_dir_all(&dir).expect("tmp dir");
    let path = dir.join("cyl.step");
    cyl.export_step(&path).expect("export");
    let text = std::fs::read_to_string(&path).expect("read step");
    let edge_curves = text.matches("EDGE_CURVE").count();
    assert_eq!(
        edge_curves, 3,
        "expected 3 shared EDGE_CURVEs (two rims + seam), got {edge_curves}"
    );
    // The column configuration (sockets both ends) must also export
    // without forked rims: 3 original edges + 2×12 socket-wall edges
    // + 2×12 socket rim edges + 2 floor... assert via the shell tag
    // instead: a forked export downgrades to OPEN_SHELL.
    assert!(
        text.contains("CLOSED_SHELL"),
        "cylinder must export as a CLOSED_SHELL"
    );
}

/// W2a pin: a big slab cutter rotated about X and subtracted from a
/// prism — the cut plane exits through the side faces (the holder's
/// laptop rest plane). Includes a volume oracle and a face-count
/// change assert (a silent no-op returns the input).
#[test]
fn rotated_slab_subtract_cuts_exact_wedge() {
    let prism = Profile::on(Plane::XY)
        .rect(40.0, 60.0)
        .extrude(50.0)
        .expect("prism");
    let tilt = -18.0_f64.to_radians();
    let slab = Profile::on(Plane::XY)
        .rect(400.0, 400.0)
        .extrude(160.0)
        .expect("slab")
        .rotate_x(tilt)
        .expect("rotate")
        .translate(0.0, 0.0, 40.0)
        .expect("lift");
    let wedge = prism.subtract(&slab).expect("wedge");
    assert_eq!(wedge.open_edge_count().expect("mesh"), 0);
    // No face-count no-op guard: a wedge legitimately keeps the
    // prism's 6 faces — the volume oracle below is the no-op check.
    let mesh = wedge
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let volume = mesh_volume(&mesh);
    let t = tilt.tan();
    let n = 4000;
    let mut expected = 0.0;
    for i in 0..n {
        let y = -30.0 + 60.0 * (i as f64 + 0.5) / n as f64;
        expected += (40.0 + y * t).clamp(0.0, 50.0) * 40.0 * (60.0 / n as f64);
    }
    assert!(
        (volume - expected).abs() / expected < 0.001,
        "wedge volume {volume:.0} vs analytic {expected:.0}"
    );
}

/// W2a companion: the HOLDER's deck — a REFLEX plan outline (tower
/// 24-gon arc + backward arms) extruded and cut by the rotated slab.
#[test]
fn rotated_slab_cuts_reflex_holder_deck() {
    let tilt = 18.0_f64.to_radians();
    let mut outline: Vec<(f64, f64)> = Vec::new();
    for k in 0..=12 {
        let a = std::f64::consts::PI * (1.0 - k as f64 / 12.0);
        outline.push((35.0 * a.cos(), 35.0 * a.sin()));
    }
    outline.extend_from_slice(&[
        (22.0, -27.0),
        (95.0, -175.0),
        (55.0, -190.0),
        (0.0, -60.0),
        (-55.0, -190.0),
        (-95.0, -175.0),
        (-22.0, -27.0),
    ]);
    let deck = Profile::on(Plane::XY)
        .polygon(&outline)
        .extrude(75.0)
        .expect("deck");
    let cutter = Profile::on(Plane::XY)
        .rect(400.0, 400.0)
        .extrude(160.0)
        .expect("cutter")
        .rotate_x(-tilt)
        .expect("rot")
        .translate(0.0, -60.0, 30.0)
        .expect("place");
    let cut = deck.subtract(&cutter).expect("rest plane cut");
    assert_eq!(cut.open_edge_count().expect("mesh"), 0);
}

/// W2a pin: BOUNDED rotated-box subtract on a convex prism — the
/// cutter's bottom and side planes cross a wall in an L whose
/// junction lies mid-face. Closed by the T-junction pre-split in
/// split_one_face: a trim curve endpoint landing mid-edge (typically
/// on an edge created by an earlier split of the same face) now
/// splits that boundary edge first; previously the whole boolean
/// hard-failed with "trim curve end point does not match any
/// vertex". The holder's stepped cuts (arms at 18°, tower flat)
/// build on this.
#[test]
fn bounded_rotated_cutter_is_watertight() {
    let prism = Profile::on(Plane::XY)
        .rect(40.0, 60.0)
        .extrude(50.0)
        .expect("prism");
    let cutter = Profile::on(Plane::XY)
        .rect(400.0, 40.0)
        .extrude(160.0)
        .expect("cutter")
        .rotate_x(-18.0_f64.to_radians())
        .expect("rot")
        .translate(0.0, -20.0, 40.0)
        .expect("place");
    let cut = prism.subtract(&cutter).expect("bounded cut");
    assert_eq!(cut.open_edge_count().expect("mesh"), 0);
    assert!(
        cut.face_count() > prism.face_count(),
        "bounded cut must add step faces"
    );
}

/// FRONTIER PIN (flips when closed): union of a rotated pre-filleted
/// box onto a tilted plane (the designed lip attach). The wall walk
/// turns with inconsistent handedness at the two ring junctions —
/// the planar-face angle rule uses a normal GUESSED from the first
/// candidate's cross product (candidate-order-dependent); replacing
/// it with the chart or oriented normal breaks seven tuned boolean
/// suites, so the fix is a walk-convention design fork, not a patch.
/// The holder ships with profile-tab lips instead.
#[test]
fn rotated_lip_union_gap_is_pinned() {
    let prism = Profile::on(Plane::XY)
        .rect(60.0, 60.0)
        .extrude(40.0)
        .expect("prism");
    let tilt = -18.0_f64.to_radians();
    let slab = Profile::on(Plane::XY)
        .rect(400.0, 400.0)
        .extrude(160.0)
        .expect("slab")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, 30.0)
        .expect("lift");
    let wedge = prism.subtract(&slab).expect("wedge");
    let lip = Profile::on(Plane::XY)
        .rect(35.0, 8.0)
        .extrude(17.0)
        .expect("lip")
        .fillet(3.0, |e: EdgeQuery| e.straight().direction(Dir::Z))
        .expect("lip fillets")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, 29.0)
        .expect("place");
    match wedge.union(&lip) {
        Ok(held) => {
            let opens = held.open_edge_count().expect("mesh");
            assert!(
                opens > 0,
                "rotated-lip-union gap has closed ({opens} opens) — flip this \
                 pin and restore the designed pre-filleted lip attach"
            );
        }
        Err(_) => {}
    }
}

/// FRONTIER PIN (flips when closed): bounded rotated cutter whose
/// CORNER EDGE pierces a wall — the cutter's bottom and side trims
/// both DEAD-END at the piercing point, so neither can split first
/// (the T-junction pre-split only splits existing boundary edges,
/// and the deferred-retry loop can't break the mutual wait). Closing
/// this needs the island chainer's open-chain sibling: chain open
/// records that meet at interior dead-ends into ONE boundary-to-
/// boundary polyline before splitting. Blocks the laptop-stand
/// holder's stepped cuts (flat tower + 18° arms + lip tabs).
#[test]
fn bounded_cutter_corner_piercing_gap_is_pinned() {
    let tilt = 18.0_f64.to_radians();
    let mut outline: Vec<(f64, f64)> = Vec::new();
    for k in 0..=12 {
        let a = std::f64::consts::PI * (1.0 - k as f64 / 12.0);
        outline.push((35.0 * a.cos(), 35.0 * a.sin()));
    }
    outline.extend_from_slice(&[
        (22.0, -27.0),
        (95.0, -175.0),
        (55.0, -190.0),
        (0.0, -60.0),
        (-55.0, -190.0),
        (-95.0, -175.0),
        (-22.0, -27.0),
    ]);
    let deck = Profile::on(Plane::XY)
        .polygon(&outline)
        .extrude(75.0)
        .expect("deck");
    let arm_cutter = Profile::on(Plane::XY)
        .rect(120.0, 130.0)
        .extrude(120.0)
        .expect("cutter")
        .rotate_x(-tilt)
        .expect("rot")
        .translate(55.0, -105.0, 42.0)
        .expect("place");
    match deck.subtract(&arm_cutter) {
        Ok(cut) => {
            let opens = cut.open_edge_count().expect("mesh");
            assert!(
                opens > 0,
                "corner-piercing bounded cut has closed ({opens} opens) — \
                 flip this pin and build the holder (task #36)"
            );
        }
        Err(_) => {}
    }
}
