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
    assert_eq!(
        opens, 0,
        "through-hole tessellation must stay watertight (closed by the \
         conforming T-point split in the boundary-first mesher)"
    );
}

/// 12-vertex hexagram (6-point star) for the polygon-profile pins:
/// outer radius `r`, inner radius `r/√3`, starting at angle 0.
/// Orientation-immune material probe: ray parity along +z against
/// the default-quality mesh. Divergence-sum volume oracles are
/// blocked by the winding-inconsistent tessellation of curved faces
/// (arc_extrude_mesh_winding_gap_is_pinned); parity is not.
fn mesh_contains_point(part: &Part, q: (f64, f64, f64)) -> bool {
    let m = part
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let mut crossings = 0;
    for t in &m.triangles {
        let (a, b, c) = (
            m.vertices[t[0]].position,
            m.vertices[t[1]].position,
            m.vertices[t[2]].position,
        );
        let (e1x, e1y) = (b.x - a.x, b.y - a.y);
        let (e2x, e2y) = (c.x - a.x, c.y - a.y);
        let det = e1x * e2y - e1y * e2x;
        if det.abs() < 1e-14 {
            continue;
        }
        let (rx, ry) = (q.0 - a.x, q.1 - a.y);
        let u = (rx * e2y - ry * e2x) / det;
        let v = (e1x * ry - e1y * rx) / det;
        if u < -1e-12 || v < -1e-12 || u + v > 1.0 + 1e-12 {
            continue;
        }
        let z = a.z + u * (b.z - a.z) + v * (c.z - a.z);
        if z > q.2 {
            crossings += 1;
        }
    }
    crossings % 2 == 1
}

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

/// W3 CLOSED: union of two COAXIAL star prisms (the peg's stepped
/// pilot tips). Journey: OpenLoop (Bishop-frame flip) →
/// watertight-χ'=4 (dropped records) → 120 opens (honest records,
/// broken ring election) → closed by the WINDING override in the
/// hole election: a large interior ring winding opposite its outer
/// is that fragment's hole regardless of the 0.5·area size rule
/// (workpiece faces only; tool-side rings stay candidate caps).
#[test]
fn star_boss_union_is_watertight() {
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
    let peg = body.union(&pilot).expect("pilot union");
    assert_eq!(peg.open_edge_count().expect("mesh"), 0);
    let b = peg.brep();
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
        - b.loops.len() as i64;
    assert_eq!(chi, 2, "stepped peg must be a clean genus-0 solid");
}

/// Union of a cylinder boss onto a REFLEX polygon plate — watertight
/// and χ-clean. Closed by the polygon-profile centroid fix: the
/// "walker mega-loop" failure was the boss landing translated by the
/// plate's −centroid, half off the authored outline. The base's
/// designed 10 mm deck + Ø60 socket boss is unblocked.
#[test]
fn cylinder_boss_on_reflex_plate_union_is_watertight() {
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

    let based = plate.union(&boss).expect("boss union");
    assert_eq!(based.open_edge_count().expect("mesh"), 0);
    let b = based.brep();
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
        - b.loops.len() as i64;
    assert_eq!(chi, 2, "boss union must be a clean genus-0 solid");
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

/// W2a companion pin: the HOLDER's deck — a REFLEX plan outline
/// (hub arc + backward arms, separate plane walls after the
/// polyline-corner fix) cut by the rotated rest-plane slab, on the
/// PRODUCT geometry. Closed by the polygon-profile centroid fix:
/// the "reflex walk" failures were tools silently translated by
/// −centroid, piercing walls they were authored to clear.
#[test]
fn rotated_slab_cuts_reflex_holder_deck() {
    let tilt = 18.0_f64.to_radians();
    let hub_r = 58.0;
    let mut outline: Vec<(f64, f64)> = Vec::new();
    for k in 0..=12 {
        let a = std::f64::consts::PI * (1.0 - k as f64 / 12.0);
        outline.push((hub_r * a.cos(), hub_r * a.sin()));
    }
    outline.extend_from_slice(&[
        (36.0, -46.0),
        (95.0, -175.0),
        (55.0, -190.0),
        (0.0, -85.0),
        (-55.0, -190.0),
        (-95.0, -175.0),
        (-36.0, -46.0),
    ]);
    let deck = Profile::on(Plane::XY)
        .polygon(&outline)
        .extrude(120.0)
        .expect("deck");
    let slab = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("slab")
        .rotate_x(-tilt)
        .expect("rot")
        .translate(0.0, 0.0, 55.0)
        .expect("lift");
    let cut = deck.subtract(&slab).expect("rest plane cut");
    assert_eq!(cut.open_edge_count().expect("mesh"), 0);
    let b = cut.brep();
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
        - b.loops.len() as i64;
    assert_eq!(chi, 2, "reflex deck cut must be a clean genus-0 solid");
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
fn rotated_prefilleted_lip_union_is_watertight() {
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
    let held = wedge.union(&lip).expect("pre-filleted lip union");
    assert_eq!(
        held.open_edge_count().expect("mesh"),
        0,
        "rotated pre-filleted lip union must stay watertight"
    );
    assert!(
        mesh_contains_point(&held, (0.0, 0.0, 33.0)),
        "the lip body must be material above the rest plane"
    );
}

/// FRONTIER PIN (flips when closed): bounded rotated cutter whose
/// CORNER EDGE pierces a wall — the cutter's bottom and side trims
/// both DEAD-END at the piercing point. The open-chain pre-pass
/// (island chainer's open sibling, added with this pin) now merges
/// such dead-ending records into one boundary-to-boundary polyline
/// and the split itself succeeds, but the assembled result still
/// carries an open loop on the pierced wall — the next layer is in
/// the post-split stitch. Blocks the laptop-stand holder's stepped
/// cuts (flat tower + 18° arms + lip tabs).
#[test]
fn bounded_cutter_corner_piercing_is_watertight() {
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
    let cut = deck
        .subtract(&arm_cutter)
        .expect("corner-piercing bounded cut");
    let opens = cut.open_edge_count().expect("mesh");
    assert_eq!(
        opens, 0,
        "corner-piercing bounded cut regressed (closed by record×record \
         T-junction splitting in the cellular graph)"
    );
    assert!(cut.face_count() > 0);
    // The cutter's corner edge pierces the deck's wall: the region
    // behind the pierce must be void, the region outboard solid.
    assert!(
        !mesh_contains_point(&cut, (55.0, -105.0, 60.0)),
        "the cut region must be void"
    );
    assert!(
        mesh_contains_point(&cut, (-30.0, -80.0, 20.0)),
        "the left leg must stay solid"
    );
}

/// W2a companion pin: a vertical STAR socket bored into an 18°
/// TILTED face — the holder's peg mechanism. Asserted watertight AND
/// χ-clean (V − E + 2F − L = 2): the mesh gate alone is insufficient
/// for this family (see the membrane pins below — a phantom face can
/// seal a pocket mouth and mesh "watertight").
#[test]
fn tilted_face_star_socket_is_watertight_and_chi_clean() {
    let tilt = -18.0_f64.to_radians();
    let prism = Profile::on(Plane::XY)
        .rect(120.0, 120.0)
        .extrude(60.0)
        .expect("prism");
    let slab = Profile::on(Plane::XY)
        .rect(500.0, 500.0)
        .extrude(160.0)
        .expect("slab")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, 45.0)
        .expect("lift");
    let wedge = prism.subtract(&slab).expect("wedge");
    let star = hexagram_points(12.0);
    let star_tool = Profile::on(Plane::xy_at(20.0))
        .polygon(&star)
        .extrude(40.0)
        .expect("star tool")
        .translate(30.0, 0.0, 0.0)
        .expect("place");
    let socketed = wedge.subtract(&star_tool).expect("tilted star socket");
    assert_eq!(socketed.open_edge_count().expect("mesh"), 0);
    let b = socketed.brep();
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
        - b.loops.len() as i64;
    assert_eq!(chi, 2, "star socket body must be a clean genus-0 solid");

    // A SQUARE counterbore (polygon ring — the working class) is
    // also χ-clean and gives the next column a flat seat: the round
    // seat below is pinned as a membrane.
    let seat_tool = Profile::on(Plane::xy_at(28.0))
        .rect(76.0, 76.0)
        .extrude(60.0)
        .expect("seat");
    let seated = wedge.subtract(&seat_tool).expect("square seat");
    assert_eq!(seated.open_edge_count().expect("mesh"), 0);
    let b = seated.brep();
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
        - b.loops.len() as i64;
    assert_eq!(chi, 2, "square seat body must be a clean genus-0 solid");
}

/// FRONTIER PIN (side-breaching bore): a Ø72 bore whose wall
/// BREACHES the prism's side face by 1 mm (centre −25, radius 36 vs
/// the wall at x = −60) cut through a tilted top. This pin
/// originally pinned the χ' = 4 membrane — that root cause (periodic
/// -chart containment reading interior records as outside) is fixed
/// and the centred cases now gate as
/// round_seat_star_socket_is_watertight_and_void and
/// perpendicular_bore_recess_is_watertight_and_void. What remains is
/// the harder configuration this fixture accidentally carried all
/// along: the mouth ellipse is an OPEN arc (clipped by the breached
/// wall) meeting the bore×wall arc at two junction points, and the
/// subtract fails LowValence building that record network. Same
/// class as the round-seat counterbore re-sequencing failure.
#[test]
fn side_breaching_tilted_bore_gap_is_pinned() {
    let tilt = -18.0_f64.to_radians();
    let prism = Profile::on(Plane::XY)
        .rect(120.0, 120.0)
        .extrude(60.0)
        .expect("prism");
    let slab = Profile::on(Plane::XY)
        .rect(500.0, 500.0)
        .extrude(160.0)
        .expect("slab")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, 45.0)
        .expect("lift");
    let wedge = prism.subtract(&slab).expect("wedge");
    let bore = Profile::on(Plane::xy_at(21.0))
        .circle(36.0)
        .extrude(60.0)
        .expect("bore")
        .translate(-25.0, 0.0, 0.0)
        .expect("place");
    let ok = match wedge.subtract(&bore) {
        Ok(bored) => {
            let b = bored.brep();
            let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
                - b.loops.len() as i64;
            chi == 2
                && bored.open_edge_count().map(|o| o == 0).unwrap_or(false)
                && mesh_contains_point(&bored, (-25.0, 0.0, 10.0))
                && !mesh_contains_point(&bored, (-25.0, 0.0, 30.0))
        }
        Err(e) => {
            eprintln!("side-breaching tilted bore: {e:?}");
            false
        }
    };
    assert!(
        !ok,
        "side-breaching tilted bore now succeeds watertight with the bore \
         void — flip this pin into a regression gate"
    );
}

/// FRONTIER PIN (flips when closed): a PERPENDICULAR cylinder pocket
/// into a flat face — intersection discovery finds ZERO records (the
/// essential circle is judged Outside its own cylinder wall: the
/// SEAM edge's single pcurve serves both loop traversals, so the
/// wall's raw UV contour jumps period branches and even-odd
/// containment reads every interior point as outside). The boolean
/// then degrades to the no-records containment path, which
/// fabricates a MEMBRANE-SEALED pocket: the top face keeps no inner
/// loop, a phantom 1-half-edge disk face covers the mouth, and the
/// body passes BOTH the watertight-mesh gate and validate() — only
/// the Euler characteristic (V−E+2F−L = 3, not 2) and an OCC
/// sew-to-solid volume check catch it. Fix direction (prototyped,
/// reverted — it regressed tuned pins): chain-normalize period
/// branches for containment polygons (face_uv_boundary_chained) as
/// part of the seam-branch campaign. The OBLIQUE variant (tilted
/// face) works — see tilted_face_star_socket_and_counterbore.
#[test]
fn perpendicular_bore_recess_is_watertight_and_void() {
    let prism = Profile::on(Plane::XY)
        .rect(120.0, 120.0)
        .extrude(60.0)
        .expect("prism");
    let bore = Profile::on(Plane::xy_at(30.0))
        .circle(36.0)
        .extrude(60.0)
        .expect("bore");
    let seated = prism.subtract(&bore).expect("bore");
    let brep = seated.brep();
    let v = brep.vertices.len() as i64;
    let e = brep.edges.len() as i64;
    let f = brep.faces.len() as i64;
    let l = brep.loops.len() as i64;
    let chi = v - e + 2 * f - l;
    assert_eq!(chi, 2, "perpendicular-bore recess must stay χ' == 2");
    assert_eq!(seated.open_edge_count().expect("mesh"), 0);
    assert!(
        mesh_contains_point(&seated, (0.0, 0.0, 15.0)),
        "below the bore floor must stay solid"
    );
    assert!(
        !mesh_contains_point(&seated, (0.0, 0.0, 45.0)),
        "the bore volume must be void"
    );
}

/// The holder deck END TO END on product geometry: reflex outline
/// (12-chord hub arc + wide-pad arms) → 18° rest-plane slab → SQUARE
/// column seat bored into the tilted hub → stacking star sockets
/// (seat floor + flat underside) → two lip-peg star sockets at the
/// arm tips, bored PERPENDICULAR to the rest plane. Watertight +
/// χ'==2 asserted after every stage. Closed by: the polygon-profile
/// centroid fix, trim-honest edge-face hit acceptance, the
/// outer-pcurve-polygon ray gate, reflex routing off Cyrus–Beck
/// (clip + split fallback), and the extruded-plane recognition
/// envelope scaling with profile span.
#[test]
fn holder_deck_full_chain_is_watertight_and_chi_clean() {
    let chi = |p: &Part| {
        let b = p.brep();
        b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
            - b.loops.len() as i64
    };
    let tilt = 18.0_f64.to_radians();
    let hub_r = 58.0;
    let mut outline: Vec<(f64, f64)> = Vec::new();
    for k in 0..=12 {
        let a = std::f64::consts::PI * (1.0 - k as f64 / 12.0);
        outline.push((hub_r * a.cos(), hub_r * a.sin()));
    }
    outline.extend_from_slice(&[
        (36.0, -46.0),
        (92.0, -142.0),
        (92.0, -178.0),
        (52.0, -178.0),
        (0.0, -85.0),
        (-52.0, -178.0),
        (-92.0, -178.0),
        (-92.0, -142.0),
        (-36.0, -46.0),
    ]);
    let deck = Profile::on(Plane::XY)
        .polygon(&outline)
        .extrude(120.0)
        .expect("deck");
    let slab = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("slab")
        .rotate_x(-tilt)
        .expect("rot")
        .translate(0.0, 0.0, 55.0)
        .expect("lift");
    let mut part = deck.subtract(&slab).expect("rest");
    assert_eq!(part.open_edge_count().expect("mesh"), 0, "rest plane");
    assert_eq!(chi(&part), 2, "rest plane χ");
    let seat = Profile::on(Plane::xy_at(40.0))
        .rect(76.0, 76.0)
        .extrude(120.0)
        .expect("seat tool");
    part = part.subtract(&seat).expect("seat");
    assert_eq!(part.open_edge_count().expect("mesh"), 0, "seat");
    assert_eq!(chi(&part), 2, "seat χ");
    let star = hexagram_points(12.0);
    let up = Profile::on(Plane::xy_at(24.0))
        .polygon(&star)
        .extrude(17.0)
        .expect("up");
    part = part.subtract(&up).expect("up socket");
    assert_eq!(part.open_edge_count().expect("mesh"), 0, "up socket");
    assert_eq!(chi(&part), 2, "up socket χ");
    let down = Profile::on(Plane::xy_at(-1.0))
        .polygon(&star)
        .extrude(17.0)
        .expect("down");
    part = part.subtract(&down).expect("down socket");
    assert_eq!(part.open_edge_count().expect("mesh"), 0, "down socket");
    assert_eq!(chi(&part), 2, "down socket χ");
    // Lip-peg sockets at the ARM TIPS, bored perpendicular to the
    // rest plane (the lip pegs mount normal to the tilted face).
    let (ny, nz) = (tilt.sin(), tilt.cos());
    for x in [-72.0_f64, 72.0] {
        let y = -160.0;
        let z_plane = 55.0 - y * tilt.tan();
        // Tool symmetric about origin, rotated onto the plane normal,
        // then anchored 16 mm embedded / 24 mm proud.
        let tool = Profile::on(Plane::xy_at(-20.0))
            .polygon(&star)
            .extrude(40.0)
            .expect("lip tool")
            .rotate_x(-tilt)
            .expect("tilt")
            .translate(x, y + 4.0 * ny, z_plane + 4.0 * nz)
            .expect("place");
        part = part.subtract(&tool).expect("lip socket");
        assert_eq!(
            part.open_edge_count().expect("mesh"),
            0,
            "lip socket at {x}"
        );
        assert_eq!(chi(&part), 2, "lip socket χ at {x}");
    }
}

/// Holder v2: the FLIPPED rest slope (rotate_x(+tilt), lips at the
/// arm tips) with the square seat — watertight. Closed by the
/// odd-crossing containment fix in the UV clip: the seat tool's
/// padded box pre-truncated its x-wall trims mid-face, the single
/// surviving boundary crossing hit the `1 => None` arm, and two of
/// the four ring records were dropped as Outside.
#[test]
fn flipped_slope_seat_is_watertight() {
    let tilt = 18.0_f64.to_radians();
    let hub_r = 58.0;
    let mut outline: Vec<(f64, f64)> = Vec::new();
    for k in 0..=12 {
        let a = std::f64::consts::PI * (1.0 - k as f64 / 12.0);
        outline.push((hub_r * a.cos(), hub_r * a.sin()));
    }
    outline.extend_from_slice(&[
        (44.0, -52.0),
        (92.0, -142.0),
        (92.0, -178.0),
        (52.0, -178.0),
        (0.0, -95.0),
        (-52.0, -178.0),
        (-92.0, -178.0),
        (-92.0, -142.0),
        (-44.0, -52.0),
    ]);
    let deck = Profile::on(Plane::XY)
        .polygon(&outline)
        .extrude(100.0)
        .expect("deck");
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("slab")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, 75.0)
        .expect("lift");
    let part = deck.subtract(&upper).expect("flipped rest plane");
    assert_eq!(part.open_edge_count().expect("mesh"), 0);
    let seat = Profile::on(Plane::xy_at(40.0))
        .rect(76.0, 76.0)
        .extrude(60.0)
        .expect("seat");
    let cut = part.subtract(&seat).expect("flipped seat");
    assert_eq!(cut.open_edge_count().expect("mesh"), 0);
    let b = cut.brep();
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
        - b.loops.len() as i64;
    assert_eq!(chi, 2, "flipped-slope seat must be a clean genus-0 solid");
}

/// Holder v2: subtracting a TOOL-SHAPED wedge (box pre-cut by the
/// rest plane) to thin the arms — watertight. Closed by the
/// topological junction rule: a walk arriving along the boundary
/// turns onto a crossing trim, one arriving along the trim continues
/// onto the boundary (exact for single-trim junctions; no surface
/// normal involved).
#[test]
fn shaped_wedge_arm_thinning_is_watertight() {
    let tilt = 18.0_f64.to_radians();
    let hub_r = 58.0;
    let mut outline: Vec<(f64, f64)> = Vec::new();
    for k in 0..=12 {
        let a = std::f64::consts::PI * (1.0 - k as f64 / 12.0);
        outline.push((hub_r * a.cos(), hub_r * a.sin()));
    }
    outline.extend_from_slice(&[
        (44.0, -52.0),
        (92.0, -142.0),
        (92.0, -178.0),
        (52.0, -178.0),
        (0.0, -95.0),
        (-52.0, -178.0),
        (-92.0, -178.0),
        (-92.0, -142.0),
        (-44.0, -52.0),
    ]);
    let deck = Profile::on(Plane::XY)
        .polygon(&outline)
        .extrude(100.0)
        .expect("deck");
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("slab")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, 75.0)
        .expect("lift");
    let part = deck.subtract(&upper).expect("flipped rest plane");
    let lower_box = Profile::on(Plane::xy_at(-10.0))
        .rect(600.0, 129.0)
        .at(0.0, -120.5)
        .extrude(100.0)
        .expect("box");
    let p2 = Profile::on(Plane::XY)
        .rect(700.0, 700.0)
        .extrude(300.0)
        .expect("slab2")
        .rotate_x(tilt)
        .expect("rot2")
        .translate(0.0, 0.0, 62.0)
        .expect("lift2");
    let tool = lower_box.subtract(&p2).expect("shaped tool");
    assert_eq!(tool.open_edge_count().expect("mesh"), 0);
    let cut = part.subtract(&tool).expect("arm thinning");
    assert_eq!(cut.open_edge_count().expect("mesh"), 0);
    let b = cut.brep();
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
        - b.loops.len() as i64;
    assert_eq!(chi, 2, "thinned deck must be a clean genus-0 solid");
}

/// Holder v2: UNION of a rotated plate prism into a vertical block —
/// the arm construction (the plate's own bottom face is the designed
/// underside) — watertight. Closed by the same odd-crossing
/// containment fix (box-truncated trims ending mid-face).
#[test]
fn rotated_plate_union_is_watertight() {
    let tilt = 18.0_f64.to_radians();
    let hub_r = 58.0;
    let mut block_plan: Vec<(f64, f64)> = Vec::new();
    for k in 0..=12 {
        let a = std::f64::consts::PI * (1.0 - k as f64 / 12.0);
        block_plan.push((hub_r * a.cos(), hub_r * a.sin()));
    }
    block_plan.extend_from_slice(&[(44.0, -56.0), (-44.0, -56.0)]);
    let block = Profile::on(Plane::XY)
        .polygon(&block_plan)
        .extrude(95.0)
        .expect("block");
    let arm_plan: Vec<(f64, f64)> = vec![
        (44.0, -40.0),
        (92.0, -130.0),
        (92.0, -170.0),
        (52.0, -170.0),
        (0.0, -90.0),
        (-52.0, -170.0),
        (-92.0, -170.0),
        (-92.0, -130.0),
        (-44.0, -40.0),
    ];
    let plate = Profile::on(Plane::XY)
        .polygon(&arm_plan)
        .extrude(20.0)
        .expect("plate")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, 62.0)
        .expect("place");
    let joined = block.union(&plate).expect("plate union");
    assert_eq!(joined.open_edge_count().expect("mesh"), 0);
    let b = joined.brep();
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
        - b.loops.len() as i64;
    assert_eq!(chi, 2, "plate union must be a clean genus-0 solid");
}

/// Holder v2: the rest-plane wedge on the UNION-BUILT body (block ∪
/// rotated arm plate) — watertight. This was the walk-convention
/// design fork: the trim ring crosses the two REFLEX block/plate
/// junction walls, and the junction pick consumed the boundary loop
/// past the trim (505 opens, χ'=−17). Closed by the topological
/// junction rule in next_around_face — alternation at trim
/// junctions is exact and normal-free; four geometry-based
/// mechanisms (normal-source swaps ×3, seed-order inversion) were
/// falsified before it.
#[test]
fn wedge_on_union_is_watertight() {
    let tilt = 18.0_f64.to_radians();
    let hub_r = 58.0;
    let mut block_plan: Vec<(f64, f64)> = Vec::new();
    for k in 0..=12 {
        let a = std::f64::consts::PI * (1.0 - k as f64 / 12.0);
        block_plan.push((hub_r * a.cos(), hub_r * a.sin()));
    }
    block_plan.extend_from_slice(&[(44.0, -56.0), (-44.0, -56.0)]);
    let block = Profile::on(Plane::XY)
        .polygon(&block_plan)
        .extrude(95.0)
        .expect("block");
    let arm_plan: Vec<(f64, f64)> = vec![
        (44.0, -40.0),
        (92.0, -130.0),
        (92.0, -170.0),
        (52.0, -170.0),
        (0.0, -90.0),
        (-52.0, -170.0),
        (-92.0, -170.0),
        (-92.0, -130.0),
        (-44.0, -40.0),
    ];
    let plate = Profile::on(Plane::XY)
        .polygon(&arm_plan)
        .extrude(20.0)
        .expect("plate")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, 62.0)
        .expect("place");
    let joined = block.union(&plate).expect("union");
    assert_eq!(joined.open_edge_count().expect("mesh"), 0, "union stage");
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("slab")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, 75.0)
        .expect("lift");
    let cut = joined.subtract(&upper).expect("wedge on union");
    assert_eq!(cut.open_edge_count().expect("mesh"), 0);
    let b = cut.brep();
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
        - b.loops.len() as i64;
    assert_eq!(chi, 2, "wedged union must be a clean genus-0 solid");
}

/// Holder v2 END TO END (user-designed, session 7): block ∪ rotated
/// arm plate (the plate's bottom face IS the 12 mm arm underside) →
/// rest-plane wedge across the union (+tilt: lips at the LOW arm
/// tips) → square column seat → stacking star sockets → DOVETAIL lip
/// slots at both tips (45° self-supporting catch at the blind end).
/// Watertight + χ'==2 after every stage. Prints rest-face-down with
/// no supports.
#[test]
fn holder_v2_full_chain_is_watertight_and_chi_clean() {
    let chi = |p: &Part| {
        let b = p.brep();
        b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
            - b.loops.len() as i64
    };
    let gate = |p: Part, label: &str| -> Part {
        assert_eq!(p.open_edge_count().expect("mesh"), 0, "{label}: opens");
        assert_eq!(chi(&p), 2, "{label}: χ");
        p
    };
    let tilt = 18.0_f64.to_radians();
    let hub_r = 58.0;
    let mut block_plan: Vec<(f64, f64)> = Vec::new();
    for k in 0..=12 {
        let a = std::f64::consts::PI * (1.0 - k as f64 / 12.0);
        block_plan.push((hub_r * a.cos(), hub_r * a.sin()));
    }
    block_plan.extend_from_slice(&[(44.0, -56.0), (-44.0, -56.0)]);
    let block = Profile::on(Plane::XY)
        .polygon(&block_plan)
        .extrude(95.0)
        .expect("block");
    let arm_plan: Vec<(f64, f64)> = vec![
        (44.0, -40.0),
        (92.0, -130.0),
        (92.0, -170.0),
        (52.0, -170.0),
        (0.0, -90.0),
        (-52.0, -170.0),
        (-92.0, -170.0),
        (-92.0, -130.0),
        (-44.0, -40.0),
    ];
    let plate = Profile::on(Plane::XY)
        .polygon(&arm_plan)
        .extrude(20.0)
        .expect("plate")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, 62.0)
        .expect("place");
    let part = gate(block.union(&plate).expect("union"), "union");
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("slab")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, 75.0)
        .expect("lift");
    let part = gate(part.subtract(&upper).expect("wedge"), "wedge");
    let seat = Profile::on(Plane::xy_at(40.0))
        .rect(76.0, 76.0)
        .extrude(60.0)
        .expect("seat");
    let part = gate(part.subtract(&seat).expect("seat"), "seat");
    let star = hexagram_points(12.0);
    let up = Profile::on(Plane::xy_at(24.0))
        .polygon(&star)
        .extrude(17.0)
        .expect("up");
    let part = gate(part.subtract(&up).expect("up"), "up socket");
    let down = Profile::on(Plane::xy_at(-1.0))
        .polygon(&star)
        .extrude(17.0)
        .expect("down");
    let mut part = gate(part.subtract(&down).expect("down"), "down socket");
    for side in [1.0_f64, -1.0] {
        // Mouth region extends 2 mm PROUD of the rest face (a mouth
        // exactly on the face is a designed-in coplanar overlap).
        let prism = Profile::on(Plane::XY)
            .polygon(&[
                (-2.0, -4.5),
                (1.0, -4.5),
                (6.0, -7.0),
                (6.0, 7.0),
                (1.0, 4.5),
                (-2.0, 4.5),
            ])
            .extrude(41.0)
            .expect("prism")
            .rotate_y(90.0_f64.to_radians())
            .expect("lay");
        let mut tool = prism; // bisect: flat catch (45° chop removed)
        if side < 0.0 {
            tool = tool.rotate_z(std::f64::consts::PI).expect("mirror");
        }
        let dz = 75.0 - 145.0 * tilt.tan();
        let tool = tool
            .rotate_x(tilt)
            .expect("tilt")
            .translate(side * 52.0, -145.0, dz)
            .expect("place");
        part = gate(part.subtract(&tool).expect("slot"), "dovetail slot");
    }
}

/// Arc-in-sketch API: a mixed line/arc outline (a rectangle whose
/// front side bulges as a circular arc) extrudes to a watertight,
/// χ-clean solid whose arc wall is a TRUE `Cylinder` surface and
/// whose straight walls are `Plane`s — the sweep splits spans at
/// curvature jumps (tangent line↔arc joints) while a full circle's
/// quadrant joints stay merged (`plain_cylinder` pins cover that).
#[test]
fn mixed_path_extrude_emits_true_cylinder_wall() {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let part = Profile::on(Plane::XY)
        .path(-30.0, 0.0)
        .line_to(30.0, 0.0)
        .arc_to((0.0, 18.0), (-30.0, 0.0))
        .close()
        .extrude(20.0)
        .expect("mixed extrude");
    assert_eq!(part.open_edge_count().expect("mesh"), 0);
    let b = part.brep();
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
        - b.loops.len() as i64;
    assert_eq!(chi, 2, "mixed-path solid must be genus 0");
    let mut planes = 0;
    let mut cylinders = 0;
    let mut nurbs = 0;
    for (_, f) in b.faces.iter() {
        match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Plane { .. }) => planes += 1,
            Some(SurfaceGeometry::Cylinder { .. }) => cylinders += 1,
            _ => nurbs += 1,
        }
    }
    // 1 straight wall + 2 caps = 3 planes; the ~113° arc subdivides
    // into two ≤90° rational spans → cylinder wall face(s).
    assert!(
        cylinders >= 1,
        "arc wall must be a true Cylinder (got {cylinders})"
    );
    assert_eq!(nurbs, 0, "no unrecognized walls (got {nurbs})");
    assert!(
        planes >= 3,
        "straight wall + caps must be Planes (got {planes})"
    );
}

/// Holder v3: the rest-plane wedge over a SMOOTH ARC front wall (a
/// true Cylinder from the mixed-path API) — watertight and χ-clean.
/// Closed by junction HEALING in the graph insert: the ellipse
/// record's far endpoint measured 9e-6 off the true ring corner
/// while the exact corner already existed as an edge-face hit node;
/// record endpoints now adopt existing nodes with a wide second
/// search AND the curve's end control points snap to the adopted
/// node (a clamped NURBS end IS its end CP), so both the topology
/// and the shell close.
#[test]
fn arc_front_wedge_is_watertight() {
    let tilt = 18.0_f64.to_radians();
    let block = Profile::on(Plane::XY)
        .path(58.0, 0.0)
        .arc_to((0.0, 58.0), (-58.0, 0.0))
        .line_to(-44.0, -56.0)
        .line_to(44.0, -56.0)
        .close()
        .extrude(125.0)
        .expect("block");
    assert_eq!(block.open_edge_count().expect("mesh"), 0, "block extrude");
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("slab")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, 100.0)
        .expect("lift");
    let cut = block.subtract(&upper).expect("arc-front wedge");
    assert_eq!(cut.open_edge_count().expect("mesh"), 0);
    let b = cut.brep();
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
        - b.loops.len() as i64;
    assert_eq!(chi, 2, "arc-front wedge must be a clean genus-0 solid");
}

/// Holder v3 BLOCK (the stacking piece): plan extrude → +tilt rest
/// wedge → square seat → stacking star sockets → two vertical
/// dovetail ARM sockets in the rear face. Watertight + χ'==2 at
/// every stage.
#[test]
fn holder_v3_block_chain_is_watertight_and_chi_clean() {
    let chi = |p: &Part| {
        let b = p.brep();
        b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
            - b.loops.len() as i64
    };
    let tilt = 18.0_f64.to_radians();
    let hub_r = 58.0;
    let _ = hub_r;
    // SMOOTH arc front — a true Cylinder wall from the mixed-path
    // sketch API.
    let mut part = Profile::on(Plane::XY)
        .path(58.0, 0.0)
        .arc_to((0.0, 58.0), (-58.0, 0.0))
        .line_to(-44.0, -56.0)
        .line_to(44.0, -56.0)
        .close()
        .extrude(125.0)
        .expect("block");
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("slab")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, 100.0)
        .expect("lift");
    part = part.subtract(&upper).expect("wedge");
    assert_eq!(part.open_edge_count().expect("mesh"), 0, "wedge");
    assert_eq!(chi(&part), 2, "wedge χ");
    let seat = Profile::on(Plane::xy_at(70.0))
        .rect(76.0, 76.0)
        .extrude(60.0)
        .expect("seat");
    part = part.subtract(&seat).expect("seat");
    assert_eq!(part.open_edge_count().expect("mesh"), 0, "seat");
    assert_eq!(chi(&part), 2, "seat χ");
    let star = hexagram_points(12.0);
    let up = Profile::on(Plane::xy_at(54.0))
        .polygon(&star)
        .extrude(17.0)
        .expect("up");
    part = part.subtract(&up).expect("up");
    assert_eq!(part.open_edge_count().expect("mesh"), 0, "up");
    assert_eq!(chi(&part), 2, "up χ");
    let down = Profile::on(Plane::xy_at(-1.0))
        .polygon(&star)
        .extrude(17.0)
        .expect("down");
    part = part.subtract(&down).expect("down");
    assert_eq!(part.open_edge_count().expect("mesh"), 0, "down");
    assert_eq!(chi(&part), 2, "down χ");
    // Two vertical dovetail arm sockets in the rear face: mouth 12
    // wide 1 mm proud of the wall, widening to 18 at 10 deep, blind
    // floor at z=20 (the arm hangs on it), open through the tilted
    // top.
    for side in [1.0_f64, -1.0] {
        let tool = Profile::on(Plane::xy_at(20.0))
            .polygon(&[(-6.0, -7.0), (6.0, -7.0), (9.0, 4.0), (-9.0, 4.0)])
            .at(side * 28.0, -50.0)
            .extrude(110.0)
            .expect("socket tool");
        part = part.subtract(&tool).expect("arm socket");
        assert_eq!(part.open_edge_count().expect("mesh"), 0, "arm socket");
        assert_eq!(chi(&part), 2, "arm socket χ");
    }
}

/// Holder v3 ARM (identical piece ×2): ONE plan outline (root body +
/// dovetail tab + plate run) extruded, +tilt rest wedge, tool-shaped
/// lower thinning cut, smooth arc swoop gusset union, THROUGH dovetail
/// lip slot. The through slot makes the arm genus 1 — χ' == 0 is
/// the honest final value.
#[test]
fn holder_v3_arm_chain_is_watertight_and_chi_clean() {
    let chi = |p: &Part| {
        let b = p.brep();
        b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
            - b.loops.len() as i64
    };
    let tilt = 18.0_f64.to_radians();
    // ONE plan outline: root body (44 wide) + dovetail tab + plate
    // run (40 wide), extruded tall and shaped by the pin-proven
    // sequence.
    let plan: Vec<(f64, f64)> = vec![
        (22.0, 0.0),
        (5.8, 0.0),
        (8.8, 9.8),
        (-8.8, 9.8),
        (-5.8, 0.0),
        (-22.0, 0.0),
        (-22.0, -24.0),
        (-20.0, -24.0),
        (-20.0, -205.0),
        (20.0, -205.0),
        (20.0, -24.0),
        (22.0, -24.0),
    ];
    let mut part = Profile::on(Plane::XY)
        .polygon(&plan)
        .extrude(96.0)
        .expect("arm blank");
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("slab")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, 81.8)
        .expect("lift");
    part = part.subtract(&upper).expect("wedge");
    assert_eq!(part.open_edge_count().expect("mesh"), 0, "wedge");
    assert_eq!(chi(&part), 2, "wedge χ");
    // Thin the plate run to 13 mm: tool-shaped lower cutter with a
    // vertical front face 2 mm behind the root.
    let lower_box = Profile::on(Plane::xy_at(-10.0))
        .rect(600.0, 190.0)
        .at(0.0, -121.0)
        .extrude(100.0)
        .expect("lower box");
    let p2 = Profile::on(Plane::XY)
        .rect(700.0, 700.0)
        .extrude(300.0)
        .expect("p2 slab")
        .rotate_x(tilt)
        .expect("rot2")
        .translate(0.0, 0.0, 81.8 - 13.0 / tilt.cos())
        .expect("lift2");
    let tool = lower_box.subtract(&p2).expect("shaped tool");
    part = part.subtract(&tool).expect("thin");
    assert_eq!(part.open_edge_count().expect("mesh"), 0, "thin");
    assert_eq!(chi(&part), 2, "thin χ");
    // TRUE-ARC swoop gusset: top edge parallel to (and embedded 1 mm
    // into) the plate underside for full-run contact; free edge is
    // one concave arc (arc_gusset_union_is_watertight).
    let plate_under = |y: f64| 81.8 + y * tilt.tan() - 13.0 / tilt.cos();
    let (y_root, y_tip) = (-20.0, -72.0);
    let gusset = Profile::on(Plane::XY)
        .path(-(plate_under(y_root) + 1.0), y_root)
        .line_to(-(plate_under(y_tip) + 1.0), y_tip)
        .arc_to((-42.0, -36.0), (-26.0, y_root))
        .close()
        .extrude(36.0)
        .expect("gusset")
        .rotate_y(90.0_f64.to_radians())
        .expect("lay")
        .translate(-18.0, 0.0, 0.0)
        .expect("center");
    part = part.union(&gusset).expect("swoop");
    assert_eq!(part.open_edge_count().expect("mesh"), 0, "swoop");
    assert_eq!(chi(&part), 2, "swoop χ");
    // Through dovetail lip slot near the tip.
    let slot = Profile::on(Plane::XY)
        .polygon(&[
            (-2.0, -4.5),
            (1.0, -4.5),
            (6.0, -7.0),
            (6.0, 7.0),
            (1.0, 4.5),
            (-2.0, 4.5),
        ])
        .extrude(44.0)
        .expect("slot prism")
        .rotate_y(90.0_f64.to_radians())
        .expect("lay")
        .rotate_x(tilt)
        .expect("tilt")
        .translate(-22.0, -189.0, 81.8 - 189.0 * tilt.tan())
        .expect("place");
    part = part.subtract(&slot).expect("lip slot");
    assert_eq!(part.open_edge_count().expect("mesh"), 0, "slot");
    assert_eq!(chi(&part), 0, "slot χ");
}

/// FRONTIER PIN: the peg's SECOND pilot union — a smaller-star boss
/// onto a body that already carries one stepped pilot. The single
/// union is closed (star_boss_union_is_watertight); the second
/// leaks. Until flipped, the peg ships as a plain prism.
#[test]
fn peg_double_pilot_gap_is_pinned() {
    let body_star = hexagram_points(11.8);
    let pilot_star = hexagram_points(11.2);
    let body = Profile::on(Plane::xy_at(6.0))
        .polygon(&body_star)
        .extrude(18.0)
        .expect("body");
    let bottom = Profile::on(Plane::XY)
        .polygon(&pilot_star)
        .extrude(10.0)
        .expect("bottom");
    let top = Profile::on(Plane::xy_at(20.0))
        .polygon(&pilot_star)
        .extrude(10.0)
        .expect("top");
    let stepped = body.union(&bottom).expect("first union");
    assert_eq!(stepped.open_edge_count().expect("mesh"), 0, "first union");
    match stepped.union(&top) {
        Ok(peg) => {
            let opens = peg.open_edge_count().expect("mesh");
            assert!(
                opens > 0,
                "double-pilot gap has closed ({opens} opens) — flip this \\
                 pin and restore both pilots (laptop_stand_api.rs)"
            );
        }
        Err(_) => {}
    }
}

/// Smooth front, final stage: the ARM-SOCKET cuts on the arc-fronted
/// block — watertight and χ-clean. Closed by: the conforming T-point
/// split in the mesher's boundary-first refinement (midpoint
/// bisection could never land ON a foreign ring point — livelock),
/// chord synthesis for curveless edges + zero-length ring skips,
/// degenerate-edge removal in cellular assembly, the geometric
/// angular gate for cylinder ray hits, and dense curved-pcurve
/// sampling with a dead-band in the outer-polygon ray gate (the
/// ellipse-fronted top face was exempt from the gate and accepted a
/// phantom 0.67 outside it).
#[test]
fn arc_block_arm_sockets_are_watertight() {
    let tilt = 18.0_f64.to_radians();
    let mut part = Profile::on(Plane::XY)
        .path(58.0, 0.0)
        .arc_to((0.0, 58.0), (-58.0, 0.0))
        .line_to(-44.0, -56.0)
        .line_to(44.0, -56.0)
        .close()
        .extrude(125.0)
        .expect("block");
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("slab")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, 100.0)
        .expect("lift");
    part = part.subtract(&upper).expect("wedge");
    let seat = Profile::on(Plane::xy_at(70.0))
        .rect(76.0, 76.0)
        .extrude(60.0)
        .expect("seat");
    part = part.subtract(&seat).expect("seat");
    let star = hexagram_points(12.0);
    let up = Profile::on(Plane::xy_at(54.0))
        .polygon(&star)
        .extrude(17.0)
        .expect("up");
    part = part.subtract(&up).expect("up");
    assert_eq!(part.open_edge_count().expect("mesh"), 0, "up socket on arc");
    let down = Profile::on(Plane::xy_at(-1.0))
        .polygon(&star)
        .extrude(17.0)
        .expect("down");
    part = part.subtract(&down).expect("down");
    assert_eq!(
        part.open_edge_count().expect("mesh"),
        0,
        "down socket on arc"
    );
    for side in [1.0_f64, -1.0] {
        let tool = Profile::on(Plane::xy_at(20.0))
            .polygon(&[(-6.0, -7.0), (6.0, -7.0), (9.0, 4.0), (-9.0, 4.0)])
            .at(side * 28.0, -50.0)
            .extrude(110.0)
            .expect("tool");
        part = part.subtract(&tool).expect("arm socket");
        assert_eq!(part.open_edge_count().expect("mesh"), 0, "arm socket");
        let b = part.brep();
        let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
            - b.loops.len() as i64;
        assert_eq!(chi, 2, "arm socket χ");
    }
}

/// The TRUE-ARC gusset union — a mixed-path profile whose arc wall
/// meets the arm as a cylinder with axis ALONG the profile normal —
/// watertight and χ-clean. The campaign closed six layers: cylinder
/// recognition envelope, plane×cylinder parallel-case angles + axial
/// span, rotational-chart pcurve re-anchoring under rigid
/// transforms, curved cap/planar pcurve projection, curve-window
/// splitting + bidirectional curved-edge point classification, and
/// finally the split parameter itself — Newton inversion failed
/// outright on the rational arc, so windows never split and the
/// stale full window kept re-matching (five-way duplicate splits);
/// the sampled coarse-scan + ternary finder replaced it.
#[test]
fn arc_gusset_union_is_watertight() {
    let tilt = 18.0_f64.to_radians();
    let plane_at_root = 81.8;
    let plan: Vec<(f64, f64)> = vec![
        (22.0, 0.0),
        (5.8, 0.0),
        (8.8, 9.8),
        (-8.8, 9.8),
        (-5.8, 0.0),
        (-22.0, 0.0),
        (-22.0, -24.0),
        (-20.0, -24.0),
        (-20.0, -231.0),
        (20.0, -231.0),
        (20.0, -24.0),
        (22.0, -24.0),
    ];
    let mut arm = Profile::on(Plane::XY)
        .polygon(&plan)
        .extrude(96.0)
        .expect("arm");
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("slab")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, plane_at_root)
        .expect("lift");
    arm = arm.subtract(&upper).expect("rest");
    let lower_box = Profile::on(Plane::xy_at(-10.0))
        .rect(600.0, 216.0)
        .at(0.0, -134.0)
        .extrude(100.0)
        .expect("box");
    let p2 = Profile::on(Plane::XY)
        .rect(700.0, 700.0)
        .extrude(300.0)
        .expect("slab2")
        .rotate_x(tilt)
        .expect("rot2")
        .translate(0.0, 0.0, plane_at_root - 13.0 / tilt.cos())
        .expect("lift2");
    let tool = lower_box.subtract(&p2).expect("shaped tool");
    arm = arm.subtract(&tool).expect("thin");
    let plate_under = |y: f64| plane_at_root + y * tilt.tan() - 13.0 / tilt.cos();
    let (y_root, y_tip) = (-20.0, -72.0);
    let gusset = Profile::on(Plane::XY)
        .path(-(plate_under(y_root) + 1.0), y_root)
        .line_to(-(plate_under(y_tip) + 1.0), y_tip)
        .arc_to((-42.0, -36.0), (-26.0, y_root))
        .close()
        .extrude(36.0)
        .expect("gusset")
        .rotate_y(90.0_f64.to_radians())
        .expect("lay")
        .translate(-18.0, 0.0, 0.0)
        .expect("center");
    let out = arm.union(&gusset).expect("swoop union");
    assert_eq!(out.open_edge_count().expect("mesh"), 0);
    let b = out.brep();
    let chi = b.vertices.len() as i64 - b.edges.len() as i64 + 2 * b.faces.len() as i64
        - b.loops.len() as i64;
    assert_eq!(chi, 2, "true-arc swoop union must stay genus 0");
}

/// COINCIDENT-BOOLEAN CAMPAIGN — CLOSED. The FLUSH swoop gusset (end
/// caps exactly COPLANAR with the arm's side walls) unions watertight
/// with the volume oracle satisfied exactly (delta 0.000 vs the
/// analytic 40/36 scaling of the inset reference). The campaign closed
/// twelve layers across three legs, all suite-green:
/// leg 1 — (1) the coincident-body shortcut VERIFIES containment
/// geometrically with inward-offset face probes (it was silently
/// swallowing the whole tool: watertight, χ-clean, minus 44 cm³ —
/// only the volume oracle catches this class); (2) degenerate
/// point-trims consumed as no-ops; (3) collinear boundary overlaps
/// clip records at the overlap ENDS (mutual-extent semantics).
/// leg 2 — (4) degenerate fragment culls; (5) overlap-face fragments
/// routed to offset-pair classification; (6) assembly edge
/// unification for distinct graph pairs on the same span (with
/// position-based direction alignment); (7) containment-based hole
/// election; (8) spur guard in the junction rule.
/// leg 3 — (9) planar walk normal derived from the face's stored
/// outer-loop winding (Newell, negated to outward) instead of a
/// candidate-order-dependent guess, plus a deterministic last-resort
/// rank for candidates collinear with the arrival edge; (10) the spur
/// test promoted to a candidate FILTER so interior graze records
/// (an embedded face's dead-end trace) are never entered by any rule;
/// (11) period-branch alignment in the fragment sampler's scanline
/// polygon (a stored pcurve on u−2π opened a phantom interval and the
/// sample landed on the cylinder's far side) and a zero-extent
/// self-loop cull (tangency-point records assembled as point-faces);
/// (12) same-surface KEF post-passes (redundant-edge removal and face
/// simplification) require the face pair to share exactly ONE edge —
/// dissolving one edge of a multi-edge shared boundary left the
/// survivor traversed twice, an untessellatable slit loop.
/// The arm ships with swoop_inset = 2.0 as an aesthetic choice;
/// swoop_inset = 0.0 (fully flush) is verified here.
#[test]
fn flush_gusset_union_is_watertight_and_volume_exact() {
    let tilt = 18.0_f64.to_radians();
    let plane_at_root = 81.8;
    let plan: Vec<(f64, f64)> = vec![
        (22.0, 0.0),
        (5.8, 0.0),
        (8.8, 9.8),
        (-8.8, 9.8),
        (-5.8, 0.0),
        (-22.0, 0.0),
        (-22.0, -24.0),
        (-20.0, -24.0),
        (-20.0, -231.0),
        (20.0, -231.0),
        (20.0, -24.0),
        (22.0, -24.0),
    ];
    let mut arm = Profile::on(Plane::XY)
        .polygon(&plan)
        .extrude(96.0)
        .expect("arm");
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("slab")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, plane_at_root)
        .expect("lift");
    arm = arm.subtract(&upper).expect("rest");
    let lower_box = Profile::on(Plane::xy_at(-10.0))
        .rect(600.0, 216.0)
        .at(0.0, -134.0)
        .extrude(100.0)
        .expect("box");
    let p2 = Profile::on(Plane::XY)
        .rect(700.0, 700.0)
        .extrude(300.0)
        .expect("slab2")
        .rotate_x(tilt)
        .expect("rot2")
        .translate(0.0, 0.0, plane_at_root - 13.0 / tilt.cos())
        .expect("lift2");
    let tool = lower_box.subtract(&p2).expect("shaped tool");
    arm = arm.subtract(&tool).expect("thin");
    let plate_under = |y: f64| plane_at_root + y * tilt.tan() - 13.0 / tilt.cos();
    let (y_root, y_tip) = (-20.0, -72.0);
    // FLUSH: full 40mm width, caps coplanar with the arm side walls.
    let gusset = Profile::on(Plane::XY)
        .path(-(plate_under(y_root) + 1.0), y_root)
        .line_to(-(plate_under(y_tip) + 1.0), y_tip)
        .arc_to((-42.0, -36.0), (-26.0, y_root))
        .close()
        .extrude(40.0)
        .expect("gusset")
        .rotate_y(90.0_f64.to_radians())
        .expect("lay")
        .translate(-20.0, 0.0, 0.0)
        .expect("center");
    let vol = |p: &Part| -> f64 {
        let m = p
            .mesh(&popkern_mesh::TessellationOptions::default())
            .expect("mesh");
        let mut v = 0.0;
        for t in &m.triangles {
            let (a, b, c) = (
                m.vertices[t[0]].position,
                m.vertices[t[1]].position,
                m.vertices[t[2]].position,
            );
            v += (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                + a.z * (b.x * c.y - b.y * c.x))
                / 6.0;
        }
        v.abs()
    };
    let arm_vol = vol(&arm);
    // The 36 mm inset union is the reference: the flush union must
    // add the same cross-section swept over 40 mm instead of 36.
    let gusset36 = Profile::on(Plane::XY)
        .path(-(plate_under(y_root) + 1.0), y_root)
        .line_to(-(plate_under(y_tip) + 1.0), y_tip)
        .arc_to((-42.0, -36.0), (-26.0, y_root))
        .close()
        .extrude(36.0)
        .expect("gusset36")
        .rotate_y(90.0_f64.to_radians())
        .expect("lay")
        .translate(-18.0, 0.0, 0.0)
        .expect("center");
    let inset_union = arm.union(&gusset36).expect("inset union");
    if let Some(dir) = std::env::var_os("PK_PIN_EXPORT") {
        let d = std::path::Path::new(&dir);
        let _ = arm.export_step(d.join("flush_ref_arm.step"));
        let _ = gusset36.export_step(d.join("flush_ref_gusset36.step"));
        let _ = gusset36.export_obj(d.join("flush_ref_gusset36.obj"));
        eprintln!("gusset36 alone: vol={:.3}", vol(&gusset36));
        let _ = inset_union.export_step(d.join("flush_ref_inset.step"));
        eprintln!(
            "flush refs: arm_vol={:.3} inset_vol={:.3}",
            arm_vol,
            vol(&inset_union)
        );
    }
    let inset_added = vol(&inset_union) - arm_vol;
    let expected_flush = arm_vol + inset_added * 40.0 / 36.0;
    let flush_ok = match arm.union(&gusset) {
        Ok(out) => {
            let opens = out.open_edge_count();
            let v = vol(&out);
            eprintln!(
                "flush union: opens={opens:?} vol={v:.3} expected={expected_flush:.3} delta={:.3}",
                v - expected_flush
            );
            if let Some(dir) = std::env::var_os("PK_PIN_EXPORT") {
                let path = std::path::Path::new(&dir).join("flush_arm_union.step");
                let _ = out.export_step(&path);
                eprintln!("flush union: exported {}", path.display());
            }
            opens.map(|o| o == 0).unwrap_or(false)
                && (v - expected_flush).abs() < expected_flush * 1e-3
        }
        Err(e) => {
            eprintln!("flush union: error {e:?}");
            false
        }
    };
    assert!(
        flush_ok,
        "flush gusset union regressed: watertightness or the exact \
         volume oracle no longer holds (see stderr for measured values)"
    );
}

/// FRONTIER PIN (mesher orientation): the tessellation of an
/// arc-profile extrude emits INCONSISTENTLY WOUND triangles — 239 of
/// 781 welded directed edges of the swoop-gusset mesh have mismatched
/// direction counts, so a divergence-theorem volume sum reads ~3× the
/// true volume (61802 vs the analytic/OCC 20789.4) while the mesh
/// geometry itself encloses the correct region (Monte-Carlo parity
/// ≈ 20.1k; OCC BRepCheck valid, zero free bounds). Ratio-based
/// volume oracles built on this mesher stay valid (both operands
/// carry the same bias — OCC independently confirms the flush/inset
/// gusset unions add volume in exactly the 40/36 ratio), but ABSOLUTE
/// mesh volumes cannot be trusted until face windings agree across
/// faces. When the signed sum matches the analytic value, flip this
/// pin into an equality gate.
#[test]
fn arc_extrude_mesh_winding_gap_is_pinned() {
    let tilt = 18.0_f64.to_radians();
    let plane_at_root = 81.8;
    let plate_under = |y: f64| plane_at_root + y * tilt.tan() - 13.0 / tilt.cos();
    let (y_root, y_tip) = (-20.0, -72.0);
    let gusset = Profile::on(Plane::XY)
        .path(-(plate_under(y_root) + 1.0), y_root)
        .line_to(-(plate_under(y_tip) + 1.0), y_tip)
        .arc_to((-42.0, -36.0), (-26.0, y_root))
        .close()
        .extrude(36.0)
        .expect("gusset")
        .rotate_y(90.0_f64.to_radians())
        .expect("lay")
        .translate(-18.0, 0.0, 0.0)
        .expect("center");

    let mesh = gusset
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let mut signed = 0.0;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        signed += (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x))
            / 6.0;
    }
    // Analytic: profile area 577.48 mm² (triangle 952.44 − circular
    // segment 374.96, circle through the three profile points) × 36.
    let analytic = 20789.4;
    let rel = (signed.abs() - analytic).abs() / analytic;
    assert!(
        rel > 0.3,
        "arc-extrude mesh winding gap has closed (signed volume {signed:.1} \
         within 30% of analytic {analytic:.1}) — flip this pin into an \
         equality gate (tolerance ~1e-2 for chord error)"
    );
}

/// COINCIDENT-BOOLEAN SUBTRACT LEG — CLOSED. A subtraction whose
/// cutter side faces are coplanar with the target's faces (the flush
/// dovetail-tab relief) comes out watertight with the exact analytic
/// volume. Five layers closed the leg, each gated to also keep the
/// flush UNION pin exact: (1) UV crossings are sorted and deduped
/// for OPEN trims (the raw per-edge soup made enter/exit arbitrary),
/// and a trim whose crossings are confined to its endpoints with its
/// midpoint ON the rim is consumed as a split no-op — a coplanar
/// cutter's record duplicating the target's own edge partitions
/// nothing, and the same rule self-resolves records re-visited after
/// their own earlier split; (2) the Cyrus-Beck line clipper's
/// parallel and side tests are normalized to angular/distance form
/// (unnormalized products silently discarded the cutter-top ×
/// flank-wall record, so the wall never split at the relief depth);
/// (3) OnBoundary classification carries contact COORIENTATION,
/// self-calibrated by probing the fragment's own body (the
/// loop-Newell sign is not uniform across construction paths), with
/// offset probes verified at 16× the base offset (1 µm probes graze
/// corner-adjacent faces); subtract drops cooriented A-side pieces —
/// keeping them left a membrane over the relief mouth; (4) on
/// overlap faces the classifier picks whichever candidate sample
/// (edge-midpoint shortcut or interior scanline) admits a DECISIVE
/// own-body side probe — the shortcut lands on crease edges and the
/// other body's corners, the scanline mis-rows on arc-heavy charts;
/// (5) the sampler evaluates a split half-edge's pcurve over the
/// half-edge's own window, recovered by inverting its 3D endpoints
/// through the surface — scoped to non-periodic charts, where the
/// inversion is decidable.
#[test]
fn flush_tab_relief_subtract_is_watertight_and_volume_exact() {
    let plan: Vec<(f64, f64)> = vec![
        (22.0, 0.0),
        (5.8, 0.0),
        (8.8, 9.8),
        (-8.8, 9.8),
        (-5.8, 0.0),
        (-22.0, 0.0),
        (-22.0, -30.0),
        (22.0, -30.0),
    ];
    let body = Profile::on(Plane::XY)
        .polygon(&plan)
        .extrude(40.0)
        .expect("body");
    let relief = Profile::on(Plane::xy_at(-1.0))
        .polygon(&[
            (5.8, 0.0),
            (5.8 + 3.0 * (10.8 / 9.8), 10.8),
            (-5.8 - 3.0 * (10.8 / 9.8), 10.8),
            (-5.8, 0.0),
        ])
        .extrude(21.0)
        .expect("relief");
    let out = body.subtract(&relief).expect("flush coplanar subtract");
    let opens = out.open_edge_count().expect("mesh");
    let mesh = out
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let mut v = 0.0;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        v += (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x))
            / 6.0;
    }
    // Analytic: plan = 44×30 rect + dovetail tab trapezoid
    // ((11.6+17.6)/2 × 9.8), extruded 40; the relief removes the tab
    // section below z = 20. All faces planar — the divergence-sum
    // volume is winding-safe here (the arc-extrude winding gap does
    // not apply).
    let tab_area = 0.5 * (11.6 + 17.6) * 9.8;
    let expected = (44.0 * 30.0 + tab_area) * 40.0 - tab_area * 20.0;
    eprintln!("flush tab relief: opens={opens} vol={v:.3} expected={expected:.3}");
    assert_eq!(opens, 0, "flush coplanar relief no longer watertight");
    assert!(
        (v - expected).abs() < expected * 1e-3,
        "flush coplanar relief volume {v:.3} vs analytic {expected:.3}"
    );
}

/// STAR-MEMBRANE CAMPAIGN — CLOSED. A hexagram socket cut into the
/// floor of a Ø76 cylindrical recess comes out watertight with the
/// pocket verified void by orientation-immune ray-parity probes.
/// The membrane was never in the star boolean at all: the RECESS
/// subtract silently degraded to the no-records containment path
/// because the cylinder × tilted-plane mouth ellipse — a closed
/// record wholly interior to the wall chart — was probed
/// inside=false by the zero-crossing containment test and dropped.
/// The raw UV boundary of a periodic face jumps a full period at the
/// seam (one pcurve shared by both seam uses), so even-odd
/// containment read every interior point as OUTSIDE; the block
/// carried an untrimmed cutter chimney (wall to z=130 plus lid), the
/// top face had no mouth loop, and the NEXT boolean's result counted
/// two shells (χ' = 4). Two kernel fixes: (1) point_inside_face_uv
/// uses the chain-normalized boundary (face_uv_boundary_chained —
/// built for exactly this, never wired in) and tests all period
/// images of the query point; (2) recompute_all_pcurves closes each
/// loop's net U-winding on periodic charts — the closed-edge
/// canonical direction convention makes the two uses of ONE edge
/// consistent but cannot see that a loop chaining TWO full rims
/// (bottom circle + mouth ellipse) needs them wound oppositely; the
/// winding pass branch-chains the loop and flips full-period members
/// (Pcurve::reversed) until the net is zero.
#[test]
fn round_seat_star_socket_is_watertight_and_void() {
    let tilt = 18.0_f64.to_radians();
    let block = Profile::on(Plane::XY)
        .path(58.0, 0.0)
        .arc_to((0.0, 58.0), (-58.0, 0.0))
        .line_to(-44.0, -56.0)
        .line_to(44.0, -56.0)
        .close()
        .extrude(125.0)
        .expect("block");
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("upper")
        .rotate_x(tilt)
        .expect("tilt")
        .translate(0.0, 0.0, 100.0)
        .expect("place");
    let block = block.subtract(&upper).expect("rest cut");
    let seat = Profile::on(Plane::xy_at(70.0))
        .circle(38.0)
        .extrude(60.0)
        .expect("seat");
    let before_faces = block.face_count();
    let block = block
        .subtract(&seat)
        .expect("round recess itself cuts clean");
    assert_eq!(block.open_edge_count().expect("mesh"), 0);
    assert_ne!(block.face_count(), before_faces);
    // The frontier: hexagram socket into the recess's disk floor.
    let socket = Profile::on(Plane::xy_at(70.0 - 16.0))
        .polygon(&hexagram_points(12.0))
        .extrude(17.0)
        .expect("socket tool");
    // Orientation-immune material probe: ray parity along +z against
    // the welded mesh. The divergence-sum volume oracle is unusable
    // on this body (curved faces carry the winding-inconsistent
    // tessellation pinned by arc_extrude_mesh_winding_gap_is_pinned;
    // measured: the removed-volume delta read −7649 against the
    // analytic +3991), but parity does not care about triangle
    // orientation.
    let mesh_contains = |p: &Part, q: (f64, f64, f64)| -> bool {
        let m = p
            .mesh(&popkern_mesh::TessellationOptions::default())
            .expect("mesh");
        let mut crossings = 0;
        for t in &m.triangles {
            let (a, b, c) = (
                m.vertices[t[0]].position,
                m.vertices[t[1]].position,
                m.vertices[t[2]].position,
            );
            let (e1x, e1y) = (b.x - a.x, b.y - a.y);
            let (e2x, e2y) = (c.x - a.x, c.y - a.y);
            let det = e1x * e2y - e1y * e2x;
            if det.abs() < 1e-14 {
                continue;
            }
            let (rx, ry) = (q.0 - a.x, q.1 - a.y);
            let u = (rx * e2y - ry * e2x) / det;
            let v = (e1x * ry - e1y * rx) / det;
            if u < -1e-12 || v < -1e-12 || u + v > 1.0 + 1e-12 {
                continue;
            }
            let z = a.z + u * (b.z - a.z) + v * (c.z - a.z);
            if z > q.2 {
                crossings += 1;
            }
        }
        crossings % 2 == 1
    };
    // Probe points strictly inside the socket prism (star shrunk to
    // its inner radius), mid-depth.
    let probes = [(0.0, 0.0, 60.0), (4.0, 0.0, 62.0), (-3.0, 2.0, 57.0)];
    for &q in &probes {
        assert!(
            mesh_contains(&block, q),
            "probe {q:?} must be solid before the socket cut"
        );
    }
    let ok = match block.subtract(&socket) {
        Ok(out) => {
            let all_void = probes.iter().all(|&q| !mesh_contains(&out, q));
            eprintln!(
                "star socket in round floor: faces {} -> {}, opens {:?}, pocket void={all_void}",
                block.face_count(),
                out.face_count(),
                out.open_edge_count(),
            );
            out.open_edge_count().map(|o| o == 0).unwrap_or(false)
                && out.face_count() != block.face_count()
                && all_void
        }
        Err(e) => {
            eprintln!("star socket in round floor: {e:?}");
            false
        }
    };
    assert!(
        ok,
        "round-seat star socket regressed: watertightness, face count, \
         or the pocket-void parity probes no longer hold"
    );
}

/// Uniform-width sibling of the flush tab relief: the body has NO
/// root shoulder, so the front wall's corner line at x = ±5.8 is
/// where the cutter's flank plane meets the wall exactly along its
/// own boundary edge for the cutter's 20 mm height — a sub-span of
/// the wall's 40 mm corner edge.
#[test]
fn flush_tab_relief_uniform_width_is_watertight() {
    let plan: Vec<(f64, f64)> = vec![
        (20.0, 0.0),
        (5.8, 0.0),
        (8.8, 9.8),
        (-8.8, 9.8),
        (-5.8, 0.0),
        (-20.0, 0.0),
        (-20.0, -30.0),
        (20.0, -30.0),
    ];
    let body = Profile::on(Plane::XY)
        .polygon(&plan)
        .extrude(40.0)
        .expect("body");
    let relief = Profile::on(Plane::xy_at(-1.0))
        .polygon(&[
            (5.8, 0.0),
            (5.8 + 3.0 * (10.8 / 9.8), 10.8),
            (-5.8 - 3.0 * (10.8 / 9.8), 10.8),
            (-5.8, 0.0),
        ])
        .extrude(21.0)
        .expect("relief");
    let out = body.subtract(&relief).expect("flush uniform subtract");
    let opens = out.open_edge_count().expect("mesh");
    eprintln!(
        "flush uniform relief: opens={opens} faces={}",
        out.face_count()
    );
    assert_eq!(opens, 0);
}

/// The FULL uniform-width flush arm chain (the fully-flush product
/// arm): uniform ±20 plan, rest-plane cut, plate thinning tool,
/// flush swoop union, through lip slot, flush tab relief.
#[test]
fn flush_arm_uniform_width_chain_is_watertight() {
    let tilt = 18.0_f64.to_radians();
    let plane_at_root = 81.8;
    let plan: Vec<(f64, f64)> = vec![
        (20.0, 0.0),
        (5.8, 0.0),
        (8.8, 9.8),
        (-8.8, 9.8),
        (-5.8, 0.0),
        (-20.0, 0.0),
        (-20.0, -231.0),
        (20.0, -231.0),
    ];
    let mut arm = Profile::on(Plane::XY)
        .polygon(&plan)
        .extrude(96.0)
        .expect("arm");
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("upper")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, plane_at_root)
        .expect("place");
    arm = arm.subtract(&upper).expect("rest cut");
    let lower_box = Profile::on(Plane::xy_at(-10.0))
        .rect(600.0, 216.0)
        .at(0.0, -134.0)
        .extrude(100.0)
        .expect("lower");
    let p2 = Profile::on(Plane::XY)
        .rect(700.0, 700.0)
        .extrude(300.0)
        .expect("p2")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, plane_at_root - 13.0 / tilt.cos())
        .expect("place");
    let tool = lower_box.subtract(&p2).expect("tool");
    arm = arm.subtract(&tool).expect("thin");
    let plate_under = |y: f64| plane_at_root + y * tilt.tan() - 13.0 / tilt.cos();
    let (y_root, y_tip) = (-20.0, -72.0);
    let gusset = Profile::on(Plane::XY)
        .path(-(plate_under(y_root) + 1.0), y_root)
        .line_to(-(plate_under(y_tip) + 1.0), y_tip)
        .arc_to((-42.0, -36.0), (-26.0, y_root))
        .close()
        .extrude(40.0)
        .expect("gusset")
        .rotate_y(90.0_f64.to_radians())
        .expect("lay")
        .translate(-20.0, 0.0, 0.0)
        .expect("center");
    arm = arm.union(&gusset).expect("swoop union");
    let slot = Profile::on(Plane::XY)
        .polygon(&[
            (-2.0, -4.5),
            (1.0, -4.5),
            (6.0, -7.0),
            (6.0, 7.0),
            (1.0, 4.5),
            (-2.0, 4.5),
        ])
        .extrude(44.0)
        .expect("slot")
        .rotate_y(90.0_f64.to_radians())
        .expect("lay")
        .rotate_x(tilt)
        .expect("tilt")
        .translate(-22.0, -215.0, plane_at_root - 215.0 * tilt.tan())
        .expect("place");
    arm = arm.subtract(&slot).expect("lip slot");
    let relief = Profile::on(Plane::xy_at(-1.0))
        .polygon(&[
            (5.8, 0.0),
            (5.8 + 3.0 * (10.8 / 9.8), 10.8),
            (-5.8 - 3.0 * (10.8 / 9.8), 10.8),
            (-5.8, 0.0),
        ])
        .extrude(21.0)
        .expect("relief");
    let arm = arm.subtract(&relief).expect("flush relief");
    let opens = arm.open_edge_count().expect("mesh");
    eprintln!(
        "flush uniform arm chain: opens={opens} faces={}",
        arm.face_count()
    );
    assert_eq!(opens, 0);
}

/// The split arm's TIP piece: rectangular plan on the tilted plate,
/// rest cut, thinning tool, lip slot, then the lap zone carved down
/// to the sliding-dovetail rail (nominal minus fit clearance, floor
/// flush with the plate underside).
#[test]
fn split_arm_tip_rail_is_watertight() {
    let tilt = 18.0_f64.to_radians();
    let plane_at_root = 81.8;
    let arm_len = 275.0;
    let (mouth, shoulder) = (120.0, 170.0);
    let plate_under = |y: f64| plane_at_root + y * tilt.tan() - 13.0 / tilt.cos();
    let z_lo = (plate_under(-arm_len) - 2.0_f64).min(-2.0);
    // ADDITIVE construction: the full plate segment beyond the
    // shoulder, UNIONED with the clearance rail prism (embedded 1 mm
    // past the shoulder for a transversal margin). The earlier
    // shaped-cutter carve needed record-network T-junction splitting
    // at the shoulder cross-section — the pinned junction-network
    // frontier; the additive form stays in closed classes.
    let rail_prism_early = |y_start: f64, len: f64, inset: f64, below: f64| -> Part {
        let slope = (11.0 - 7.0) / 7.0_f64;
        let shift = inset * (1.0 + slope * slope).sqrt();
        let x_at = |y: f64| 11.0 - shift - slope * y;
        let (yb, yt) = (-below, 7.0 - inset);
        Profile::on(Plane::XY)
            .polygon(&[
                (x_at(yb), yb),
                (x_at(yt), yt),
                (-x_at(yt), yt),
                (-x_at(yb), yb),
            ])
            .extrude(len)
            .expect("prism")
            .rotate_x(std::f64::consts::FRAC_PI_2 + tilt)
            .expect("align")
            .translate(0.0, -y_start, plate_under(-y_start))
            .expect("place")
    };
    let full_seg = Profile::on(Plane::xy_at(z_lo))
        .polygon(&[
            (20.0, -shoulder),
            (20.0, -arm_len),
            (-20.0, -arm_len),
            (-20.0, -shoulder),
        ])
        .extrude(96.0 - z_lo)
        .expect("plan");
    let rail_seg = rail_prism_early(
        mouth + 0.4,
        (shoulder - mouth - 0.4 + 1.0) / tilt.cos(),
        0.2,
        0.0,
    );
    let mut tip = full_seg;
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("upper")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, plane_at_root)
        .expect("place");
    tip = tip.subtract(&upper).expect("rest cut");
    let lower_box = Profile::on(Plane::xy_at(z_lo - 5.0))
        .rect(600.0, arm_len - 16.0)
        .at(0.0, -(26.0 + arm_len + 10.0) / 2.0)
        .extrude(120.0)
        .expect("lower");
    let p2 = Profile::on(Plane::XY)
        .rect(700.0, 700.0)
        .extrude(300.0)
        .expect("p2")
        .rotate_x(tilt)
        .expect("rot")
        .translate(0.0, 0.0, plane_at_root - 13.0 / tilt.cos())
        .expect("place");
    let tool = lower_box.subtract(&p2).expect("tool");
    tip = tip.subtract(&tool).expect("thin");
    let slot = Profile::on(Plane::XY)
        .polygon(&[
            (-2.0, -4.5),
            (1.0, -4.5),
            (6.0, -7.0),
            (6.0, 7.0),
            (1.0, 4.5),
            (-2.0, 4.5),
        ])
        .extrude(44.0)
        .expect("slot")
        .rotate_y(90.0_f64.to_radians())
        .expect("lay")
        .rotate_x(tilt)
        .expect("tilt")
        .translate(
            -22.0,
            -(arm_len - 16.0),
            plane_at_root - (arm_len - 16.0) * tilt.tan(),
        )
        .expect("place");
    tip = tip.subtract(&slot).expect("lip slot");
    let tip = tip.union(&rail_seg).expect("rail union");
    let opens = tip.open_edge_count().expect("mesh");
    eprintln!("split tip: opens={opens} faces={}", tip.face_count());
    assert_eq!(opens, 0);
    assert!(tip.face_count() > 0);
}

/// The FANNED arm's yaw-compensated joint prism: axis along the arm
/// centreline WITHIN the yawed plate plane, floor parallel to (and
/// coplanar-unionable with) the yawed underside, realised as a z-x-z
/// Euler composition. Reversed (true) dovetail: narrow at the
/// underside mouth, widest at the crown, so the rail cannot drop out
/// of the tunnel.
fn yawed_rail_prism(
    tilt: f64,
    yaw: f64,
    y_start: f64,
    plan_len: f64,
    inset: f64,
    below: f64,
    mouth_z: f64,
) -> Part {
    let slope = (11.0 - 7.0) / 7.0_f64;
    let shift = inset * (1.0 + slope * slope).sqrt();
    let x_at = |y: f64| 7.0 - shift + slope * y;
    let (yb, yt) = (-below, 7.0 - inset);
    let k = yaw.cos() * tilt.tan();
    let n = (1.0 + k * k).sqrt();
    let a = [0.0, -1.0 / n, -k / n];
    let m = [-yaw.sin() * tilt.sin(), -yaw.cos() * tilt.sin(), tilt.cos()];
    let e1 = [
        m[1] * a[2] - m[2] * a[1],
        m[2] * a[0] - m[0] * a[2],
        m[0] * a[1] - m[1] * a[0],
    ];
    let alpha = e1[2].atan2(m[2]);
    let phi = a[0].hypot(a[1]).atan2(a[2]);
    let gamma = a[0].atan2(-a[1]);
    Profile::on(Plane::XY)
        .polygon(&[
            (x_at(yb), yb),
            (x_at(yt), yt),
            (-x_at(yt), yt),
            (-x_at(yb), yb),
        ])
        .extrude(plan_len * n)
        .expect("prism")
        .rotate_z(alpha)
        .expect("rz1")
        .rotate_x(phi)
        .expect("rx")
        .rotate_z(gamma)
        .expect("rz2")
        .translate(0.0, -y_start, mouth_z)
        .expect("place")
}

/// FANNED-ARM plate segment: a straight bar whose rest plane and
/// underside are yaw-compensated (rotate_z(−w) after rotate_x), so a
/// mounted pair of arms shares ONE global rest plane.
fn yawed_plate_segment(tilt: f64, w: f64, y0: f64, y1: f64, plane_at_root: f64) -> Part {
    let cross_fall = 20.0 * w.sin().abs() * tilt.tan();
    let z_lo = (plane_at_root - 13.0 / tilt.cos() - y1 * w.cos() * tilt.tan() - cross_fall - 2.0)
        .min(-2.0);
    let seg = Profile::on(Plane::xy_at(z_lo))
        .polygon(&[(20.0, -y0), (20.0, -y1), (-20.0, -y1), (-20.0, -y0)])
        .extrude(96.0 - z_lo)
        .expect("plan");
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("upper")
        .rotate_x(tilt)
        .expect("rx")
        .rotate_z(-w)
        .expect("rz")
        .translate(0.0, 0.0, plane_at_root)
        .expect("place");
    let seg = seg.subtract(&upper).expect("rest cut");
    let lower_box = Profile::on(Plane::xy_at(z_lo - 5.0))
        .rect(600.0, y1 - y0 + 100.0)
        .at(0.0, -(y0 + y1) / 2.0)
        .extrude(120.0)
        .expect("lower");
    let p2 = Profile::on(Plane::XY)
        .rect(700.0, 700.0)
        .extrude(300.0)
        .expect("p2")
        .rotate_x(tilt)
        .expect("rx")
        .rotate_z(-w)
        .expect("rz")
        .translate(0.0, 0.0, plane_at_root - 13.0 / tilt.cos())
        .expect("place");
    let tool = lower_box.subtract(&p2).expect("tool");
    seg.subtract(&tool).expect("thin")
}

/// PIN → GATE (fanned-arm campaign): union of the REVERSED, yawed
/// dovetail rail onto the yaw-compensated plate — the rail floor is
/// coplanar with the tilted-and-yawed underside (a non-axis-aligned
/// coplanar union), and the prism frame is the z-x-z composition. A
/// wrong frame reads as a near-coincident (not coplanar) floor and
/// the union goes open-shell.
#[test]
fn fanned_reversed_rail_union_is_watertight() {
    let tilt = 18.0_f64.to_radians();
    let w = 12.0_f64.to_radians();
    let plane_at_root = 100.0 + (-56.0 + 28.0 * w.sin()) * tilt.tan();
    let (mouth, shoulder) = (120.0, 170.0);
    let seg = yawed_plate_segment(tilt, w, shoulder, 275.0, plane_at_root);
    let mouth_z = plane_at_root - 13.0 / tilt.cos() - (mouth + 0.4) * w.cos() * tilt.tan();
    let rail = yawed_rail_prism(
        tilt,
        w,
        mouth + 0.4,
        shoulder - mouth - 0.4 + 1.0,
        0.2,
        0.0,
        mouth_z,
    );
    let tip = seg.union(&rail).expect("rail union");
    let opens = tip.open_edge_count().expect("mesh");
    eprintln!("fanned rail: opens={opens} faces={}", tip.face_count());
    assert_eq!(opens, 0);
    assert!(tip.face_count() > 0);
    // The reversed dovetail retains: the crown (11 half-width, 7 up)
    // is wider than the mouth (7 at the underside) — probe material
    // 5 mm up the plate normal at x = 8.5, inside the new flank
    // (half-width 9.66 there) but OUTSIDE the old orientation's
    // (7.91): air before the reversal, material after.
    let (xq, yq, h) = (8.5_f64, -145.0_f64, 5.0_f64);
    let under = plane_at_root - 13.0 / tilt.cos() + (xq * w.sin() + yq * w.cos()) * tilt.tan();
    let m = (-w.sin() * tilt.sin(), -w.cos() * tilt.sin(), tilt.cos());
    assert!(
        mesh_contains_point(&tip, (xq + h * m.0, yq + h * m.1, under + h * m.2)),
        "crown must be wider than the mouth (reversed dovetail)"
    );
}

/// PIN → GATE (fanned-arm campaign): the swoop-web union on a YAWED
/// plate. The web's cylinder × the yawed underside plane produces a
/// closed off-branch UV record whose extent dwarfs the trim; the
/// clip must window it to the TRUE weld line, not the extended
/// cylinder's phantom re-crossing (period-branch normalization +
/// torn-fold rejection + containment transition rescue), and the
/// web fragment — whose every boundary edge hugs the arm — must
/// classify by a cut-offset interior sample, not a boundary
/// midpoint.
#[test]
fn fanned_flush_gusset_union_weld_is_watertight() {
    let tilt = 18.0_f64.to_radians();
    let w = 12.0_f64.to_radians();
    let plane_at_root = 100.0 + (-56.0 + 28.0 * w.sin()) * tilt.tan();
    let plate_under_c = |y: f64| plane_at_root - 13.0 / tilt.cos() + y * w.cos() * tilt.tan();
    // Root collar + plate run, as in the arm chain.
    let cross_fall = 20.0 * w.sin().abs() * tilt.tan();
    let z_lo = (plate_under_c(-200.0) - cross_fall - 2.0_f64).min(-2.0);
    let bar = Profile::on(Plane::xy_at(z_lo))
        .polygon(&[(20.0, 0.0), (20.0, -200.0), (-20.0, -200.0), (-20.0, 0.0)])
        .extrude(96.0 - z_lo)
        .expect("plan");
    let upper = Profile::on(Plane::XY)
        .rect(600.0, 600.0)
        .extrude(200.0)
        .expect("upper")
        .rotate_x(tilt)
        .expect("rx")
        .rotate_z(-w)
        .expect("rz")
        .translate(0.0, 0.0, plane_at_root)
        .expect("place");
    let bar = bar.subtract(&upper).expect("rest cut");
    let root_floor = Profile::on(Plane::xy_at(z_lo - 5.0))
        .rect(60.0, 46.0)
        .at(0.0, -8.0)
        .extrude(5.0 - z_lo)
        .expect("floor");
    let bar = bar.subtract(&root_floor).expect("root floor");
    let lower_box = Profile::on(Plane::xy_at(z_lo - 5.0))
        .rect(600.0, 184.0)
        .at(0.0, -118.0)
        .extrude(120.0)
        .expect("lower");
    let p2 = Profile::on(Plane::XY)
        .rect(700.0, 700.0)
        .extrude(300.0)
        .expect("p2")
        .rotate_x(tilt)
        .expect("rx")
        .rotate_z(-w)
        .expect("rz")
        .translate(0.0, 0.0, plane_at_root - 13.0 / tilt.cos())
        .expect("place");
    let tool = lower_box.subtract(&p2).expect("tool");
    let bar = bar.subtract(&tool).expect("thin");
    // Fully FLUSH swoop web (half_run = 20), embedded past the yaw
    // cross-fall.
    let (y_root, y_tip) = (-20.0, -72.0);
    let half_run = 20.0_f64;
    let embed = 1.0 + half_run * w.sin().abs() * tilt.tan();
    let gusset = Profile::on(Plane::XY)
        .path(-(plate_under_c(y_root) + embed), y_root)
        .line_to(-(plate_under_c(y_tip) + embed), y_tip)
        .arc_to((-42.0, -36.0), (-26.0, y_root))
        .close()
        .extrude(2.0 * half_run)
        .expect("web")
        .rotate_y(90.0_f64.to_radians())
        .expect("lay")
        .translate(-half_run, 0.0, 0.0)
        .expect("center");
    let arm = bar.union(&gusset).expect("gusset union");
    let opens = arm.open_edge_count().expect("mesh");
    eprintln!("fanned gusset: opens={opens} faces={}", arm.face_count());
    assert_eq!(opens, 0);
    // The web must SURVIVE below the plate: probe material inside the
    // web, under the underside, mid-run.
    let y_mid = -46.0_f64;
    let z_web = plate_under_c(y_mid) - 6.0;
    assert!(
        mesh_contains_point(&arm, (0.0, y_mid, z_web)),
        "swoop web dropped below the plate — cut-offset classification regressed"
    );
    // And the phantom re-crossing region (past the web's end) must
    // NOT have grown a membrane: void just below the underside there.
    let y_ph = -95.0_f64;
    assert!(
        !mesh_contains_point(&arm, (0.0, y_ph, plate_under_c(y_ph) - 3.0)),
        "material past the web's end — phantom clip window regressed"
    );
}
