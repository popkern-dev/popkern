//! STEP export of every part a general body holds (PK-1114): a sheet
//! glued beside a solid, a lamina glued inside one, and a lamina, wire
//! and acorn floating in a solid beside a wire and an acorn outside it.
//! Each writes as the items `export_general` names (the solid's
//! `MANIFOLD_SOLID_BREP`, a `SHELL_BASED_SURFACE_MODEL`, a
//! `GEOMETRIC_CURVE_SET`), OpenCASCADE reads the solid's volume and
//! every face, edge and vertex from the file (an edge two items share
//! once per item, as it reads its own compounds), and popkern reads it back
//! under `StepNonManifold::General` as the same general body: the same
//! counts, each shell in the region it was in, the same region volumes,
//! and valid.

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind, RegionKind};
use popkern_io::{StepExportOptions, StepNonManifold, export_step};

#[path = "step_general_forms_support/harness.rs"]
mod harness;
use harness::{assert_close, bodies, census, document, occ_read, solid_volumes};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube() -> BRep {
    BRep::unit_cube(tol())
        .and_then(|c| c.transform(&Transform3::identity()))
        .expect("cube")
}

fn line(p: Point3, q: Point3) -> CurveOnEdge {
    CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: p,
            direction: q - p,
        },
        t_start: 0.0,
        t_end: 1.0,
    }
}

/// A wire edge from `p` to `q` in the outside of `body`.
fn wire(body: &mut BRep, p: Point3, q: Point3) {
    let outside = body.outside;
    let a = body.make_acorn(outside, p).expect("acorn").vertex;
    let b = body.make_acorn(outside, q).expect("acorn").vertex;
    let made = body.make_wire_edge(outside, a, b).expect("wire");
    body.edges[made.edge].curve = Some(line(p, q));
}

/// A one-face planar sheet in the void over the convex polygon `points`.
fn polygon(points: &[(f64, f64, f64)]) -> BRep {
    let mut body = BRep::new(tol());
    let outside = body.outside;
    let pts: Vec<Point3> = points
        .iter()
        .map(|&(x, y, z)| Point3::new(x, y, z))
        .collect();
    let v: Vec<_> = pts
        .iter()
        .map(|&p| body.make_acorn(outside, p).expect("acorn").vertex)
        .collect();
    let n = pts.len();
    let fins: Vec<_> = (0..n)
        .map(|i| {
            let made = body
                .make_wire_edge(outside, v[i], v[(i + 1) % n])
                .expect("wire");
            body.edges[made.edge].curve = Some(line(pts[i], pts[(i + 1) % n]));
            made.fin
        })
        .collect();
    let u: Vec3 = pts[1] - pts[0];
    let normal = u.cross(&(pts[2] - pts[0])).normalized().expect("plane");
    let u = u.normalized().expect("axis");
    let plane = SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin: pts[0],
            normal,
            u_axis: u,
            v_axis: normal.cross(&u),
        },
        orientation: true,
    };
    let made = body
        .make_face_from_wire_cycle(&fins, Some(plane))
        .expect("face");
    body.recompute_face_pcurves(made.face).expect("pcurves");
    body
}

/// A cube with a sheet glued on its edge along the z axis, the sheet in
/// the plane x = y reaching to `(c, c)`: beside the cube for negative
/// `c`, a lamina inside it for positive.
fn cube_with_fin(c: f64) -> BRep {
    let mut body = cube();
    body.append_body(&polygon(&[
        (0.0, 0.0, 0.0),
        (c, c, 0.3),
        (c, c, 0.7),
        (0.0, 0.0, 1.0),
    ]));
    let glued = body.glue_coincident().expect("glue");
    assert_eq!((glued.vertices, glued.edges), (2, 1));
    body
}

/// A cube holding a floating lamina, a wire and an acorn, with a wire
/// and an acorn outside it.
fn cube_holding_free_parts() -> BRep {
    let mut body = cube();
    body.append_body(&polygon(&[
        (0.25, 0.25, 0.4),
        (0.75, 0.25, 0.4),
        (0.75, 0.75, 0.4),
        (0.25, 0.75, 0.4),
    ]));
    wire(
        &mut body,
        Point3::new(0.2, 0.2, 0.8),
        Point3::new(0.8, 0.7, 0.8),
    );
    wire(
        &mut body,
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(3.0, 0.5, 0.0),
    );
    let outside = body.outside;
    body.make_acorn(outside, Point3::new(0.5, 0.5, 0.6))
        .expect("acorn");
    body.make_acorn(outside, Point3::new(-1.0, 0.5, 0.5))
        .expect("acorn");
    body.fill_enclosed_cells();
    body
}

/// Three sheets glued on one hinge edge (the z axis): a general body of
/// sheets alone, its edge a ring of three.
fn three_sheets_on_a_hinge() -> BRep {
    let mut body = BRep::new(tol());
    for (x, y) in [(1.0, 0.0), (0.0, 1.0), (-1.0, -1.0)] {
        body.append_body(&polygon(&[
            (0.0, 0.0, 0.0),
            (x, y, 0.0),
            (x, y, 1.0),
            (0.0, 0.0, 1.0),
        ]));
    }
    body.glue_coincident().expect("glue");
    body
}

/// Two wires meeting at a vertex (one a quarter circle) and an acorn
/// apart: no faces at all.
fn wires_and_an_acorn() -> BRep {
    use popkern_core::geometry::curve::CurveEval;
    let mut body = BRep::new(tol());
    let outside = body.outside;
    let arc = CurveOnEdge {
        geometry: CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        },
        t_start: 0.0,
        t_end: std::f64::consts::FRAC_PI_2,
    };
    let (o, p, q) = (
        Point3::new(0.0, 0.0, 0.0),
        arc.geometry.point_at(arc.t_start),
        arc.geometry.point_at(arc.t_end),
    );
    let mut acorn = |at| body.make_acorn(outside, at).expect("acorn").vertex;
    let (a, b, c) = (acorn(o), acorn(p), acorn(q));
    let ab = body.make_wire_edge(outside, a, b).expect("wire");
    body.edges[ab.edge].curve = Some(line(o, p));
    let bc = body.make_wire_edge(outside, b, c).expect("wire");
    body.edges[bc.edge].curve = Some(arc);
    body.make_acorn(outside, Point3::new(0.0, 0.0, 2.0))
        .expect("acorn");
    body
}

/// A cube's surface as a closed sheet bounding a void.
fn closed_sheet_around_a_void() -> BRep {
    let mut body = cube();
    let solid = body.solid_regions()[0];
    body.kill_face_region(solid).expect("kill");
    let face = body.faces.keys().next().expect("face");
    let shell = body.face_back_shell(face).expect("shell");
    body.make_face_region(
        shell,
        popkern_core::topology::FaceSide::Back,
        RegionKind::Void,
    )
    .expect("void");
    body
}

/// Bodies without solids round trip under both options: a general body
/// of three sheets on a hinge (one `OPEN_SHELL` using its edge three
/// times), a closed sheet bounding a void (a `CLOSED_SHELL` in the
/// surface model, its inside a void again, not a solid), and wires with
/// an acorn (a curve set: a line, a quarter circle trimmed at its end
/// points, a point).
#[test]
fn sheet_bodies_without_solids_round_trip() {
    for (name, body, kind) in [
        (
            "three_sheets_on_a_hinge",
            three_sheets_on_a_hinge(),
            BodyKind::General,
        ),
        (
            "closed_sheet_around_a_void",
            closed_sheet_around_a_void(),
            BodyKind::Sheet,
        ),
        (
            "wires_and_an_acorn",
            wires_and_an_acorn(),
            BodyKind::General,
        ),
    ] {
        assert_eq!(body.body_kind(), kind, "{name}");
        let findings = body.validate_all();
        assert!(findings.is_empty(), "{name}: {findings:?}");
        let want = census(&body);
        let written =
            export_step(&document(&body, name), &StepExportOptions::default()).expect("export");
        let text = String::from_utf8_lossy(&written);
        assert!(
            !text.contains("MANIFOLD_SOLID_BREP("),
            "{name}: no solid written"
        );
        for option in [StepNonManifold::General, StepNonManifold::Pieces] {
            if kind == BodyKind::General && option == StepNonManifold::Pieces {
                continue;
            }
            let again = bodies(&written, option, name);
            assert_eq!(again.len(), 1, "{name} {option:?}: one body");
            assert_eq!(again[0].body_kind(), kind, "{name} {option:?}: re-read");
            assert_eq!(census(&again[0]), want, "{name} {option:?}: re-read census");
            let findings = again[0].validate_all();
            assert!(findings.is_empty(), "{name} {option:?}: {findings:?}");
        }
        if let Some((occ_volumes, occ_counts)) = occ_read(&written, name) {
            assert!(occ_volumes.is_empty(), "{name}: no solid");
            assert_eq!(occ_counts[0], body.faces.len() as u64, "{name}: faces");
        }
    }
}

#[test]
fn general_bodies_with_sheets_wires_and_acorns_round_trip() {
    for (name, body, items) in [
        (
            "sheet_beside_a_cube",
            cube_with_fin(-0.5),
            ["MANIFOLD_SOLID_BREP", "SHELL_BASED_SURFACE_MODEL", ""],
        ),
        (
            "lamina_in_a_cube",
            cube_with_fin(0.5),
            ["MANIFOLD_SOLID_BREP", "SHELL_BASED_SURFACE_MODEL", ""],
        ),
        (
            "free_parts_in_and_beside_a_cube",
            cube_holding_free_parts(),
            [
                "MANIFOLD_SOLID_BREP",
                "SHELL_BASED_SURFACE_MODEL",
                "GEOMETRIC_CURVE_SET",
            ],
        ),
    ] {
        assert_eq!(body.body_kind(), BodyKind::General, "{name}");
        let findings = body.validate_all();
        assert!(findings.is_empty(), "{name}: {findings:?}");
        let want = census(&body);
        let volumes = solid_volumes(&body);
        assert_close(name, &volumes, &[1.0]);

        let written =
            export_step(&document(&body, name), &StepExportOptions::default()).expect("export");
        let text = String::from_utf8_lossy(&written);
        for item in items.iter().filter(|i| !i.is_empty()) {
            assert!(text.contains(&format!("{item}(")), "{name}: writes {item}");
        }
        assert!(
            !text.contains("ADVANCED_BREP_SHAPE_REPRESENTATION("),
            "{name}: a plain SHAPE_REPRESENTATION holds the mixed items"
        );

        let again = bodies(&written, StepNonManifold::General, name);
        assert_eq!(again.len(), 1, "{name}: one body");
        assert_eq!(again[0].body_kind(), BodyKind::General, "{name}: re-read");
        assert_eq!(census(&again[0]), want, "{name}: re-read census");
        assert_close(name, &solid_volumes(&again[0]), &volumes);
        let findings = again[0].validate_all();
        assert!(findings.is_empty(), "{name}: re-read {findings:?}");

        // The default option keeps the items as bodies of their own.
        let pieces = bodies(&written, StepNonManifold::Pieces, name);
        assert!(pieces.len() > 1, "{name}: pieces");

        if let Some((occ_volumes, occ_counts)) = occ_read(&written, name) {
            assert_close(name, &occ_volumes, &volumes);
            // OpenCASCADE reads each item on its own: an edge or vertex
            // two items share (the sheet's edge on the cube) reads once
            // per item, as OpenCASCADE writes such a compound itself.
            let ours = [
                body.faces.len() as u64,
                body.edges.len() as u64,
                body.vertices.len() as u64,
            ];
            assert_eq!(occ_counts[0], ours[0], "{name}: OpenCASCADE's faces");
            assert!(
                occ_counts[1] >= ours[1] && occ_counts[2] >= ours[2],
                "{name}: OpenCASCADE's edges and vertices {occ_counts:?} vs {ours:?}"
            );
        }
    }
}

/// A cube standing on its corner on another's top face: the face holds
/// the corner as a vertex loop (PK-1116), written as a `VERTEX_LOOP`
/// bound sharing the corner's `VERTEX_POINT`. It re-reads under the
/// general option as the same general body (the vertex loop again,
/// glued to the corner), and OpenCASCADE reads both solids.
#[test]
fn a_vertex_loop_round_trips() {
    let mut body = cube();
    let diagonal = Vec3::new(1.0, 1.0, 1.0);
    let axis = diagonal.cross(&Vec3::new(0.0, 0.0, 1.0));
    let angle = (1.0 / 3.0_f64.sqrt()).acos();
    let turn =
        Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), axis, angle).expect("turn");
    let standing = BRep::unit_cube(tol())
        .and_then(|c| c.transform(&Transform3::translation(0.5, 0.5, 1.0).compose(&turn)))
        .expect("standing cube");
    let map = body.append_body(&standing);
    let tip = Point3::new(0.5, 0.5, 1.0);
    let corner = map
        .vertices
        .values()
        .copied()
        .find(|&v| body.vertices[v].position.distance_to(&tip) < 1e-9)
        .expect("corner");
    let top = body
        .faces
        .keys()
        .find(|&f| {
            body.face_point_trim_containment(f, &tip)
                .is_some_and(|c| c.distance < 1e-9)
        })
        .expect("top face");
    let made = body.make_vertex_loop(top, tip).expect("vertex loop");
    body.glue_vertices(made.vertex, corner).expect("glue");
    assert_eq!(body.body_kind(), BodyKind::General);
    let findings = body.validate_all();
    assert!(findings.is_empty(), "{findings:?}");

    let name = "cube_on_a_corner";
    let written =
        export_step(&document(&body, name), &StepExportOptions::default()).expect("export");
    assert!(String::from_utf8_lossy(&written).contains("VERTEX_LOOP("));
    let again = bodies(&written, StepNonManifold::General, name);
    assert_eq!(again.len(), 1, "one body");
    assert_eq!(again[0].vertex_loops().len(), 1, "the vertex loop again");
    assert_eq!(census(&again[0]), census(&body), "re-read census");
    assert_close(name, &solid_volumes(&again[0]), &[1.0, 1.0]);
    let findings = again[0].validate_all();
    assert!(findings.is_empty(), "re-read {findings:?}");
    if let Some((occ_volumes, _)) = occ_read(&written, name) {
        assert_close(name, &occ_volumes, &[1.0, 1.0]);
    }
}
