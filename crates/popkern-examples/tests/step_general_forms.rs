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

use std::sync::Arc;

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind, RegionKind};
use popkern_io::step::diagnostics::Severity;
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepImportOptions, StepNonManifold, StepProductNode, UnitContext,
    export_step, import_step,
};

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

fn document(brep: &BRep, name: &str) -> StepDocument {
    StepDocument {
        products: vec![StepProductNode {
            name: name.into(),
            id: name.into(),
            children: Vec::new(),
            shapes: vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids: vec![brep.clone()],
                context: RepresentationContext::default(),
                non_brep_items: Vec::new(),
            })],
            metadata: ProductMetadata::default(),
        }],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    }
}

/// Every body the file holds under `non_manifold`, with no error
/// diagnostics.
fn bodies(data: &[u8], non_manifold: StepNonManifold, name: &str) -> Vec<BRep> {
    let imported = import_step(
        data,
        &StepImportOptions {
            non_manifold,
            ..Default::default()
        },
    )
    .expect("import");
    let errors: Vec<_> = imported
        .diagnostics
        .entries()
        .iter()
        .filter(|e| e.severity == Severity::Error)
        .collect();
    assert!(errors.is_empty(), "{name}: {errors:?}");
    imported
        .document
        .products
        .iter()
        .flat_map(|p| &p.shapes)
        .flat_map(|s| s.solids.iter().cloned())
        .collect()
}

/// What the test compares: entity counts, then for each region kind
/// (outside included) how many shells hold faces, wires and acorns.
/// Vertex, edge and face counts, then one row per shell: its region
/// kind and how many faces, wire edges and acorns it holds.
type Census = ((usize, usize, usize), Vec<(String, usize, usize, usize)>);

fn census(body: &BRep) -> Census {
    let counts = (body.vertices.len(), body.edges.len(), body.faces.len());
    let mut rows: Vec<(String, usize, usize, usize)> = body
        .shells
        .iter()
        .map(|(h, s)| {
            let kind = if s.region == body.outside {
                "outside".to_string()
            } else {
                format!("{:?}", body.region_kind(s.region))
            };
            let faces = body
                .faces
                .values()
                .filter(|f| f.back == h || f.front == h)
                .count();
            (
                kind,
                faces,
                s.wire_edges.len(),
                usize::from(s.acorn.is_some()),
            )
        })
        .collect();
    rows.sort();
    (counts, rows)
}

fn solid_volumes(body: &BRep) -> Vec<f64> {
    let mut v: Vec<f64> = popkern_measure::region_volumes(body, &body.tolerance)
        .expect("volumes")
        .into_iter()
        .filter(|(r, _)| body.region_kind(*r) == Some(RegionKind::Solid))
        .map(|(_, v)| v)
        .collect();
    v.sort_by(f64::total_cmp);
    v
}

/// OpenCASCADE's solid volumes and distinct face, edge and vertex
/// counts from `data`; `None` without the OCP venv.
fn occ_read(data: &[u8], name: &str) -> Option<(Vec<f64>, [u64; 3])> {
    let home = std::path::PathBuf::from(std::env::var("HOME").unwrap_or_default());
    let python = home.join("venvs/ocp/bin/python");
    let oracle = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tools/occ-oracle/occ_solids.py");
    if !python.exists() || !oracle.exists() {
        return None;
    }
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("step_general_forms");
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join(format!("{name}.step"));
    std::fs::write(&path, data).expect("write");
    let out = std::process::Command::new(&python)
        .arg(&oracle)
        .arg(&path)
        .output()
        .expect("oracle");
    let report: serde_json::Value = serde_json::from_slice(&out.stdout).expect("json");
    let mut volumes: Vec<f64> = report["solids"]
        .as_array()
        .expect("solids")
        .iter()
        .map(|s| s["volume"].as_f64().expect("volume"))
        .collect();
    volumes.sort_by(f64::total_cmp);
    let count = |k: &str| report["counts"][k].as_u64().expect("count");
    Some((volumes, [count("faces"), count("edges"), count("vertices")]))
}

fn assert_close(name: &str, ours: &[f64], want: &[f64]) {
    assert_eq!(ours.len(), want.len(), "{name}: {ours:?} vs {want:?}");
    for (o, w) in ours.iter().zip(want) {
        assert!(
            (o - w).abs() < 1e-6 * w.abs().max(1.0),
            "{name}: {o} vs {w}"
        );
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
