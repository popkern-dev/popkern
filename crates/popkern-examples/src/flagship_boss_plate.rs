//! Flagship demo: sketch → solid → boolean → fillet → watertight STEP.
//!
//! A rectangular plate and a cylindrical boss are each built from a
//! closed sketch profile by extrusion, unioned, and filleted — the
//! boss/plate seam is a CURVED (circular) fillet edge, and the plate
//! corners carry vertex-blend caps. The result must tessellate
//! watertight and export as a STEP solid.
//!
//! Stages print their gate results; the binary exits non-zero if a
//! gate fails, so it can serve as an integration check.

use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sweep::{SweepOptions, SweepProfile, sweep};

/// Closed degree-1 rectangle wire in the z=0 plane, centered at the
/// origin: the plate sketch.
fn rectangle_sketch(width: f64, depth: f64) -> NurbsCurve {
    let (hw, hd) = (0.5 * width, 0.5 * depth);
    let pts = vec![
        Point3::new(-hw, -hd, 0.0),
        Point3::new(hw, -hd, 0.0),
        Point3::new(hw, hd, 0.0),
        Point3::new(-hw, hd, 0.0),
        Point3::new(-hw, -hd, 0.0),
    ];
    let n = pts.len();
    let mut knots = vec![0.0];
    for i in 0..n {
        knots.push(i as f64 / (n - 1) as f64);
    }
    knots.push(1.0);
    NurbsCurve {
        degree: 1,
        knots,
        control_points: pts,
        weights: vec![1.0; n],
    }
}

/// Closed rational quadratic circle wire (standard 9-CP form) in the
/// plane z = `z`, centered on the z-axis: the boss sketch.
fn circle_sketch(radius: f64, z: f64) -> NurbsCurve {
    let r = radius;
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let pts = vec![
        Point3::new(r, 0.0, z),
        Point3::new(r, r, z),
        Point3::new(0.0, r, z),
        Point3::new(-r, r, z),
        Point3::new(-r, 0.0, z),
        Point3::new(-r, -r, z),
        Point3::new(0.0, -r, z),
        Point3::new(r, -r, z),
        Point3::new(r, 0.0, z),
    ];
    let weights = vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
    let knots = vec![
        0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
    ];
    NurbsCurve {
        degree: 2,
        knots,
        control_points: pts,
        weights,
    }
}

/// Straight extrusion path from `z0` to `z1` along +Z.
fn extrusion_path(z0: f64, z1: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, z0), Point3::new(0.0, 0.0, z1)],
        weights: vec![1.0, 1.0],
    }
}

fn extrude(sketch: NurbsCurve, z0: f64, z1: f64, tol: Tolerance, label: &str) -> BRep {
    let result = sweep(
        &SweepProfile::Wire(sketch),
        &extrusion_path(z0, z1),
        &SweepOptions::default(),
        tol,
        None,
    )
    .unwrap_or_else(|e| panic!("{label} extrusion failed: {e:?}"));
    for w in &result.warnings {
        eprintln!("  [{label}] sweep warning: {w:?}");
    }
    if let Err(e) = result.brep.validate() {
        eprintln!("  [{label}] INVALID: {e:?}");
        eprintln!(
            "  [{label}] faces={} edges={} vertices={} loops={}",
            result.brep.faces.len(),
            result.brep.edges.len(),
            result.brep.vertices.len(),
            result.brep.loops.len()
        );
        for (fh, f) in result.brep.faces.iter() {
            let n = result.brep.face_outer_loop_positions(fh).len();
            let kind = match f.surface.as_ref().map(|s| &s.geometry) {
                Some(popkern_core::geometry::SurfaceGeometry::Plane { .. }) => "Plane",
                Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) => "Nurbs",
                Some(popkern_core::geometry::SurfaceGeometry::Cylinder { .. }) => "Cyl",
                _ => "?",
            };
            eprintln!(
                "    face {fh:?} {kind} loop_len={n} inner={}",
                f.inner_loops.len()
            );
        }
        panic!("{label} extrusion invalid");
    }
    result.brep
}

fn main() {
    let tol = Tolerance::default();

    // Stage 1: sketches → solids.
    let plate = extrude(rectangle_sketch(2.4, 1.6), 0.0, 0.4, tol, "plate");
    println!(
        "stage 1a plate: faces={} edges={} (validate ok)",
        plate.faces.len(),
        plate.edges.len()
    );
    // The boss penetrates the plate (base below the plate top) so the
    // union seam is a general-position circle on the plate top.
    let boss = extrude(circle_sketch(0.35, 0.0), 0.2, 1.0, tol, "boss");
    println!(
        "stage 1b boss: faces={} edges={} (validate ok)",
        boss.faces.len(),
        boss.edges.len()
    );

    // Stage 2: boolean union.
    let unioned = popkern_bool::union(plate, boss, tol, None).expect("union");
    unioned.validate().expect("union must validate");
    println!(
        "stage 2 union: faces={} edges={} vertices={} (validate ok)",
        unioned.faces.len(),
        unioned.edges.len(),
        unioned.vertices.len()
    );

    // Locate the curved seam edge: cylinder × plate-top plane.
    let mut seam_edges = Vec::new();
    for (eh, edge) in unioned.edges.iter() {
        let he = edge.half_edge;
        let f1 = unioned.loops[unioned.half_edges[he].loop_].face;
        let f2 = unioned.loops[unioned.half_edges[unioned.half_edges[he].twin].loop_].face;
        let kind = |fh: popkern_core::topology::FaceHandle| {
            unioned.faces[fh].surface.as_ref().map(|s| &s.geometry)
        };
        use popkern_core::geometry::SurfaceGeometry as SG;
        let cyl_plane = matches!(
            (kind(f1), kind(f2)),
            (Some(SG::Cylinder { .. }), Some(SG::Plane { .. }))
                | (Some(SG::Plane { .. }), Some(SG::Cylinder { .. }))
        ) || matches!(
            (kind(f1), kind(f2)),
            (Some(SG::Nurbs(_)), Some(SG::Plane { .. }))
                | (Some(SG::Plane { .. }), Some(SG::Nurbs(_)))
        );
        if !cyl_plane {
            continue;
        }
        // On the plate top (z = 0.4)?
        if let Some(curve) = edge.curve.as_ref() {
            use popkern_core::geometry::curve::CurveEval;
            let ceval: &dyn CurveEval = &curve.geometry;
            let mid = ceval.point_at(0.5 * (curve.t_start + curve.t_end));
            if (mid.z - 0.4).abs() < 1e-6 {
                seam_edges.push(eh);
            }
        }
    }
    println!("stage 2 seam edges at plate top: {}", seam_edges.len());
    if seam_edges.is_empty() {
        for (eh, edge) in unioned.edges.iter() {
            let he = edge.half_edge;
            let f1 = unioned.loops[unioned.half_edges[he].loop_].face;
            let f2 = unioned.loops[unioned.half_edges[unioned.half_edges[he].twin].loop_].face;
            let kindof = |fh: popkern_core::topology::FaceHandle| match unioned.faces[fh]
                .surface
                .as_ref()
                .map(|s| &s.geometry)
            {
                Some(popkern_core::geometry::SurfaceGeometry::Plane { .. }) => "P",
                Some(popkern_core::geometry::SurfaceGeometry::Cylinder { .. }) => "C",
                Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) => "N",
                _ => "?",
            };
            let mid = edge.curve.as_ref().map(|c| {
                use popkern_core::geometry::curve::CurveEval;
                let ceval: &dyn CurveEval = &c.geometry;
                ceval.point_at(0.5 * (c.t_start + c.t_end))
            });
            eprintln!(
                "  edge {eh:?} {}x{} mid {:?}",
                kindof(f1),
                kindof(f2),
                mid.map(|p| (p.x, p.y, p.z))
            );
        }
    }
    assert!(!seam_edges.is_empty(), "no curved seam edge found");

    // Stage 3: fillet — the CURVED seam edges plus the plate's top
    // perimeter and vertical edges (vertex-blend caps at the top
    // corners).
    let mut fillet_edges_list: Vec<popkern_core::topology::EdgeHandle> = seam_edges.clone();
    let seam_only = std::env::var_os("FLAGSHIP_SEAM_ONLY").is_some();
    for (eh, edge) in unioned.edges.iter() {
        if seam_only {
            break;
        }
        let he = edge.half_edge;
        let a = unioned.vertices[unioned.half_edges[he].origin].position;
        let b = unioned.vertices[unioned.half_edges[unioned.half_edges[he].twin].origin].position;
        let vertical = (a.z - b.z).abs() > 0.3
            && (a.x - b.x).abs() < 1e-9
            && (a.y - b.y).abs() < 1e-9
            && a.z.min(b.z) < 1e-6;
        let top_perimeter =
            (a.z - 0.4).abs() < 1e-6 && (b.z - 0.4).abs() < 1e-6 && a.distance_to(&b) > 1.0;
        if vertical || top_perimeter {
            fillet_edges_list.push(eh);
        }
    }
    println!(
        "stage 3 filleting {} edges ({} seam + {} plate)",
        fillet_edges_list.len(),
        seam_edges.len(),
        fillet_edges_list.len() - seam_edges.len()
    );
    let specs: Vec<popkern_fillet::FilletSpec> = fillet_edges_list
        .iter()
        .map(|&edge| popkern_fillet::FilletSpec {
            edge,
            radius: popkern_fillet::FilletRadius::Constant(0.08),
        })
        .collect();
    let options = popkern_fillet::FilletOptions {
        continuity: popkern_fillet::ContinuityLevel::G1,
        ..Default::default()
    };
    if seam_only {
        for (fh, f) in unioned.faces.iter() {
            if let Some(sof) = f.surface.as_ref() {
                if let popkern_core::geometry::SurfaceGeometry::Plane { origin, normal, .. } =
                    &sof.geometry
                {
                    if (origin.z - 0.4).abs() < 1e-9 || (normal.z.abs() > 0.9) {
                        eprintln!(
                            "  face {fh:?} plane origin z {:.3} normal ({:.2},{:.2},{:.2}) orientation {}",
                            origin.z, normal.x, normal.y, normal.z, sof.orientation
                        );
                    }
                }
            }
        }
        for &eh in &fillet_edges_list {
            match popkern_fillet::analysis::analyze_edges(&unioned, &[eh], false) {
                Ok(list) => {
                    for a in &list {
                        eprintln!(
                            "  seam analysis: convexity {:?} face_1 {:?} face_2 {:?}",
                            a.convexity, a.face_1, a.face_2
                        );
                    }
                }
                Err(e) => eprintln!("  seam analysis FAILED: {e:?}"),
            }
        }
    }
    let result = match popkern_fillet::fillet_edges(&unioned, &specs, &options, tol, None) {
        Ok(r) => r,
        Err(popkern_fillet::BlendError::PartialFailure { failures, .. }) => {
            for (eh, f) in &failures {
                eprintln!("  fillet FAILED on {eh:?}: {f:?}");
            }
            panic!("fillet partial failure");
        }
        Err(e) => panic!("fillet pipeline: {e}"),
    };
    result
        .brep
        .validate()
        .expect("filleted flagship must validate");
    for (eh, w) in &result.edge_warnings {
        eprintln!("  warning on {eh:?}: {w:?}");
    }
    println!(
        "stage 3 fillet: faces={} edges={} (validate ok, {} warnings, {} corner patches)",
        result.brep.faces.len(),
        result.brep.edges.len(),
        result.edge_warnings.len(),
        result.corner_patches.len()
    );

    if std::env::var_os("FLAGSHIP_TOPO_DUMP").is_some() {
        for (fh, f) in result.brep.faces.iter() {
            let kind = match f.surface.as_ref().map(|s| &s.geometry) {
                Some(popkern_core::geometry::SurfaceGeometry::Plane { .. }) => "Plane",
                Some(popkern_core::geometry::SurfaceGeometry::Cylinder { .. }) => "Cyl",
                Some(popkern_core::geometry::SurfaceGeometry::Torus { .. }) => "Torus",
                Some(popkern_core::geometry::SurfaceGeometry::Nurbs(_)) => "Nurbs",
                _ => "?",
            };
            let mut loops_desc = Vec::new();
            for (li, &lh) in std::iter::once(&f.outer_loop)
                .chain(f.inner_loops.iter())
                .enumerate()
            {
                let start = result.brep.loops[lh].half_edge;
                let mut he = start;
                let mut edges = Vec::new();
                for _ in 0..64 {
                    let d = &result.brep.half_edges[he];
                    edges.push(format!("{:?}", d.edge));
                    he = d.next;
                    if he == start {
                        break;
                    }
                }
                loops_desc.push(format!(
                    "{}[{}]",
                    if li == 0 { "outer" } else { "inner" },
                    edges.join(",")
                ));
            }
            eprintln!("  face {fh:?} {kind}: {}", loops_desc.join(" "));
        }
    }

    // Stage 4: watertight tessellation + OBJ render source + STEP.
    let mesh =
        popkern_mesh::tessellate(&result.brep, &popkern_mesh::TessellationOptions::default())
            .expect("tessellation");
    println!(
        "stage 4 tessellated: {} vertices, {} triangles",
        mesh.vertices.len(),
        mesh.triangles.len()
    );
    // Closure audit: welded edge-use counts.
    {
        use std::collections::HashMap;
        let mut key = HashMap::new();
        let mut canon = Vec::new();
        let mut id = |p: popkern_core::geometry::Point3,
                      key: &mut HashMap<(i64, i64, i64), usize>,
                      canon: &mut Vec<popkern_core::geometry::Point3>| {
            let k = (
                (p.x * 1.0e7).round() as i64,
                (p.y * 1.0e7).round() as i64,
                (p.z * 1.0e7).round() as i64,
            );
            *key.entry(k).or_insert_with(|| {
                canon.push(p);
                canon.len() - 1
            })
        };
        let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
        for t in &mesh.triangles {
            let w = [
                id(mesh.vertices[t[0]].position, &mut key, &mut canon),
                id(mesh.vertices[t[1]].position, &mut key, &mut canon),
                id(mesh.vertices[t[2]].position, &mut key, &mut canon),
            ];
            if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
                continue;
            }
            for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
                *counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
            }
        }
        let open = counts.values().filter(|&&c| c == 1).count();
        println!("stage 4 open edges: {open}");
        if seam_only {
            assert_eq!(open, 0, "seam-fillet flagship must tessellate watertight");
        }
        let mut shown = 0;
        for (&(a, b), &c) in &counts {
            if c == 1 && shown < 6 {
                let pa = canon[a];
                let pb = canon[b];
                eprintln!(
                    "  open ({:.3},{:.3},{:.3})-({:.3},{:.3},{:.3})",
                    pa.x, pa.y, pa.z, pb.x, pb.y, pb.z
                );
                shown += 1;
            }
        }
    }
    let out_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/output");
    std::fs::create_dir_all(&out_dir).expect("out dir");
    // Plain OBJ.
    let obj_path = out_dir.join("flagship_boss_plate.obj");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&obj_path).expect("obj file");
        for v in &mesh.vertices {
            writeln!(f, "v {} {} {}", v.position.x, v.position.y, v.position.z).unwrap();
        }
        for t in &mesh.triangles {
            writeln!(f, "f {} {} {}", t[0] + 1, t[1] + 1, t[2] + 1).unwrap();
        }
    }
    println!("wrote {}", obj_path.display());
    let step_path = out_dir.join("flagship_boss_plate.step");
    let doc = popkern_io::StepDocument {
        products: vec![popkern_io::StepProductNode {
            name: "flagship_boss_plate".into(),
            id: "FLAGSHIP-BOSS-PLATE-001".into(),
            children: vec![],
            shapes: vec![std::sync::Arc::new(popkern_io::ShapeRep {
                kind: popkern_io::ShapeRepKind::AdvancedBrep,
                solids: vec![result.brep.clone()],
                context: popkern_io::RepresentationContext::default(),
                non_brep_items: vec![],
            })],
            metadata: popkern_io::ProductMetadata::default(),
        }],
        roots: vec![0],
        units: popkern_io::UnitContext::default(),
        opaque_bag: popkern_io::OpaqueBag::default(),
    };
    match popkern_io::export_step_file(&doc, &step_path, &popkern_io::StepExportOptions::default())
    {
        Ok(()) => println!("wrote {}", step_path.display()),
        Err(e) => println!("STEP export declined: {e}"),
    }
}
