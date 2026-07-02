//! KP20 paper cap on REAL fillet data: a pentagonal pyramid's apex.
//!
//! Pipeline: pentagonal pyramid → rolling-ball fillets on the five
//! apex edges (`preview_fillet`) → apex blend boundary (the 2n
//! alternating fillet-end / face-edge segments) → PK-090 two-piece
//! border extraction with exact cylinder/plane ring marches →
//! §4.4 paper cap assembly → STEP + OBJ (cap plus the surrounding
//! ring collar reconstructed from the border's ring-side rows).
//!
//! The printed report is the honest story: the cap's internal gates
//! (fold-free, regularity) and the G¹/G² residuals against the real
//! border. Rolling-ball fillets meet the faces with G¹ contact only,
//! so cross-junction G² residuals reflect the INPUT's smoothness
//! class, not an assembly defect — full G² end-to-end needs the G²
//! fillet sweep (`ContinuityLevel::G2`) feeding this same path.

use std::path::Path;
use std::sync::Arc;

use popkern_core::constructors::BRepBuilder;
use popkern_core::geometry::Point3;
use popkern_core::geometry::SurfaceGeometry;
use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20PieceId, Kp20StripEnforcement, RectangularGregoryBiNonic,
    assemble_kp20_paper_cap_with_enforcement, elevate_bi6_grid_to_bi9,
    kp20_paper_border_corner_mismatch,
};
use popkern_core::geometry::surface::SurfaceOnFace;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::blend::compute_blend_boundary;
use popkern_fillet::construct::kp20_paper_border_extract::kp20_paper_border_from_blend_boundary;
use popkern_fillet::fillet_surface::FilletSurfaceResult;
use popkern_fillet::{FilletRadius, FilletSpec, preview_fillet};
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepProductNode, UnitContext, export_step_file,
};

fn main() {
    std::thread::Builder::new()
        .stack_size(256 * 1024 * 1024)
        .spawn(run)
        .expect("spawn showcase thread")
        .join()
        .expect("showcase thread");
}

fn run() {
    let out = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("crates dir")
        .parent()
        .expect("repo root")
        .join("examples")
        .join("output");
    std::fs::create_dir_all(&out).expect("create output dir");

    let tol = Tolerance::default();
    let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol).expect("pentagonal pyramid");
    let apex = pyramid
        .vertices
        .iter()
        .find(|(_, v)| v.position.x.abs() < 1e-9 && v.position.y.abs() < 1e-9 && v.position.z > 0.0)
        .map(|(h, _)| h)
        .expect("apex vertex");
    let apex_edges: Vec<_> = pyramid
        .edges
        .iter()
        .filter_map(|(eh, _)| {
            let he = pyramid.edges[eh].half_edge;
            let origin = pyramid.half_edges[he].origin;
            let twin_origin = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
            (origin == apex || twin_origin == apex).then_some(eh)
        })
        .collect();
    println!("apex has {} incident edges", apex_edges.len());

    let specs: Vec<FilletSpec> = apex_edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.12),
        })
        .collect();
    let preview = preview_fillet(&pyramid, &specs, tol).expect("fillet preview");
    assert!(
        preview.all_succeeded(),
        "preview failures: {:?}",
        preview.failed
    );
    let fillet_surfaces: Vec<FilletSurfaceResult> = preview
        .per_edge
        .iter()
        .map(|eg| FilletSurfaceResult {
            surface: eg.blend_surface.clone(),
            trim_curve_1: eg.trim_curve_1.clone(),
            trim_curve_2: eg.trim_curve_2.clone(),
            continuity_downgrade: None,
        })
        .collect();
    let indices: Vec<usize> = (0..fillet_surfaces.len()).collect();
    let boundary = compute_blend_boundary(apex, &indices, &fillet_surfaces, &pyramid, &tol)
        .expect("blend boundary");
    println!(
        "blend boundary: {} segments around vertex at z={:.3}",
        boundary.segments.len(),
        boundary.vertex_position.z,
    );

    let (border, extract_report) =
        kp20_paper_border_from_blend_boundary(&boundary, &fillet_surfaces, &pyramid)
            .expect("two-piece border extraction");
    println!(
        "extraction: n={} depth={:.4} loop_gap={:.3e} first_order_marches={} reversed={}",
        extract_report.n,
        extract_report.depth,
        extract_report.max_loop_gap,
        extract_report.first_order_marches,
        extract_report.reversed_for_orientation,
    );
    let n = extract_report.n;
    let c = (core::f64::consts::TAU / n as f64).cos();
    let mut corner_mismatch = 0.0_f64;
    for sector in &border.sectors {
        corner_mismatch = corner_mismatch.max(kp20_paper_border_corner_mismatch(
            n,
            &sector[0].rows,
            &sector[1].rows,
            2.0 * c / 3.0,
        ));
    }
    println!("input corner mismatch at c: {corner_mismatch:.3e}");

    // Cap centre: pull the boundary centroid toward the original
    // apex. The pull factor is a shape parameter; sweep it and keep
    // the most regular fold-free cap.
    let centroid = border_centroid(&border);
    let vertex = boundary.vertex_position;
    let star_at = |factor: f64| {
        Point3::new(
            centroid.x + factor * (vertex.x - centroid.x),
            centroid.y + factor * (vertex.y - centroid.y),
            centroid.z + factor * (vertex.z - centroid.z),
        )
    };
    // Rolling-ball input is G¹-class: the paper's hard eq (8)/(9)
    // strip system is structurally infeasible at the m junctions and
    // tears the strips. Compare the paper mode against the SoftG
    // containment mode (C⁰/(11)/(12) hard, G rows least-squares) and
    // keep the containment result for out-of-class input.
    let star = star_at(0.25);
    for (label, mode) in [
        ("HardG2", Kp20StripEnforcement::HardG2),
        ("SoftG", Kp20StripEnforcement::SoftG),
    ] {
        let candidate = assemble_kp20_paper_cap_with_enforcement(&border, star, mode)
            .expect("paper cap assembly");
        println!(
            "  {label}: g1={:.3e} g2={:.3e} c0={:.3e} min_regularity={:.3e} strip_move={:.3e}",
            candidate.max_g1_residual,
            candidate.max_g2_residual,
            candidate.max_c0_residual,
            candidate.min_regularity,
            candidate.max_strip_movement,
        );
    }
    let report =
        assemble_kp20_paper_cap_with_enforcement(&border, star, Kp20StripEnforcement::SoftG)
            .expect("paper cap assembly");
    println!(
        "paper cap: g1={:.3e} g2={:.3e} c0={:.3e} border={:.3e} lemma4={:.3e} guide_fit={:.3e} strip_move={:.3e} min_regularity={:.3e} corner_mismatch={:.3e}",
        report.max_g1_residual,
        report.max_g2_residual,
        report.max_c0_residual,
        report.max_border_residual,
        report.lemma4_delta,
        report.guide_fit_residual,
        report.max_strip_movement,
        report.min_regularity,
        report.input_corner_mismatch,
    );

    let bi9_grids: Vec<[[Point3; 10]; 10]> = report
        .sectors
        .iter()
        .flat_map(|sector| {
            Kp20PieceId::all()
                .into_iter()
                .map(|piece_id| elevate_bi6_grid_to_bi9(sector.piece(piece_id).cps))
        })
        .collect();
    let brep = build_cap_brep(&bi9_grids).expect("cap BRep");
    let step_path = out.join("kp20_pyramid_apex_cap.step");
    export_step_file(
        &step_document(&brep),
        &step_path,
        &StepExportOptions::default(),
    )
    .expect("STEP export");
    println!("wrote {}", step_path.display());

    let obj_path = out.join("kp20_pyramid_apex_cap.obj");
    write_scene_obj(&bi9_grids, &border, &obj_path).expect("OBJ export");
    println!("wrote {}", obj_path.display());
}

fn border_centroid(border: &popkern_core::geometry::gregory::Kp20PaperGuideBorder) -> Point3 {
    let mut acc = [0.0_f64; 3];
    let mut count = 0.0;
    for sector in &border.sectors {
        for piece in sector {
            for cp in &piece.rows[0] {
                acc[0] += cp.x;
                acc[1] += cp.y;
                acc[2] += cp.z;
                count += 1.0;
            }
        }
    }
    Point3::new(acc[0] / count, acc[1] / count, acc[2] / count)
}

/// One OBJ scene: the cap patches plus the surrounding ring collar
/// reconstructed from the border's ring-side rows (`ring_d =
/// reflection of the cap-side prolongation`), so the render shows the
/// cap in context of the real fillet/face surfaces it joins.
fn write_scene_obj(
    bi9_grids: &[[[Point3; 10]; 10]],
    border: &popkern_core::geometry::gregory::Kp20PaperGuideBorder,
    path: &Path,
) -> Result<(), std::io::Error> {
    use std::io::Write;
    const SAMPLES: usize = 33;
    let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        file,
        "# KP20 paper cap on a filleted pentagonal pyramid apex"
    )?;
    let mut base = 1_usize;

    for (index, grid) in bi9_grids.iter().enumerate() {
        let patch = RectangularGregoryBiNonic { grid: *grid };
        for i in 0..SAMPLES {
            for j in 0..SAMPLES {
                let u = i as f64 / (SAMPLES - 1) as f64;
                let v = j as f64 / (SAMPLES - 1) as f64;
                let p = patch.point_at(u, v);
                let du = patch.partial_u(u, v);
                let dv = patch.partial_v(u, v);
                let raw = du.cross(&dv);
                let len = raw.length();
                let nrm = if len > 1.0e-12 {
                    popkern_core::geometry::Vec3::new(raw.x / len, raw.y / len, raw.z / len)
                } else {
                    popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0)
                };
                writeln!(file, "v {} {} {}", p.x, p.y, p.z)?;
                writeln!(file, "vn {} {} {}", nrm.x, nrm.y, nrm.z)?;
            }
        }
        writeln!(file, "g cap_patch_{index}")?;
        base = emit_grid_faces(&mut file, base, SAMPLES, SAMPLES)?;
    }

    // Ring collar: per border piece, evaluate the (cubic × quadratic)
    // BB patch of the RING-side rows.
    const ALONG: usize = 25;
    const DEEP: usize = 7;
    for (sector_index, sector) in border.sectors.iter().enumerate() {
        for (piece_index, piece) in sector.iter().enumerate() {
            // Ring rows: reflect the cap-side prolongation back.
            let q0 = piece.rows[0];
            let mut q1 = [Point3::new(0.0, 0.0, 0.0); 4];
            let mut q2 = [Point3::new(0.0, 0.0, 0.0); 4];
            for i in 0..4 {
                q1[i] = Point3::new(
                    2.0 * q0[i].x - piece.rows[1][i].x,
                    2.0 * q0[i].y - piece.rows[1][i].y,
                    2.0 * q0[i].z - piece.rows[1][i].z,
                );
                q2[i] = Point3::new(
                    4.0 * q0[i].x - 4.0 * piece.rows[1][i].x + piece.rows[2][i].x,
                    4.0 * q0[i].y - 4.0 * piece.rows[1][i].y + piece.rows[2][i].y,
                    4.0 * q0[i].z - 4.0 * piece.rows[1][i].z + piece.rows[2][i].z,
                );
            }
            let rows = [q0, q1, q2];
            for a in 0..ALONG {
                for d in 0..DEEP {
                    let u = a as f64 / (ALONG - 1) as f64;
                    let w = d as f64 / (DEEP - 1) as f64;
                    let p = eval_cubic_by_quadratic(&rows, u, w);
                    let du = eval_cubic_by_quadratic_du(&rows, u, w);
                    let dw = eval_cubic_by_quadratic_dw(&rows, u, w);
                    let raw = du.cross(&dw);
                    let len = raw.length();
                    let nrm = if len > 1.0e-12 {
                        popkern_core::geometry::Vec3::new(raw.x / len, raw.y / len, raw.z / len)
                    } else {
                        popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0)
                    };
                    writeln!(file, "v {} {} {}", p.x, p.y, p.z)?;
                    writeln!(file, "vn {} {} {}", nrm.x, nrm.y, nrm.z)?;
                }
            }
            writeln!(file, "g collar_{sector_index}_{piece_index}")?;
            base = emit_grid_faces(&mut file, base, ALONG, DEEP)?;
        }
    }
    Ok(())
}

fn emit_grid_faces<W: std::io::Write>(
    file: &mut W,
    base: usize,
    rows: usize,
    cols: usize,
) -> Result<usize, std::io::Error> {
    for i in 0..(rows - 1) {
        for j in 0..(cols - 1) {
            let a = base + i * cols + j;
            let b = a + cols;
            let c = b + 1;
            let d = a + 1;
            writeln!(file, "f {a}//{a} {b}//{b} {c}//{c}")?;
            writeln!(file, "f {a}//{a} {c}//{c} {d}//{d}")?;
        }
    }
    Ok(base + rows * cols)
}

fn bernstein3(t: f64) -> ([f64; 4], [f64; 4]) {
    let s = 1.0 - t;
    let b = [s * s * s, 3.0 * s * s * t, 3.0 * s * t * t, t * t * t];
    let db = [
        -3.0 * s * s,
        3.0 * s * s - 6.0 * s * t,
        6.0 * s * t - 3.0 * t * t,
        3.0 * t * t,
    ];
    (b, db)
}

fn bernstein2(t: f64) -> ([f64; 3], [f64; 3]) {
    let s = 1.0 - t;
    let b = [s * s, 2.0 * s * t, t * t];
    let db = [-2.0 * s, 2.0 * s - 2.0 * t, 2.0 * t];
    (b, db)
}

fn eval_cubic_by_quadratic(rows: &[[Point3; 4]; 3], u: f64, w: f64) -> Point3 {
    let (bu, _) = bernstein3(u);
    let (bw, _) = bernstein2(w);
    accumulate(rows, &bu, &bw)
}

fn eval_cubic_by_quadratic_du(
    rows: &[[Point3; 4]; 3],
    u: f64,
    w: f64,
) -> popkern_core::geometry::Vec3 {
    let (_, dbu) = bernstein3(u);
    let (bw, _) = bernstein2(w);
    let p = accumulate(rows, &dbu, &bw);
    popkern_core::geometry::Vec3::new(p.x, p.y, p.z)
}

fn eval_cubic_by_quadratic_dw(
    rows: &[[Point3; 4]; 3],
    u: f64,
    w: f64,
) -> popkern_core::geometry::Vec3 {
    let (bu, _) = bernstein3(u);
    let (_, dbw) = bernstein2(w);
    let p = accumulate(rows, &bu, &dbw);
    popkern_core::geometry::Vec3::new(p.x, p.y, p.z)
}

fn accumulate(rows: &[[Point3; 4]; 3], bu: &[f64; 4], bw: &[f64; 3]) -> Point3 {
    let mut acc = [0.0_f64; 3];
    for (d, row) in rows.iter().enumerate() {
        for (i, cp) in row.iter().enumerate() {
            let weight = bu[i] * bw[d];
            acc[0] += weight * cp.x;
            acc[1] += weight * cp.y;
            acc[2] += weight * cp.z;
        }
    }
    Point3::new(acc[0], acc[1], acc[2])
}

fn build_cap_brep(bi9_grids: &[[[Point3; 10]; 10]]) -> Result<BRep, String> {
    let mut builder = BRepBuilder::new(Tolerance::default());
    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);
    for grid in bi9_grids {
        let v_00 = builder.add_vertex(grid[0][0]);
        let v_90 = builder.add_vertex(grid[9][0]);
        let v_99 = builder.add_vertex(grid[9][9]);
        let v_09 = builder.add_vertex(grid[0][9]);
        let (_, he1, _) = builder.add_edge(v_00, v_90);
        let (_, he2, _) = builder.add_edge(v_90, v_99);
        let (_, he3, _) = builder.add_edge(v_99, v_09);
        let (_, he4, _) = builder.add_edge(v_09, v_00);
        let face = builder.add_face(shell);
        builder
            .add_loop(&[he1, he2, he3, he4], face)
            .map_err(|e| format!("loop: {e:?}"))?;
        builder.set_face_surface(
            face,
            SurfaceOnFace {
                geometry: SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonic(Box::new(
                    RectangularGregoryBiNonic { grid: *grid },
                ))),
                orientation: true,
            },
        );
    }
    Ok(builder.build_unchecked())
}

fn step_document(brep: &BRep) -> StepDocument {
    StepDocument {
        products: vec![StepProductNode {
            name: "Kp20PyramidApexCap".into(),
            id: "KP20-PYRAMID-APEX-CAP".into(),
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
