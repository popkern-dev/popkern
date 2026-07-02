//! KP20 paper cap on REAL fillet data: a pentagonal pyramid's apex.
//!
//! Pipeline: pentagonal pyramid → G² fillet sweeps on the five apex
//! edges (`preview_fillet_with_continuity`, degree-5 curvature-
//! matched cross-sections) → apex blend boundary (2n alternating
//! fillet-end / face-edge segments) → PK-090 two-piece border
//! extraction with surface-projected ring marches → §4.4 paper cap
//! assembly → STEP for the cap, plus a WHOLE-OBJECT render: the full
//! `fillet_edges` BRep (G² options) tessellated with its own apex
//! cap faces removed and the paper cap merged in their place.
//!
//! The printed report is the honest story: the extraction evidence,
//! the input corner mismatch at `c`, and the cap's residuals against
//! the real border.

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
use popkern_fillet::{
    ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
    preview_fillet_with_continuity,
};
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

    let radius: f64 = std::env::var("KP20_APEX_RADIUS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0.12);
    println!("fillet radius {radius}");
    let specs: Vec<FilletSpec> = apex_edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        })
        .collect();
    let preview = preview_fillet_with_continuity(&pyramid, &specs, tol, ContinuityLevel::G2)
        .expect("fillet preview");
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
        "blend boundary: {} segments around vertex at z={:.3}, bridge_fallbacks={:?}, sphere_fit={}",
        boundary.segments.len(),
        boundary.vertex_position.z,
        boundary.bridge_fallbacks,
        boundary.sphere_fit_substituted,
    );

    for (k, segment) in boundary.segments.iter().enumerate() {
        use popkern_core::geometry::curve::CurveEval;
        let (label, curve) = match segment {
            popkern_fillet::blend::BoundarySegment::FilletEnd {
                fillet_index,
                curve,
            } => (format!("FilletEnd[{fillet_index}]"), curve),
            popkern_fillet::blend::BoundarySegment::FaceEdge { curve, .. } => {
                ("FaceEdge".to_string(), curve)
            }
        };
        let (t0, t1) = curve.parameter_range();
        let a = curve.point_at(t0);
        let b = curve.point_at(t1);
        println!(
            "  raw seg {k:2} {label:14} ({:6.3},{:6.3}) -> ({:6.3},{:6.3})  atan2 {:6.1} -> {:6.1}",
            a.x,
            a.y,
            b.x,
            b.y,
            a.y.atan2(a.x).to_degrees(),
            b.y.atan2(b.x).to_degrees(),
        );
    }
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
    // Calibration probe: the characteristic border's own cross-to-rim
    // leg ratios at the m and c ends (the scale the tau/b machinery
    // was built around).
    {
        let chi = popkern_core::geometry::gregory::cc_lifted_ring_border_pieces(n, |_, _| 0.0)
            .expect("chi");
        let [near_c, near_m] = &chi[0];
        let leg = |a: Point3, b: Point3| {
            ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt()
        };
        // near_c runs [m^0 -> c]: its m end is column 0; near_m runs
        // [c -> m^1]: its m end is column 3.
        for (label, rows, col) in [
            ("chi m-end", near_c, 0_usize),
            ("chi c-end (near_c)", near_c, 3),
            ("chi c-start (near_m)", near_m, 0),
        ] {
            let rim_leg = if col == 0 {
                leg(rows[0][1], rows[0][0])
            } else {
                leg(rows[0][3], rows[0][2])
            };
            let cross_leg = leg(rows[1][col], rows[0][col]);
            let rim_vec = if col == 0 {
                [
                    rows[0][1].x - rows[0][0].x,
                    rows[0][1].y - rows[0][0].y,
                    rows[0][1].z - rows[0][0].z,
                ]
            } else {
                [
                    rows[0][3].x - rows[0][2].x,
                    rows[0][3].y - rows[0][2].y,
                    rows[0][3].z - rows[0][2].z,
                ]
            };
            let cross_vec = [
                rows[1][col].x - rows[0][col].x,
                rows[1][col].y - rows[0][col].y,
                rows[1][col].z - rows[0][col].z,
            ];
            let dot =
                rim_vec[0] * cross_vec[0] + rim_vec[1] * cross_vec[1] + rim_vec[2] * cross_vec[2];
            let angle = (dot / (rim_leg * cross_leg))
                .clamp(-1.0, 1.0)
                .acos()
                .to_degrees();
            println!(
                "  {label}: rim_leg={rim_leg:.4} cross_leg={cross_leg:.4} ratio={:.3} angle_to_rim={angle:.1} deg",
                cross_leg / rim_leg,
            );
        }
    }
    // PROBE: raw geometry of sector 0's two pieces.
    for (label, rows) in [
        ("piece A (FilletEnd)", &border.sectors[0][0].rows),
        ("piece B (FaceEdge)", &border.sectors[0][1].rows),
    ] {
        let r = &rows[0];
        println!(
            "  {label} row0: ({:.3},{:.3},{:.3}) ({:.3},{:.3},{:.3}) ({:.3},{:.3},{:.3}) ({:.3},{:.3},{:.3})",
            r[0].x,
            r[0].y,
            r[0].z,
            r[1].x,
            r[1].y,
            r[1].z,
            r[2].x,
            r[2].y,
            r[2].z,
            r[3].x,
            r[3].y,
            r[3].z,
        );
    }
    // Rim turning angles at the two junction types (m = sector
    // boundaries, c = in-sector corners), from the extracted row-0
    // cubics. The characteristic layout expects ~100-115 deg turns
    // at c; near-zero turns mean the shape class fights g∘χ̃.
    for (label, junction) in [("m", 0_usize), ("c", 1_usize)] {
        let mut worst: f64 = 0.0;
        let mut sum = 0.0;
        for (k, sector) in border.sectors.iter().enumerate() {
            let (incoming, outgoing) = if junction == 1 {
                (sector[0].rows[0], sector[1].rows[0])
            } else {
                let previous = &border.sectors[(k + n - 1) % n];
                (previous[1].rows[0], sector[0].rows[0])
            };
            let ta = [
                incoming[3].x - incoming[2].x,
                incoming[3].y - incoming[2].y,
                incoming[3].z - incoming[2].z,
            ];
            let tb = [
                outgoing[1].x - outgoing[0].x,
                outgoing[1].y - outgoing[0].y,
                outgoing[1].z - outgoing[0].z,
            ];
            let na = (ta[0] * ta[0] + ta[1] * ta[1] + ta[2] * ta[2]).sqrt();
            let nb = (tb[0] * tb[0] + tb[1] * tb[1] + tb[2] * tb[2]).sqrt();
            let cosv = (ta[0] * tb[0] + ta[1] * tb[1] + ta[2] * tb[2]) / (na * nb);
            let angle = cosv.clamp(-1.0, 1.0).acos().to_degrees();
            worst = worst.max(angle);
            sum += angle;
        }
        println!(
            "  input rim turn at {label}: mean {:.1} deg, worst {worst:.1} deg",
            sum / n as f64,
        );
    }

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
    // Note: normalize-assemble-denormalize was probed and falsified —
    // the assembly is affine-clean; residuals are genuinely relative
    // to the input's smoothness class.
    let report =
        assemble_kp20_paper_cap_with_enforcement(&border, star, Kp20StripEnforcement::SoftG)
            .expect("paper cap assembly");

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

    // Whole object: run the real fillet pipeline (G² options), drop
    // its own apex cap faces from the tessellation, and merge the
    // paper cap into the hole.
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..FilletOptions::default()
    };
    let filleted = fillet_edges(&pyramid, &specs, &options, tol, None).expect("fillet_edges");
    println!(
        "fillet_edges: {} faces, {} corner patches, {} warnings",
        filleted.brep.faces.len(),
        filleted.corner_patches.len(),
        filleted.edge_warnings.len(),
    );
    let excluded: Vec<_> = filleted
        .corner_patches
        .iter()
        .map(|patch| patch.face)
        .collect();
    let mesh = popkern_mesh::tessellate(
        &filleted.brep,
        &popkern_mesh::TessellationOptions {
            chord_height: 5.0e-4,
            min_subdivisions: 12,
            max_subdivisions: 64,
        },
    )
    .expect("tessellate filleted BRep");
    let whole_path = out.join("kp20_pyramid_apex_whole.obj");
    write_whole_object_obj(&mesh, &excluded, &bi9_grids, &whole_path).expect("whole OBJ");
    println!("wrote {}", whole_path.display());
}

/// Whole-object OBJ: the tessellated filleted BRep minus the
/// pipeline's own apex cap faces, plus the paper cap patches.
fn write_whole_object_obj(
    mesh: &popkern_mesh::TriMesh,
    excluded_faces: &[popkern_core::topology::FaceHandle],
    bi9_grids: &[[[Point3; 10]; 10]],
    path: &Path,
) -> Result<(), std::io::Error> {
    use std::io::Write;
    let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        file,
        "# filleted pentagonal pyramid with the KP20 paper cap"
    )?;
    for v in &mesh.vertices {
        writeln!(file, "v {} {} {}", v.position.x, v.position.y, v.position.z)?;
        writeln!(file, "vn {} {} {}", v.normal.x, v.normal.y, v.normal.z)?;
    }
    writeln!(file, "g filleted_brep")?;
    for (triangle, face) in mesh.triangles.iter().zip(&mesh.face_map) {
        if excluded_faces.contains(face) {
            continue;
        }
        let [a, b, c] = [triangle[0] + 1, triangle[1] + 1, triangle[2] + 1];
        writeln!(file, "f {a}//{a} {b}//{b} {c}//{c}")?;
    }
    let mut base = mesh.vertices.len() + 1;
    const SAMPLES: usize = 33;
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
        writeln!(file, "g paper_cap_{index}")?;
        base = emit_grid_faces(&mut file, base, SAMPLES, SAMPLES)?;
    }
    Ok(())
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
