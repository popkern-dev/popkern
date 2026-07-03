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
    assemble_kp20_paper_cap_conditioned, assemble_kp20_paper_cap_with_enforcement,
    elevate_bi6_grid_to_bi9, kp20_paper_border_corner_mismatch,
};
use popkern_core::geometry::surface::SurfaceOnFace;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::blend::compute_blend_boundary;
use popkern_fillet::construct::kp20_paper_border_extract::{
    kp20_paper_border_from_blend_boundary, kp20_paper_border_from_fillet_setback,
};
use popkern_fillet::fillet_surface::FilletSurfaceResult;
use popkern_fillet::{ContinuityLevel, FilletRadius, FilletSpec, preview_fillet_with_continuity};
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
        .unwrap_or(0.16);
    println!("fillet radius {radius}");
    let specs: Vec<FilletSpec> = apex_edges
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        })
        .collect();
    let continuity = if std::env::var("KP20_G1").map(|v| v == "1").unwrap_or(false) {
        ContinuityLevel::G1
    } else {
        ContinuityLevel::G2
    };
    let preview =
        preview_fillet_with_continuity(&pyramid, &specs, tol, continuity).expect("fillet preview");
    assert!(
        preview.all_succeeded(),
        "preview failures: {:?}",
        preview.failed
    );
    for (i, eg) in preview.per_edge.iter().enumerate() {
        let variant = match &eg.blend_surface {
            popkern_core::geometry::SurfaceGeometry::Cylinder { .. } => "Cylinder",
            popkern_core::geometry::SurfaceGeometry::Nurbs(_) => "Nurbs",
            popkern_core::geometry::SurfaceGeometry::Torus { .. } => "Torus",
            _ => "other",
        };
        println!("  fillet[{i}] surface variant: {variant}");
    }
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
    // Side check: a true inscribed fillet's mid cross-section lies
    // INSIDE the original wedge (negative height under both adjacent
    // face planes). Positive = the sweep bulges outward.
    {
        use popkern_core::geometry::surface::SurfaceEval;
        let interior = Point3::new(0.0, 0.0, 0.35);
        let planes: Vec<(Point3, popkern_core::geometry::Vec3)> = pyramid
            .faces
            .iter()
            .filter_map(|(_, face)| {
                let sof = face.surface.as_ref()?;
                if let popkern_core::geometry::SurfaceGeometry::Plane { normal, origin, .. } =
                    &sof.geometry
                {
                    if normal.z.abs() > 0.9 {
                        return None;
                    }
                    let n = if normal.dot(&(*origin - interior)) >= 0.0 {
                        *normal
                    } else {
                        -*normal
                    };
                    Some((*origin, n))
                } else {
                    None
                }
            })
            .collect();
        for (k, fillet) in fillet_surfaces.iter().enumerate() {
            let mid = fillet.surface.point_at(0.5, 0.5);
            let mut heights: Vec<f64> = planes.iter().map(|(o, n)| (mid - *o).dot(n)).collect();
            heights.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
            // Tangency check: the sweep normal at the contact lines
            // must match the face normals (G1 flush).
            let mut worst_tangency = 0.0_f64;
            for &v_edge in &[0.0_f64, 1.0] {
                for step in 0..5 {
                    let u = 0.1 + 0.2 * step as f64;
                    let p_edge = fillet.surface.point_at(u, v_edge);
                    let n_sweep = fillet.surface.normal_at(u, v_edge);
                    // Nearest face plane normal.
                    let (_, n_face) = planes
                        .iter()
                        .min_by(|a, b| {
                            let da = ((p_edge - a.0).dot(&a.1)).abs();
                            let db = ((p_edge - b.0).dot(&b.1)).abs();
                            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                        })
                        .copied()
                        .unwrap();
                    let align = n_sweep.dot(&n_face).abs();
                    worst_tangency = worst_tangency.max(1.0 - align);
                }
            }
            // Finite-difference normal at the contact (trust position
            // eval over normal_at).
            {
                let e = 1.0e-4;
                let p0 = fillet.surface.point_at(0.5, 0.0);
                let du = (fillet.surface.point_at(0.5 + e, 0.0)
                    - fillet.surface.point_at(0.5 - e, 0.0))
                    * (0.5 / e);
                let dv = (fillet.surface.point_at(0.5, e) - p0) * (1.0 / e);
                let n_fd = du.cross(&dv).normalized().unwrap_or(du);
                let (_, n_face) = planes
                    .iter()
                    .min_by(|a, b| {
                        let da = ((p0 - a.0).dot(&a.1)).abs();
                        let db = ((p0 - b.0).dot(&b.1)).abs();
                        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                    })
                    .copied()
                    .unwrap();
                let n_eval = fillet.surface.normal_at(0.5, 0.0);
                println!(
                    "  fillet[{k}] FD-normal vs face: dev {:.2e}; normal_at vs FD: dev {:.2e}",
                    1.0 - n_fd.dot(&n_face).abs(),
                    1.0 - n_eval.normalized().unwrap_or(n_eval).dot(&n_fd).abs(),
                );
            }
            let d_v0 =
                trim_curve_distance(&fillet.trim_curve_1, fillet.surface.point_at(0.5, 0.0)).min(
                    trim_curve_distance(&fillet.trim_curve_2, fillet.surface.point_at(0.5, 0.0)),
                );
            let d_v1 =
                trim_curve_distance(&fillet.trim_curve_1, fillet.surface.point_at(0.5, 1.0)).min(
                    trim_curve_distance(&fillet.trim_curve_2, fillet.surface.point_at(0.5, 1.0)),
                );
            println!(
                "  fillet[{k}] mid height: {:.4}, {:.4}; tangency dev: {:.2e}; v-edge-to-trim: {:.4}, {:.4}",
                heights[0], heights[1], worst_tangency, d_v0, d_v1,
            );
        }
    }
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
    let use_setback = std::env::var("KP20_SETBACK")
        .map(|v| v == "1")
        .unwrap_or(true);
    let (border, extract_report) = if use_setback {
        kp20_paper_border_from_fillet_setback(&fillet_surfaces, boundary.vertex_position)
            .expect("setback border extraction")
    } else {
        kp20_paper_border_from_blend_boundary(&boundary, &fillet_surfaces, &pyramid)
            .expect("two-piece border extraction")
    };
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
    for (label, rows) in [
        ("piece A", &border.sectors[0][0].rows),
        ("piece B", &border.sectors[0][1].rows),
    ] {
        for col in [0_usize, 3] {
            let leg = Point3::new(
                rows[1][col].x - rows[0][col].x,
                rows[1][col].y - rows[0][col].y,
                rows[1][col].z - rows[0][col].z,
            );
            println!(
                "  {label} cross leg at col {col}: ({:.5}, {:.5}, {:.5}) len {:.5}",
                leg.x,
                leg.y,
                leg.z,
                (leg.x * leg.x + leg.y * leg.y + leg.z * leg.z).sqrt(),
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
    let star_factor: f64 = std::env::var("KP20_STAR_FACTOR")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0.5);
    let star = star_at(star_factor);
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
    // The setback border satisfies the paper's premises well enough
    // for the paper-faithful hard-constraint mode.
    let fairing: f64 = std::env::var("KP20_GUIDE_FAIRING")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(0.3);
    let report =
        assemble_kp20_paper_cap_conditioned(&border, star, Kp20StripEnforcement::HardG2, fairing)
            .expect("paper cap assembly");

    {
        use popkern_core::geometry::gregory::Kp20PieceId;
        let p4 = report.sectors[0].piece(Kp20PieceId::NearMidpointLeft);
        for (label, cp) in [
            ("P4[0][6] (m)", p4.cps[0][6]),
            ("P4[1][6] (rim+1)", p4.cps[1][6]),
            ("P4[0][5] (cross+1)", p4.cps[0][5]),
            ("P4[0][4] (cross+2)", p4.cps[0][4]),
            ("P4[1][5] (diag)", p4.cps[1][5]),
        ] {
            println!("  {label}: ({:.5}, {:.5}, {:.5})", cp.x, cp.y, cp.z);
        }
    }
    println!(
        "final: fairing={fairing:.1e} guide_fit={:.4e} g1={:.3e} g2={:.3e} min_reg={:.3e}",
        report.guide_fit_residual,
        report.max_g1_residual,
        report.max_g2_residual,
        report.min_regularity,
    );
    // Worst-sample measurement: for each sector corner c, the signed
    // height of the nearby cap surface over ITS OWN face plane
    // (negative = dips below -> wedge interpenetration). Samples are
    // limited to a small lateral disc around c.
    {
        use popkern_core::geometry::gregory::Kp20PieceId;
        let planes: Vec<(Point3, popkern_core::geometry::Vec3)> = pyramid
            .faces
            .iter()
            .filter_map(|(_, face)| {
                let surface = face.surface.as_ref()?;
                if let popkern_core::geometry::SurfaceGeometry::Plane { normal, origin, .. } =
                    &surface.geometry
                {
                    if normal.z.abs() > 0.9 {
                        return None;
                    }
                    let outward_n = if normal.dot(&(*origin - Point3::new(0.0, 0.0, 0.35))) >= 0.0 {
                        *normal
                    } else {
                        -*normal
                    };
                    Some((*origin, outward_n))
                } else {
                    None
                }
            })
            .collect();
        let mut worst = f64::INFINITY;
        let mut worst_sector = 0usize;
        for (s_idx, sector) in report.sectors.iter().enumerate() {
            let piece = sector.piece(Kp20PieceId::NearMidpointRight);
            let c_point = piece.cps[6][6];
            // The corner's own face: nearest plane.
            let (origin, normal) = planes
                .iter()
                .copied()
                .min_by(|a, b| {
                    let da = ((c_point - a.0).dot(&a.1)).abs();
                    let db = ((c_point - b.0).dot(&b.1)).abs();
                    da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
                })
                .unwrap();
            for a in 0..13 {
                for b in 0..13 {
                    let u = 0.5 + 0.5 * a as f64 / 12.0;
                    let v = 0.5 + 0.5 * b as f64 / 12.0;
                    let mut acc = [0.0f64; 3];
                    for i in 0..7 {
                        for j in 0..7 {
                            let w = bern6(u, i) * bern6(v, j);
                            acc[0] += w * piece.cps[i][j].x;
                            acc[1] += w * piece.cps[i][j].y;
                            acc[2] += w * piece.cps[i][j].z;
                        }
                    }
                    let pnt = Point3::new(acc[0], acc[1], acc[2]);
                    if pnt.distance_to(&c_point) > 0.15 {
                        continue;
                    }
                    let h = (pnt - origin).dot(&normal);
                    if h < worst {
                        worst = h;
                        worst_sector = s_idx;
                    }
                }
            }
        }
        println!(
            "cap min signed height over own face plane near c: {worst:.5} (sector {worst_sector})",
        );
        // CP-level dip map for sector 0's corner piece.
        let piece = report.sectors[0].piece(Kp20PieceId::NearMidpointRight);
        let c_point = piece.cps[6][6];
        let (origin, normal) = planes
            .iter()
            .copied()
            .min_by(|a, b| {
                let da = ((c_point - a.0).dot(&a.1)).abs();
                let db = ((c_point - b.0).dot(&b.1)).abs();
                da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap();
        for i in (3..7).rev() {
            let mut row = String::new();
            for j in 3..7 {
                let h = (piece.cps[i][j] - origin).dot(&normal);
                row.push_str(&format!("{h:8.4} "));
            }
            println!("  P3 heights i={i}: {row}");
        }
    }
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

    // Whole object from parts, honestly trimmed to the setback
    // boundary: fillet sweeps sampled from their base ends up to OUR
    // setback arcs, faces trimmed at the contact lines up to the
    // corner points, the base pentagon, and the paper cap.
    let whole_path = out.join("kp20_pyramid_apex_whole.obj");
    write_parts_whole_obj(
        &fillet_surfaces,
        &border,
        &pyramid,
        &bi9_grids,
        boundary.vertex_position,
        &extract_report.trim_corners,
        &whole_path,
    )
    .expect("whole OBJ");
    println!("wrote {}", whole_path.display());
}

/// Whole-object OBJ from parts: fillet strips over their sweep
/// rectangles clipped at the setback arcs (the cap rim), lateral
/// faces as fans over their trimmed boundary polylines, the base
/// pentagon, and the paper cap.
#[allow(clippy::too_many_arguments)]
fn write_parts_whole_obj(
    fillet_surfaces: &[popkern_fillet::fillet_surface::FilletSurfaceResult],
    border: &popkern_core::geometry::gregory::Kp20PaperGuideBorder,
    pyramid: &BRep,
    bi9_grids: &[[[Point3; 10]; 10]],
    vertex: Point3,
    trim_corners: &[Point3],
    path: &Path,
) -> Result<(), std::io::Error> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    use std::io::Write;
    let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        file,
        "# filleted pentagonal pyramid + KP20 paper cap (from parts)"
    )?;
    let mut base_index = 1_usize;
    // Cap sample cloud for mutual trimming of the face wedges (the
    // BRep-level version is face-vs-cap surface trimming; here the
    // wedge triangles covered by cap material are dropped).
    let mut cap_cloud: Vec<Point3> = Vec::new();
    {
        const CS: usize = 17;
        for grid in bi9_grids {
            let patch = RectangularGregoryBiNonic { grid: *grid };
            for i in 0..CS {
                for j in 0..CS {
                    let u = i as f64 / (CS - 1) as f64;
                    let v = j as f64 / (CS - 1) as f64;
                    cap_cloud.push(patch.point_at(u, v));
                }
            }
        }
    }
    let under_cap = |p: Point3, n: popkern_core::geometry::Vec3| -> bool {
        // Covered when some cap sample sits above p along the plane
        // normal within a small lateral radius.
        for q in &cap_cloud {
            let d = *q - p;
            let h = d.dot(&n);
            if h > 1.0e-4 {
                let lateral = d - n * h;
                if lateral.length() < 0.02 {
                    return true;
                }
            }
        }
        false
    };
    // Popkern surface normals point INTO the solid; orient every
    // emitted normal outward from an interior point so viewers shade
    // front faces.
    let interior = Point3::new(0.0, 0.0, 0.35);
    let outward = |p: Point3, n: popkern_core::geometry::Vec3| {
        if n.dot(&(p - interior)) >= 0.0 { n } else { -n }
    };

    let n = fillet_surfaces.len();
    let corners: Vec<Point3> = border
        .sectors
        .iter()
        .map(|sector| sector[0].rows[0][3])
        .collect();
    for (k, fillet) in fillet_surfaces.iter().enumerate() {
        // Clip line = the ACTUAL cap rim on this sweep: project every
        // border rim onto the sweep, keep on-surface samples, and
        // follow them as u(v).
        let mut clip: Vec<(f64, f64)> = Vec::new();
        for sector in &border.sectors {
            for piece in sector {
                let row = &piece.rows[0];
                for step in 0..=12 {
                    let t = step as f64 / 12.0;
                    let mt = 1.0 - t;
                    let pnt = Point3::new(
                        mt * mt * mt * row[0].x
                            + 3.0 * mt * mt * t * row[1].x
                            + 3.0 * mt * t * t * row[2].x
                            + t * t * t * row[3].x,
                        mt * mt * mt * row[0].y
                            + 3.0 * mt * mt * t * row[1].y
                            + 3.0 * mt * t * t * row[2].y
                            + t * t * t * row[3].y,
                        mt * mt * mt * row[0].z
                            + 3.0 * mt * mt * t * row[1].z
                            + 3.0 * mt * t * t * row[2].z
                            + t * t * t * row[3].z,
                    );
                    if let Some((u, v)) = project_point_to_uv(&fillet.surface, &pnt)
                        && fillet.surface.point_at(u, v).distance_to(&pnt) < 2.0e-2
                    {
                        clip.push((v, u));
                    }
                }
            }
        }
        if clip.len() < 4 {
            continue;
        }
        clip.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        let u_arc_at = |v: f64| -> f64 {
            if v <= clip[0].0 {
                let (v0, u0) = clip[0];
                let (v1, u1) = clip[1.min(clip.len() - 1)];
                let span = (v1 - v0).max(1.0e-12);
                return u0 + (u1 - u0) * (v - v0) / span;
            }
            for w in clip.windows(2) {
                if v <= w[1].0 {
                    let span = (w[1].0 - w[0].0).max(1.0e-12);
                    return w[0].1 + (w[1].1 - w[0].1) * (v - w[0].0) / span;
                }
            }
            let (v0, u0) = clip[clip.len() - 2];
            let (v1, u1) = clip[clip.len() - 1];
            let span = (v1 - v0).max(1.0e-12);
            u1 + (u1 - u0) * (v - v1) / span
        };
        // Base end of the sweep = the u extreme farthest from the
        // blend region.
        let probe_v = 0.5;
        let far = corners.iter().fold(Point3::new(0.0, 0.0, 0.0), |acc, c| {
            Point3::new(
                acc.x + c.x / n as f64,
                acc.y + c.y / n as f64,
                acc.z + c.z / n as f64,
            )
        });
        let d0 = fillet.surface.point_at(0.0, probe_v).distance_to(&far);
        let d1 = fillet.surface.point_at(1.0, probe_v).distance_to(&far);
        let u_base = if d0 > d1 { 0.0 } else { 1.0 };
        const AL: usize = 33;
        const AC: usize = 25;
        let mut vote = 0.0_f64;
        for j in 0..AC {
            let v = j as f64 / (AC - 1) as f64;
            let u = 0.5 * (u_arc_at(v) + u_base);
            let p = fillet.surface.point_at(u, v);
            vote += fillet.surface.normal_at(u, v).dot(&(p - interior));
        }
        let flip = vote < 0.0;
        let mut grid_points = vec![Point3::new(0.0, 0.0, 0.0); AL * AC];
        let mut chosen_normal = popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0);
        for i in 0..AL {
            for j in 0..AC {
                let v = j as f64 / (AC - 1) as f64;
                let u_arc = u_arc_at(v);
                let u = u_arc + (u_base - u_arc) * (i as f64 / (AL - 1) as f64);
                let p = fillet.surface.point_at(u, v);
                let raw = fillet.surface.normal_at(u, v);
                let nrm = if flip { -raw } else { raw };
                grid_points[i * AC + j] = p;
                if i == AL / 2 && j == AC / 2 {
                    chosen_normal = nrm;
                }
                writeln!(file, "v {} {} {}", p.x, p.y, p.z)?;
                writeln!(file, "vn {} {} {}", nrm.x, nrm.y, nrm.z)?;
            }
        }
        let ti = grid_points[(AL / 2 + 1) * AC + AC / 2] - grid_points[(AL / 2) * AC + AC / 2];
        let tj = grid_points[(AL / 2) * AC + AC / 2 + 1] - grid_points[(AL / 2) * AC + AC / 2];
        let wind_flip = ti.cross(&tj).dot(&chosen_normal) < 0.0;
        writeln!(file, "g fillet_{k}")?;
        base_index = emit_grid_faces_wound(&mut file, base_index, AL, AC, wind_flip)?;
    }

    for (_face_handle, face) in pyramid.faces.iter() {
        let Some(surface) = face.surface.as_ref() else {
            continue;
        };
        let popkern_core::geometry::SurfaceGeometry::Plane { normal, origin, .. } =
            &surface.geometry
        else {
            continue;
        };
        if normal.z.abs() > 0.9 {
            continue;
        }
        let plane_distance = |p: Point3| ((p - *origin).dot(normal)).abs();
        let mut on_plane: Vec<&popkern_core::geometry::NurbsCurve> = Vec::new();
        for fillet in fillet_surfaces {
            for curve in [&fillet.trim_curve_1, &fillet.trim_curve_2] {
                let (t0, t1) = curve.parameter_range();
                let mid = curve.point_at(0.5 * (t0 + t1));
                if plane_distance(mid) < 1.0e-6 {
                    on_plane.push(curve);
                }
            }
        }
        if on_plane.len() != 2 {
            continue;
        }
        // The face is bounded by its contact lines up to their TRUE
        // intersection, exported by the extractor (independent of the
        // cap's pulled-back corner); the mutual trim against the cap
        // handles the overlap.
        let corner = trim_corners
            .iter()
            .copied()
            .min_by(|a, b| {
                plane_distance(*a)
                    .partial_cmp(&plane_distance(*b))
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(vertex);
        let base_points: Vec<Point3> = pyramid
            .vertices
            .iter()
            .map(|(_, v)| v.position)
            .filter(|p| p.z < 1.0e-6 && plane_distance(*p) < 1.0e-6)
            .collect();
        if base_points.len() != 2 {
            continue;
        }
        let d00 = trim_curve_distance(on_plane[0], base_points[0]);
        let d01 = trim_curve_distance(on_plane[0], base_points[1]);
        let (left_base, right_base) = if d00 < d01 {
            (base_points[0], base_points[1])
        } else {
            (base_points[1], base_points[0])
        };
        // Fan the wedge from the corner point itself, with each trim
        // sampled strictly from just below the corner down to its
        // base end — the previous centroid fan let closest-param
        // imprecision run the polyline past the trim crossing,
        // spraying inverted slivers over the cap.
        let refined_closest = |curve: &popkern_core::geometry::NurbsCurve, target: Point3| -> f64 {
            let (t0, t1) = curve.parameter_range();
            let mut best = (f64::INFINITY, t0);
            for i in 0..=512 {
                let t = t0 + (t1 - t0) * i as f64 / 512.0;
                let d = curve.point_at(t).distance_to(&target);
                if d < best.0 {
                    best = (d, t);
                }
            }
            best.1
        };
        let side_samples =
            |curve: &popkern_core::geometry::NurbsCurve, base_point: Point3| -> Vec<Point3> {
                let t_c = refined_closest(curve, corner);
                let t_b = refined_closest(curve, base_point);
                // Strictly inside (skip the corner itself; the fan apex
                // supplies it exactly).
                (1..=16)
                    .map(|i| {
                        let t = t_c + (t_b - t_c) * i as f64 / 16.0;
                        curve.point_at(t)
                    })
                    .collect()
            };
        // Boundary from just-below-c on the left trim, down to the
        // base edge, back up the right trim.
        let (left_trim, right_trim) = {
            let a_near_left = trim_curve_distance(on_plane[0], left_base);
            let b_near_left = trim_curve_distance(on_plane[1], left_base);
            if a_near_left < b_near_left {
                (on_plane[0], on_plane[1])
            } else {
                (on_plane[1], on_plane[0])
            }
        };
        let mut polyline: Vec<Point3> = side_samples(left_trim, left_base);
        polyline.push(right_base);
        let mut up = side_samples(right_trim, right_base);
        up.reverse();
        polyline.extend(up);
        let face_normal = outward(corner, *normal);
        writeln!(file, "v {} {} {}", corner.x, corner.y, corner.z)?;
        writeln!(
            file,
            "vn {} {} {}",
            face_normal.x, face_normal.y, face_normal.z
        )?;
        let apex_idx = base_index;
        for p in &polyline {
            writeln!(file, "v {} {} {}", p.x, p.y, p.z)?;
            writeln!(
                file,
                "vn {} {} {}",
                face_normal.x, face_normal.y, face_normal.z
            )?;
        }
        writeln!(file, "g face")?;
        let count = polyline.len();
        let e0 = polyline[0] - corner;
        let e1 = polyline[1] - corner;
        let wind_flip = e0.cross(&e1).dot(&face_normal) < 0.0;
        // Mutual trim: drop fan triangles whose off-apex vertices are
        // covered by cap material (the wedge tip beyond the cap
        // intersection is interior to the blend solid).
        let covered: Vec<bool> = polyline
            .iter()
            .map(|p| under_cap(*p, face_normal))
            .collect();
        let apex_covered = under_cap(corner, face_normal);
        for i in 0..count - 1 {
            if apex_covered && (covered[i] || covered[i + 1]) {
                continue;
            }
            if covered[i] && covered[i + 1] {
                continue;
            }
            let a = apex_idx;
            let b = apex_idx + 1 + i;
            let c = apex_idx + 1 + i + 1;
            if wind_flip {
                writeln!(file, "f {a}//{a} {c}//{c} {b}//{b}")?;
            } else {
                writeln!(file, "f {a}//{a} {b}//{b} {c}//{c}")?;
            }
        }
        base_index = apex_idx + 1 + count;
    }

    let base_ring: Vec<Point3> = pyramid
        .vertices
        .iter()
        .map(|(_, v)| v.position)
        .filter(|p| p.z < 1.0e-6)
        .collect();
    if base_ring.len() >= 3 {
        let mut ring = base_ring.clone();
        ring.sort_by(|a, b| {
            a.y.atan2(a.x)
                .partial_cmp(&b.y.atan2(b.x))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        writeln!(file, "v 0 0 0")?;
        writeln!(file, "vn 0 0 -1")?;
        let center_idx = base_index;
        for p in &ring {
            writeln!(file, "v {} {} {}", p.x, p.y, p.z)?;
            writeln!(file, "vn 0 0 -1")?;
        }
        writeln!(file, "g base")?;
        let count = ring.len();
        for i in 0..count {
            let a = center_idx;
            let b = center_idx + 1 + i;
            let c = center_idx + 1 + ((i + 1) % count);
            writeln!(file, "f {a}//{a} {c}//{c} {b}//{b}")?;
        }
        base_index = center_idx + 1 + count;
    }

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
                    outward(
                        p,
                        popkern_core::geometry::Vec3::new(raw.x / len, raw.y / len, raw.z / len),
                    )
                } else {
                    popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0)
                };
                writeln!(file, "v {} {} {}", p.x, p.y, p.z)?;
                writeln!(file, "vn {} {} {}", nrm.x, nrm.y, nrm.z)?;
            }
        }
        writeln!(file, "g paper_cap_{index}")?;
        // The emitted cap normals are outward-flipped du×dv; winding
        // follows the same flip decision at the patch center.
        let mid = 0.5;
        let patch = RectangularGregoryBiNonic { grid: *grid };
        let du = patch.partial_u(mid, mid);
        let dv = patch.partial_v(mid, mid);
        let raw = du.cross(&dv);
        let p_mid = patch.point_at(mid, mid);
        let wind_flip = raw.dot(&(p_mid - interior)) < 0.0;
        base_index = emit_grid_faces_wound(&mut file, base_index, SAMPLES, SAMPLES, wind_flip)?;
    }
    Ok(())
}

fn trim_curve_distance(curve: &popkern_core::geometry::NurbsCurve, target: Point3) -> f64 {
    use popkern_core::geometry::curve::CurveEval;
    let (t0, t1) = curve.parameter_range();
    let mut best = f64::INFINITY;
    for i in 0..=32 {
        let t = t0 + (t1 - t0) * i as f64 / 32.0;
        best = best.min(curve.point_at(t).distance_to(&target));
    }
    best
}

fn bern6(t: f64, i: usize) -> f64 {
    let c = [1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0];
    c[i] * t.powi(i as i32) * (1.0 - t).powi(6 - i as i32)
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
    emit_grid_faces_wound(file, base, rows, cols, false)
}

/// Grid triangulation with optional winding flip so the front face
/// matches the emitted normals (viewers backface-shade otherwise).
fn emit_grid_faces_wound<W: std::io::Write>(
    file: &mut W,
    base: usize,
    rows: usize,
    cols: usize,
    flip: bool,
) -> Result<usize, std::io::Error> {
    for i in 0..(rows - 1) {
        for j in 0..(cols - 1) {
            let a = base + i * cols + j;
            let b = a + cols;
            let c = b + 1;
            let d = a + 1;
            if flip {
                writeln!(file, "f {a}//{a} {c}//{c} {b}//{b}")?;
                writeln!(file, "f {a}//{a} {d}//{d} {c}//{c}")?;
            } else {
                writeln!(file, "f {a}//{a} {b}//{b} {c}//{c}")?;
                writeln!(file, "f {a}//{a} {c}//{c} {d}//{d}")?;
            }
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
