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
            // ADJACENT planes = the two whose contact-line distance
            // is smallest; measure the mid point against those only.
            let c0 = fillet.surface.point_at(0.5, 0.0);
            let c1 = fillet.surface.point_at(0.5, 1.0);
            let mut by_contact: Vec<(f64, f64)> = planes
                .iter()
                .map(|(o, n)| {
                    let d = ((c0 - *o).dot(n)).abs().min(((c1 - *o).dot(n)).abs());
                    (d, (mid - *o).dot(n))
                })
                .collect();
            by_contact.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            let heights: Vec<f64> = by_contact.iter().take(2).map(|x| x.1).collect();
            // Contact side: each contact's height vs the OTHER
            // adjacent plane (negative = on-material, ball inside).
            let mut c0_planes: Vec<(f64, f64)> = planes
                .iter()
                .map(|(o, n)| (((c0 - *o).dot(n)).abs(), (c0 - *o).dot(n)))
                .collect();
            c0_planes.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
            println!(
                "  fillet[{k}] contact0: on-plane residual {:.1e}, height vs other plane {:+.4}",
                c0_planes[0].1, c0_planes[1].1,
            );
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
    // Curvature distribution along one mid-edge cross-section (v-dir).
    {
        use popkern_core::geometry::surface::SurfaceEval;
        let f0 = &fillet_surfaces[0].surface;
        let mut line = String::new();
        for step in 0..=10 {
            let v = 0.02 + 0.96 * step as f64 / 10.0;
            let e = 1.0e-4;
            let p_m = f0.point_at(0.5, v - e);
            let p_0 = f0.point_at(0.5, v);
            let p_p = f0.point_at(0.5, v + e);
            let d1 = (p_p - p_m) * (0.5 / e);
            let d2 = (p_p - p_0 - (p_0 - p_m)) * (1.0 / (e * e));
            let speed2 = d1.dot(&d1);
            let kappa = d1.cross(&d2).length() / speed2.powf(1.5).max(1.0e-30);
            line.push_str(&format!("{kappa:6.2} "));
        }
        println!(
            "  cross-section curvature v=0..1: {line} (circular arc would be constant {:.2})",
            1.0 / 0.16
        );
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
    _vertex: Point3,
    _trim_corners: &[Point3],
    path: &Path,
) -> Result<(), std::io::Error> {
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    use std::io::Write;

    // Watertight-by-construction: every pairwise seam (strip↔cap,
    // strip↔face, strip↔base, face↔base) reuses the SAME boundary
    // point sequence on both sides, so welding coincident vertices
    // yields a closed surface. The remaining self-intersection class
    // (cap grazing the face planes near the setback corners) is a
    // recorded geometry item (PK-090), not a hole.
    let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        file,
        "# filleted pentagonal pyramid + KP20 paper cap (watertight from parts)"
    )?;
    let mut base_index = 1_usize;
    let interior = Point3::new(0.0, 0.0, 0.35);
    let outward = |p: Point3, n: popkern_core::geometry::Vec3| {
        if n.dot(&(p - interior)) >= 0.0 { n } else { -n }
    };

    const AL: usize = 33; // rows: rim → base plane
    const PIECE_SAMPLES: usize = 33; // must equal the cap grid sampling
    const AC: usize = 2 * PIECE_SAMPLES - 1; // columns across both rim pieces

    let eval_rim = |row: &[Point3; 4], t: f64| -> Point3 {
        let mt = 1.0 - t;
        Point3::new(
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
        )
    };

    // ----- strips: one grid per fillet, top row = the ACTUAL rim -----
    // Each fillet owns two border pieces (the two halves of its
    // setback arc). Identify them by projection residual, order by
    // sweep parameter, and use their row-0 cubics verbatim as the
    // strip's top row so cap and strip share rim vertices exactly.
    let mut strip_edge_columns: Vec<[Vec<Point3>; 2]> = Vec::new();
    let mut strip_bottom_rows: Vec<Vec<Point3>> = Vec::new();
    let mut strip_rim_halves: Vec<Vec<Point3>> = Vec::new();
    for (k, fillet) in fillet_surfaces.iter().enumerate() {
        let mut owned: Vec<(f64, f64, [Point3; 4])> = Vec::new();
        for sector in &border.sectors {
            for piece in sector {
                let row = [
                    piece.rows[0][0],
                    piece.rows[0][1],
                    piece.rows[0][2],
                    piece.rows[0][3],
                ];
                let mut residual = 0.0;
                let mut mean_v = 0.0;
                let mut ok = true;
                for step in 0..5 {
                    let t = 0.1 + 0.2 * step as f64;
                    let pnt = eval_rim(&row, t);
                    if let Some((u, v)) = project_point_to_uv(&fillet.surface, &pnt) {
                        residual += fillet.surface.point_at(u, v).distance_to(&pnt);
                        mean_v += v;
                    } else {
                        ok = false;
                    }
                }
                if ok {
                    owned.push((residual / 5.0, mean_v / 5.0, row));
                }
            }
        }
        // The fillet's rim = the two halves of ONE arc: pick the
        // lowest-residual pair that CHAINS (first piece ends exactly
        // where the second starts, at the arc midpoint m). Residual
        // thresholding alone can grab a neighbouring arc's half near
        // the corners.
        // The fillet's own two arc halves project onto its sweep at
        // ~machine residual; every other piece is mostly off-surface.
        // Take the two lowest-residual pieces and orient by chaining
        // (they share the arc midpoint m).
        owned.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        if owned.len() < 2 {
            eprintln!("  [whole-obj] fillet {k}: rim pieces not found — skipped");
            strip_edge_columns.push([Vec::new(), Vec::new()]);
            strip_bottom_rows.push(Vec::new());
            continue;
        }
        owned.truncate(2);
        let (row_first, row_second) = if owned[0].2[3].distance_to(&owned[1].2[0]) < 1.0e-9 {
            (owned[0].2, owned[1].2)
        } else if owned[1].2[3].distance_to(&owned[0].2[0]) < 1.0e-9 {
            (owned[1].2, owned[0].2)
        } else {
            eprintln!("  [whole-obj] fillet {k}: rim halves do not chain — skipped");
            strip_edge_columns.push([Vec::new(), Vec::new()]);
            strip_bottom_rows.push(Vec::new());
            continue;
        };
        let owned: [(f64, [Point3; 4]); 2] = [(0.0, row_first), (0.0, row_second)];

        // Base end of the sweep = u extreme farther from the rim.
        let rim_mid = eval_rim(&owned[0].1, 1.0);
        let d0 = fillet.surface.point_at(0.0, 0.5).distance_to(&rim_mid);
        let d1 = fillet.surface.point_at(1.0, 0.5).distance_to(&rim_mid);
        let u_base = if d0 > d1 { 0.0 } else { 1.0 };

        // Columns: top point from the rim cubics (sampled at the SAME
        // parameters the cap grids use), interior from the sweep,
        // bottom snapped to the base plane z = 0.
        let mut grid = vec![Point3::new(0.0, 0.0, 0.0); AL * AC];
        let mut degenerate = false;
        for j in 0..AC {
            let (piece, t) = if j < PIECE_SAMPLES {
                (0_usize, j as f64 / (PIECE_SAMPLES - 1) as f64)
            } else {
                (
                    1_usize,
                    (j - (PIECE_SAMPLES - 1)) as f64 / (PIECE_SAMPLES - 1) as f64,
                )
            };
            let top = eval_rim(&owned[piece].1, t);
            let Some((u_top, v)) = project_point_to_uv(&fillet.surface, &top) else {
                degenerate = true;
                break;
            };
            // Base-plane crossing along this column (bisection).
            let z_top = fillet.surface.point_at(u_top, v).z;
            let z_end = fillet.surface.point_at(u_base, v).z;
            let u_bottom = if (z_end > 0.0) == (z_top > 0.0) {
                u_base
            } else {
                let (mut lo, mut hi) = (u_top, u_base);
                for _ in 0..48 {
                    let mid = 0.5 * (lo + hi);
                    if (fillet.surface.point_at(mid, v).z > 0.0) == (z_top > 0.0) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                0.5 * (lo + hi)
            };
            for i in 0..AL {
                let idx = i * AC + j;
                if i == 0 {
                    grid[idx] = top;
                } else {
                    let frac = i as f64 / (AL - 1) as f64;
                    let u = u_top + (u_bottom - u_top) * frac;
                    let mut pnt = fillet.surface.point_at(u, v);
                    if i == AL - 1 {
                        pnt = Point3::new(pnt.x, pnt.y, 0.0);
                    }
                    grid[idx] = pnt;
                }
            }
        }
        if degenerate {
            strip_edge_columns.push([Vec::new(), Vec::new()]);
            strip_bottom_rows.push(Vec::new());
            continue;
        }

        // Normal orientation: majority vote at mid-band.
        let mut vote = 0.0_f64;
        for j in (0..AC).step_by(4) {
            let p = grid[(AL / 2) * AC + j];
            if let Some((u, v)) = project_point_to_uv(&fillet.surface, &p) {
                vote += fillet.surface.normal_at(u, v).dot(&(p - interior));
            }
        }
        let flip = vote < 0.0;
        for i in 0..AL {
            for j in 0..AC {
                let p = grid[i * AC + j];
                let nrm = match project_point_to_uv(&fillet.surface, &p) {
                    Some((u, v)) => {
                        let raw = fillet.surface.normal_at(u, v);
                        if flip { -raw } else { raw }
                    }
                    None => outward(p, popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0)),
                };
                writeln!(file, "v {} {} {}", p.x, p.y, p.z)?;
                writeln!(file, "vn {} {} {}", nrm.x, nrm.y, nrm.z)?;
            }
        }
        let center = grid[(AL / 2) * AC + AC / 2];
        let ti = grid[(AL / 2 + 1) * AC + AC / 2] - center;
        let tj = grid[(AL / 2) * AC + AC / 2 + 1] - center;
        let chosen = match project_point_to_uv(&fillet.surface, &center) {
            Some((u, v)) => {
                let raw = fillet.surface.normal_at(u, v);
                if flip { -raw } else { raw }
            }
            None => popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
        };
        let wind_flip = ti.cross(&tj).dot(&chosen) < 0.0;
        writeln!(file, "g fillet_{k}")?;
        base_index = emit_grid_faces_wound(&mut file, base_index, AL, AC, wind_flip)?;

        strip_rim_halves.push(grid[..PIECE_SAMPLES].to_vec());
        strip_rim_halves.push(grid[PIECE_SAMPLES - 1..AC].to_vec());
        let col0: Vec<Point3> = (0..AL).map(|i| grid[i * AC]).collect();
        let col1: Vec<Point3> = (0..AL).map(|i| grid[i * AC + AC - 1]).collect();
        strip_edge_columns.push([col0, col1]);
        strip_bottom_rows.push(grid[(AL - 1) * AC..].to_vec());
    }

    // ----- lateral faces: loop = strip edge ∪ base segment ∪ strip edge -----
    // The two adjacent strips' facing edge columns share their top
    // point (the rim corner c), so the face boundary closes exactly.
    let mut base_face_segments: Vec<Vec<Point3>> = Vec::new();
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
            continue; // base face handled separately
        }
        let plane_distance = |p: Point3| ((p - *origin).dot(normal)).abs();
        // The two closest strip edge columns bound this face.
        let mut candidates: Vec<(f64, usize, usize)> = Vec::new();
        for (k, cols) in strip_edge_columns.iter().enumerate() {
            for (side, col) in cols.iter().enumerate() {
                if col.is_empty() {
                    continue;
                }
                let mid = col[col.len() / 2];
                candidates.push((plane_distance(mid), k, side));
            }
        }
        candidates.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
        if candidates.len() < 2 {
            continue;
        }
        let col_a = &strip_edge_columns[candidates[0].1][candidates[0].2];
        let col_b = &strip_edge_columns[candidates[1].1][candidates[1].2];
        // Shared base segment between the two column bottoms.
        let seg: Vec<Point3> = (0..=8)
            .map(|i| {
                let t = i as f64 / 8.0;
                let a = col_a[AL - 1];
                let b = col_b[AL - 1];
                Point3::new(a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t, 0.0)
            })
            .collect();
        // Boundary loop: down column A, along the base, up column B.
        // Tops coincide (the rim corner c) — drop B's top.
        let mut boundary: Vec<Point3> = col_a.clone();
        boundary.extend(seg.iter().skip(1).take(7).copied());
        let mut up_b: Vec<Point3> = col_b.clone();
        up_b.reverse();
        boundary.extend(up_b.iter().take(AL - 1).copied());
        base_face_segments.push(seg);

        let face_centroid = {
            let mut acc = [0.0_f64; 3];
            for p in &boundary {
                acc[0] += p.x;
                acc[1] += p.y;
                acc[2] += p.z;
            }
            let count = boundary.len() as f64;
            Point3::new(acc[0] / count, acc[1] / count, acc[2] / count)
        };
        let face_normal = outward(face_centroid, *normal);
        writeln!(
            file,
            "v {} {} {}",
            face_centroid.x, face_centroid.y, face_centroid.z
        )?;
        writeln!(
            file,
            "vn {} {} {}",
            face_normal.x, face_normal.y, face_normal.z
        )?;
        let center_idx = base_index;
        for p in &boundary {
            writeln!(file, "v {} {} {}", p.x, p.y, p.z)?;
            writeln!(
                file,
                "vn {} {} {}",
                face_normal.x, face_normal.y, face_normal.z
            )?;
        }
        writeln!(file, "g face")?;
        let count = boundary.len();
        let e0 = boundary[0] - face_centroid;
        let e1 = boundary[1] - face_centroid;
        let wind_flip = e0.cross(&e1).dot(&face_normal) < 0.0;
        for i in 0..count {
            let a = center_idx;
            let b = center_idx + 1 + i;
            let c = center_idx + 1 + ((i + 1) % count);
            if wind_flip {
                writeln!(file, "f {a}//{a} {c}//{c} {b}//{b}")?;
            } else {
                writeln!(file, "f {a}//{a} {b}//{b} {c}//{c}")?;
            }
        }
        base_index = center_idx + 1 + count;
    }

    // ----- base face: fillet end arcs + the SAME face base segments -----
    let mut ring: Vec<Point3> = Vec::new();
    for row in &strip_bottom_rows {
        ring.extend(row.iter().copied());
    }
    for seg in &base_face_segments {
        ring.extend(seg.iter().skip(1).take(7).copied());
    }
    if ring.len() >= 3 {
        ring.sort_by(|a, b| {
            a.y.atan2(a.x)
                .partial_cmp(&b.y.atan2(b.x))
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        ring.dedup_by(|a, b| a.distance_to(b) < 1.0e-12);
        writeln!(file, "v 0 0 0")?;
        writeln!(file, "vn 0 0 -1")?;
        let center_idx = base_index;
        for p in &ring {
            writeln!(file, "v {} {} 0", p.x, p.y)?;
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

    // ----- cap: bi-9 piece grids with SIDE-PAIRED seams -----
    // The assembled cap's residual C⁰ mismatch (the c0 gate: piece
    // seams and border-vs-rim, ~5e-3 average with larger local
    // excursions) is a recorded core item. For a closed render mesh,
    // whole grid SIDES are paired index-to-index: sides matching a
    // strip rim row take the rim verbatim (authoritative), cap-cap
    // seam pairs are averaged pointwise, and shared corners are
    // clustered exactly. Interior samples are untouched, so the edit
    // is bounded by the c0 metric itself.
    const SAMPLES: usize = PIECE_SAMPLES;
    const SEAM_TOL: f64 = 4.0e-2;
    let mut cap_grids_pts: Vec<Vec<Point3>> = Vec::with_capacity(bi9_grids.len());
    for grid in bi9_grids {
        let patch = RectangularGregoryBiNonic { grid: *grid };
        let mut pts = vec![Point3::new(0.0, 0.0, 0.0); SAMPLES * SAMPLES];
        for i in 0..SAMPLES {
            for j in 0..SAMPLES {
                let u = i as f64 / (SAMPLES - 1) as f64;
                let v = j as f64 / (SAMPLES - 1) as f64;
                pts[i * SAMPLES + j] = patch.point_at(u, v);
            }
        }
        cap_grids_pts.push(pts);
    }
    // Grid side accessors: 4 sides, each SAMPLES indices into a grid.
    let side_indices = |side: usize| -> Vec<usize> {
        (0..SAMPLES)
            .map(|t| match side {
                0 => t,                           // i = 0 row
                1 => (SAMPLES - 1) * SAMPLES + t, // i = max row
                2 => t * SAMPLES,                 // j = 0 col
                _ => t * SAMPLES + SAMPLES - 1,   // j = max col
            })
            .collect()
    };
    let side_points = |pts: &[Point3], side: usize| -> Vec<Point3> {
        side_indices(side).iter().map(|&ix| pts[ix]).collect()
    };
    let match_score = |a: &[Point3], b: &[Point3]| -> (f64, bool) {
        // (mean distance, reversed?) — orientation from three probes.
        let probe = |rev: bool| -> f64 {
            [0, SAMPLES / 2, SAMPLES - 1]
                .iter()
                .map(|&t| {
                    let o = if rev { SAMPLES - 1 - t } else { t };
                    a[t].distance_to(&b[o])
                })
                .sum::<f64>()
        };
        let reversed = probe(true) < probe(false);
        let mut acc = 0.0;
        for t in 0..SAMPLES {
            let other = if reversed { b[SAMPLES - 1 - t] } else { b[t] };
            acc += a[t].distance_to(&other);
        }
        (acc / SAMPLES as f64, reversed)
    };
    // Rim sides from the strips (two per strip), authoritative.
    let rim_sides: Vec<Vec<Point3>> = strip_rim_halves.clone();
    // 1. Snap cap sides onto rim sides — globally greedy and UNIQUE
    //    (each rim half claims its single best cap side; first-under-
    //    threshold matching grabbed wrong sides at this cap's scale).
    let n_grids = cap_grids_pts.len();
    let mut rim_matched = vec![[false; 4]; n_grids];
    let mut rim_candidates: Vec<(f64, usize, usize, usize, bool)> = Vec::new();
    for (r, rim) in rim_sides.iter().enumerate() {
        if rim.is_empty() {
            continue;
        }
        for g in 0..n_grids {
            for side in 0..4 {
                let (score, reversed) = match_score(&side_points(&cap_grids_pts[g], side), rim);
                if score < SEAM_TOL {
                    rim_candidates.push((score, r, g, side, reversed));
                }
            }
        }
    }
    rim_candidates.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut rim_used = vec![false; rim_sides.len()];
    for (_, r, g, side, reversed) in rim_candidates {
        if rim_used[r] || rim_matched[g][side] {
            continue;
        }
        let idxs = side_indices(side);
        let rim = &rim_sides[r];
        for (t, &ix) in idxs.iter().enumerate() {
            // match_score compared side (a) against rim (b): a[t] pairs
            // with rim[SAMPLES-1-t] when reversed.
            let src = if reversed {
                rim[SAMPLES - 1 - t]
            } else {
                rim[t]
            };
            cap_grids_pts[g][ix] = src;
        }
        rim_used[r] = true;
        rim_matched[g][side] = true;
    }
    if std::env::var("KP20_SEAM_TRACE").is_ok() {
        eprintln!(
            "  [seam] rims used {}/{}",
            rim_used.iter().filter(|u| **u).count(),
            rim_sides.len()
        );
    }
    // Cap piece seams meet with T-junctions (adjacent pieces share
    // half-length sides), so side pairing cannot close them; the
    // residual sliver gaps (bounded by the c0 gate) are closed by the
    // mesh-finish weld/fill pass in tools/render-apex-showcase.sh.
    // The BRep-level fix — driving assembly c0 to machine zero — is
    // the recorded core item.
    for (index, grid) in bi9_grids.iter().enumerate() {
        let patch = RectangularGregoryBiNonic { grid: *grid };
        let pts = &cap_grids_pts[index];
        for i in 0..SAMPLES {
            for j in 0..SAMPLES {
                let u = i as f64 / (SAMPLES - 1) as f64;
                let v = j as f64 / (SAMPLES - 1) as f64;
                let p = pts[i * SAMPLES + j];
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
        let mid = 0.5;
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
