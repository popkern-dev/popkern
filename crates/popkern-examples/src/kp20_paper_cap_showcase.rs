//! KP20 paper-faithful G² cap showcase.
//!
//! Builds the §4.4 paper cap (`assemble_kp20_paper_cap`) on lifted
//! Catmull-Clark ring borders — a bowl, a cone, and a volcano rim —
//! and exports each as STEP alongside the printed acceptance report
//! (raw eq (8)/(9) G¹/G² residuals, C⁰/border/apex gates).
//!
//! Input validity matters twice here. Chord-type borders with rim
//! corners at the spokes are not Hermite data of any C² ring and are
//! declined outright. More subtly, splitting a SMOOTH cubic at its
//! parameter midpoint manufactures two border pieces that are not
//! corner-compatible at the sector corner `c` (the two halves' depth
//! rows do not form one consistent corner jet), and the eq (14)
//! corner bands then average conflicting data — measured as a rim
//! tangent kink of ~110-150° at every `c`. The paper's native inputs
//! are two genuine surface pieces per sector meeting at `c`
//! (production fillet boundaries arrive as 2n segments); the lifted
//! characteristic-net CC ring used here is the reference-grade
//! corner-compatible generator.
//!
//! ```bash
//! cargo run --bin kp20_paper_cap_showcase
//! f3d --output examples/output/kp20_paper_cap_saddle6.png \
//!     --resolution 1024,768 --edges --up +Z \
//!     examples/output/kp20_paper_cap_saddle6.step
//! ```

use std::path::Path;
use std::sync::Arc;

use popkern_core::constructors::BRepBuilder;
use popkern_core::geometry::Point3;
use popkern_core::geometry::SurfaceGeometry;
use popkern_core::geometry::gregory::{
    GregoryPatch, Kp20PaperGuideBorder, Kp20PaperGuideBorderPiece, Kp20PieceId,
    RectangularGregoryBiNonic, assemble_kp20_paper_cap, cc_lifted_ring_border_pieces,
    elevate_bi6_grid_to_bi9,
};
use popkern_core::geometry::surface::SurfaceOnFace;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_io::{
    OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind, StepDocument,
    StepExportOptions, StepProductNode, UnitContext, export_step_file,
};

fn main() {
    // The paper assembly's dense solves use deep debug-profile frames;
    // run on a roomy stack like the geometry test harnesses do.
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

    for fixture in [Fixture::bowl5(), Fixture::cone6(), Fixture::volcano8()] {
        let report =
            assemble_kp20_paper_cap(&fixture.border, fixture.star).expect("paper cap assembly");
        println!(
            "{name}: n={n} g1={g1:.3e} g2={g2:.3e} c0={c0:.3e} border={b:.3e} lemma4={l4:.3e} \
             guide_fit={gf:.3e} strip_move={sm:.3e} min_regularity={mr:.3e}",
            name = fixture.name,
            n = report.n,
            g1 = report.max_g1_residual,
            g2 = report.max_g2_residual,
            c0 = report.max_c0_residual,
            b = report.max_border_residual,
            l4 = report.lemma4_delta,
            gf = report.guide_fit_residual,
            sm = report.max_strip_movement,
            mr = report.min_regularity,
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
        let step_path = out.join(format!("kp20_paper_cap_{}.step", fixture.name));
        export_step_file(
            &step_document(&brep, fixture.name),
            &step_path,
            &StepExportOptions::default(),
        )
        .expect("STEP export");
        println!("  wrote {}", step_path.display());

        // Untrimmed direct tessellation of the patch grids for smooth
        // rendering: viewers re-mesh STEP coarsely, and BRep-based
        // tessellation clips against the loop's chord edges, notching
        // the rim where the true boundary bulges past the chord.
        let obj_path = out.join(format!("kp20_paper_cap_{}.obj", fixture.name));
        write_patch_grid_obj(&bi9_grids, &obj_path).expect("OBJ export");
        println!("  wrote {}", obj_path.display());
        rim_probe(&report.sectors);
    }
}

struct Fixture {
    name: &'static str,
    border: Kp20PaperGuideBorder,
    star: Point3,
}

impl Fixture {
    fn lifted(name: &'static str, n: usize, star_z: f64, lift: impl Fn(f64, f64) -> f64) -> Self {
        let pieces = cc_lifted_ring_border_pieces(n, lift).expect("lifted CC ring");
        let border = Kp20PaperGuideBorder {
            sectors: pieces
                .into_iter()
                .map(|[near_c, near_m]| {
                    [
                        Kp20PaperGuideBorderPiece::from_bi3(near_c),
                        Kp20PaperGuideBorderPiece::from_bi3(near_m),
                    ]
                })
                .collect(),
        };
        Fixture {
            name,
            border,
            star: Point3::new(0.0, 0.0, star_z),
        }
    }

    /// Paraboloid bowl.
    fn bowl5() -> Self {
        Self::lifted("bowl5", 5, 0.0, |x, y| 0.45 * (x * x + y * y))
    }

    /// Near-conical rise toward the apex.
    fn cone6() -> Self {
        Self::lifted("cone6", 6, 0.55, |x, y| {
            0.55 * (1.0 - (x * x + y * y).sqrt()).max(0.0)
        })
    }

    /// Volcano rim: high ring collapsing toward a low apex.
    fn volcano8() -> Self {
        Self::lifted("volcano8", 8, 0.05, |x, y| {
            let r2 = x * x + y * y;
            0.9 * r2 * (2.0 - r2).max(0.0) * 0.5
        })
    }
}

/// Diagnostic: tangent turning angle of the rim curve at each
/// piece-corner along the cap boundary.
fn rim_probe(sectors: &[popkern_core::geometry::gregory::Kp20MacroPatch]) {
    use popkern_core::geometry::gregory::Kp20PieceId;
    let mut rim_rows: Vec<[Point3; 7]> = Vec::new();
    for sector in sectors {
        let p4 = sector.piece(Kp20PieceId::NearMidpointLeft);
        let p3 = sector.piece(Kp20PieceId::NearMidpointRight);
        let p2 = sector.piece(Kp20PieceId::NearOriginRight);
        rim_rows.push(core::array::from_fn(|i| p4.cps[i][6]));
        rim_rows.push(core::array::from_fn(|i| p3.cps[i][6]));
        rim_rows.push(core::array::from_fn(|j| p3.cps[6][j]));
        rim_rows.push(core::array::from_fn(|j| p2.cps[6][j]));
    }
    let mut worst_gap = 0.0_f64;
    let mut worst_angle = 0.0_f64;
    for k in 0..rim_rows.len() {
        let a = rim_rows[k];
        let b = rim_rows[(k + 1) % rim_rows.len()];
        let d = |p: Point3, q: Point3| {
            ((p.x - q.x).powi(2) + (p.y - q.y).powi(2) + (p.z - q.z).powi(2)).sqrt()
        };
        let candidates = [d(a[6], b[0]), d(a[6], b[6]), d(a[0], b[0]), d(a[0], b[6])];
        let gap = candidates.iter().cloned().fold(f64::INFINITY, f64::min);
        worst_gap = worst_gap.max(gap);
        let best = candidates
            .iter()
            .enumerate()
            .min_by(|x, y| x.1.partial_cmp(y.1).unwrap())
            .map(|(i, _)| i)
            .unwrap_or(0);
        let (ta, tb) = match best {
            0 => (sub(a[6], a[5]), sub(b[1], b[0])),
            1 => (sub(a[6], a[5]), sub(b[5], b[6])),
            2 => (sub(a[0], a[1]), sub(b[1], b[0])),
            _ => (sub(a[0], a[1]), sub(b[5], b[6])),
        };
        let na = norm(ta);
        let nb = norm(tb);
        if na > 1.0e-14 && nb > 1.0e-14 {
            let cosv = (ta.0 * tb.0 + ta.1 * tb.1 + ta.2 * tb.2) / (na * nb);
            let angle = cosv.clamp(-1.0, 1.0).acos().to_degrees();
            if angle > 5.0 {
                println!(
                    "  rim corner {k} (row {} -> {}): turn={angle:.1} deg legs |ta|={na:.3e} |tb|={nb:.3e} best={best}",
                    k % 4,
                    (k + 1) % 4,
                );
            }
            if angle > worst_angle {
                worst_angle = angle;
            }
        } else {
            println!("  rim corner {k}: DEGENERATE tangent (|ta|={na:.3e} |tb|={nb:.3e})");
        }
    }
    println!(
        "  rim probe: worst corner gap={worst_gap:.3e}, worst tangent turn={worst_angle:.2} deg"
    );
}

fn sub(a: Point3, b: Point3) -> (f64, f64, f64) {
    (a.x - b.x, a.y - b.y, a.z - b.z)
}

fn norm(v: (f64, f64, f64)) -> f64 {
    (v.0 * v.0 + v.1 * v.1 + v.2 * v.2).sqrt()
}

/// Sample each bi-9 patch on a dense uniform grid and write one OBJ
/// with analytic per-vertex normals (no trims — the cap surface
/// exactly as constructed).
fn write_patch_grid_obj(
    bi9_grids: &[[[Point3; 10]; 10]],
    path: &Path,
) -> Result<(), std::io::Error> {
    use std::io::Write;
    const SAMPLES: usize = 33;
    let mut file = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(file, "# KP20 paper cap — direct patch tessellation")?;
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
                let n = if len > 1.0e-12 {
                    popkern_core::geometry::Vec3::new(raw.x / len, raw.y / len, raw.z / len)
                } else {
                    popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0)
                };
                writeln!(file, "v {} {} {}", p.x, p.y, p.z)?;
                writeln!(file, "vn {} {} {}", n.x, n.y, n.z)?;
            }
        }
        writeln!(file, "g patch_{index}")?;
        for i in 0..(SAMPLES - 1) {
            for j in 0..(SAMPLES - 1) {
                let a = base + i * SAMPLES + j;
                let b = a + SAMPLES;
                let c = b + 1;
                let d = a + 1;
                writeln!(file, "f {a}//{a} {b}//{b} {c}//{c}")?;
                writeln!(file, "f {a}//{a} {c}//{c} {d}//{d}")?;
            }
        }
        base += SAMPLES * SAMPLES;
    }
    Ok(())
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

fn step_document(brep: &BRep, label: &str) -> StepDocument {
    StepDocument {
        products: vec![StepProductNode {
            name: format!("Kp20PaperCap_{label}"),
            id: format!("KP20-PAPER-CAP-{}", label.to_uppercase()),
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
