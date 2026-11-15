//! PK-1189 (GitHub issue #1): a boss near a plate's edge, its root and
//! the plate's top and vertical edges rounded together (G1, 0.08). The
//! boss root fillet's foot circle (r 0.43) runs past the plate edge
//! round's rail (x 1.12): it used to serve with the plate top's foot
//! loop crossing its outer loop, validate-clean and 52 open mesh edges.
//! The root fillet now rolls over the round: its ball keeps the boss
//! wall and rolls onto the round where the plate top ends
//! (`construct::strip_overflow`), with G1 arcs or G2 sections.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// The issue's script, the boss moved to `(bx, by)`, the model turned
/// `turn` degrees about z, rounds of `radius` with `continuity`.
fn script(bx: f64, by: f64, turn: f64, radius: f64, continuity: &str) -> String {
    format!(
        r#"let m = model();
let plate = m.profile(xy()).rect(2.4, 1.6).sketch("plate_sketch").extrude(blind(0.4));
let boss = m.tool(|m| m.profile(xy_at(0.2)).circle(0.35).sketch("boss_sketch")
             .extrude(blind(0.8))).translate({bx}, {by}, 0);
plate = plate.union(boss);
plate = plate.rotate(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), {turn});
plate.fillet(
    edges().between("cylinder", "plane").at_height(0.4)
        .or(edges().straight().at_height(0.4))
        .or(edges().straight().direction("z")),
    {radius}, #{{ continuity: "{continuity}", label: "rounds" }});
m.assert_watertight();
m"#
    )
}

fn serve(source: &str) -> Result<BRep, String> {
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host
        .run(source, "boss_near_plate_edge")
        .map_err(|e| e.to_string())?;
    let mut model = run.model.ok_or("no model")?;
    Ok(model.to_parts_all().map_err(|e| e.to_string())?[0]
        .brep()
        .clone())
}

/// Faces, volume, area; the body validated (self-intersection included)
/// and meshed watertight at the display chord.
fn measure(b: &BRep) -> (usize, f64, f64) {
    let tol = Tolerance::default();
    let errors = b.validate_all();
    assert!(errors.is_empty(), "validate_all: {errors:?}");
    let mesh = popkern_mesh::tessellate(b, &popkern_mesh::TessellationOptions::display()).unwrap();
    assert_eq!(
        popkern_mesh::welded_open_edges(&mesh),
        0,
        "display mesh open"
    );
    let volume = popkern_measure::body_volume(b, &tol).unwrap();
    let area = b
        .faces
        .keys()
        .map(|f| popkern_measure::face_area(b, f, &tol).unwrap())
        .sum();
    (b.faces.len(), volume, area)
}

#[test]
fn the_issue_script_serves_one_watertight_body() {
    let source = r#"// A boss plate: plate, boss, union, rounds.
let m = model();
let plate = m.profile(xy()).rect(2.4, 1.6).sketch("plate_sketch")
             .extrude(blind(0.4));
let boss = m.tool(|m| m.profile(xy_at(0.2)).circle(0.35).sketch("boss_sketch")
             .extrude(blind(0.8))).translate(0.7, 0.1, 0);
plate = plate.union(boss);
plate.fillet(
    edges().between("cylinder", "plane").at_height(0.4)
        .or(edges().straight().at_height(0.4))
        .or(edges().straight().direction("z")),
    0.08, #{ continuity: "g1", label: "rounds" });
m.assert_watertight();
print(m.feature("rounds").describe());
m
"#;
    let b = serve(source).expect("the issue's script serves");
    let (faces, _, _) = measure(&b);
    // The 21 faces of the boss plate and the overflow.
    assert_eq!(faces, 22);
}

#[test]
fn the_overflow_reads_the_rolled_over_volume() {
    // Below the onset (the foot short of the rail) the boss's place does
    // not move the volume; past it the root fillet rolling over the round
    // reads what OCC builds rounding the plate first and the root after
    // (its ball rolls onto the round there too): 1.757237231.
    let (_, below, _) = measure(&serve(&script(0.69, 0.1, 0.0, 0.08, "g1")).unwrap());
    let (_, apart, _) = measure(&serve(&script(0.6, 0.1, 0.0, 0.08, "g1")).unwrap());
    assert!((below - apart).abs() < 1e-9, "{below} vs {apart}");
    assert!((apart - 1.757239352).abs() < 1e-8, "{apart}");
    let (_, rolled, _) = measure(&serve(&script(0.7, 0.1, 0.0, 0.08, "g1")).unwrap());
    assert!((rolled - 1.757237231).abs() < 5e-9, "{rolled}");
}

#[test]
fn mirrored_placements_read_alike() {
    // The same overflow against the right, left and back edges, and with
    // the model turned (23 degrees too, since PK-1191): one topology, one
    // volume and area.
    for group in [
        vec![
            (0.7, 0.1, 0.0),
            (-0.7, 0.1, 0.0),
            (0.2, 0.3, 0.0),
            (0.7, 0.1, 90.0),
            (0.7, 0.1, 180.0),
            (0.7, 0.1, 23.0),
        ],
        vec![(0.72, 0.1, 0.0), (-0.72, -0.2, 0.0), (0.2, 0.32, 0.0)],
    ] {
        let reads: Vec<(usize, f64, f64)> = group
            .iter()
            .map(|&(x, y, t)| measure(&serve(&script(x, y, t, 0.08, "g1")).unwrap()))
            .collect();
        for (k, r) in reads.iter().enumerate() {
            assert_eq!(r.0, reads[0].0, "{:?} faces", group[k]);
            assert!(
                (r.1 - reads[0].1).abs() < 2e-9,
                "{:?} volume {} vs {}",
                group[k],
                r.1,
                reads[0].1
            );
            assert!(
                (r.2 - reads[0].2).abs() < 2e-8,
                "{:?} area {} vs {}",
                group[k],
                r.2,
                reads[0].2
            );
        }
    }
}

#[test]
fn one_ring_rolls_over_two_rounds() {
    // Near a corner the foot overflows both the x and the y edge round;
    // between them it cuts the plate top's corner off as its own face.
    let a = measure(&serve(&script(0.7, 0.3, 0.0, 0.08, "g1")).unwrap());
    let b = measure(&serve(&script(0.7, -0.3, 0.0, 0.08, "g1")).unwrap());
    assert_eq!(a.0, 25);
    assert_eq!(a.0, b.0);
    assert!(
        (a.1 - b.1).abs() < 2e-9 && (a.2 - b.2).abs() < 2e-8,
        "{a:?} vs {b:?}"
    );
}

#[test]
fn the_overflow_lies_on_the_rolling_balls_envelope() {
    // The ball's centre keeps 0.43 from the boss axis and 0.16 from the
    // round's axis (x 1.12, z 0.32): every point of the overflow face lies
    // 0.08 from that curve. A 17-section skin strayed 6.7e-6 at x 0.8.
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::surface::SurfaceEval;
    for bx in [0.7, 0.75, 0.8] {
        let b = serve(&script(bx, 0.1, 0.0, 0.08, "g1")).unwrap();
        let spine = |th: f64| {
            let (x, y) = (bx + 0.43 * th.cos(), 0.1 + 0.43 * th.sin());
            let dx = x - 1.12;
            let z = if dx > 0.0 {
                0.32 + (0.0256 - dx * dx).sqrt()
            } else {
                0.48
            };
            Point3::new(x, y, z)
        };
        let to_spine = |p: &Point3| {
            let n = 4000;
            let k = (0..=n)
                .min_by(|&a, &c| {
                    let d = |k: usize| spine(-0.8 + 1.6 * k as f64 / n as f64).distance_to(p);
                    d(a).total_cmp(&d(c))
                })
                .unwrap();
            let (mut lo, mut hi) = (
                -0.8 + 1.6 * (k as f64 - 1.0) / n as f64,
                -0.8 + 1.6 * (k as f64 + 1.0) / n as f64,
            );
            for _ in 0..80 {
                let (m1, m2) = (lo + (hi - lo) / 3.0, hi - (hi - lo) / 3.0);
                if spine(m1).distance_to(p) < spine(m2).distance_to(p) {
                    hi = m2
                } else {
                    lo = m1
                }
            }
            spine(0.5 * (lo + hi)).distance_to(p)
        };
        let overflow: Vec<_> = b
            .faces
            .values()
            .filter_map(|f| f.surface.as_ref())
            .filter(|s| matches!(s.geometry, SurfaceGeometry::Nurbs(_)))
            .collect();
        assert_eq!(overflow.len(), 1, "one overflow face at x {bx}");
        let SurfaceGeometry::Nurbs(n) = &overflow[0].geometry else {
            unreachable!()
        };
        let ((u0, u1), (v0, v1)) = n.parameter_range();
        for i in 0..=24 {
            for j in 0..=12 {
                let p = overflow[0].geometry.point_at(
                    u0 + (u1 - u0) * i as f64 / 24.0,
                    v0 + (v1 - v0) * j as f64 / 12.0,
                );
                let off = (to_spine(&p) - 0.08).abs();
                assert!(
                    off < 1e-7,
                    "x {bx}: overflow point {p:?} {off:.3e} off the envelope"
                );
            }
        }
    }
}

#[test]
fn other_radii_and_offsets_roll_over() {
    for (bx, radius) in [(0.66, 0.12), (0.78, 0.08)] {
        let (faces, _, _) = measure(&serve(&script(bx, 0.1, 0.0, radius, "g1")).unwrap());
        assert_eq!(faces, 22, "x {bx} r {radius}");
    }
}

#[test]
fn a_g2_overflow_rolls_over_alike_everywhere() {
    // With G2 rounds the root strip is swept G2 sections; its overflow is
    // G2 sections too, read against the plate top at its ends (where it
    // meets the strip) and against the round beyond, the round's
    // curvature ramped in over a fifth of it (PK-1192). Mirrored and
    // turned placements read one body.
    let reads: Vec<(usize, f64, f64)> = [(0.7, 0.1, 0.0), (-0.7, 0.1, 0.0), (0.7, 0.1, 23.0)]
        .iter()
        .map(|&(x, y, t)| measure(&serve(&script(x, y, t, 0.08, "g2")).unwrap()))
        .collect();
    for r in &reads {
        assert_eq!(r.0, reads[0].0, "faces {reads:?}");
        assert!((r.1 - reads[0].1).abs() < 2e-9, "volume {reads:?}");
        assert!((r.2 - reads[0].2).abs() < 2e-8, "area {reads:?}");
    }
}
