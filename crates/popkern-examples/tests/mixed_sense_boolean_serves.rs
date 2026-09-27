//! PK-478 acceptance: a MIXED material-sense operand pair — a
//! folded-inward body cut by a fresh (folded-outward) tool — serves
//! watertight with the correct volume. Blend rebuilds used to fold
//! bodies inward (a single tip chamfer flipped a whole shank — the
//! closed-rim surgery's contact-circle axis sign, fixed in PK-564),
//! so every blend-then-boolean chain was such a pair; inward-folded
//! operands still arrive from other producers (primitive shells,
//! imports), so the pair is now built explicitly by mirroring the
//! chamfered body's material side. Before the sense-aware entry
//! normalization the merge wrong-twinned seam-laddered wall loops on
//! exactly this class (the PK-475 seam-ladder gap and the falsified
//! repair routes). Validator-clean alone is NOT accepted here: the measured
//! failure modes include validator-clean volume lies, so the mesh
//! census + exact volume oracle are the gate.

use popkern_core::validation::MaterialSense;
use popkern_part::prelude::*;

/// Chamfered puck (mirrored to fold inward) minus a
/// through-slot prism (folded-outward fresh extrusion) that PIERCES
/// the periodic wall on both sides — the wall's fragment loops then
/// mix shared (SSI) edges with operand-own edges, the class the
/// mixed-sense defect wrong-twinned.
///
/// Exact volume: puck pi*r^2*h, minus the chamfer ring (Pappus: a
/// 0.4 x 0.4 right-triangle cross-section, centroid at r - 0.4/3),
/// minus the slot's in-puck volume A*2 where A is the width-1 chord
/// strip area of the r=2 circle. The slot spans z in [0,2], clear of
/// the chamfer band (z >= 3.6), and its caps sit at z=-1 and z=2 so
/// no face is coincident with the puck's caps — the two cuts do not
/// interact and neither coincidence class is triggered.
#[test]
fn chamfered_body_cut_by_fresh_tool_serves_with_exact_volume() {
    let chamfered = Profile::on(Plane::XY)
        .circle(2.0)
        .extrude(4.0)
        .expect("puck")
        .chamfer(0.4, |e| e.circular().at_height(4.0))
        .expect("chamfered puck rim");
    let tool = Profile::on(Plane::XY)
        .rect(5.0, 1.0)
        .extrude(3.0)
        .expect("slot tool")
        .translate(0.0, 0.0, -1.0)
        .expect("slot placement");

    // The pair must actually be mixed-sense or the fixture guards
    // nothing (the sense conventions are measured, not assumed).
    let sense_of = |b: &popkern_core::topology::BRep| {
        b.faces
            .iter()
            .next()
            .and_then(|(fh, _)| b.component_material_sense(fh))
    };
    // The chamfer serves outward (PK-564); mirror it for the pair.
    assert_eq!(
        sense_of(chamfered.brep()),
        Some(MaterialSense::FoldedOutward),
        "the chamfer keeps its body folded outward"
    );
    let mut mirrored = chamfered.brep().clone();
    mirrored.reverse_body_orientation();
    let chamfered = Part::from_brep(mirrored, chamfered.tolerance(), "mirrored chamfered puck")
        .expect("mirrored puck validates");
    assert_eq!(
        sense_of(chamfered.brep()),
        Some(MaterialSense::FoldedInward),
        "mirrored operand folds inward — the fixture premise"
    );
    assert_eq!(
        sense_of(tool.brep()),
        Some(MaterialSense::FoldedOutward),
        "fresh extrusion folds outward — the fixture premise"
    );

    let keyed = chamfered.subtract(&tool).expect("mixed-sense subtract");
    let brep = keyed.brep();
    brep.validate().expect("valid result");

    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .expect("result meshes");
    // Weld vertices by quantized position (per-face meshes share no
    // indices), then census edge use counts and the signed volume.
    let mut key: std::collections::HashMap<(i64, i64, i64), usize> =
        std::collections::HashMap::new();
    let mut welded: Vec<usize> = Vec::with_capacity(mesh.vertices.len());
    for v in &mesh.vertices {
        let p = v.position;
        let q = (
            (p.x * 1.0e7).round() as i64,
            (p.y * 1.0e7).round() as i64,
            (p.z * 1.0e7).round() as i64,
        );
        let n = key.len();
        welded.push(*key.entry(q).or_insert(n));
    }
    let mut boundary: std::collections::HashMap<(usize, usize), i64> =
        std::collections::HashMap::new();
    let mut volume6 = 0.0_f64;
    for t in &mesh.triangles {
        let [a, b, c] = *t;
        let (wa, wb, wc) = (welded[a], welded[b], welded[c]);
        for (i, j) in [(wa, wb), (wb, wc), (wc, wa)] {
            if i == j {
                continue;
            }
            let entry = boundary.entry((i.min(j), i.max(j))).or_default();
            *entry += if i < j { 1 } else { -1 };
        }
        let (pa, pb, pc) = (
            mesh.vertices[a].position,
            mesh.vertices[b].position,
            mesh.vertices[c].position,
        );
        volume6 += pa.x * (pb.y * pc.z - pc.y * pb.z) - pa.y * (pb.x * pc.z - pc.x * pb.z)
            + pa.z * (pb.x * pc.y - pc.x * pb.y);
    }
    let open = boundary.values().filter(|&&v| v != 0).count();
    assert_eq!(open, 0, "{open} open mesh edges on the mixed-sense cut");

    // Visual verification artifact (geometry-affecting change rule).
    let out_dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR"));
    std::fs::create_dir_all(out_dir).expect("output dir");
    let camera = popkern_viz::camera::Camera {
        eye: glam::Vec3::new(6.0, -5.0, 5.0),
        target: glam::Vec3::new(0.0, 0.0, 2.0),
        up: glam::Vec3::new(0.0, 0.0, 1.0),
        fov_y: 45.0_f32.to_radians(),
        aspect: 800.0 / 600.0,
        near: 0.1,
        far: 100.0,
    };
    let png = out_dir.join("pk478_mixed_sense_cut.png");
    popkern_viz::headless::render_to_png(brep, &camera, 800, 600, &png, true).expect("render");
    println!("[render] {}", png.display());

    let volume = (volume6 / 6.0).abs();
    let puck = std::f64::consts::PI * 2.0 * 2.0 * 4.0;
    let chamfer_ring = 2.0 * std::f64::consts::PI * (2.0 - 0.4 / 3.0) * (0.5 * 0.4 * 0.4);
    // Chord-strip area for |y| <= 0.5 in the r=2 disc:
    // A = 2 * [ (y/2)sqrt(4-y^2) + 2 asin(y/2) ] over [-1/2, 1/2].
    let f = |y: f64| 0.5 * y * (4.0 - y * y).sqrt() + 2.0 * (y / 2.0).asin();
    let slot = 2.0 * (f(0.5) - f(-0.5)) * 2.0;
    let expected = puck - chamfer_ring - slot;
    let rel = (volume - expected).abs() / expected;
    assert!(
        rel < 0.01,
        "mixed-sense volume {volume:.4} vs expected {expected:.4} ({:.2}% off)",
        rel * 100.0
    );
}
