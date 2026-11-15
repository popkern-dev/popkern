//! PK-1190: a user's 54 mm portafilter dosing funnel (measured from
//! their mesh), built as popkern scripts: a revolved flared body, three
//! raised collar segments with round ends on its ring, the segments' tops
//! rounded and the ring's inner edge chamfered. Three of these declined.
//! - The slot-shaped collar united with the ring declined non-manifold.
//!   Each round end's bottom arc lies along the ring top, and its circle
//!   record crosses both the end cylinder's seam and the plane chart's
//!   circle start: the clip took overlap ends at the loop's own
//!   parameter end for crossings, then kept the half circle off the face
//!   (`clip_overlap_ends`, `clip_wrap_complement`).
//! - The inner chamfer on the ring's slightly conical bore (a 1.4 degree
//!   taper at radius 28) offset half its stations up off the face:
//!   aiming at the boundary centroid reads a shallow cone wrong
//!   (`chamfer_offset_rows`, `InteriorHint::ChartSigned`).
//! - The collar tops' rounds, which follow once the collar unites.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

const PI: f64 = std::f64::consts::PI;
const RI: f64 = 31.45;
const RO: f64 = 33.14;
const H: f64 = 4.6;

/// The funnel's revolved body and, per `collar`, its three segments
/// (`"slot"` profiles with round ends, `"caps"` flat-ended sectors with
/// end cylinders, `"round"` flat-ended sectors whose vertical end edges
/// are filleted to a full round), then `finish` (`"tops"`, `"chamfer"`).
fn funnel(collar: &str, finish: &[&str]) -> String {
    let mut s = String::from(
        r#"let m = model();
let pi = 3.141592653589793;
let body = m.profile(xz())
    .path(34.5, -11.6).line_to(36.3, -11.6).line_to(36.5, -11.1).line_to(36.43, -10.7)
    .line_to(32.1, -1.4).line_to(31.7, 0.0).line_to(31.65, 0.9).line_to(31.92, 2.2)
    .line_to(32.73, 3.9).line_to(33.05, 4.9).line_to(33.14, 5.9).line_to(33.14, 7.0)
    .line_to(27.95, 7.0).line_to(27.98, 5.8).line_to(28.2, 3.6).line_to(28.6, 1.8)
    .line_to(29.05, 0.3).line_to(29.4, -0.6).line_to(30.15, -2.25).close()
    .sketch("funnel_section")
    .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
let ri = 31.45;
let ro = 33.14;
let rc = 0.5 * (ri + ro);
let w = 0.5 * (ro - ri);
fn pt(r, a) { [r * cos(a), r * sin(a)] }
for k in 0..3 {
    let c = (57.0 + 120.0 * k) * pi / 180.0;
    let a0 = c - 30.0 * pi / 180.0;
    let a1 = c + 30.0 * pi / 180.0;
    let o0 = pt(ro, a0); let o1 = pt(ro, a1); let om = pt(ro, c);
    let i0 = pt(ri, a0); let i1 = pt(ri, a1); let im = pt(ri, c);
"#,
    );
    match collar {
        "slot" => s.push_str(
            r#"    let e1 = [rc * cos(a1) - w * sin(a1), rc * sin(a1) + w * cos(a1)];
    let e0 = [rc * cos(a0) + w * sin(a0), rc * sin(a0) - w * cos(a0)];
    body = body.union(m.profile(xy_at(7.0)).path(o0[0], o0[1]).arc_to(om, o1).arc_to(e1, i1)
        .arc_to(im, i0).arc_to(e0, o0).close().sketch("collar_" + k).extrude(blind(4.6)));
}
"#,
        ),
        "caps" => s.push_str(
            r#"    body = body.union(m.profile(xy_at(7.0)).path(o0[0], o0[1]).arc_to(om, o1).line_to(i1[0], i1[1])
        .arc_to(im, i0).close().sketch("collar_" + k).extrude(blind(4.6)));
    for a in [a0, a1] {
        body = body.union(m.tool(|m| m.profile(xy_at(7.0)).circle(w).sketch("cap").extrude(blind(4.6)))
            .translate(rc * cos(a), rc * sin(a), 0));
    }
}
"#,
        ),
        _ => s.push_str(
            r#"    body = body.union(m.profile(xy_at(7.0)).path(o0[0], o0[1]).arc_to(om, o1).line_to(i1[0], i1[1])
        .arc_to(im, i0).close().sketch("collar_" + k).extrude(blind(4.6)));
}
body.fillet(edges().straight().direction("z"), w, #{ continuity: "g1", label: "collar_ends" });
"#,
        ),
    }
    if finish.contains(&"tops") {
        s.push_str(
            "body.fillet(edges().at_height(11.6), 0.5, #{ continuity: \"g1\", label: \"collar_tops\" });\n",
        );
    }
    if finish.contains(&"chamfer") {
        s.push_str(
            "body.chamfer(edges().between(\"cone\", \"plane\").at_height(7.0), 0.45, #{ label: \"inner_chamfer\" });\n",
        );
    }
    s.push_str("m.assert_watertight();\nm\n");
    s
}

fn serve(source: &str) -> BRep {
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host
        .run(source, "portafilter_funnel")
        .unwrap_or_else(|e| panic!("{e}"));
    run.model.expect("model").to_parts_all().expect("parts")[0]
        .brep()
        .clone()
}

/// The body validated (self-intersection included) and watertight at
/// the display chord; its volume.
fn volume(b: &BRep) -> f64 {
    let errors = b.validate_all();
    assert!(errors.is_empty(), "validate_all: {errors:?}");
    let mesh = popkern_mesh::tessellate(b, &popkern_mesh::TessellationOptions::display()).unwrap();
    assert_eq!(
        popkern_mesh::welded_open_edges(&mesh),
        0,
        "display mesh open"
    );
    popkern_measure::body_volume(b, &Tolerance::default()).unwrap()
}

/// The revolved body alone.
fn body_only() -> BRep {
    let source = funnel("slot", &[]);
    let cut = source.find("let ri = 31.45;").unwrap();
    serve(&format!("{}m\n", &source[..cut]))
}

#[test]
fn slot_collars_unite_with_the_ring() {
    // The collars stand on the ring top, touching it along their bottom
    // faces only: the union's volume is the body's and three slot prisms'
    // (an annular sector of 60 degrees and a disc of half the wall).
    let base = volume(&body_only());
    let united = volume(&serve(&funnel("slot", &[])));
    let w = 0.5 * (RO - RI);
    let slot = (PI / 3.0 / 2.0 * (RO * RO - RI * RI) + PI * w * w) * H;
    let want = base + 3.0 * slot;
    assert!(
        ((united - want) / want).abs() < 1e-7,
        "united {united} want {want} (body {base} + 3 x {slot})"
    );
}

#[test]
fn every_collar_end_rounds_and_tops_serve() {
    // The three ways of making the collar's round ends, each with its
    // tops rounded (0.5) and the inner edge chamfered, serve one valid
    // watertight body. Slots and end caps are one shape. Filleting the
    // flat ends to a full round keeps the round inside the sector, its
    // arcs Δ = asin(w / rc) shorter at each end: per segment it loses
    // Δ (ro² - ri²) H of collar and keeps the quarter-round rounds' cut
    // (area A, centroids 0.2234 r in from the walls) along those arcs.
    let reads: Vec<f64> = ["slot", "caps", "round"]
        .iter()
        .map(|collar| volume(&serve(&funnel(collar, &["tops", "chamfer"]))))
        .collect();
    let (slot, caps, round) = (reads[0], reads[1], reads[2]);
    assert!(((caps - slot) / slot).abs() < 1e-7, "{reads:?}");
    let (w, rc, r) = (0.5 * (RO - RI), 0.5 * (RO + RI), 0.5);
    let delta = (w / rc).asin();
    let a = (1.0 - PI / 4.0) * r * r;
    let c = r * (10.0 - 3.0 * PI) / (3.0 * (4.0 - PI));
    let per_radian = (RO * RO - RI * RI) / 2.0 * H - a * ((RO - c) + (RI + c));
    let want = slot - 3.0 * 2.0 * delta * per_radian;
    assert!(
        ((round - want) / want).abs() < 1e-7,
        "round {round} want {want} ({reads:?})"
    );
}

#[test]
fn the_inner_chamfer_removes_its_exact_ring() {
    // On a ring whose bore tapers 0.03 over 1.2 (1.4 degrees at radius
    // 28), the 0.45 chamfer of the bore's top edge removes the triangle
    // between the edge, 0.45 in along the top and 0.45 down the bore,
    // revolved: Pappus, 2 pi times its centroid radius times its area.
    let ring = |chamfer: bool| {
        let mut s = String::from(
            r#"let m = model();
let body = m.profile(xz())
    .path(27.98, 5.8).line_to(33.0, 5.8).line_to(33.0, 7.0).line_to(27.95, 7.0).close()
    .sketch("ring").revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
"#,
        );
        if chamfer {
            s.push_str("body.chamfer(edges().between(\"cone\", \"plane\").at_height(7.0), 0.45, #{ label: \"inner_chamfer\" });\n");
        }
        s.push_str("m.assert_watertight();\nm\n");
        volume(&serve(&s))
    };
    let removed = ring(false) - ring(true);
    let corner = (27.95, 7.0);
    let along_top = (27.95 + 0.45, 7.0);
    let (dr, dz) = (0.03_f64, -1.2_f64);
    let len = dr.hypot(dz);
    let down_bore = (27.95 + 0.45 * dr / len, 7.0 + 0.45 * dz / len);
    let pts = [corner, along_top, down_bore];
    let area = 0.5
        * ((pts[1].0 - pts[0].0) * (pts[2].1 - pts[0].1)
            - (pts[2].0 - pts[0].0) * (pts[1].1 - pts[0].1))
            .abs();
    let r_bar = (pts[0].0 + pts[1].0 + pts[2].0) / 3.0;
    let want = 2.0 * PI * r_bar * area;
    assert!(
        ((removed - want) / want).abs() < 1e-6,
        "removed {removed} want {want}"
    );
}
