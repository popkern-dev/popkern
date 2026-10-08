//! Closed tangent-chain fillets (PK-733): a rounded rectangle's rim —
//! four lines and four arcs joined on one tangent line at every
//! junction — fillets as one chain. The per-edge pipeline read each
//! junction as a corner and walked its anchor `r` along the tangent arc
//! (off the junction section); the chain is now built whole, the strips
//! meeting on the exact section at every junction.
//!
//! Volumes are analytic: the fillet removes (convex) or adds (concave)
//! the section `(1 − π/4) r²` along the rim — straight runs by length,
//! corner arcs by Pappus with the section centroid `r(10 − 3π)/(3(4 − π))`
//! inside the corner radius.

use popkern_rhai::ScriptHost;

fn served_volume(source: &str) -> f64 {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(source, "chain.rhai")
        .unwrap_or_else(|e| panic!("{e}\n{}", e.payload));
    let mut model = outcome.model.expect("model");
    let part = model.to_part().expect("part");
    let brep = part.brep();
    popkern_measure::body_volume(brep, &brep.tolerance).expect("volume")
}

/// Section area and its centroid offset from the sharp corner.
fn section(r: f64) -> (f64, f64) {
    let pi = std::f64::consts::PI;
    (
        (1.0 - pi / 4.0) * r * r,
        r * (10.0 - 3.0 * pi) / (3.0 * (4.0 - pi)),
    )
}

fn rim_change(w: f64, h: f64, corner: f64, r: f64) -> f64 {
    let (a, x) = section(r);
    let straight = 2.0 * (w - 2.0 * corner) + 2.0 * (h - 2.0 * corner);
    a * (straight + 2.0 * std::f64::consts::PI * (corner - x))
}

fn rounded_area(w: f64, h: f64, corner: f64) -> f64 {
    w * h - (4.0 - std::f64::consts::PI) * corner * corner
}

fn assert_rel(got: f64, want: f64, what: &str) {
    let rel = (got - want).abs() / want;
    assert!(
        rel < 1e-4,
        "{what}: volume {got} vs analytic {want} (rel {rel:.2e})"
    );
}

#[test]
fn rounded_box_top_rim_fillets_as_one_chain() {
    for r in [1.0, 3.0, 7.0] {
        let source = format!(
            r#"
let m = model();
let b = m.profile(xy()).rect(60.0, 40.0).sketch("block").extrude(blind(20.0));
b = b.fillet(edges().straight().direction("z"), 8.0, #{{ label: "corners" }});
b = b.fillet(edges().at_height(20.0), {r:?}, #{{ label: "top_rim" }});
m.assert_watertight();
m
"#
        );
        let want = rounded_area(60.0, 40.0, 8.0) * 20.0 - rim_change(60.0, 40.0, 8.0, r);
        assert_rel(served_volume(&source), want, &format!("convex rim r = {r}"));
    }
}

#[test]
fn rounded_pocket_floor_rim_fillets_as_one_concave_chain() {
    for r in [1.0, 2.0, 5.0] {
        let source = format!(
            r#"
let m = model();
let b = m.profile(xy()).rect(80.0, 60.0).sketch("block").extrude(blind(20.0));
let pocket = m.tool(|t| {{
    let p = t.profile(xy_at(8.0)).rect(50.0, 30.0).sketch("pocket").extrude(blind(20.0));
    p.fillet(edges().straight().direction("z"), 6.0, #{{ label: "pocket_corners" }});
}});
b = b.subtract(pocket, #{{ label: "pocket" }});
b = b.fillet(edges().at_height(8.0), {r:?}, #{{ label: "floor_rim" }});
m.assert_watertight();
m
"#
        );
        let pocket = rounded_area(50.0, 30.0, 6.0) * 12.0;
        let want = 80.0 * 60.0 * 20.0 - pocket + rim_change(50.0, 30.0, 6.0, r);
        assert_rel(
            served_volume(&source),
            want,
            &format!("concave rim r = {r}"),
        );
    }
}

/// The r40 × r30 cylinder tee's saddle, split into three G1 edges at
/// the two cylinders' seam crossings, fillets as one closed chain
/// (PK-729/PK-734): every strip row is exact (rows on their own solved
/// ball centres) and the walker refines its stations until the rails
/// hold the supports. No closed form for the volume; the concave
/// fillet's added section scales as r² along a fixed saddle, so
/// (V − V₀) / r² must agree across radii.
#[test]
fn cylinder_tee_saddle_fillets_as_one_chain() {
    // r30 and r20 branches into the r40 main (the r20 saddle's fitted
    // edges meet at its seam junctions with 1.2e-6 tangent noise).
    for (rb, radii) in [(30.0, vec![2.0, 6.0, 10.0]), (20.0, vec![2.0, 6.0])] {
        let tee = |fillet: Option<f64>| {
            let line = fillet
                .map(|r| {
                    format!(
                        r#"b = b.fillet(edges().between("cylinder", "cylinder"), {r:?}, #{{ label: "saddle" }});"#
                    )
                })
                .unwrap_or_default();
            format!(
                r#"
let m = model();
let main = m.profile(xy()).circle(40.0).sketch("main").extrude(blind(120.0));
let branch = m.tool(|t| t.profile(yz()).circle({rb:?}).at(0.0, 60.0).sketch("branch").extrude(blind(100.0)));
let b = main.union(branch, #{{ label: "tee" }});
{line}
m.assert_watertight();
m
"#
            )
        };
        let v0 = served_volume(&tee(None));
        let ratios: Vec<f64> = radii
            .iter()
            .map(|&r| (served_volume(&tee(Some(r))) - v0) / (r * r))
            .collect();
        let (lo, hi) = ratios
            .iter()
            .fold((f64::MAX, f64::MIN), |(a, b), &x| (a.min(x), b.max(x)));
        assert!(
            lo > 0.0,
            "rb {rb}: a concave fillet adds material: {ratios:?}"
        );
        assert!(
            (hi - lo) / hi < 0.03,
            "rb {rb}: section scaling inconsistent: {ratios:?}"
        );
    }
}
