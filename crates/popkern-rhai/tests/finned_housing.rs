//! Finned housings (PK-730): rectangular fins unioned onto a cylinder
//! housing, their straight root edges filleted. The housing face wraps
//! round each fin, so every fillet end meets the support at a REFLEX
//! corner (the fin's bottom/top root vertex) and settles on a terminal
//! disc in the fin's cap plane.
//!
//! Volumes are analytic. A fin of half-width `w` on a housing of radius
//! `R` adds its cross-section outside the circle; each root fillet of
//! radius `r` adds the section bounded by the circle, the fin side and
//! the rolling-ball arc (centre at `R + r` from the axis, `w + r` from
//! the fin's mid-plane). `body_volume` meshes curved faces (its chord
//! error on the housing is ~20, above one fillet's ~12), so each case
//! asserts the fillet's volume CHANGE against the same body unfilleted
//! — the housing's chord error cancels — and the unfilleted body
//! against its analytic volume at the mesh's resolution. The OCC
//! oracle measured both rings exact (985485.40 analytic) on STEP.

use popkern_rhai::ScriptHost;

const R: f64 = 55.0;
const W: f64 = 1.5;
const TIP: f64 = 63.0;
const ROOT: f64 = 53.0;
const Z0: f64 = 20.0;
const H: f64 = 60.0;

fn served_volume(source: &str) -> f64 {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(source, "fins.rhai")
        .unwrap_or_else(|e| panic!("{e}\n{}", e.payload));
    let mut model = outcome.model.expect("model");
    let part = model.to_part().expect("part");
    let brep = part.brep();
    popkern_measure::body_volume(brep, &brep.tolerance).expect("volume")
}

/// Fin cross-section outside the housing circle.
fn fin_area() -> f64 {
    // ∫_{-w}^{w} (TIP − √(R² − y²)) dy
    let arc = W * (R * R - W * W).sqrt() + R * R * (W / R).asin();
    2.0 * W * TIP - arc
}

/// Concave root section between the circle, the fin side and the arc.
fn root_section(r: f64) -> f64 {
    let c = (((R + r).powi(2) - (W + r).powi(2)).sqrt(), W + r);
    let t1 = (c.0 * R / (R + r), c.1 * R / (R + r));
    let t2 = (c.0, W);
    let p = ((R * R - W * W).sqrt(), W);
    let quad = [p, t2, c, t1];
    let shoelace: f64 = (0..4)
        .map(|i| {
            let (a, b) = (quad[i], quad[(i + 1) % 4]);
            a.0 * b.1 - b.0 * a.1
        })
        .sum::<f64>()
        .abs()
        * 0.5;
    let sweep = ((t2.1 - c.1).atan2(t2.0 - c.0) - (t1.1 - c.1).atan2(t1.0 - c.0)).abs();
    let chord = ((t1.0 - p.0).powi(2) + (t1.1 - p.1).powi(2)).sqrt();
    let phi = 2.0 * (chord / (2.0 * R)).asin();
    shoelace - 0.5 * r * r * sweep - 0.5 * R * R * (phi - phi.sin())
}

fn housing() -> f64 {
    std::f64::consts::PI * R * R * 100.0
}

fn assert_rel(got: f64, want: f64, rel_tol: f64, what: &str) {
    let rel = (got - want).abs() / want.abs();
    assert!(
        rel < rel_tol,
        "{what}: volume {got} vs analytic {want} (diff {:.3}, rel {rel:.2e})",
        got - want
    );
}

/// Fillet volume change: served body minus the same script with the
/// fillet line removed.
fn fillet_change(source: &str) -> (f64, f64) {
    let unfilleted: String = source
        .lines()
        .filter(|l| !l.contains(".fillet("))
        .collect::<Vec<_>>()
        .join("\n");
    let base = served_volume(&unfilleted);
    (base, served_volume(source) - base)
}

/// Corners of a fin at polar angle `a` (degrees), anticlockwise.
fn fin_corners(a: f64) -> [[f64; 2]; 4] {
    let (s, c) = a.to_radians().sin_cos();
    [
        [ROOT * c + W * s, ROOT * s - W * c],
        [TIP * c + W * s, TIP * s - W * c],
        [TIP * c - W * s, TIP * s + W * c],
        [ROOT * c - W * s, ROOT * s + W * c],
    ]
}

/// Root vertex of the fin side at `+w` (`side = 1`) or `−w` (`−1`), at `z`.
fn root_vertex(a: f64, side: f64, z: f64) -> String {
    let (s, c) = a.to_radians().sin_cos();
    let x = (R * R - W * W).sqrt();
    let y = side * W;
    format!("[{}, {}, {z}]", x * c - y * s, x * s + y * c)
}

fn single_fin(a: f64, sides: &[f64]) -> String {
    let pts: Vec<String> = fin_corners(a)
        .iter()
        .map(|p| format!("[{}, {}]", p[0], p[1]))
        .collect();
    let selector = sides
        .iter()
        .map(|&s| {
            format!(
                "edges().straight().direction(\"z\").at_vertex({})",
                root_vertex(a, s, Z0)
            )
        })
        .reduce(|acc, e| format!("{acc}.or({e})"))
        .unwrap_or_default();
    format!(
        r#"
let m = model();
let h = m.profile(xy()).circle({R:?}).sketch("housing").extrude(blind(100.0));
let fin = m.tool(|t| t.profile(xy_at({Z0:?})).polygon([{}]).sketch("fin").extrude(blind({H:?})));
h = h.union(fin, #{{ label: "fin" }});
h = h.fillet({selector}, 1.0, #{{ label: "root" }});
m.assert_watertight();
m
"#,
        pts.join(", ")
    )
}

/// The ring is laid out in the script (rhai trig), exactly as a user
/// writes it: the corner bits decide whether the union serves on the
/// cellular pipeline, and 16 seam-aligned fins take the graph-merge
/// fallback whose closed top rim once shipped a direction lie that
/// turned the seam fin inside out under the fillet (PK-730).
fn fin_ring(fins: usize, phase: f64, r: f64) -> String {
    format!(
        r#"
let pi = 3.141592653589793;
let m = model();
let h = m.profile(xy()).circle({R:?}).sketch("housing").extrude(blind(100.0));
let pts = [];
for k in 0..{fins} {{
    let a = 2.0 * pi * (k + {phase:?}) / {fins};
    let ca = cos(a); let sa = sin(a);
    pts.push([{ROOT:?} * ca + {W:?} * sa, {ROOT:?} * sa - {W:?} * ca]);
    pts.push([{TIP:?} * ca + {W:?} * sa, {TIP:?} * sa - {W:?} * ca]);
    pts.push([{TIP:?} * ca - {W:?} * sa, {TIP:?} * sa + {W:?} * ca]);
    pts.push([{ROOT:?} * ca - {W:?} * sa, {ROOT:?} * sa + {W:?} * ca]);
}}
let ring = m.tool(|t| t.profile(xy_at({Z0:?})).polygon(pts).sketch("fins").extrude(blind({H:?})));
h = h.union(ring, #{{ label: "fins" }});
h = h.fillet(edges().straight().direction("z").between("plane", "cylinder"), {r}, #{{ label: "fin_roots" }});
m.assert_watertight();
m
"#
    )
}

#[test]
fn a_fin_root_fillets_on_and_off_the_housing_seam() {
    let fin = housing() + fin_area() * H;
    for a in [0.0, 30.0] {
        for sides in [&[1.0][..], &[1.0, -1.0][..]] {
            let (base, change) = fillet_change(&single_fin(a, sides));
            let what = format!("{} root edge(s), fin at {a}°", sides.len());
            assert_rel(base, fin, 1e-4, &what);
            assert_rel(
                change,
                sides.len() as f64 * root_section(1.0) * H,
                1e-2,
                &what,
            );
        }
    }
}

#[test]
fn finned_housings_fillet_every_root_edge() {
    // 24 fins on and off the seam serve cellular; 16 seam-aligned
    // fins serve through the graph-merge fallback.
    for (fins, phase) in [(24, 0.0), (24, 0.5), (16, 0.0)] {
        let (base, change) = fillet_change(&fin_ring(fins, phase, 1.0));
        let what = format!("{fins} fins, phase {phase:?} pitch");
        let n = fins as f64;
        assert_rel(base, housing() + n * fin_area() * H, 1e-4, &what);
        assert_rel(change, 2.0 * n * root_section(1.0) * H, 1e-2, &what);
    }
}

/// The union alone, 16 and 24 seam-aligned fins: it serves on the
/// CELLULAR pipeline (no recorded boolean fallback), every closed rim's
/// canonical half-edge traverses its curve forward, and every face's
/// flag agrees with its loop winding. The graph-merge assembler once
/// created a closed edge from whichever side the first fragment needed
/// while the curve belonged to the cap: the top rim lied, the cap
/// served inward and the wall's loop wound 4π (PK-730). The fallback
/// itself came from Bezier clipping bounding a root by the coefficient
/// polygon instead of its hull — one fin side × circle crossing was
/// clipped away (PK-739).
#[test]
fn seam_fin_ring_unions_serve_cellular_with_honest_rims() {
    for fins in [16, 24] {
        let source: String = fin_ring(fins, 0.0, 1.0)
            .lines()
            .filter(|l| !l.contains(".fillet("))
            .collect::<Vec<_>>()
            .join("\n");
        let mut host = ScriptHost::new();
        let outcome = host
            .run(&source, "fins.rhai")
            .unwrap_or_else(|e| panic!("{e}\n{}", e.payload));
        let mut model = outcome.model.expect("model");
        let fallbacks: Vec<String> = model
            .report()
            .fallbacks()
            .map(|(f, r)| format!("{}: {r:?}", f.label))
            .collect();
        assert!(fallbacks.is_empty(), "{fins} fins fell back: {fallbacks:?}");
        let part = model.to_part().expect("part");
        let brep = part.brep();
        let violated = brep.validate_canonical_forward_all();
        assert!(
            violated.is_empty(),
            "{fins} fins: canonical-forward {violated:?}"
        );
        let dishonest = brep.dishonest_winding_faces();
        assert!(
            dishonest.is_empty(),
            "{fins} fins: dishonest faces {dishonest:?}"
        );
    }
}

/// A fin's root ARC (cap plane ⊥ housing axis × housing): a torus strip
/// whose ends meet the housing at reflex corners and settle on terminal
/// discs in the fin's OBLIQUE side planes, bounded by the fitted section
/// `torus ∩ side plane` (PK-740). Its volume (0.644) sits under the
/// housing mesh's chord noise, so this pins the structure — one torus
/// strip, two planar discs, every edge on its faces; the OCC oracle
/// measured the STEP exact against `∬ 2ρ·asin(w/ρ) dA` (951773.65).
#[test]
fn a_fin_root_arc_fillets_with_oblique_terminal_discs() {
    let (s, c) = (30.0_f64.to_radians().sin(), 30.0_f64.to_radians().cos());
    let x = (R * R - W * W).sqrt();
    let corner = format!("[{}, {}, {Z0:?}]", x * c - W * s, x * s + W * c);
    let filleted = single_fin(30.0, &[1.0]).replace(
        &format!(
            "edges().straight().direction(\"z\").at_vertex({})",
            root_vertex(30.0, 1.0, Z0)
        ),
        &format!("edges().circular().at_vertex({corner}).at_height({Z0:?})"),
    );
    assert!(filleted.contains("circular()"), "selector rewrite");
    let unfilleted: String = filleted
        .lines()
        .filter(|l| !l.contains(".fillet("))
        .collect::<Vec<_>>()
        .join("\n");
    let brep_of = |source: &str| {
        let mut host = ScriptHost::new();
        let outcome = host
            .run(source, "arc.rhai")
            .unwrap_or_else(|e| panic!("{e}\n{}", e.payload));
        let mut model = outcome.model.expect("model");
        model.to_part().expect("part").brep().clone()
    };
    let (before, after) = (brep_of(&unfilleted), brep_of(&filleted));
    assert_eq!(
        after.faces.len(),
        before.faces.len() + 3,
        "strip + two discs"
    );
    let torus = after.faces.values().filter(|f| {
        f.surface.as_ref().is_some_and(|s| {
            matches!(
                s.geometry,
                popkern_core::geometry::surface::SurfaceGeometry::Torus { .. }
            )
        })
    });
    assert_eq!(torus.count(), 1, "one torus strip");
    assert!(after.validate().is_ok(), "{:?}", after.validate());
    assert!(
        after.validate_curve_surface().is_ok(),
        "{:?}",
        after.validate_curve_surface()
    );
}
